#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>

#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/window/blackman_harris_3term.hpp>
#include <qwqdsp/window/blackman_harris.hpp>

namespace pitch_shift_rt {

/**
 * @brief 实时移调器使用的相位重建算法
 */
enum class Algorithm {
    Pghi,
    PghiFlux,
    PghiSuperFlux,
    PhaseLock,
    TransientVocoder,
};

namespace detail {

static inline float wrapToPi(float phase) noexcept {
    constexpr float kTwoPi = 2.0f * std::numbers::pi_v<float>;
    return phase - kTwoPi * std::round(phase / kTwoPi);
}

} // namespace detail

/**
 * @brief 使用固定容量缓冲的实时频域移调器
 *
 * 处理器采用 2048 点 STFT 和八倍重叠。可变分析 hop 完成时间拉伸，
 * 固定合成 hop 保证 COLA，随后对连续 OLA 信号流式重采样完成移调。
 * 五种模式只切换相位传播和瞬态处理策略。初始化完成后音频处理路径
 * 不进行动态内存分配。
 */
class RealtimePitchShifter {
public:
    static constexpr std::size_t kFrameSize = 2048;
    static constexpr std::size_t kHopSize = 256;
    static constexpr std::size_t kNumBins = kFrameSize / 2 + 1;

    /**
     * @brief 初始化 FFT、窗函数和处理状态
     * @param[in] sample_rate 输入输出采样率，单位 Hz
     */
    void init(float sample_rate) {
        sample_rate_ = sample_rate;
        fft_.Init(kFrameSize);
        qwqdsp_window::BlackmanHarris::Window(window_, true);
        buildSuperFluxFilterBank();
        reset();
    }

    /**
     * @brief 设置音高移动量
     * @param[in] semitones 半音数，正值升高音高，负值降低音高
     */
    void setPitchShift(float semitones) {
        target_ratio_ = std::exp2(std::clamp(semitones, -12.0f, 12.0f) / 12.0f);
        if (first_frame_) {
            current_ratio_ = target_ratio_;
        }
    }

    /**
     * @brief 切换相位重建算法
     * @param[in] algorithm 要使用的算法
     */
    void setAlgorithm(Algorithm algorithm) {
        if (algorithm_ == algorithm) {
            return;
        }
        algorithm_ = algorithm;
        resetPhaseState();
    }

    /**
     * @brief 清空流式缓冲和全部算法状态
     */
    void reset() {
        input_ring_.fill(0.0f);
        synthesis_ring_.fill(0.0f);
        input_write_ = 0;
        collected_samples_ = 0;
        samples_since_frame_ = 0;
        samples_to_frame_ = 0.0f;
        analysis_started_ = false;
        synthesis_frame_start_ = 0;
        synthesis_available_ = 0;
        synthesis_cleared_until_ = 0;
        resample_position_ = 0.0;
        resampler_started_ = false;
        current_ratio_ = target_ratio_;
        resetPhaseState();
    }

    /**
     * @brief 处理一段单声道音频
     * @param[in] input 单声道浮点输入缓冲
     * @param[out] output 单声道浮点输出缓冲
     * @param[in] frame_count 缓冲中的采样数
     */
    void process(const float* input, float* output, std::size_t frame_count) {
        if (output == nullptr) {
            return;
        }
        if (input == nullptr) {
            std::fill_n(output, frame_count, 0.0f);
            return;
        }

        for (std::size_t i = 0; i < frame_count; ++i) {
            input_ring_[input_write_] = input[i];
            input_write_ = (input_write_ + 1) % kFrameSize;
            collected_samples_ = std::min(collected_samples_ + 1, kFrameSize);
            ++samples_since_frame_;

            if (!analysis_started_ && collected_samples_ == kFrameSize) {
                processFrame(kHopSize);
                analysis_started_ = true;
                samples_since_frame_ = 0;
                samples_to_frame_ = static_cast<float>(kHopSize) / current_ratio_;
            }
            else if (analysis_started_) {
                samples_to_frame_ -= 1.0f;
                if (samples_to_frame_ <= 0.0f) {
                    processFrame(std::max<std::size_t>(1, samples_since_frame_));
                    samples_since_frame_ = 0;
                    samples_to_frame_ += static_cast<float>(kHopSize) / current_ratio_;
                }
            }
            output[i] = renderOutputSample();
        }
    }

    /**
     * @brief 返回零移调时的稳态处理延迟
     * @return 零移调延迟采样数
     */
    static constexpr std::size_t latencySamples() noexcept {
        return kFrameSize - 1;
    }
private:
    static constexpr std::size_t kSynthesisRingSize = kFrameSize * 16;
    static constexpr std::size_t kFluxHistorySize = 32;
    static constexpr std::size_t kHeapCapacity = kNumBins * 2;
    static constexpr std::size_t kMaxSuperFluxBands = 256;
    static constexpr std::size_t kMaxSuperFluxWeights = kNumBins * 3;
    // 八倍重叠下周期三项 Blackman-Harris 窗平方和为常数。
    static constexpr float kOlaGain = 0.4073837f;

    struct HeapItem {
        float magnitude = 0.0f;
        std::size_t bin = 0;
        bool is_previous = false;
    };

    void resetPhaseState() {
        previous_re_.fill(0.0f);
        previous_im_.fill(0.0f);
        previous_analysis_phase_.fill(0.0f);
        previous_synthesis_phase_.fill(0.0f);
        previous_magnitude_.fill(0.0f);
        flux_history_.fill(0.0f);
        flux_history_count_ = 0;
        flux_history_write_ = 0;
        flux_cooldown_ = 0;
        transient_mean_ = 0.0f;
        transient_variance_ = 0.0f;
        transient_frames_ = 0;
        transient_cooldown_ = 0;
        transient_has_mean_ = false;
        superflux_current_.fill(0.0f);
        superflux_previous_max_.fill(0.0f);
        first_frame_ = true;
    }

    void processFrame(std::size_t analysis_hop) {
        current_ratio_ += 0.25f * (target_ratio_ - current_ratio_);

        const bool uses_pghi = algorithm_ == Algorithm::Pghi || algorithm_ == Algorithm::PghiFlux
                            || algorithm_ == Algorithm::PghiSuperFlux;
        if (uses_pghi) {
            constexpr std::size_t kHalfFrame = kFrameSize / 2;
            for (std::size_t i = 0; i < kHalfFrame; ++i) {
                fft_input_[i] = input_ring_[(input_write_ + kHalfFrame + i) % kFrameSize] * window_[kHalfFrame + i];
                fft_input_[kHalfFrame + i] = input_ring_[(input_write_ + i) % kFrameSize] * window_[i];
            }
        }
        else {
            for (std::size_t i = 0; i < kFrameSize; ++i) {
                fft_input_[i] = input_ring_[(input_write_ + i) % kFrameSize] * window_[i];
            }
        }
        fft_.FFT(fft_input_.data(), fft_output_.data());

        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            current_re_[bin] = fft_output_[2 * bin];
            current_im_[bin] = fft_output_[2 * bin + 1];
            current_magnitude_[bin] = std::hypot(current_re_[bin], current_im_[bin]);
            current_analysis_phase_[bin] = std::atan2(current_im_[bin], current_re_[bin]);
        }

        bool reset_phase = first_frame_;
        if (algorithm_ == Algorithm::PghiSuperFlux) {
            reset_phase = detectSuperFluxTransient();
        }
        else if (!first_frame_) {
            if (algorithm_ == Algorithm::PghiFlux) {
                reset_phase = detectFluxTransient();
            }
            else if (algorithm_ == Algorithm::TransientVocoder) {
                reset_phase = detectVocoderTransient();
            }
        }

        if (reset_phase) {
            current_synthesis_phase_ = current_analysis_phase_;
        }
        else if (algorithm_ == Algorithm::Pghi || algorithm_ == Algorithm::PghiFlux
                 || algorithm_ == Algorithm::PghiSuperFlux) {
            propagatePghi(analysis_hop);
        }
        else {
            propagateVocoder(analysis_hop);
            if (!reset_phase) {
                lockPhase();
            }
        }

        synthesizeSpectrum();
        fft_.IFFT(fft_output_.data(), fft_input_.data());
        for (std::size_t i = 0; i < kFrameSize; ++i) {
            const std::size_t index = static_cast<std::size_t>((synthesis_frame_start_ + i) % kSynthesisRingSize);
            const std::size_t source = uses_pghi ? (i + kFrameSize / 2) % kFrameSize : i;
            synthesis_ring_[index] += fft_input_[source] * window_[i] * kOlaGain;
        }
        // 新合成帧在 hop 边界处的窗值为零，因此该边界样本已完整可读。
        synthesis_available_ = synthesis_frame_start_ + kHopSize + 1;
        synthesis_frame_start_ += kHopSize;

        previous_re_ = current_re_;
        previous_im_ = current_im_;
        previous_analysis_phase_ = current_analysis_phase_;
        previous_synthesis_phase_ = current_synthesis_phase_;
        previous_magnitude_ = current_magnitude_;
        first_frame_ = false;
    }

    void propagateVocoder(std::size_t analysis_hop) noexcept {
        const float analysis_scale =
            2.0f * std::numbers::pi_v<float> * static_cast<float>(analysis_hop) / static_cast<float>(kFrameSize);
        constexpr float kSynthesisScale =
            2.0f * std::numbers::pi_v<float> * static_cast<float>(kHopSize) / static_cast<float>(kFrameSize);
        const float stretch_ratio = static_cast<float>(kHopSize) / static_cast<float>(analysis_hop);
        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            const float expected = analysis_scale * static_cast<float>(bin);
            const float residual =
                detail::wrapToPi(current_analysis_phase_[bin] - previous_analysis_phase_[bin] - expected);
            current_synthesis_phase_[bin] = detail::wrapToPi(
                previous_synthesis_phase_[bin] + kSynthesisScale * static_cast<float>(bin) + stretch_ratio * residual);
        }
    }

    void propagatePghi(std::size_t analysis_hop) {
        processed_bins_.fill(0);
        heap_size_ = 0;
        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            heapPush({previous_magnitude_[bin], bin, true});
        }

        const float analysis_scale =
            2.0f * std::numbers::pi_v<float> * static_cast<float>(analysis_hop) / static_cast<float>(kFrameSize);
        constexpr float kSynthesisScale =
            2.0f * std::numbers::pi_v<float> * static_cast<float>(kHopSize) / static_cast<float>(kFrameSize);
        const float stretch_ratio = static_cast<float>(kHopSize) / static_cast<float>(analysis_hop);
        std::size_t processed_count = 0;
        while (processed_count < kNumBins && heap_size_ > 0) {
            const HeapItem item = heapPop();
            const std::size_t bin = item.bin;
            if (item.is_previous) {
                if (processed_bins_[bin] != 0) {
                    continue;
                }
                const float expected = analysis_scale * static_cast<float>(bin);
                const float residual =
                    detail::wrapToPi(current_analysis_phase_[bin] - previous_analysis_phase_[bin] - expected);
                current_synthesis_phase_[bin] =
                    detail::wrapToPi(previous_synthesis_phase_[bin] + kSynthesisScale * static_cast<float>(bin)
                                     + stretch_ratio * residual);
                processed_bins_[bin] = 1;
                ++processed_count;
                heapPush({current_magnitude_[bin], bin, false});
                continue;
            }

            if (bin + 1 < kNumBins && processed_bins_[bin + 1] == 0) {
                const float delta = detail::wrapToPi(current_analysis_phase_[bin + 1] - current_analysis_phase_[bin]);
                current_synthesis_phase_[bin + 1] =
                    detail::wrapToPi(current_synthesis_phase_[bin] + stretch_ratio * delta);
                processed_bins_[bin + 1] = 1;
                ++processed_count;
                heapPush({current_magnitude_[bin + 1], bin + 1, false});
            }
            if (bin > 0 && processed_bins_[bin - 1] == 0) {
                const float delta = detail::wrapToPi(current_analysis_phase_[bin - 1] - current_analysis_phase_[bin]);
                current_synthesis_phase_[bin - 1] =
                    detail::wrapToPi(current_synthesis_phase_[bin] + stretch_ratio * delta);
                processed_bins_[bin - 1] = 1;
                ++processed_count;
                heapPush({current_magnitude_[bin - 1], bin - 1, false});
            }
        }
    }

    void lockPhase() {
        peak_count_ = 0;
        float maximum = 0.0f;
        for (float magnitude : current_magnitude_) {
            maximum = std::max(maximum, magnitude);
        }
        const float floor = maximum * 0.015f;
        for (std::size_t bin = 1; bin + 1 < kNumBins; ++bin) {
            if (current_magnitude_[bin] >= floor && current_magnitude_[bin] >= current_magnitude_[bin - 1]
                && current_magnitude_[bin] > current_magnitude_[bin + 1]) {
                peak_bins_[peak_count_++] = bin;
            }
        }
        if (peak_count_ == 0) {
            peak_bins_[peak_count_++] = static_cast<std::size_t>(
                std::max_element(current_magnitude_.begin(), current_magnitude_.end()) - current_magnitude_.begin());
        }

        for (std::size_t peak_index = 0; peak_index < peak_count_; ++peak_index) {
            const std::size_t peak = peak_bins_[peak_index];
            const std::size_t start = peak_index == 0 ? 0 : (peak_bins_[peak_index - 1] + peak) / 2 + 1;
            const std::size_t end =
                peak_index + 1 == peak_count_ ? kNumBins - 1 : (peak + peak_bins_[peak_index + 1]) / 2;
            const float offset = current_synthesis_phase_[peak] - current_analysis_phase_[peak];
            const float lock_floor = std::max(1.0e-10f, current_magnitude_[peak] * 0.03f);
            for (std::size_t bin = start; bin <= end; ++bin) {
                if (current_magnitude_[bin] >= lock_floor) {
                    current_synthesis_phase_[bin] = detail::wrapToPi(current_analysis_phase_[bin] + offset);
                }
            }
        }
    }

    bool detectFluxTransient() {
        float positive_difference = 0.0f;
        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            positive_difference += std::max(0.0f, current_magnitude_[bin] - previous_magnitude_[bin]);
        }
        const float flux = positive_difference / static_cast<float>(kNumBins);

        float baseline = 0.0f;
        if (flux_history_count_ > 0) {
            flux_scratch_ = flux_history_;
            const auto middle = flux_scratch_.begin() + static_cast<std::ptrdiff_t>(flux_history_count_ / 2);
            std::nth_element(flux_scratch_.begin(), middle,
                             flux_scratch_.begin() + static_cast<std::ptrdiff_t>(flux_history_count_));
            baseline = *middle;
        }
        const bool transient = flux_cooldown_ == 0 && flux > std::max(4.0f * baseline, 0.02f);
        flux_cooldown_ = transient ? 6 : (flux_cooldown_ > 0 ? flux_cooldown_ - 1 : 0);
        if (!transient) {
            flux_history_[flux_history_write_] = flux;
            flux_history_write_ = (flux_history_write_ + 1) % kFluxHistorySize;
            flux_history_count_ = std::min(flux_history_count_ + 1, kFluxHistorySize);
        }
        return transient;
    }

    bool detectSuperFluxTransient() {
        if (superflux_band_count_ == 0) {
            return false;
        }

        float flux = 0.0f;
        constexpr float kReferenceFrameSize = 4096.0f;
        const float scale = kReferenceFrameSize / static_cast<float>(kFrameSize);
        for (std::size_t band = 0; band < superflux_band_count_; ++band) {
            float magnitude = 0.0f;
            for (std::size_t i = 0; i < superflux_band_sizes_[band]; ++i) {
                const std::size_t bin = superflux_band_starts_[band] + i;
                magnitude += current_magnitude_[bin] * superflux_weights_[superflux_band_offsets_[band] + i];
            }
            superflux_current_[band] = std::log10(magnitude * scale + 1.0f);
            if (!first_frame_) {
                flux += std::max(0.0f, superflux_current_[band] - superflux_previous_max_[band]);
            }
        }

        for (std::size_t band = 0; band < superflux_band_count_; ++band) {
            float maximum = superflux_current_[band];
            if (band > 0) {
                maximum = std::max(maximum, superflux_current_[band - 1]);
            }
            if (band + 1 < superflux_band_count_) {
                maximum = std::max(maximum, superflux_current_[band + 1]);
            }
            superflux_previous_max_[band] = maximum;
        }
        return !first_frame_ && flux / static_cast<float>(superflux_band_count_) > 0.03f;
    }

    bool detectVocoderTransient() {
        float flux = 0.0f;
        float energy = 0.0f;
        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            const float weight = 0.5f + 0.5f * static_cast<float>(bin) / static_cast<float>(kNumBins - 1);
            const float current_log = std::log1p(current_magnitude_[bin]);
            flux += std::max(0.0f, current_log - std::log1p(previous_magnitude_[bin]));
            energy += weight * current_log;
        }
        const float normalized_flux = flux / std::max(energy, 1.0e-9f);
        const float deviation = std::sqrt(transient_variance_);
        const bool transient = transient_frames_ > 4 && transient_cooldown_ == 0
                            && normalized_flux > transient_mean_ + 1.5f * std::max(0.07f, deviation)
                            && normalized_flux > transient_mean_ * 1.35f;

        const float alpha = transient ? 0.3f : 0.12f;
        if (!transient_has_mean_) {
            transient_mean_ = normalized_flux;
            transient_has_mean_ = true;
        }
        else {
            const float delta = normalized_flux - transient_mean_;
            transient_mean_ += alpha * delta;
            transient_variance_ = (1.0f - alpha) * (transient_variance_ + alpha * delta * delta);
        }
        transient_cooldown_ = transient ? 1 : (transient_cooldown_ > 0 ? transient_cooldown_ - 1 : 0);
        ++transient_frames_;
        return transient;
    }

    void synthesizeSpectrum() {
        fft_output_.fill(0.0f);
        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            fft_output_[2 * bin] = current_magnitude_[bin] * std::cos(current_synthesis_phase_[bin]);
            fft_output_[2 * bin + 1] = current_magnitude_[bin] * std::sin(current_synthesis_phase_[bin]);
        }
        fft_output_[1] = 0.0f;
        fft_output_[2 * (kNumBins - 1) + 1] = 0.0f;
    }

    float renderOutputSample() {
        constexpr std::uint64_t kResampleStart = kFrameSize - kHopSize;
        if (!resampler_started_) {
            if (synthesis_available_ <= kResampleStart + 1) {
                return 0.0f;
            }
            resample_position_ = static_cast<double>(kResampleStart);
            while (synthesis_cleared_until_ < kResampleStart) {
                synthesis_ring_[static_cast<std::size_t>(synthesis_cleared_until_ % kSynthesisRingSize)] = 0.0f;
                ++synthesis_cleared_until_;
            }
            resampler_started_ = true;
        }

        const auto lower = static_cast<std::uint64_t>(resample_position_);
        if (lower + 1 >= synthesis_available_) {
            return 0.0f;
        }
        const float fraction = static_cast<float>(resample_position_ - static_cast<double>(lower));
        const float first = synthesis_ring_[static_cast<std::size_t>(lower % kSynthesisRingSize)];
        const float second = synthesis_ring_[static_cast<std::size_t>((lower + 1) % kSynthesisRingSize)];
        const float output = first + fraction * (second - first);

        resample_position_ += static_cast<double>(current_ratio_);
        const auto clear_until = static_cast<std::uint64_t>(resample_position_);
        while (synthesis_cleared_until_ < clear_until) {
            synthesis_ring_[static_cast<std::size_t>(synthesis_cleared_until_ % kSynthesisRingSize)] = 0.0f;
            ++synthesis_cleared_until_;
        }
        return output;
    }

    void buildSuperFluxFilterBank() {
        constexpr float kMinimumFrequency = 27.5f;
        constexpr float kMaximumFrequency = 16000.0f;
        constexpr int kBandsPerOctave = 24;
        std::array<std::size_t, kMaxSuperFluxBands + 2> centers{};
        std::size_t center_count = 0;
        const float bin_hz = sample_rate_ / static_cast<float>(kFrameSize);
        const float maximum_frequency = std::min(kMaximumFrequency, sample_rate_ * 0.5f * 0.999f);

        for (int index = 0; center_count < centers.size(); ++index) {
            const float frequency =
                kMinimumFrequency * std::exp2(static_cast<float>(index) / static_cast<float>(kBandsPerOctave));
            if (frequency > maximum_frequency) {
                break;
            }
            const std::size_t bin =
                std::min<std::size_t>(static_cast<std::size_t>(std::lround(frequency / bin_hz)), kNumBins - 1);
            if (center_count == 0 || bin > centers[center_count - 1]) {
                centers[center_count++] = bin;
            }
        }

        superflux_band_count_ = 0;
        superflux_weight_count_ = 0;
        for (std::size_t center = 1; center + 1 < center_count; ++center) {
            const std::size_t low = centers[center - 1];
            const std::size_t middle = centers[center];
            const std::size_t high = centers[center + 1];
            const std::size_t size = high - low + 1;
            if (low >= middle || middle >= high || superflux_band_count_ == kMaxSuperFluxBands
                || superflux_weight_count_ + size > kMaxSuperFluxWeights) {
                continue;
            }

            const std::size_t band = superflux_band_count_++;
            superflux_band_starts_[band] = low;
            superflux_band_sizes_[band] = size;
            superflux_band_offsets_[band] = superflux_weight_count_;
            for (std::size_t bin = low; bin <= high; ++bin) {
                const float weight = bin <= middle ? static_cast<float>(bin - low) / static_cast<float>(middle - low)
                                                   : static_cast<float>(high - bin) / static_cast<float>(high - middle);
                superflux_weights_[superflux_weight_count_++] = weight;
            }
        }
    }

    void heapPush(HeapItem item) noexcept {
        if (heap_size_ == kHeapCapacity) {
            return;
        }
        std::size_t index = heap_size_++;
        heap_[index] = item;
        while (index > 0) {
            const std::size_t parent = (index - 1) / 2;
            if (heap_[parent].magnitude >= heap_[index].magnitude) {
                break;
            }
            std::swap(heap_[parent], heap_[index]);
            index = parent;
        }
    }

    HeapItem heapPop() noexcept {
        const HeapItem result = heap_[0];
        heap_[0] = heap_[--heap_size_];
        std::size_t index = 0;
        while (true) {
            const std::size_t left = index * 2 + 1;
            if (left >= heap_size_) {
                break;
            }
            const std::size_t right = left + 1;
            const std::size_t largest =
                right < heap_size_ && heap_[right].magnitude > heap_[left].magnitude ? right : left;
            if (heap_[index].magnitude >= heap_[largest].magnitude) {
                break;
            }
            std::swap(heap_[index], heap_[largest]);
            index = largest;
        }
        return result;
    }

    float sample_rate_ = 48000.0f;
    float target_ratio_ = 1.0f;
    float current_ratio_ = 1.0f;
    Algorithm algorithm_ = Algorithm::Pghi;
    bool first_frame_ = true;

    qwqdsp_spectral::RealFFT fft_;
    std::array<float, kFrameSize> window_{};
    std::array<float, kFrameSize> input_ring_{};
    std::array<float, kSynthesisRingSize> synthesis_ring_{};
    std::array<float, kFrameSize> fft_input_{};
    std::array<float, kFrameSize + 2> fft_output_{};

    std::array<float, kNumBins> previous_re_{};
    std::array<float, kNumBins> previous_im_{};
    std::array<float, kNumBins> previous_analysis_phase_{};
    std::array<float, kNumBins> previous_synthesis_phase_{};
    std::array<float, kNumBins> previous_magnitude_{};
    std::array<float, kNumBins> current_re_{};
    std::array<float, kNumBins> current_im_{};
    std::array<float, kNumBins> current_analysis_phase_{};
    std::array<float, kNumBins> current_synthesis_phase_{};
    std::array<float, kNumBins> current_magnitude_{};

    std::array<std::uint8_t, kNumBins> processed_bins_{};
    std::array<HeapItem, kHeapCapacity> heap_{};
    std::size_t heap_size_ = 0;
    std::array<std::size_t, kNumBins> peak_bins_{};
    std::size_t peak_count_ = 0;

    std::array<float, kFluxHistorySize> flux_history_{};
    std::array<float, kFluxHistorySize> flux_scratch_{};
    std::size_t flux_history_count_ = 0;
    std::size_t flux_history_write_ = 0;
    std::size_t flux_cooldown_ = 0;
    float transient_mean_ = 0.0f;
    float transient_variance_ = 0.0f;
    std::size_t transient_frames_ = 0;
    std::size_t transient_cooldown_ = 0;
    bool transient_has_mean_ = false;
    std::array<std::size_t, kMaxSuperFluxBands> superflux_band_starts_{};
    std::array<std::size_t, kMaxSuperFluxBands> superflux_band_sizes_{};
    std::array<std::size_t, kMaxSuperFluxBands> superflux_band_offsets_{};
    std::array<float, kMaxSuperFluxWeights> superflux_weights_{};
    std::array<float, kMaxSuperFluxBands> superflux_current_{};
    std::array<float, kMaxSuperFluxBands> superflux_previous_max_{};
    std::size_t superflux_band_count_ = 0;
    std::size_t superflux_weight_count_ = 0;

    std::size_t input_write_ = 0;
    std::size_t collected_samples_ = 0;
    std::size_t samples_since_frame_ = 0;
    float samples_to_frame_ = 0.0f;
    bool analysis_started_ = false;
    std::uint64_t synthesis_frame_start_ = 0;
    std::uint64_t synthesis_available_ = 0;
    std::uint64_t synthesis_cleared_until_ = 0;
    double resample_position_ = 0.0;
    bool resampler_started_ = false;
};

} // namespace pitch_shift_rt
