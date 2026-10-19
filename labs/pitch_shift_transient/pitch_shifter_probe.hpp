#pragma once

// example/gui/audiofx/pitch_shift_rt/pitch_shifter_rt.hpp 的实验副本，
// 仅用于在 labs 中观察 DSPark 触发并对照不同瞬态处理策略；库内原文件保持不变。

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <vector>

#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/window/blackman_harris.hpp>
#include <qwqdsp/window/blackman_harris_3term.hpp>
#include <qwqdsp/window/hann.hpp>

namespace pitch_shift_rt {

enum class WindowType {
    Hann,
    BlackmanHarrisThreeTerm,
    BlackmanHarris
};

namespace detail {
template <WindowType kType>
struct WindowTraits;
template <>
struct WindowTraits<WindowType::Hann> {
    static constexpr std::size_t kOverlap = 4;
    static constexpr float kOlaGain = 2.0f / 3.0f;
    using Type = qwqdsp_window::Hann;
};
template <>
struct WindowTraits<WindowType::BlackmanHarrisThreeTerm> {
    static constexpr std::size_t kOverlap = 6;
    static constexpr float kOlaGain = 0.5431781f;
    using Type = qwqdsp_window::BlackmanHarrisThreeTerm;
};
template <>
struct WindowTraits<WindowType::BlackmanHarris> {
    static constexpr std::size_t kOverlap = 8;
    static constexpr float kOlaGain = 0.4845649f;
    using Type = qwqdsp_window::BlackmanHarris;
};
constexpr std::size_t nextPowerOfTwo(std::size_t value) noexcept {
    std::size_t result = 1;
    while (result < value)
        result <<= 1;
    return result;
}
} // namespace detail

/**
 * @brief 实时移调器使用的相位重建算法
 */
enum class Algorithm {
    Pghi,
    PghiDSPark,
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
 * 处理器采用加窗 STFT 和固定合成 hop 完成时间拉伸，
 * 随后对连续 OLA 信号流式重采样完成移调。
 * 四种模式只切换相位传播和瞬态处理策略。初始化完成后音频处理路径
 * 不进行动态内存分配。
 */
template <WindowType kWindowType = WindowType::Hann, std::size_t kHopSize = 512>
class RealtimePitchShifter {
public:
    static constexpr std::size_t kWindowSize = kHopSize * detail::WindowTraits<kWindowType>::kOverlap;
    static constexpr std::size_t kFrameSize = detail::nextPowerOfTwo(kWindowSize * 2);
    static constexpr std::size_t kNumBins = kFrameSize / 2 + 1;

    /**
     * @brief 初始化 FFT、窗函数和处理状态
     * @param[in] sample_rate 输入输出采样率，单位 Hz
     */
    void init(float sample_rate) {
        sample_rate_ = sample_rate;
        fft_.Init(kFrameSize);
        detector_fft_.Init(kDetectorFrameSize);
        detail::WindowTraits<kWindowType>::Type::Window(window_, true);
        qwqdsp_window::Hann::Window(detector_window_, true);
        buildDetectorFilterBank();
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
            input_write_ = (input_write_ + 1) % kWindowSize;
            collected_samples_ = std::min(collected_samples_ + 1, kWindowSize);
            ++samples_since_frame_;
            ++input_position_;

            if (!analysis_started_ && collected_samples_ == kWindowSize) {
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
        return kWindowSize - 1;
    }

    enum class ResetMode {
        kFull,
        kNone,
        kWeighted,
        kPeakSeeded,
    };

    struct ProbeConfig {
        ResetMode reset_mode = ResetMode::kFull;
        float min_flux = 0.0f;
        float flux_peak_ratio = 0.0f;
        float flux_peak_decay = 1.0f;
        std::size_t min_rising_bands = 0;
        float min_max_band_rise = 0.0f;
        bool require_low_band_rise = false;
        float low_band_limit_hz = 2000.0f;
        float weighted_floor = 0.0f;
        float magnitude_baseline_alpha = 0.1f;
        float peak_rise_floor_ratio = 0.20f;
        float peak_magnitude_floor_ratio = 0.01f;
    };

    void setProbeConfig(const ProbeConfig& config) noexcept {
        config_ = config;
    }

    struct TransientEvent {
        std::size_t frame_index = 0;
        std::size_t analysis_hop = 0;
        double time_seconds = 0.0;
        float flux = 0.0f;
        float median = 0.0f;
        float threshold = 0.0f;
        float phase_jump = 0.0f;
        float phase_jump_low = 0.0f;
        float phase_jump_mid = 0.0f;
        float phase_jump_high = 0.0f;
        std::size_t rising_band_count = 0;
        float max_band_rise = 0.0f;
        float rise_f0 = 0.0f;
        float rise_f1 = 0.0f;
    };

    struct DetectorFrame {
        std::size_t frame_index = 0;
        double time_seconds = 0.0;
        std::size_t analysis_hop = 0;
        float flux = 0.0f;
        float flux_peak = 0.0f;
        float median = 0.0f;
        float threshold = 0.0f;
        bool rise = false;
        bool guard = true;
        bool transient = false;
        std::size_t cooldown_before = 0;
    };

    const std::vector<TransientEvent>& getTransientEvents() const noexcept {
        return transient_events_;
    }

    const std::vector<DetectorFrame>& getDetectorTrace() const noexcept {
        return detector_trace_;
    }

    const std::vector<float>& getAppliedJumps() const noexcept {
        return applied_jumps_;
    }

    const std::vector<std::size_t>& getAppliedBinCounts() const noexcept {
        return applied_bin_counts_;
    }

    void clearTransientEvents() {
        transient_events_.clear();
        detector_trace_.clear();
        applied_jumps_.clear();
        applied_bin_counts_.clear();
    }
private:
    static constexpr std::size_t kSynthesisRingSize = kFrameSize * 16;
    static constexpr std::size_t kHeapCapacity = kNumBins * 2;
    // 瞬态检测器使用独立的固定几何，与主链窗函数和帧长无关。
    static constexpr std::size_t kDetectorWindowSize = kWindowSize < 2048 ? kWindowSize : 2048;
    static constexpr std::size_t kDetectorFrameSize = detail::nextPowerOfTwo(kDetectorWindowSize * 2);
    static constexpr std::size_t kDetectorNumBins = kDetectorFrameSize / 2 + 1;
    static constexpr std::size_t kMaxDetectorBands = 256;
    static constexpr std::size_t kMaxDetectorWeights = kDetectorNumBins * 3;
    static constexpr std::size_t kDsFluxHistorySize = 12;
    static constexpr float kDsFluxDelta = 0.02f;
    static constexpr float kDsFluxCooldownMs = 64.0f;
    // 八倍重叠下周期三项 Blackman-Harris 窗平方和为常数。
    static constexpr float kOlaGain = detail::WindowTraits<kWindowType>::kOlaGain;

    struct HeapItem {
        float magnitude = 0.0f;
        std::size_t bin = 0;
        bool is_previous = false;
    };

    void resetPhaseState() {
        previous_analysis_phase_.fill(0.0f);
        previous_synthesis_phase_.fill(0.0f);
        previous_magnitude_.fill(0.0f);
        transient_mean_ = 0.0f;
        transient_variance_ = 0.0f;
        transient_frames_ = 0;
        transient_cooldown_ = 0;
        transient_has_mean_ = false;
        dsflux_current_.fill(0.0f);
        dsflux_previous_max_.fill(0.0f);
        dsflux_history_.fill(0.0f);
        dsflux_scratch_.fill(0.0f);
        dsflux_history_write_ = 0;
        dsflux_history_count_ = 0;
        dsflux_cooldown_ = 0;
        first_frame_ = true;
    }

    void processFrame(std::size_t analysis_hop) {
        ++frame_index_;
        current_ratio_ += 0.25f * (target_ratio_ - current_ratio_);

        const bool uses_pghi = algorithm_ == Algorithm::Pghi || algorithm_ == Algorithm::PghiDSPark;
        fft_input_.fill(0.0f);
        if (uses_pghi) {
            constexpr std::size_t kHalfWindow = kWindowSize / 2;
            for (std::size_t i = 0; i < kHalfWindow; ++i) {
                fft_input_[i] = input_ring_[(input_write_ + kHalfWindow + i) % kWindowSize] * window_[kHalfWindow + i];
            }
            const std::size_t pad = kFrameSize - kHalfWindow;
            for (std::size_t i = 0; i < kHalfWindow; ++i) {
                fft_input_[pad + i] = input_ring_[(input_write_ + i) % kWindowSize] * window_[i];
            }
        }
        else {
            for (std::size_t i = 0; i < kWindowSize; ++i) {
                fft_input_[i] = input_ring_[(input_write_ + i) % kWindowSize] * window_[i];
            }
        }
        fft_.FFT(fft_input_.data(), fft_output_.data());

        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            current_re_[bin] = fft_output_[2 * bin];
            current_im_[bin] = fft_output_[2 * bin + 1];
            current_magnitude_[bin] = std::hypot(current_re_[bin], current_im_[bin]);
            current_analysis_phase_[bin] = std::atan2(current_im_[bin], current_re_[bin]);
        }

        if (first_frame_) {
            magnitude_baseline_ = current_magnitude_;
        }

        if (algorithm_ == Algorithm::PghiDSPark) {
            processDetectorSpectrum();
        }

        bool reset_phase = first_frame_;
        if (!first_frame_) {
            if (algorithm_ == Algorithm::PghiDSPark) {
                reset_phase = detectDsParkTransient(analysis_hop);
            }
            else if (algorithm_ == Algorithm::TransientVocoder) {
                reset_phase = detectVocoderTransient();
            }
        }

        if (first_frame_) {
            current_synthesis_phase_ = current_analysis_phase_;
        }
        else if (!reset_phase) {
            if (uses_pghi) {
                propagatePghi(analysis_hop);
            }
            else {
                propagateVocoder(analysis_hop);
                lockPhase();
            }
        }
        else if (uses_pghi) {
            if (config_.reset_mode == ResetMode::kFull) {
                applied_jumps_.push_back(measurePhaseDeltaFrom(previous_synthesis_phase_));
                applied_bin_counts_.push_back(kNumBins);
                current_synthesis_phase_ = current_analysis_phase_;
            }
            else if (config_.reset_mode == ResetMode::kWeighted) {
                propagatePghi(analysis_hop);
                applyWeightedReset();
            }
            else if (config_.reset_mode == ResetMode::kPeakSeeded) {
                propagatePghi(analysis_hop, true);
            }
            else {
                propagatePghi(analysis_hop);
                applied_jumps_.push_back(0.0f);
                applied_bin_counts_.push_back(0);
            }
        }
        else {
            if (config_.reset_mode == ResetMode::kFull) {
                current_synthesis_phase_ = current_analysis_phase_;
            }
            else {
                propagateVocoder(analysis_hop);
                lockPhase();
            }
        }

        synthesizeSpectrum();
        fft_.IFFT(fft_output_.data(), fft_input_.data());
        constexpr std::size_t kHalfWindow = kWindowSize / 2;
        for (std::size_t i = 0; i < kHalfWindow; ++i) {
            const std::size_t index = static_cast<std::size_t>((synthesis_frame_start_ + i) % kSynthesisRingSize);
            const std::size_t source = uses_pghi ? kFrameSize - kHalfWindow + i : i;
            synthesis_ring_[index] += fft_input_[source] * window_[i] * kOlaGain;
        }
        for (std::size_t i = 0; i < kHalfWindow; ++i) {
            const std::size_t index =
                static_cast<std::size_t>((synthesis_frame_start_ + kHalfWindow + i) % kSynthesisRingSize);
            const std::size_t source = uses_pghi ? i : kHalfWindow + i;
            synthesis_ring_[index] += fft_input_[source] * window_[kHalfWindow + i] * kOlaGain;
        }
        // 新合成帧在 hop 边界处的窗值为零，因此该边界样本已完整可读。
        synthesis_available_ = synthesis_frame_start_ + kHopSize + 1;
        synthesis_frame_start_ += kHopSize;

        if (!reset_phase && !first_frame_) {
            const float alpha = config_.magnitude_baseline_alpha;
            for (std::size_t bin = 0; bin < kNumBins; ++bin) {
                magnitude_baseline_[bin] =
                    (1.0f - alpha) * magnitude_baseline_[bin] + alpha * current_magnitude_[bin];
            }
        }

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

    void propagatePghi(std::size_t analysis_hop, bool seed_transient_peaks = false) {
        processed_bins_.fill(0);
        heap_size_ = 0;
        std::size_t processed_count = 0;
        if (seed_transient_peaks) {
            float max_rise = 0.0f;
            float max_magnitude = 0.0f;
            float weight_total = 0.0f;
            for (std::size_t bin = 0; bin < kNumBins; ++bin) {
                max_rise = std::max(max_rise, current_magnitude_[bin] - magnitude_baseline_[bin]);
                max_magnitude = std::max(max_magnitude, current_magnitude_[bin]);
                weight_total += current_magnitude_[bin];
            }
            std::size_t seeded_bins = 0;
            float applied_sum = 0.0f;
            if (max_rise > 0.0f) {
                const float rise_floor = max_rise * config_.peak_rise_floor_ratio;
                const float magnitude_floor = max_magnitude * config_.peak_magnitude_floor_ratio;
                for (std::size_t bin = 1; bin + 1 < kNumBins; ++bin) {
                    const float rise = current_magnitude_[bin] - magnitude_baseline_[bin];
                    if (rise < rise_floor || current_magnitude_[bin] < magnitude_floor) {
                        continue;
                    }
                    if (current_magnitude_[bin] < current_magnitude_[bin - 1]
                        || current_magnitude_[bin] <= current_magnitude_[bin + 1]) {
                        continue;
                    }
                    current_synthesis_phase_[bin] = current_analysis_phase_[bin];
                    processed_bins_[bin] = 1;
                    ++processed_count;
                    ++seeded_bins;
                    heapPush({current_magnitude_[bin], bin, false});
                    applied_sum +=
                        current_magnitude_[bin]
                        * std::fabs(detail::wrapToPi(current_analysis_phase_[bin] - previous_synthesis_phase_[bin]));
                }
            }
            applied_jumps_.push_back(weight_total > 0.0f ? applied_sum / weight_total : 0.0f);
            applied_bin_counts_.push_back(seeded_bins);
        }

        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            heapPush({previous_magnitude_[bin], bin, true});
        }

        const float analysis_scale =
            2.0f * std::numbers::pi_v<float> * static_cast<float>(analysis_hop) / static_cast<float>(kFrameSize);
        constexpr float kSynthesisScale =
            2.0f * std::numbers::pi_v<float> * static_cast<float>(kHopSize) / static_cast<float>(kFrameSize);
        const float stretch_ratio = static_cast<float>(kHopSize) / static_cast<float>(analysis_hop);
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

    void processDetectorSpectrum() {
        detector_input_.fill(0.0f);
        constexpr std::size_t kHalfWindow = kDetectorWindowSize / 2;
        for (std::size_t i = 0; i < kHalfWindow; ++i) {
            detector_input_[i] = input_ring_[(input_write_ + kHalfWindow + i) % kWindowSize]
                               * detector_window_[kHalfWindow + i];
        }
        constexpr std::size_t kPad = kDetectorFrameSize - kHalfWindow;
        for (std::size_t i = 0; i < kHalfWindow; ++i) {
            detector_input_[kPad + i] = input_ring_[(input_write_ + i) % kWindowSize] * detector_window_[i];
        }
        detector_fft_.FFT(detector_input_.data(), detector_output_.data());
        for (std::size_t bin = 0; bin < kDetectorNumBins; ++bin) {
            detector_magnitude_[bin] = std::hypot(detector_output_[2 * bin], detector_output_[2 * bin + 1]);
        }
    }

    struct BandStats {
        float max_rise = 0.0f;
        std::size_t rising_count = 0;
        float f0 = 0.0f;
        float f1 = 0.0f;
        bool has_low_rise = false;
    };

    BandStats computeBandStats(float low_limit_hz) const {
        BandStats stats;
        for (std::size_t band = 0; band < detector_band_count_; ++band) {
            stats.max_rise =
                std::max(stats.max_rise, std::max(0.0f, dsflux_current_[band] - dsflux_previous_max_[band]));
        }
        const float rise_floor = stats.max_rise * 0.25f;
        std::size_t first_band = detector_band_count_;
        std::size_t last_band = 0;
        const float hz_per_bin = sample_rate_ / static_cast<float>(kDetectorFrameSize);
        for (std::size_t band = 0; band < detector_band_count_; ++band) {
            const float rise = std::max(0.0f, dsflux_current_[band] - dsflux_previous_max_[band]);
            if (rise > 0.0f && rise >= rise_floor) {
                ++stats.rising_count;
                first_band = std::min(first_band, band);
                last_band = std::max(last_band, band);
                const float band_hz = static_cast<float>(detector_band_starts_[band]) * hz_per_bin;
                if (band_hz <= low_limit_hz) {
                    stats.has_low_rise = true;
                }
            }
        }
        if (stats.rising_count > 0) {
            stats.f0 = static_cast<float>(detector_band_starts_[first_band]) * hz_per_bin;
            stats.f1 = static_cast<float>(detector_band_starts_[last_band] + detector_band_sizes_[last_band] - 1)
                     * hz_per_bin;
        }
        return stats;
    }

    float measurePhaseDeltaFrom(const std::array<float, kNumBins>& from) const {
        float weighted_sum = 0.0f;
        float weight_total = 0.0f;
        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            const float magnitude = current_magnitude_[bin];
            weighted_sum += magnitude
                          * std::fabs(detail::wrapToPi(current_analysis_phase_[bin] - from[bin]));
            weight_total += magnitude;
        }
        return weight_total > 0.0f ? weighted_sum / weight_total : 0.0f;
    }

    void applyWeightedReset() {
        float max_rise = 0.0f;
        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            max_rise = std::max(max_rise, current_magnitude_[bin] - previous_magnitude_[bin]);
        }
        if (max_rise <= 0.0f) {
            applied_jumps_.push_back(0.0f);
            applied_bin_counts_.push_back(0);
            return;
        }
        float applied_sum = 0.0f;
        float weight_total = 0.0f;
        std::size_t touched_bins = 0;
        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            const float rise = current_magnitude_[bin] - previous_magnitude_[bin];
            float weight = std::clamp(rise / max_rise, 0.0f, 1.0f);
            weight = std::max(weight, config_.weighted_floor);
            if (weight > 0.0f) {
                ++touched_bins;
            }
            const float delta =
                detail::wrapToPi(current_analysis_phase_[bin] - current_synthesis_phase_[bin]);
            current_synthesis_phase_[bin] =
                detail::wrapToPi(current_synthesis_phase_[bin] + weight * delta);
            applied_sum += current_magnitude_[bin] * std::fabs(weight * delta);
            weight_total += current_magnitude_[bin];
        }
        applied_jumps_.push_back(weight_total > 0.0f ? applied_sum / weight_total : 0.0f);
        applied_bin_counts_.push_back(touched_bins);
    }

    void measurePhaseJump(TransientEvent& event) const {
        std::array<float, 3> weighted_sums{};
        std::array<float, 3> weights{};
        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            const float magnitude = current_magnitude_[bin];
            const float delta = std::fabs(
                detail::wrapToPi(current_analysis_phase_[bin] - previous_synthesis_phase_[bin]));
            const float hz = static_cast<float>(bin) * sample_rate_ / static_cast<float>(kFrameSize);
            const std::size_t band = hz < 500.0f ? 0 : (hz < 4000.0f ? 1 : 2);
            weighted_sums[band] += magnitude * delta;
            weights[band] += magnitude;
        }
        const auto mean = [](float sum, float weight) { return weight > 0.0f ? sum / weight : 0.0f; };
        event.phase_jump_low = mean(weighted_sums[0], weights[0]);
        event.phase_jump_mid = mean(weighted_sums[1], weights[1]);
        event.phase_jump_high = mean(weighted_sums[2], weights[2]);
        event.phase_jump = mean(weighted_sums[0] + weighted_sums[1] + weighted_sums[2],
                                weights[0] + weights[1] + weights[2]);
    }

    void recordTransientEvent(std::size_t analysis_hop, float flux, float median, float threshold,
                              const BandStats& stats) {
        TransientEvent event;
        event.frame_index = frame_index_;
        event.analysis_hop = analysis_hop;
        const std::uint64_t center =
            input_position_ > kWindowSize / 2 ? input_position_ - kWindowSize / 2 : 0;
        event.time_seconds = static_cast<double>(center) / static_cast<double>(sample_rate_);
        event.flux = flux;
        event.median = median;
        event.threshold = threshold;
        event.max_band_rise = stats.max_rise;
        event.rising_band_count = stats.rising_count;
        event.rise_f0 = stats.f0;
        event.rise_f1 = stats.f1;
        measurePhaseJump(event);
        transient_events_.push_back(event);
    }

    /**
     * @brief 用 Hann 窗独立分析谱的对数频带正向通量是否超过中值基准。
     */
    bool detectDsParkTransient(std::size_t analysis_hop) {
        if (detector_band_count_ == 0) {
            return false;
        }

        float flux = 0.0f;
        constexpr float kReferenceFrameSize = 2048.0f;
        const float scale = kReferenceFrameSize / static_cast<float>(kDetectorFrameSize);
        for (std::size_t band = 0; band < detector_band_count_; ++band) {
            float magnitude = 0.0f;
            for (std::size_t i = 0; i < detector_band_sizes_[band]; ++i) {
                const std::size_t bin = detector_band_starts_[band] + i;
                magnitude += detector_magnitude_[bin] * detector_weights_[detector_band_offsets_[band] + i];
            }
            dsflux_current_[band] = std::log10(magnitude * scale + 1.0f);
            flux += std::max(0.0f, dsflux_current_[band] - dsflux_previous_max_[band]);
        }
        flux /= static_cast<float>(detector_band_count_);
        const float peak_decay =
            std::pow(config_.flux_peak_decay, static_cast<float>(analysis_hop) / static_cast<float>(kHopSize));
        flux_peak_ = std::max(flux, flux_peak_ * peak_decay);

        float median = 0.0f;
        const float scaled_delta =
            kDsFluxDelta * static_cast<float>(analysis_hop) / static_cast<float>(kHopSize);
        bool rise = false;
        if (dsflux_history_count_ >= kDsFluxHistorySize) {
            dsflux_scratch_ = dsflux_history_;
            for (std::size_t i = 1; i < kDsFluxHistorySize; ++i) {
                const float value = dsflux_scratch_[i];
                std::size_t j = i;
                while (j > 0 && dsflux_scratch_[j - 1] > value) {
                    dsflux_scratch_[j] = dsflux_scratch_[j - 1];
                    --j;
                }
                dsflux_scratch_[j] = value;
            }
            // 阈值随分析 hop 缩放，保证不同变调量下的时间灵敏度一致。
            median = dsflux_scratch_[kDsFluxHistorySize / 2];
            rise = (flux - median) >= scaled_delta;
        }
        dsflux_history_[dsflux_history_write_] = flux;
        dsflux_history_write_ = (dsflux_history_write_ + 1) % kDsFluxHistorySize;
        dsflux_history_count_ = std::min(dsflux_history_count_ + 1, kDsFluxHistorySize);

        const bool candidate = dsflux_cooldown_ == 0 && rise;
        BandStats stats;
        const bool need_stats = candidate || config_.min_rising_bands > 0 || config_.min_max_band_rise > 0.0f
                             || config_.require_low_band_rise;
        if (need_stats) {
            stats = computeBandStats(config_.low_band_limit_hz);
        }
        bool guard = true;
        if (config_.min_flux > 0.0f && flux < config_.min_flux) {
            guard = false;
        }
        if (config_.flux_peak_ratio > 0.0f && flux < config_.flux_peak_ratio * flux_peak_) {
            guard = false;
        }
        if (config_.min_rising_bands > 0 && stats.rising_count < config_.min_rising_bands) {
            guard = false;
        }
        if (config_.min_max_band_rise > 0.0f && stats.max_rise < config_.min_max_band_rise) {
            guard = false;
        }
        if (config_.require_low_band_rise && !stats.has_low_rise) {
            guard = false;
        }
        const bool transient = candidate && guard;
        if (transient) {
            recordTransientEvent(analysis_hop, flux, median, scaled_delta, stats);
        }

        for (std::size_t band = 0; band < detector_band_count_; ++band) {
            float maximum = dsflux_current_[band];
            if (band > 0) {
                maximum = std::max(maximum, dsflux_current_[band - 1]);
            }
            if (band + 1 < detector_band_count_) {
                maximum = std::max(maximum, dsflux_current_[band + 1]);
            }
            dsflux_previous_max_[band] = maximum;
        }

        DetectorFrame trace;
        trace.frame_index = frame_index_;
        trace.analysis_hop = analysis_hop;
        const std::uint64_t center =
            input_position_ > kWindowSize / 2 ? input_position_ - kWindowSize / 2 : 0;
        trace.time_seconds = static_cast<double>(center) / static_cast<double>(sample_rate_);
        trace.flux = flux;
        trace.flux_peak = flux_peak_;
        trace.median = median;
        trace.threshold = scaled_delta;
        trace.rise = rise;
        trace.guard = guard;
        trace.transient = transient;
        trace.cooldown_before = dsflux_cooldown_;
        detector_trace_.push_back(trace);

        dsflux_cooldown_ = transient ? cooldownFramesForHop(kDsFluxCooldownMs, analysis_hop)
                                     : (dsflux_cooldown_ > 0 ? dsflux_cooldown_ - 1 : 0);
        return transient;
    }

    std::size_t cooldownFramesForHop(float milliseconds, std::size_t analysis_hop) const noexcept {
        const float samples = milliseconds * sample_rate_ / 1000.0f;
        const std::size_t hop = analysis_hop > 0 ? analysis_hop : 1;
        const std::size_t frames = static_cast<std::size_t>(samples / static_cast<float>(hop) + 0.5f);
        return frames < 1 ? 1 : frames;
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

    void buildDetectorFilterBank() {
        constexpr float kMinimumFrequency = 27.5f;
        constexpr float kMaximumFrequency = 16000.0f;
        constexpr int kBandsPerOctave = 24;
        std::array<std::size_t, kMaxDetectorBands + 2> centers{};
        std::size_t center_count = 0;
        const float bin_hz = sample_rate_ / static_cast<float>(kDetectorFrameSize);
        const float maximum_frequency = std::min(kMaximumFrequency, sample_rate_ * 0.5f * 0.999f);

        for (int index = 0; center_count < centers.size(); ++index) {
            const float frequency =
                kMinimumFrequency * std::exp2(static_cast<float>(index) / static_cast<float>(kBandsPerOctave));
            if (frequency > maximum_frequency) {
                break;
            }
            const std::size_t bin =
                std::min<std::size_t>(static_cast<std::size_t>(std::lround(frequency / bin_hz)), kDetectorNumBins - 1);
            if (center_count == 0 || bin > centers[center_count - 1]) {
                centers[center_count++] = bin;
            }
        }

        detector_band_count_ = 0;
        detector_weight_count_ = 0;
        for (std::size_t center = 1; center + 1 < center_count; ++center) {
            const std::size_t low = centers[center - 1];
            const std::size_t middle = centers[center];
            const std::size_t high = centers[center + 1];
            const std::size_t size = high - low + 1;
            if (low >= middle || middle >= high || detector_band_count_ == kMaxDetectorBands
                || detector_weight_count_ + size > kMaxDetectorWeights) {
                continue;
            }

            const std::size_t band = detector_band_count_++;
            detector_band_starts_[band] = low;
            detector_band_sizes_[band] = size;
            detector_band_offsets_[band] = detector_weight_count_;
            for (std::size_t bin = low; bin <= high; ++bin) {
                const float weight = bin <= middle ? static_cast<float>(bin - low) / static_cast<float>(middle - low)
                                                   : static_cast<float>(high - bin) / static_cast<float>(high - middle);
                detector_weights_[detector_weight_count_++] = weight;
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
    qwqdsp_spectral::RealFFT detector_fft_;
    std::array<float, kWindowSize> window_{};
    std::array<float, kWindowSize> input_ring_{};
    std::array<float, kSynthesisRingSize> synthesis_ring_{};
    std::array<float, kFrameSize> fft_input_{};
    std::array<float, kFrameSize + 2> fft_output_{};
    std::array<float, kDetectorWindowSize> detector_window_{};
    std::array<float, kDetectorFrameSize> detector_input_{};
    std::array<float, kDetectorFrameSize + 2> detector_output_{};
    std::array<float, kDetectorNumBins> detector_magnitude_{};

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

    float transient_mean_ = 0.0f;
    float transient_variance_ = 0.0f;
    std::size_t transient_frames_ = 0;
    std::size_t transient_cooldown_ = 0;
    bool transient_has_mean_ = false;
    std::array<std::size_t, kMaxDetectorBands> detector_band_starts_{};
    std::array<std::size_t, kMaxDetectorBands> detector_band_sizes_{};
    std::array<std::size_t, kMaxDetectorBands> detector_band_offsets_{};
    std::array<float, kMaxDetectorWeights> detector_weights_{};
    std::array<float, kMaxDetectorBands> dsflux_current_{};
    std::array<float, kMaxDetectorBands> dsflux_previous_max_{};
    std::array<float, kDsFluxHistorySize> dsflux_history_{};
    std::array<float, kDsFluxHistorySize> dsflux_scratch_{};
    std::size_t detector_band_count_ = 0;
    std::size_t detector_weight_count_ = 0;
    std::size_t dsflux_history_count_ = 0;
    std::size_t dsflux_history_write_ = 0;
    std::size_t dsflux_cooldown_ = 0;

    ProbeConfig config_{};
    float flux_peak_ = 0.0f;
    std::array<float, kNumBins> magnitude_baseline_{};

    std::size_t frame_index_ = 0;
    std::uint64_t input_position_ = 0;
    std::vector<TransientEvent> transient_events_;
    std::vector<DetectorFrame> detector_trace_;
    std::vector<float> applied_jumps_;
    std::vector<std::size_t> applied_bin_counts_;

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
