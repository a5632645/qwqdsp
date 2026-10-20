/// @brief 纯 PGHI phase vocoder（无瞬态检测）+ 同帧 FFT/YIN 音高检测

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numbers>
#include <span>
#include <vector>

#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/window/hann.hpp>

#include "autotune_common.hpp"

// PGHI 相位传播与时间拉伸/重采样结构参考
// example/gui/audiofx/pitch_shift_rt/pitch_shifter_rt.hpp；本文件只保留纯 PGHI 路径，
// 不包含任何瞬态检测。
//
// PGHI 的帧为了 OLA 做了半窗重排；直接用它的功率谱算 YIN 会让自相关与前缀能量
// 不一致，周期出现明显偏置。因此 YIN 对同一个 Hann 分析帧（同样的窗长/hop）再取
// 一次自然顺序的 FFT，用其功率谱做自相关和 CMNDF。

namespace autotune_lab {

constexpr std::size_t kPvHopSize = 512;
constexpr std::size_t kPvWindowSize = kPvHopSize * 4;
constexpr std::size_t kPvFrameSize = kPvWindowSize * 2;
constexpr std::size_t kPvNumBins = kPvFrameSize / 2 + 1;
constexpr std::size_t kPvSynthesisRingSize = kPvFrameSize * 16;
constexpr std::size_t kPvHeapCapacity = kPvNumBins * 2;
constexpr float kPvOlaGain = 2.0f / 3.0f;

class PghiAutotune {
public:
    static constexpr std::size_t LatencySamples() noexcept {
        return kPvWindowSize - 1;
    }

    void Init(float sample_rate) {
        sample_rate_ = sample_rate;
        fft_.Init(kPvFrameSize);
        yin_fft_.Init(kPvFrameSize);
        qwqdsp_window::Hann::Window(window_, true);
        Reset();
    }

    std::vector<float> Process(std::span<const float> input) {
        std::size_t const latency = LatencySamples();
        std::vector<float> out(input.size() + latency, 0.0f);
        if (!input.empty()) {
            ProcessBlock(input.data(), out.data(), input.size());
        }
        std::vector<float> zeros(latency, 0.0f);
        ProcessBlock(zeros.data(), out.data() + input.size(), latency);
        return std::vector<float>(out.begin() + static_cast<std::ptrdiff_t>(latency), out.end());
    }

private:
    struct HeapItem {
        float magnitude = 0.0f;
        std::size_t bin = 0;
        bool is_previous = false;
    };

    static float WrapToPi(float phase) noexcept {
        constexpr float kTwoPi = 2.0f * std::numbers::pi_v<float>;
        return phase - kTwoPi * std::round(phase / kTwoPi);
    }

    void Reset() {
        input_ring_.fill(0.0f);
        synthesis_ring_.fill(0.0f);
        previous_analysis_phase_.fill(0.0f);
        previous_synthesis_phase_.fill(0.0f);
        previous_magnitude_.fill(0.0f);
        current_re_.fill(0.0f);
        current_im_.fill(0.0f);
        current_analysis_phase_.fill(0.0f);
        current_synthesis_phase_.fill(0.0f);
        current_magnitude_.fill(0.0f);

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
        first_frame_ = true;
        target_ratio_ = 1.0f;
        current_ratio_ = 1.0f;
    }

    void ProcessBlock(const float* input, float* output, std::size_t count) {
        for (std::size_t i = 0; i < count; ++i) {
            input_ring_[input_write_] = input[i];
            input_write_ = (input_write_ + 1) % kPvWindowSize;
            collected_samples_ = std::min(collected_samples_ + 1, kPvWindowSize);
            ++samples_since_frame_;

            if (!analysis_started_ && collected_samples_ == kPvWindowSize) {
                ProcessFrame(kPvHopSize);
                analysis_started_ = true;
                samples_since_frame_ = 0;
                samples_to_frame_ = static_cast<float>(kPvHopSize) / current_ratio_;
            }
            else if (analysis_started_) {
                samples_to_frame_ -= 1.0f;
                if (samples_to_frame_ <= 0.0f) {
                    ProcessFrame(std::max<std::size_t>(1, samples_since_frame_));
                    samples_since_frame_ = 0;
                    samples_to_frame_ += static_cast<float>(kPvHopSize) / current_ratio_;
                }
            }
            output[i] = RenderOutputSample();
        }
    }

    void ProcessFrame(std::size_t analysis_hop) {
        fft_input_.fill(0.0f);
        constexpr std::size_t kHalfWindow = kPvWindowSize / 2;
        for (std::size_t i = 0; i < kHalfWindow; ++i) {
            fft_input_[i] = input_ring_[(input_write_ + kHalfWindow + i) % kPvWindowSize] * window_[kHalfWindow + i];
        }
        constexpr std::size_t kPad = kPvFrameSize - kHalfWindow;
        for (std::size_t i = 0; i < kHalfWindow; ++i) {
            fft_input_[kPad + i] = input_ring_[(input_write_ + i) % kPvWindowSize] * window_[i];
        }
        fft_.FFT(fft_input_.data(), fft_output_.data());

        for (std::size_t bin = 0; bin < kPvNumBins; ++bin) {
            current_re_[bin] = fft_output_[2 * bin];
            current_im_[bin] = fft_output_[2 * bin + 1];
            current_magnitude_[bin] = std::hypot(current_re_[bin], current_im_[bin]);
            current_analysis_phase_[bin] = std::atan2(current_im_[bin], current_re_[bin]);
        }

        float const f0 = DetectPitch();
        float target = 1.0f;
        if (f0 >= kMinF0 && f0 <= kMaxF0) {
            float const midi = 69.0f + 12.0f * std::log2(f0 / 440.0f);
            int const nearest = FindNearestNote(midi, kScaleMask);
            float const target_midi = midi + kCorrection * (static_cast<float>(nearest) - midi);
            target = std::exp2((target_midi - midi) / 12.0f);
            float const max_ratio = std::exp2(kMaxShiftSemitones / 12.0f);
            target = std::clamp(target, 1.0f / max_ratio, max_ratio);
        }
        target_ratio_ = target;
        if (first_frame_) {
            current_ratio_ = target_ratio_;
        }
        else {
            float const alpha = 1.0f - std::exp(-static_cast<float>(analysis_hop)
                                                / (sample_rate_ * kRetuneSeconds));
            current_ratio_ += alpha * (target_ratio_ - current_ratio_);
        }

        if (first_frame_) {
            current_synthesis_phase_ = current_analysis_phase_;
        }
        else {
            PropagatePghi(analysis_hop);
        }

        SynthesizeSpectrum();
        fft_.IFFT(fft_output_.data(), fft_input_.data());
        for (std::size_t i = 0; i < kHalfWindow; ++i) {
            const std::size_t index = (synthesis_frame_start_ + i) % kPvSynthesisRingSize;
            const std::size_t source = kPvFrameSize - kHalfWindow + i;
            synthesis_ring_[index] += fft_input_[source] * window_[i] * kPvOlaGain;
        }
        for (std::size_t i = 0; i < kHalfWindow; ++i) {
            const std::size_t index = (synthesis_frame_start_ + kHalfWindow + i) % kPvSynthesisRingSize;
            const std::size_t source = i;
            synthesis_ring_[index] += fft_input_[source] * window_[kHalfWindow + i] * kPvOlaGain;
        }
        synthesis_available_ = synthesis_frame_start_ + kPvHopSize + 1;
        synthesis_frame_start_ += kPvHopSize;

        previous_analysis_phase_ = current_analysis_phase_;
        previous_synthesis_phase_ = current_synthesis_phase_;
        previous_magnitude_ = current_magnitude_;
        first_frame_ = false;
    }

    /// @brief 对同一 Hann 帧做 YIN；PGHI 的半窗重排会偏置周期检测，因此这里用自然顺序的帧
    float DetectPitch() {
        yin_input_.fill(0.0f);
        for (std::size_t i = 0; i < kPvWindowSize; ++i) {
            yin_input_[i] = input_ring_[(input_write_ + i) % kPvWindowSize] * window_[i];
        }
        yin_fft_.FFT(yin_input_.data(), yin_spectrum_.data());

        acf_in_.fill(0.0f);
        acf_in_[0] = yin_spectrum_[0] * yin_spectrum_[0] + yin_spectrum_[1] * yin_spectrum_[1];
        acf_in_[1] = 0.0f;
        for (std::size_t bin = 1; bin + 1 < kPvNumBins; ++bin) {
            float const re = yin_spectrum_[2 * bin];
            float const im = yin_spectrum_[2 * bin + 1];
            acf_in_[2 * bin] = re * re + im * im;
            acf_in_[2 * bin + 1] = 0.0f;
        }
        float const nyquist_re = yin_spectrum_[2 * (kPvNumBins - 1)];
        float const nyquist_im = yin_spectrum_[2 * (kPvNumBins - 1) + 1];
        acf_in_[kPvFrameSize] = nyquist_re * nyquist_re + nyquist_im * nyquist_im;
        acf_in_[kPvFrameSize + 1] = 0.0f;
        yin_fft_.IFFT(acf_in_.data(), yin_acf_.data());

        double energy = 0.0;
        prefix_[0] = 0.0f;
        for (std::size_t i = 0; i < kPvFrameSize; ++i) {
            energy += static_cast<double>(yin_input_[i]) * static_cast<double>(yin_input_[i]);
            prefix_[i + 1] = static_cast<float>(energy);
        }
        if (energy <= 1.0e-12 || yin_acf_[0] <= 1.0e-12f) {
            return 0.0f;
        }
        float const rms = static_cast<float>(std::sqrt(energy / static_cast<double>(kPvFrameSize)));
        if (rms < kSilenceRms) {
            return 0.0f;
        }
        float const scale = static_cast<float>(energy) / yin_acf_[0];

        int const min_lag = std::max(2, static_cast<int>(std::lround(sample_rate_ / kMaxF0)));
        int const max_lag = std::min(static_cast<int>(kPvWindowSize / 2),
                                     static_cast<int>(std::lround(sample_rate_ / kMinF0)));
        if (max_lag <= min_lag + 1) {
            return 0.0f;
        }

        diff_[0] = 0.0f;
        cmndf_[0] = 1.0f;
        float running = 0.0f;
        for (int tau = 1; tau <= max_lag; ++tau) {
            float const e0 = prefix_[kPvFrameSize - static_cast<std::size_t>(tau)];
            float const e1 = prefix_[kPvFrameSize] - prefix_[static_cast<std::size_t>(tau)];
            float const d = std::max(0.0f, e0 + e1 - 2.0f * yin_acf_[static_cast<std::size_t>(tau)] * scale);
            diff_[static_cast<std::size_t>(tau)] = d;
            running += d;
            cmndf_[static_cast<std::size_t>(tau)] = (running > 0.0f)
                ? d * static_cast<float>(tau) / running
                : 1.0f;
        }

        constexpr float kThreshold = 0.3f;
        int best = -1;
        for (int tau = min_lag; tau < max_lag; ++tau) {
            float const value = cmndf_[static_cast<std::size_t>(tau)];
            if (value < kThreshold && value < cmndf_[static_cast<std::size_t>(tau + 1)]) {
                best = tau;
                break;
            }
        }
        if (best < 0) {
            float best_value = std::numeric_limits<float>::infinity();
            for (int tau = min_lag; tau <= max_lag; ++tau) {
                float const value = cmndf_[static_cast<std::size_t>(tau)];
                if (value < best_value) {
                    best_value = value;
                    best = tau;
                }
            }
            if (best < 0 || best_value > kThreshold) {
                return 0.0f;
            }
        }

        float period = static_cast<float>(best);
        if (best > 1 && best < max_lag) {
            float const y0 = cmndf_[static_cast<std::size_t>(best - 1)];
            float const y1 = cmndf_[static_cast<std::size_t>(best)];
            float const y2 = cmndf_[static_cast<std::size_t>(best + 1)];
            float const denom = y0 - 2.0f * y1 + y2;
            if (std::abs(denom) > 1.0e-8f) {
                period += 0.5f * (y0 - y2) / denom;
                period = std::clamp(period, static_cast<float>(best) - 0.5f, static_cast<float>(best) + 0.5f);
            }
        }
        float const f0 = sample_rate_ / period;
        return (f0 >= kMinF0 && f0 <= kMaxF0) ? f0 : 0.0f;
    }

    void PropagatePghi(std::size_t analysis_hop) {
        processed_bins_.fill(0);
        heap_size_ = 0;
        std::size_t processed_count = 0;
        for (std::size_t bin = 0; bin < kPvNumBins; ++bin) {
            HeapPush({previous_magnitude_[bin], bin, true});
        }

        float const analysis_scale = 2.0f * std::numbers::pi_v<float>
                                   * static_cast<float>(analysis_hop) / static_cast<float>(kPvFrameSize);
        constexpr float kSynthesisScale = 2.0f * std::numbers::pi_v<float>
                                        * static_cast<float>(kPvHopSize) / static_cast<float>(kPvFrameSize);
        float const stretch_ratio = static_cast<float>(kPvHopSize) / static_cast<float>(analysis_hop);

        while (processed_count < kPvNumBins && heap_size_ > 0) {
            HeapItem const item = HeapPop();
            std::size_t const bin = item.bin;
            if (item.is_previous) {
                if (processed_bins_[bin] != 0) {
                    continue;
                }
                float const expected = analysis_scale * static_cast<float>(bin);
                float const residual = WrapToPi(current_analysis_phase_[bin] - previous_analysis_phase_[bin] - expected);
                current_synthesis_phase_[bin] = WrapToPi(previous_synthesis_phase_[bin]
                                                         + kSynthesisScale * static_cast<float>(bin)
                                                         + stretch_ratio * residual);
                processed_bins_[bin] = 1;
                ++processed_count;
                HeapPush({current_magnitude_[bin], bin, false});
                continue;
            }

            if (bin + 1 < kPvNumBins && processed_bins_[bin + 1] == 0) {
                float const delta = WrapToPi(current_analysis_phase_[bin + 1] - current_analysis_phase_[bin]);
                current_synthesis_phase_[bin + 1] = WrapToPi(current_synthesis_phase_[bin] + stretch_ratio * delta);
                processed_bins_[bin + 1] = 1;
                ++processed_count;
                HeapPush({current_magnitude_[bin + 1], bin + 1, false});
            }
            if (bin > 0 && processed_bins_[bin - 1] == 0) {
                float const delta = WrapToPi(current_analysis_phase_[bin - 1] - current_analysis_phase_[bin]);
                current_synthesis_phase_[bin - 1] = WrapToPi(current_synthesis_phase_[bin] + stretch_ratio * delta);
                processed_bins_[bin - 1] = 1;
                ++processed_count;
                HeapPush({current_magnitude_[bin - 1], bin - 1, false});
            }
        }
    }

    void SynthesizeSpectrum() {
        fft_output_.fill(0.0f);
        for (std::size_t bin = 0; bin < kPvNumBins; ++bin) {
            fft_output_[2 * bin] = current_magnitude_[bin] * std::cos(current_synthesis_phase_[bin]);
            fft_output_[2 * bin + 1] = current_magnitude_[bin] * std::sin(current_synthesis_phase_[bin]);
        }
        fft_output_[1] = 0.0f;
        fft_output_[2 * (kPvNumBins - 1) + 1] = 0.0f;
    }

    float RenderOutputSample() {
        constexpr std::uint64_t kResampleStart = kPvFrameSize - kPvHopSize;
        if (!resampler_started_) {
            if (synthesis_available_ <= kResampleStart + 1) {
                return 0.0f;
            }
            resample_position_ = static_cast<double>(kResampleStart);
            while (synthesis_cleared_until_ < kResampleStart) {
                synthesis_ring_[static_cast<std::size_t>(synthesis_cleared_until_ % kPvSynthesisRingSize)] = 0.0f;
                ++synthesis_cleared_until_;
            }
            resampler_started_ = true;
        }

        auto const lower = static_cast<std::uint64_t>(resample_position_);
        if (lower + 1 >= synthesis_available_) {
            return 0.0f;
        }
        float const fraction = static_cast<float>(resample_position_ - static_cast<double>(lower));
        float const first = synthesis_ring_[static_cast<std::size_t>(lower % kPvSynthesisRingSize)];
        float const second = synthesis_ring_[static_cast<std::size_t>((lower + 1) % kPvSynthesisRingSize)];
        float const output = first + fraction * (second - first);

        resample_position_ += static_cast<double>(current_ratio_);
        auto const clear_until = static_cast<std::uint64_t>(resample_position_);
        while (synthesis_cleared_until_ < clear_until) {
            synthesis_ring_[static_cast<std::size_t>(synthesis_cleared_until_ % kPvSynthesisRingSize)] = 0.0f;
            ++synthesis_cleared_until_;
        }
        return output;
    }

    void HeapPush(HeapItem item) {
        if (heap_size_ == kPvHeapCapacity) {
            return;
        }
        std::size_t index = heap_size_++;
        heap_[index] = item;
        while (index > 0) {
            std::size_t const parent = (index - 1) / 2;
            if (heap_[parent].magnitude >= heap_[index].magnitude) {
                break;
            }
            std::swap(heap_[parent], heap_[index]);
            index = parent;
        }
    }

    HeapItem HeapPop() {
        HeapItem const result = heap_[0];
        heap_[0] = heap_[--heap_size_];
        std::size_t index = 0;
        while (true) {
            std::size_t const left = index * 2 + 1;
            if (left >= heap_size_) {
                break;
            }
            std::size_t const right = left + 1;
            std::size_t const largest =
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
    bool first_frame_ = true;

    qwqdsp_spectral::RealFFT fft_;
    qwqdsp_spectral::RealFFT yin_fft_;
    std::array<float, kPvWindowSize> window_{};
    std::array<float, kPvWindowSize> input_ring_{};
    std::array<float, kPvSynthesisRingSize> synthesis_ring_{};
    std::array<float, kPvFrameSize> fft_input_{};
    std::array<float, kPvFrameSize + 2> fft_output_{};
    std::array<float, kPvFrameSize> yin_input_{};
    std::array<float, kPvFrameSize + 2> yin_spectrum_{};
    std::array<float, kPvFrameSize + 2> acf_in_{};
    std::array<float, kPvFrameSize> yin_acf_{};
    std::array<float, kPvFrameSize + 1> prefix_{};
    std::array<float, kPvWindowSize / 2 + 2> diff_{};
    std::array<float, kPvWindowSize / 2 + 2> cmndf_{};

    std::array<float, kPvNumBins> previous_analysis_phase_{};
    std::array<float, kPvNumBins> previous_synthesis_phase_{};
    std::array<float, kPvNumBins> previous_magnitude_{};
    std::array<float, kPvNumBins> current_re_{};
    std::array<float, kPvNumBins> current_im_{};
    std::array<float, kPvNumBins> current_analysis_phase_{};
    std::array<float, kPvNumBins> current_synthesis_phase_{};
    std::array<float, kPvNumBins> current_magnitude_{};

    std::array<std::uint8_t, kPvNumBins> processed_bins_{};
    std::array<HeapItem, kPvHeapCapacity> heap_{};
    std::size_t heap_size_ = 0;

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

inline bool RunPghi(std::span<const float> input, int fs, std::string const& name) {
    auto pv = std::make_unique<PghiAutotune>();
    pv->Init(static_cast<float>(fs));
    std::vector<float> output = pv->Process(input);

    std::cout << std::format("  pghi: Hann {} / hop {} / FFT {}, latency {} samples ({:.1f} ms)\n",
                             kPvWindowSize, kPvHopSize, kPvFrameSize,
                             PghiAutotune::LatencySamples(),
                             1000.0 * static_cast<double>(PghiAutotune::LatencySamples()) / static_cast<double>(fs));
    return SaveAndReport(input, output, fs, name);
}

} // namespace autotune_lab
