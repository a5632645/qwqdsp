#pragma once
#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

#include "oneset_common.hpp"

namespace qwqdsp_test {

// ------------------------------------------------------------
// BaselineTransientDetector
// ------------------------------------------------------------
/**
 * @brief 基线瞬态检测器 — 即当前 transient_vocoder 的检测逻辑。
 *
 * 从 pitch_shift/transient_vocoder.hpp 的 DetectTransient 抽离：
 * 逐 bin log1p 幅度谱通量（正增量累加）除以能量加权归一化，
 * EMA 跟踪均值/方差，双阈值触发：
 *   norm_flux > mean + threshold·max(0.07, std) 且 norm_flux > mean·1.35
 * 帧数 > 4 且冷却期（cooldown）结束后才能触发。
 *
 * 作为对照基线：新算法（能量法/谱通量/复合域/SuperFlux）应在此基础上
 * 有更准确的检测或更少的误报。
 */
class BaselineTransientDetector {
public:
    static constexpr size_t kDefaultFrameSize = 2048;
    static constexpr float kDefaultThreshold = 1.5f;

    /** 设置分析帧长度 */
    void SetFrameSize(size_t n) noexcept {
        frame_size_ = n;
    }

    /** 设置分析 hop */
    void SetHopSize(size_t n) noexcept {
        hop_ = n;
    }

    /** 设置瞬态检测阈值（越大越少触发） */
    void SetThreshold(float t) noexcept {
        threshold_ = t;
    }

    /** 检测输入信号中的瞬态 onset */
    OnsetResult Detect(std::span<const float> input, float /*sample_rate*/) {
        const size_t La = frame_size_;
        const size_t hop = hop_;

        if (La == 0)
            return {};

        const size_t bins = La / 2 + 1;

        // 帧数：分析窗口落在输入内
        const size_t N = (input.size() > La) ? (input.size() - La) / hop + 1 : 0;

        std::vector<float> window(La);
        detail::MakePeriodicHann(window);

        qwqdsp_spectral::RealFFT fft;
        fft.Init(La);
        std::vector<float> fft_in(La);
        std::vector<float> fft_out(La + 2);

        std::vector<float> mag(bins);
        std::vector<float> prev_mag(bins, 0.0f);

        float flux_mean = 0.0f;
        float flux_var = 0.0f;
        bool has_mean = false;
        size_t frames = 0;
        size_t cooldown = 0;
        bool first = true;

        OnsetResult result;
        result.odf.reserve(N);

        const auto update_flux_stats = [&](float value, float alpha) {
            if (!has_mean) {
                flux_mean = value;
                flux_var = 0.0f;
                has_mean = true;
                return;
            }
            const float delta = value - flux_mean;
            flux_mean += alpha * delta;
            flux_var = (1.0f - alpha) * (flux_var + alpha * delta * delta);
        };

        for (size_t i = 0; i < N; ++i) {
            const float* frame = input.data() + i * hop;
            for (size_t j = 0; j < La; ++j)
                fft_in[j] = frame[j] * window[j];

            fft.FFT(fft_in.data(), fft_out.data());

            for (size_t j = 0; j < bins; ++j) {
                const float re = fft_out[2 * j];
                const float im = fft_out[2 * j + 1];
                mag[j] = std::sqrt(re * re + im * im);
            }

            // ---- 瞬态检测 ----
            bool is_transient = false;
            float odf_value = 0.0f;
            if (!first) {
                float flux = 0.0f;
                float energy = 0.0f;
                for (size_t j = 0; j < bins; ++j) {
                    const float weight = 0.5f + 0.5f * static_cast<float>(j) / std::max<size_t>(1, bins - 1);
                    const float d = std::log1p(mag[j]) - std::log1p(prev_mag[j]);
                    if (d > 0.0f)
                        flux += d;
                    energy += weight * std::log1p(mag[j]);
                }
                const float norm_flux = (energy > 1e-10f) ? flux / energy : 0.0f;
                const float mean = has_mean ? flux_mean : norm_flux;
                const float std = std::sqrt(flux_var);
                // Std floor 0.07：稳定拍频 normFlux ≤ ~0.07，真实起始 ≥ ~0.19
                is_transient = frames > 4 && cooldown == 0 && norm_flux > mean + threshold_ * std::max(0.07f, std)
                               && norm_flux > mean * 1.35f;
                update_flux_stats(norm_flux, is_transient ? 0.3f : 0.12f);
                cooldown = is_transient ? 1 : (cooldown > 0 ? cooldown - 1 : 0);
                odf_value = norm_flux;
            }

            prev_mag = mag;
            first = false;
            ++frames;

            result.odf.push_back(odf_value);
            if (is_transient)
                result.onset_frames.push_back(i);
        }

        for (size_t f : result.onset_frames)
            result.onset_samples.push_back(f * hop);

        return result;
    }

private:
    size_t frame_size_ = kDefaultFrameSize;
    size_t hop_ = kDefaultFrameSize / 4;
    float threshold_ = kDefaultThreshold;
};

} // namespace qwqdsp_test
