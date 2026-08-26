#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

#include "aubio_peak_pick.hpp"
#include "oneset_common.hpp"

namespace qwqdsp_test {

// ------------------------------------------------------------
// HpssPercussiveOnsetDetector
// ------------------------------------------------------------
/**
 * @brief HPSS 打击乐分量 onset 检测器。
 *
 * 移植自 ref_repos/librosa-master/src/decompose.cpp 的 hpss
 * （Harmonic-Percussive Source Separation）：
 *
 * 1. 对整段信号做 STFT，得到二维幅度谱（帧 × bin）
 * 2. 沿时间轴中值滤波得到谐波分量（harmonic），沿频率轴中值滤波
 *    得到打击乐分量（percussive）
 * 3. Wiener 软掩码分离：mask_perc = softmask(perc, harm·margin, power)
 * 4. 打击乐分量做谱通量 ODF（半波整流差分）
 *
 * 离线批处理（需整段二维谱图），对打击乐起始最敏感。中值滤波用
 * std::nth_element（参考 librosa median_filter_2d，reflect 填充）。
 */
class HpssPercussiveOnsetDetector {
public:
    static constexpr size_t kDefaultFrameSize = 2048;
    /** 中值滤波核大小（librosa 默认 31） */
    static constexpr int kDefaultKernelSize = 31;
    /** 软掩码幂（librosa 默认 2.0） */
    static constexpr float kDefaultPower = 2.0f;

    /** 设置分析帧长度 */
    void SetFrameSize(size_t n) noexcept {
        frame_size_ = n;
    }

    /** 设置分析 hop */
    void SetHopSize(size_t n) noexcept {
        hop_ = n;
    }

    /** 检测输入信号中的 onset */
    OnsetResult Detect(std::span<const float> input, float sample_rate) {
        const size_t La = frame_size_;
        const size_t hop = hop_;

        if (La == 0)
            return {};

        const size_t bins = La / 2 + 1;
        const size_t N = (input.size() > La) ? (input.size() - La) / hop + 1 : 0;
        if (N == 0)
            return {};

        // ---- STFT 得到二维幅度谱（帧 × bin，行主序） ----
        detail::StftFrame stft;
        stft.Init(La);

        std::vector<float> spectrogram(N * bins, 0.0f);
        for (size_t i = 0; i < N; ++i) {
            const float* frame = input.data() + i * hop;
            stft.Analyze(frame);
            std::copy(stft.mag_.begin(), stft.mag_.end(), spectrogram.begin() + static_cast<std::ptrdiff_t>(i * bins));
        }

        // ---- 中值滤波（时间轴→谐波，频率轴→打击乐） ----
        const int kernel = std::max(1, kDefaultKernelSize);
        const int half = kernel / 2;

        // 谐波：沿时间轴（每 bin 跨帧中值滤波）
        std::vector<float> harmonic(N * bins, 0.0f);
        for (size_t k = 0; k < bins; ++k) {
            for (size_t i = 0; i < N; ++i) {
                std::vector<float> window;
                window.reserve(static_cast<size_t>(kernel));
                for (int j = -half; j <= half; ++j) {
                    // reflect 填充
                    int idx = static_cast<int>(i) + j;
                    if (idx < 0)
                        idx = -idx;
                    if (idx >= static_cast<int>(N))
                        idx = 2 * static_cast<int>(N) - 1 - idx;
                    window.push_back(spectrogram[static_cast<size_t>(idx) * bins + k]);
                }
                std::nth_element(window.begin(), window.begin() + window.size() / 2, window.end());
                harmonic[i * bins + k] = window[window.size() / 2];
            }
        }

        // 打击乐：沿频率轴（每帧跨 bin 中值滤波）
        std::vector<float> percussive(N * bins, 0.0f);
        for (size_t i = 0; i < N; ++i) {
            for (size_t k = 0; k < bins; ++k) {
                std::vector<float> window;
                window.reserve(static_cast<size_t>(kernel));
                for (int j = -half; j <= half; ++j) {
                    int idx = static_cast<int>(k) + j;
                    if (idx < 0)
                        idx = -idx;
                    if (idx >= static_cast<int>(bins))
                        idx = 2 * static_cast<int>(bins) - 1 - idx;
                    window.push_back(spectrogram[i * bins + static_cast<size_t>(idx)]);
                }
                std::nth_element(window.begin(), window.begin() + window.size() / 2, window.end());
                percussive[i * bins + k] = window[window.size() / 2];
            }
        }

        // ---- 软掩码（Wiener）：mask_perc = perc^p / (harm^p + perc^p) ----
        std::vector<float> perc_mask(N * bins, 0.0f);
        for (size_t i = 0; i < N * bins; ++i) {
            const float h = harmonic[i];
            const float p = percussive[i];
            const float hp = std::pow(h, kDefaultPower);
            const float pp = std::pow(p, kDefaultPower);
            perc_mask[i] = (hp + pp > 1e-12f) ? pp / (hp + pp) : 0.0f;
        }

        // ---- 打击乐分量 + 谱通量 ODF ----
        OnsetResult result;
        result.odf.reserve(N);
        for (size_t i = 0; i < N; ++i) {
            float perc_flux = 0.0f;
            for (size_t k = 0; k < bins; ++k) {
                const float cur = spectrogram[i * bins + k] * perc_mask[i * bins + k];
                const float prev = (i > 0) ? spectrogram[(i - 1) * bins + k] * perc_mask[(i - 1) * bins + k] : 0.0f;
                const float d = cur - prev;
                if (d > 0.0f)
                    perc_flux += d;
            }
            result.odf.push_back(perc_flux);
        }

        detail::AubioPeakPicker picker;
        picker.SetThreshold(0.3f); // 与 aubio specflux 同量级
        result.onset_frames = picker.Pick(result.odf, hop, sample_rate);
        for (size_t f : result.onset_frames)
            result.onset_samples.push_back(f * hop);

        return result;
    }

private:
    size_t frame_size_ = kDefaultFrameSize;
    size_t hop_ = kDefaultFrameSize / 4;
};

} // namespace qwqdsp_test
