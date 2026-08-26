#pragma once
#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

#include "aubio_peak_pick.hpp"
#include "oneset_common.hpp"

namespace qwqdsp_test {

// ------------------------------------------------------------
// SpectralDifferenceOnsetDetector
// ------------------------------------------------------------
/**
 * @brief 谱差（Spectral Difference）onset 检测器。
 *
 * 移植自 ref_repos/aubio/src/spectral/specdesc.c:182-208 的
 * aubio_specdesc_specdiff。取平方幅度差的绝对值：
 *
 *   D(t) = Σ_k √| X_t(k)² − X_{t-1}(k)² |
 *
 * 与正谱通量（只响应上升）不同，幅度下降也产生响应（平方差绝对值）。
 * 低幅 bin（<阈值）置 0 抑制噪声。峰值拾取用 aubio 自适应阈值。
 */
class SpectralDifferenceOnsetDetector {
public:
    static constexpr size_t kDefaultFrameSize = 2048;
    /** 低幅 bin 抑制阈值（aubio 默认 0.03 附近） */
    static constexpr float kDefaultMagThreshold = 0.03f;

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

        detail::StftFrame stft;
        stft.Init(La);

        OnsetResult result;
        result.odf.reserve(N);

        for (size_t i = 0; i < N; ++i) {
            const float* frame = input.data() + i * hop;
            stft.Analyze(frame);

            // 首帧无上一帧幅度，warm-up 置 0
            if (i == 0) {
                result.odf.push_back(0.0f);
                continue;
            }

            float odf = 0.0f;
            for (size_t k = 0; k < bins; ++k) {
                if (stft.mag_[k] < kDefaultMagThreshold)
                    continue; // 低幅 bin 抑制
                const float cur2 = stft.mag_[k] * stft.mag_[k];
                const float prev2 = stft.prev_mag_[k] * stft.prev_mag_[k];
                odf += std::sqrt(std::abs(cur2 - prev2));
            }

            result.odf.push_back(odf);
        }

        detail::AubioPeakPicker picker;
        picker.SetThreshold(0.3f); // aubio specdiff 默认
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
