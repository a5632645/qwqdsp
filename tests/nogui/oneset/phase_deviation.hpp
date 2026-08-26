#pragma once
#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

#include "aubio_peak_pick.hpp"
#include "oneset_common.hpp"

namespace qwqdsp_test {

// ------------------------------------------------------------
// PhaseDeviationOnsetDetector
// ------------------------------------------------------------
/**
 * @brief 相位偏差（Phase Deviation）onset 检测器。
 *
 * 移植自 ref_repos/aubio/src/spectral/specdesc.c:135-165 的
 * aubio_specdesc_phase。用两帧相位历史线性外推预测当前相位，
 * 偏差（unwrap 后）反映相位突变：
 *
 *   D(t) = mean_k | unwrap(φt(k) − 2·φt−1(k) + φt−2(k)) |
 *
 * 幅度低于阈值（默认 0.03）的 bin 置 0（抑制噪声相位）。
 * 峰值拾取用 aubio 自适应阈值。
 */
class PhaseDeviationOnsetDetector {
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

            // 首帧无相位历史，warm-up 置 0
            if (i < 2) {
                result.odf.push_back(0.0f);
                continue;
            }

            float odf = 0.0f;
            size_t count = 0;
            for (size_t k = 0; k < bins; ++k) {
                if (stft.mag_[k] < kDefaultMagThreshold)
                    continue; // 低幅 bin 抑制
                // 相位偏差：unwrap(φt − 2φt−1 + φt−2)
                const float dev = detail::PrincArg(stft.phase_[k] - 2.0f * stft.prev_phase_[k] + stft.prev_phase2_[k]);
                odf += std::abs(dev);
                ++count;
            }
            odf /= (count > 0) ? static_cast<float>(count) : 1.0f;

            result.odf.push_back(odf);
        }

        detail::AubioPeakPicker picker;
        picker.SetThreshold(0.3f); // aubio phase 默认（未覆盖）
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
