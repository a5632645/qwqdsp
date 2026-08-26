#pragma once
#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

#include "aubio_peak_pick.hpp"
#include "oneset_common.hpp"

namespace qwqdsp_test {

// ------------------------------------------------------------
// WeightedPhaseDeviationOnsetDetector
// ------------------------------------------------------------
/**
 * @brief 加权相位偏差（Weighted Phase Deviation）onset 检测器。
 *
 * 移植自 ref_repos/aubio/src/spectral/specdesc.c:166-180 的
 * aubio_specdesc_wphase。先算相位偏差，再按当前幅度加权：
 *
 *   D(t) = mean_k |X_t(k)| · | unwrap(φt(k) − 2·φt−1(k) + φt−2(k)) |
 *
 * 幅度加权让强分量的相位突变贡献更大，弱分量（可能是噪声）影响小。
 * 峰值拾取用 aubio 自适应阈值。
 */
class WeightedPhaseDeviationOnsetDetector {
public:
    static constexpr size_t kDefaultFrameSize = 2048;

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
                if (stft.mag_[k] < 1e-6f)
                    continue;
                // 相位偏差 × 当前幅度
                const float dev = detail::PrincArg(stft.phase_[k] - 2.0f * stft.prev_phase_[k] + stft.prev_phase2_[k]);
                odf += stft.mag_[k] * std::abs(dev);
                ++count;
            }
            odf /= (count > 0) ? static_cast<float>(count) : 1.0f;

            result.odf.push_back(odf);
        }

        detail::AubioPeakPicker picker;
        picker.SetThreshold(0.3f); // aubio wphase 默认
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
