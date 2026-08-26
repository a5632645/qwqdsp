#pragma once
#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

#include "aubio_peak_pick.hpp"
#include "oneset_common.hpp"

namespace qwqdsp_test {

// ------------------------------------------------------------
// HfcOnsetDetector
// ------------------------------------------------------------
/**
 * @brief 高频含量（High-Frequency Content）onset 检测器。
 *
 * 移植自 ref_repos/aubio/src/spectral/specdesc.c:101-110 的
 * aubio_specdesc_hfc。ODF 为幅度按 bin 序号加权的和：
 *
 *   D(t) = Σ_{k=0}^{bins-1} (k+1) · |X_t(k)|
 *
 * 高频 bin 权重大，对打击乐/瞬态的宽带高频突增敏感，实现极简
 * （只需当前帧幅度谱）。峰值拾取用 aubio 自适应阈值。
 */
class HfcOnsetDetector {
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

            // HFC：Σ (k+1)·|X(k)|
            float hfc = 0.0f;
            for (size_t k = 0; k < bins; ++k)
                hfc += static_cast<float>(k + 1) * stft.mag_[k];

            result.odf.push_back(hfc);
        }

        detail::AubioPeakPicker picker;
        picker.SetThreshold(0.058f); // aubio hfc 默认
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
