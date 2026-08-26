#pragma once
#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

#include "oneset_common.hpp"

namespace qwqdsp_test {

// ------------------------------------------------------------
// ComplexDomainDetector
// ------------------------------------------------------------
/**
 * @brief 复合域 onset 检测器（Dixon 2006，DAFx）。
 *
 * 在复平面上比较实际频谱与预测频谱：用上一帧与上上帧的相位外推
 * 预测当前帧相位（target = princArg(2·prev_phase − prev_phase2)），
 * 保持上一帧幅度，计算 |X − Xhat|；仅当当前幅度 ≥ 上一帧幅度时累加
 * （半波整流，抑制衰减段）。归一化到 bin 数后用自适应峰值拾取。
 *
 * 相比纯谱通量，复合域同时利用相位与幅度信息，对相位相干信号的
 * 起始更敏感。
 */
class ComplexDomainDetector {
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

    /** 设置峰值拾取 delta（默认 0.05，复合域 ODF 尺度更小） */
    void SetDelta(float d) noexcept {
        delta_ = d;
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
            if (i == 0) {
                result.odf.push_back(0.0f);
                continue;
            }

            float odf = 0.0f;
            for (size_t j = 0; j < bins; ++j) {
                // 相位预测（两帧线性外推）
                const float target = detail::PrincArg(2.0f * stft.prev_phase_[j] - stft.prev_phase2_[j]);
                const float pm = stft.prev_mag_[j];
                const float cm = stft.mag_[j];

                // X 与预测 Xhat 之差（幅度变化时相位保持不变）
                const float re = cm * std::cos(stft.phase_[j]) - pm * std::cos(target);
                const float im = cm * std::sin(stft.phase_[j]) - pm * std::sin(target);

                // 半波整流：仅幅度增加时累加
                if (cm >= pm)
                    odf += std::sqrt(re * re + im * im);
            }
            odf /= static_cast<float>(bins);

            result.odf.push_back(odf);
        }

        result.onset_frames = detail::PeakPick(result.odf, hop, sample_rate, delta_);
        for (size_t f : result.onset_frames)
            result.onset_samples.push_back(f * hop);

        return result;
    }

private:
    size_t frame_size_ = kDefaultFrameSize;
    size_t hop_ = kDefaultFrameSize / 4;
    float delta_ = 0.05f;
};

} // namespace qwqdsp_test
