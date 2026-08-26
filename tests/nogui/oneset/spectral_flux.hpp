#pragma once
#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

#include "oneset_common.hpp"

namespace qwqdsp_test {

// ------------------------------------------------------------
// SpectralFluxDetector
// ------------------------------------------------------------
/**
 * @brief 谱通量 onset 检测器（librosa onset_strength / DSPark SpectralFlux）。
 *
 * ODF = 线性幅度谱对前一帧的半波整流差分，除以 bin 数归一化：
 *   flux = Σ max(0, mag[k] - prev_mag[k]) / numBins
 * 可选频率方向最大滤波（vibrato/tremolo 抑制，librosa max_size=3）。
 * 最后用自适应峰值拾取（detail::PeakPick）得到 onset 位置。
 *
 * 比基线（log1p 能量加权版）更接近标准定义：线性幅度差分对
 * 中低频幅度突变更直接。
 */
class SpectralFluxDetector {
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

    /** 设置峰值拾取 delta（默认 0.07，librosa 默认） */
    void SetDelta(float d) noexcept {
        delta_ = d;
    }

    /** 启用频率最大滤波（vibrato 抑制，默认开） */
    void SetVibratoSuppress(bool on) noexcept {
        vibrato_suppress_ = on;
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

            // 首帧无上一帧参考，warm-up 置 0（否则首帧差分虚高误报）
            if (i == 0) {
                result.odf.push_back(0.0f);
                continue;
            }

            // 频率最大滤波（librosa max_size=3 沿频率轴）
            std::vector<float> mag_ref(bins);
            if (vibrato_suppress_) {
                for (size_t j = 0; j < bins; ++j) {
                    float mx = stft.prev_mag_[j];
                    if (j > 0)
                        mx = std::max(mx, stft.prev_mag_[j - 1]);
                    if (j + 1 < bins)
                        mx = std::max(mx, stft.prev_mag_[j + 1]);
                    mag_ref[j] = mx;
                }
            }
            else {
                mag_ref = stft.prev_mag_;
            }

            // 半波整流差分
            float flux = 0.0f;
            for (size_t j = 0; j < bins; ++j) {
                const float d = stft.mag_[j] - mag_ref[j];
                if (d > 0.0f)
                    flux += d;
            }
            flux /= static_cast<float>(bins);

            result.odf.push_back(flux);
        }

        result.onset_frames = detail::PeakPick(result.odf, hop, sample_rate, delta_);
        for (size_t f : result.onset_frames)
            result.onset_samples.push_back(f * hop);

        return result;
    }

private:
    size_t frame_size_ = kDefaultFrameSize;
    size_t hop_ = kDefaultFrameSize / 4;
    float delta_ = 0.07f;
    bool vibrato_suppress_ = true;
};

} // namespace qwqdsp_test
