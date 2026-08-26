#pragma once
#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

#include "aubio_peak_pick.hpp"
#include "oneset_common.hpp"

namespace qwqdsp_test {

// ------------------------------------------------------------
// KlDivergenceOnsetDetector
// ------------------------------------------------------------
/**
 * @brief KL 散度（Kullback-Leibler）onset 检测器。
 *
 * 移植自 ref_repos/aubio/src/spectral/specdesc.c:210-223 的
 * aubio_specdesc_kl。以幅度谱相对上一帧的 KL 散度为 ODF：
 *
 *   D(t) = Σ_k X_t(k) · log(1 + X_t(k) / (X_{t-1}(k) + 0.1))
 *
 * 常数 0.1 为稳定项（aubio 原值）。NaN 结果置 0。
 * 峰值拾取用 aubio 自适应阈值。
 */
class KlDivergenceOnsetDetector {
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

            // 首帧无上一帧幅度，warm-up 置 0
            if (i == 0) {
                result.odf.push_back(0.0f);
                continue;
            }

            float odf = 0.0f;
            for (size_t k = 0; k < bins; ++k) {
                const float cur = stft.mag_[k];
                const float prev = stft.prev_mag_[k];
                const float term = cur * std::log1p(cur / (prev + 0.1f));
                if (std::isfinite(term))
                    odf += term;
            }

            result.odf.push_back(odf);
        }

        detail::AubioPeakPicker picker;
        picker.SetThreshold(0.35f); // aubio kl 默认
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
