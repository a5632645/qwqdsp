#pragma once
#include <algorithm>
#include <cmath>
#include <deque>
#include <span>
#include <vector>

namespace qwqdsp_test {
namespace detail {

// ------------------------------------------------------------
// AubioPeakPicker
// ------------------------------------------------------------
/**
 * @brief aubio 自适应峰值拾取（aubio_peakpicker_do）。
 *
 * 移植自 ref_repos/aubio/src/onset/peakpicker.c:86-123。与 librosa 风格
 * 的 PeakPick（移动均值 + delta）不同：对 ODF 历史做 biquad 低通平滑，
 * 自适应阈值 T = median(历史) + mean(历史)·threshold，阈值化后做三点
 * 局部峰判断，命中时用二次插值细化峰位置。
 *
 * 默认参数：threshold=0.1、win_post=5、win_pre=1（与 aubio 一致），
 * biquad 为 butter(2, 0.34) 低通。
 */
class AubioPeakPicker {
public:
    static constexpr float kThreshold = 0.1f;
    static constexpr int kWinPost = 5;
    static constexpr int kWinPre = 1;
    static constexpr int kHistoryLen = kWinPost + kWinPre + 1; // 7
    /** butter(2, 0.34) 低通系数（aubio peakpicker.c:181-184） */
    static constexpr float kBiquadB0 = 0.15998789f;
    static constexpr float kBiquadB1 = 0.31997577f;
    static constexpr float kBiquadB2 = 0.15998789f;
    static constexpr float kBiquadA1 = 0.23484048f;
    static constexpr float kBiquadA2 = 0.0f;

    /** 设置阈值（越大越少触发） */
    void SetThreshold(float t) noexcept {
        threshold_ = t;
    }

    /**
     * @brief 对 ODF 序列做自适应峰值拾取。
     *
     * @param odf      onset 强度序列（每帧一个值）
     * @param hop      帧 hop（采样数）
     * @param sample_rate 采样率（用于 min_ioi 换算）
     * @param min_ioi_ms 最小 onset 间隔（毫秒，aubio 默认 50）
     * @return onset 帧下标（升序）
     */
    std::vector<size_t> Pick(std::span<const float> odf, size_t hop = 1, float sample_rate = 1.0f,
                             float min_ioi_ms = 50.0f) {
        std::vector<size_t> onsets;
        if (odf.empty())
            return onsets;

        // 最小间隔换算成帧数
        const size_t min_ioi_frames = std::max<size_t>(1, static_cast<size_t>(
            std::lround(min_ioi_ms * 0.001f * sample_rate / static_cast<float>(std::max<size_t>(1, hop)))));
        size_t last_onset = SIZE_MAX;

        // 阈值化后的延迟序列（每帧一个，对应 aubio 的 thresholded）
        std::vector<float> thresholded(odf.size(), 0.0f);
        std::deque<float> history; // 最近 kHistoryLen 个平滑 ODF
        float bq_state0 = 0.0f, bq_state1 = 0.0f; // biquad 直接型 II 状态

        for (size_t i = 0; i < odf.size(); ++i) {
            const float x = odf[i];

            // biquad 低通（butter(2,0.34)）：y = b0·x + s0; s0 = b1·x - a1·y + s1; s1 = b2·x - a2·y
            const float y = kBiquadB0 * x + bq_state0;
            bq_state0 = kBiquadB1 * x - kBiquadA1 * y + bq_state1;
            bq_state1 = kBiquadB2 * x - kBiquadA2 * y;

            history.push_back(y);
            if (history.size() > static_cast<size_t>(kHistoryLen))
                history.pop_front();

            if (history.size() < static_cast<size_t>(kHistoryLen))
                continue; // 历史未满，不判定（aubio 缓冲未满时阈值无意义）

            // 历史 mean 和 median
            std::vector<float> hist(history.begin(), history.end());
            float mean = 0.0f;
            for (float v : hist)
                mean += v;
            mean /= static_cast<float>(hist.size());

            std::vector<float> sorted = hist;
            std::nth_element(sorted.begin(), sorted.begin() + sorted.size() / 2, sorted.end());
            const float median = sorted[sorted.size() / 2];

            // 阈值化：取延迟样本（win_post 前），T = median + mean·threshold
            // aubio: thresholded[0] = onset_proc[win_post] - median - mean*threshold
            if (i >= static_cast<size_t>(kWinPost)) {
                const float proc = odf[i - static_cast<size_t>(kWinPost)];
                thresholded[i] = proc - median - mean * threshold_;
            }
        }

        // 三点局部峰判断（fvec_peakpick）：thresholded[i-1] 是局部最大且 > 0
        for (size_t i = 1; i + 1 < thresholded.size(); ++i) {
            const float t1 = thresholded[i - 1];
            const float t2 = thresholded[i];
            const float t3 = thresholded[i + 1];
            if (t2 >= t1 && t2 > t3 && t2 > 0.0f) {
                // 二次插值细化峰位置（fvec_quadratic_peak_pos）
                size_t peak = i;
                const float denom = t1 - 2.0f * t2 + t3;
                if (std::abs(denom) > 1e-10f) {
                    const float offset = 0.5f * (t1 - t3) / denom;
                    if (offset > 0.5f && i + 1 < thresholded.size())
                        peak = i + 1;
                    else if (offset < -0.5f && i > 0)
                        peak = i - 1;
                }
                if (last_onset == SIZE_MAX || peak - last_onset >= min_ioi_frames) {
                    onsets.push_back(peak);
                    last_onset = peak;
                }
            }
        }

        return onsets;
    }

private:
    float threshold_ = kThreshold;
};

} // namespace detail
} // namespace qwqdsp_test
