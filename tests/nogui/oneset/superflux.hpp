#pragma once
#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

#include "oneset_common.hpp"

namespace qwqdsp_test {

// ------------------------------------------------------------
// SuperFluxDetector
// ------------------------------------------------------------
/**
 * @brief SuperFlux onset 检测器（Boeck & Widmer, DAFx 2013）。
 *
 * 移植自 DSPark Analysis/OnsetDetector.h 的 SuperFlux ODF：
 * 1. 对数频率三角滤波器组（四分之一音分辨率，27.5 Hz–16 kHz）
 *    得到 band 幅度；band 幅度先乘帧不变量 2048/fftSize 再 log10(x+1)。
 * 2. 频率方向 3 邻域最大滤波（对上一帧），抑制 vibrato/tremolo 误报。
 * 3. 对前帧（mu=1）做半波整流谱通量，除以 band 数归一化。
 *
 * 可选 Stowell-Plumbley 自适应白化（默认关）：逐 bin 除以运行峰值，
 * 白化后 band 不再做帧不变量缩放。
 *
 * 相比纯谱通量，对数滤波器组对低音区更敏感，最大滤波抑制颤音。
 */
class SuperFluxDetector {
public:
    static constexpr size_t kDefaultFrameSize = 2048;
    /** 滤波器组低边缘 (Hz) */
    static constexpr float kFMin = 27.5f;
    /** 滤波器组高边缘 (Hz) */
    static constexpr float kFMaxHz = 16000.0f;
    /** 滤波器组分辨率：每倍频程 band 数（四分之一音） */
    static constexpr int kBandsPerOctave = 24;
    /** ODF 帧不变量参考帧长 */
    static constexpr double kOdfRefFrame = 2048.0;
    /** 自适应白化衰减 */
    static constexpr float kWhitenDecay = 0.9995f;
    static constexpr float kWhitenFloor = 1e-4f;

    /** 设置分析帧长度 */
    void SetFrameSize(size_t n) noexcept {
        frame_size_ = n;
    }

    /** 设置分析 hop */
    void SetHopSize(size_t n) noexcept {
        hop_ = n;
    }

    /** 设置峰值拾取 delta（默认 0.03，DSPark SuperFlux 默认） */
    void SetDelta(float d) noexcept {
        delta_ = d;
    }

    /** 启用自适应白化（默认关） */
    void SetAdaptiveWhitening(bool on) noexcept {
        whitening_ = on;
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

        // ---- 滤波器组 ----
        const float bin_hz = static_cast<float>(sample_rate) / static_cast<float>(La);
        const float f_max = std::min(kFMaxHz, static_cast<float>(sample_rate) * 0.5f * 0.999f);

        std::vector<int> centres;
        for (int i = 0;; ++i) {
            const float f = kFMin * std::pow(2.0f, static_cast<float>(i) / static_cast<float>(kBandsPerOctave));
            if (f > f_max)
                break;
            int bin = static_cast<int>(std::lround(f / bin_hz));
            bin = std::clamp(bin, 0, static_cast<int>(bins) - 1);
            if (centres.empty() || bin > centres.back())
                centres.push_back(bin);
        }

        // 三角滤波器：相邻三个中心 (lo, ce, hi) 构成一个 band
        struct Band {
            int start;
            std::vector<float> weights;
        };
        std::vector<Band> bands;
        for (size_t j = 1; j + 1 < centres.size(); ++j) {
            const int lo = centres[j - 1];
            const int ce = centres[j];
            const int hi = centres[j + 1];
            if (!(lo < ce && ce < hi))
                continue;
            Band b;
            b.start = lo;
            b.weights.resize(static_cast<size_t>(hi - lo + 1));
            for (int k = lo; k <= hi; ++k) {
                b.weights[static_cast<size_t>(k - lo)] =
                    (k <= ce) ? static_cast<float>(k - lo) / static_cast<float>(ce - lo)
                              : static_cast<float>(hi - k) / static_cast<float>(hi - ce);
            }
            bands.push_back(std::move(b));
        }
        if (bands.empty())
            return {};

        const size_t num_bands = bands.size();
        const float odf_scale = static_cast<float>(kOdfRefFrame / static_cast<double>(La));

        std::vector<float> band_cur(num_bands, 0.0f);
        std::vector<float> band_max_prev(num_bands, 0.0f);
        std::vector<float> whiten_peak(bins, 0.0f);

        OnsetResult result;
        result.odf.reserve(N);

        for (size_t i = 0; i < N; ++i) {
            const float* frame = input.data() + i * hop;
            stft.Analyze(frame);

            // ---- 自适应白化（可选） ----
            std::vector<float> mag_used(bins);
            if (whitening_) {
                for (size_t j = 0; j < bins; ++j) {
                    float& pk = whiten_peak[j];
                    pk = std::max({stft.mag_[j], kWhitenFloor, pk * kWhitenDecay});
                    mag_used[j] = stft.mag_[j] / pk;
                }
            }
            else {
                mag_used = stft.mag_;
            }

            // ---- 滤波器组 + log 压缩 ----
            const float scale = whitening_ ? 1.0f : odf_scale;
            for (size_t b = 0; b < num_bands; ++b) {
                float acc = 0.0f;
                const Band& band = bands[b];
                for (size_t w = 0; w < band.weights.size(); ++w) {
                    const size_t k = static_cast<size_t>(band.start) + w;
                    if (k < bins)
                        acc += mag_used[k] * band.weights[w];
                }
                band_cur[b] = std::log10(acc * scale + 1.0f);
            }

            // ---- SuperFlux 通量：对前帧最大滤波参考做半波整流差分 ----
            // 首帧无前帧参考，warm-up 置 0（否则首帧差分虚高误报）
            float flux = 0.0f;
            if (i > 0) {
                for (size_t b = 0; b < num_bands; ++b) {
                    const float d = band_cur[b] - band_max_prev[b];
                    if (d > 0.0f)
                        flux += d;
                }
                flux /= static_cast<float>(num_bands);
            }

            // 旋转：当前 → 前帧，并对前帧做频率最大滤波
            std::vector<float> band_prev = band_cur;
            for (size_t b = 0; b < num_bands; ++b) {
                float mx = band_prev[b];
                if (b > 0)
                    mx = std::max(mx, band_prev[b - 1]);
                if (b + 1 < num_bands)
                    mx = std::max(mx, band_prev[b + 1]);
                band_max_prev[b] = mx;
            }

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
    float delta_ = 0.03f;
    bool whitening_ = false;
};

} // namespace qwqdsp_test
