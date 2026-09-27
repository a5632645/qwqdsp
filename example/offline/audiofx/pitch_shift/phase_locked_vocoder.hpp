#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <qwqdsp/spectral/real_fft.hpp>
#include <span>
#include <vector>

namespace qwqdsp_test {
namespace detail {

/** 将相位折叠至 [-π, π) */
static inline float WrapToPi(float x) noexcept {
    return x - 2.0f * std::numbers::pi_v<float> * std::round(x / (2.0f * std::numbers::pi_v<float>));
}

/** periodic Hann 窗：w[i] = 0.5 * (1 - cos(2πi/N)) */
static inline void MakePeriodicHann(std::span<float> w) noexcept {
    const size_t n = w.size();
    for (size_t i = 0; i < n; ++i) {
        w[i] = 0.5f * (1.0f - std::cos(2.0f * std::numbers::pi_v<float> * static_cast<float>(i) / static_cast<float>(n)));
    }
}

/** 线性重采样，ratio > 1 变长，< 1 变短 */
static inline std::vector<float> LinearResample(std::span<const float> in, float ratio) {
    const size_t out_len = static_cast<size_t>(std::round(static_cast<float>(in.size()) / ratio));
    if (out_len < 2)
        return {};
    std::vector<float> out(out_len);
    const float step = static_cast<float>(in.size() - 1) / static_cast<float>(out_len - 1);
    for (size_t i = 0; i < out_len; ++i) {
        const float pos = step * static_cast<float>(i);
        const size_t idx = static_cast<size_t>(pos);
        const float frac = pos - static_cast<float>(idx);
        const size_t nxt = std::min(idx + 1, in.size() - 1);
        out[i] = in[idx] + frac * (in[nxt] - in[idx]);
    }
    return out;
}

// ------------------------------------------------------------
// LockPhase — Identity Phase Locking
// ------------------------------------------------------------
/**
 * @brief 相位锁定（Laroche & Dolson 1999 identity phase locking）。
 *
 * 移植自 @audio/spectral-pvoc 的 lockPhase：检测谱峰后，把每个峰
 * 相邻区域（按相邻峰中点划分）内的非峰 bin 相位设为
 * "该 bin 分析相位 + 峰处传播相位偏移"，恢复谐波垂直相位相干，
 * 消除普通声码器的 "phasiness"。
 *
 * @param mag  当前帧幅度（size = half + 1）
 * @param prop 传播后的合成相位（in/out，size = half + 1）
 * @param ana  当前帧分析相位（size = half + 1）
 * @param half 最大 bin 序号（fft_size / 2）
 */
static inline void LockPhase(const float* mag, float* prop, const float* ana, size_t half) {
    // --- 峰值检测（peakMask） ---
    std::vector<uint8_t> peaks(half + 1, 0);
    if (half <= 1) {
        peaks[0] = 1;
        if (half == 1)
            peaks[1] = 1;
    }
    else {
        float max_mag = 0.0f;
        for (size_t k = 0; k <= half; ++k)
            max_mag = std::max(max_mag, mag[k]);

        const float min_mag = std::max(1e-8f, max_mag * 0.015f);
        const float min_prom = std::max(1e-9f, max_mag * 0.003f);
        int last_peak = -2;
        float last_peak_mag = 0.0f;

        for (size_t k = 1; k < half; ++k) {
            const float value = mag[k];
            if (value < min_mag || value < mag[k - 1] || value < mag[k + 1])
                continue;

            float shoulder = mag[k - 1];
            shoulder = std::max(shoulder, mag[k + 1]);
            if (k > 1)
                shoulder = std::max(shoulder, mag[k - 2]);
            if (k + 2 <= half)
                shoulder = std::max(shoulder, mag[k + 2]);
            if (value - shoulder < min_prom && value < max_mag * 0.1f)
                continue;

            if (static_cast<int>(k) - last_peak <= 1) {
                if (value > last_peak_mag) {
                    peaks[last_peak] = 0;
                    peaks[k] = 1;
                    last_peak = static_cast<int>(k);
                    last_peak_mag = value;
                }
                continue;
            }
            peaks[k] = 1;
            last_peak = static_cast<int>(k);
            last_peak_mag = value;
        }

        // 无峰时取最大幅度 bin
        bool found = false;
        for (size_t k = 0; k <= half; ++k) {
            if (peaks[k]) {
                found = true;
                break;
            }
        }
        if (!found) {
            size_t best = 0;
            for (size_t k = 1; k <= half; ++k) {
                if (mag[k] > mag[best])
                    best = k;
            }
            peaks[best] = 1;
        }
    }

    // --- 按峰划分区域并锁定 ---
    std::vector<size_t> peak_bins;
    peak_bins.reserve(half + 1);
    for (size_t k = 0; k <= half; ++k) {
        if (peaks[k])
            peak_bins.push_back(k);
    }
    if (peak_bins.empty())
        return;

    for (size_t i = 0; i < peak_bins.size(); ++i) {
        const size_t pk = peak_bins[i];
        const size_t start = (i == 0) ? 0 : ((peak_bins[i - 1] + pk) / 2) + 1;
        const size_t end = (i == peak_bins.size() - 1) ? half : (pk + peak_bins[i + 1]) / 2;
        const float delta = prop[pk] - ana[pk];
        const float lock_floor = std::max(1e-10f, mag[pk] * 0.03f);
        for (size_t k = start; k <= end; ++k) {
            if (k == pk || mag[k] < lock_floor)
                continue;
            prop[k] = ana[k] + delta;
        }
    }
}

} // namespace detail

// ------------------------------------------------------------
// PhaseLockedVocoder
// ------------------------------------------------------------
/**
 * @brief 相位锁定声码器（phase-locked vocoder）。
 *
 * 移植自 ref_repos/stretch 的 stretch-pvoc-lock（Laroche & Dolson 1999）：
 * 标准相位传播之后，用 identity phase locking 把非峰 bin 锁到最近谱峰
 * 的旋转上，恢复谐波相位相干，消除普通声码器的 "phasiness"。
 *
 * 时间拉伸通过 ana_hop / syn_hop 比率控制，音高移动通过输出重采样实现
 * （与 phase_vocoder2 相同的解耦方式）。
 */
class PhaseLockedVocoder {
public:
    static constexpr size_t kDefaultFrameSize = 2048;
    static constexpr size_t kDefaultHopSize = 512;

    /** 设置分析帧长度 */
    void SetFrameSize(size_t n) noexcept {
        frame_size_ = n;
    }

    /** 设置合成 hop（也决定默认分析 hop） */
    void SetHopSize(size_t n) noexcept {
        hop_size_ = n;
    }

    /** 设置时间拉伸比 */
    void SetTimeStretch(float kt) noexcept {
        kt_ = kt;
    }

    /** 设置音高比例（>1 升高，<1 降低） */
    void SetPitchShift(float kp) noexcept {
        kp_ = kp;
    }

    /** 清理相位传播状态 */
    void Reset() noexcept {
        prev_.clear();
        syn_prev_.clear();
        first_ = true;
    }

    /** 处理输入信号，返回拉伸/移调后的输出 */
    std::vector<float> Process(std::span<const float> input) {
        const float kt = kt_;
        const float kp = kp_;
        const size_t La = frame_size_;
        const size_t syn_hop = hop_size_;

        // 总拉伸比 = kt * kp（之后重采样还原 kp 部分）
        const float total_stretch = kt * kp;
        const size_t ana_hop = std::max<size_t>(1, static_cast<size_t>(std::round(static_cast<float>(syn_hop) / std::max(total_stretch, 1e-6f))));

        if (La == 0)
            return {};

        const size_t fft_size = La;
        const size_t bins = fft_size / 2 + 1;
        const float freq_per_bin = 2.0f * std::numbers::pi_v<float> / static_cast<float>(fft_size);

        // 前后各 La 零填充（stretch 的 stftBatch 行为：frameStart 从 -N 开始）
        std::vector<float> padded(2 * La + input.size(), 0.0f);
        std::copy(input.begin(), input.end(), padded.begin() + static_cast<std::ptrdiff_t>(La));

        // 帧数：分析窗口整体落在 padded 内
        const size_t N = (padded.size() > La) ? (padded.size() - La) / ana_hop + 1 : 0;
        if (N == 0)
            return {};

        // 窗函数
        std::vector<float> window(La);
        detail::MakePeriodicHann(window);

        // 输出（win² 归一化 OLA）
        const size_t out_len = (N - 1) * syn_hop + La;
        std::vector<float> out(out_len, 0.0f);
        std::vector<float> norm(out_len, 0.0f);

        qwqdsp_spectral::RealFFT fft;
        fft.Init(fft_size);
        std::vector<float> fft_in(fft_size);
        std::vector<float> fft_out(fft_size + 2);

        std::vector<float> mag(bins);
        std::vector<float> ana_phase(bins);
        std::vector<float> p(bins);

        if (prev_.size() != bins) {
            prev_.assign(bins, 0.0f);
            syn_prev_.assign(bins, 0.0f);
            first_ = true;
        }

        for (size_t i = 0; i < N; ++i) {
            // ---- 分析 ----
            const float* frame = padded.data() + static_cast<std::ptrdiff_t>(i * ana_hop);
            for (size_t j = 0; j < La; ++j)
                fft_in[j] = frame[j] * window[j];

            fft.FFT(fft_in.data(), fft_out.data());

            for (size_t j = 0; j < bins; ++j) {
                const float re = fft_out[2 * j];
                const float im = fft_out[2 * j + 1];
                mag[j] = std::sqrt(re * re + im * im);
                ana_phase[j] = std::atan2(im, re);
            }

            // ---- 相位传播 ----
            if (first_) {
                p = ana_phase;
                first_ = false;
            }
            else {
                for (size_t j = 0; j < bins; ++j) {
                    const float dp = detail::WrapToPi(ana_phase[j] - prev_[j] - static_cast<float>(j) * freq_per_bin * static_cast<float>(ana_hop));
                    p[j] = syn_prev_[j] + (static_cast<float>(j) * freq_per_bin + dp / static_cast<float>(ana_hop)) * static_cast<float>(syn_hop);
                }
                detail::LockPhase(mag.data(), p.data(), ana_phase.data(), bins - 1);
            }

            prev_ = ana_phase;
            syn_prev_ = p;

            // ---- 综合（由幅度与锁定相位重建） ----
            for (size_t j = 0; j < bins; ++j) {
                fft_out[2 * j] = mag[j] * std::cos(p[j]);
                fft_out[2 * j + 1] = mag[j] * std::sin(p[j]);
            }

            fft.IFFT(fft_out.data(), fft_in.data());

            // OLA：加窗 + win² 归一化累加
            const size_t dst = i * syn_hop;
            for (size_t j = 0; j < La; ++j) {
                out[dst + j] += fft_in[j] * window[j];
                norm[dst + j] += window[j] * window[j];
            }
        }

        // ---- 归一化 + 输出偏移 ----
        // 输出偏移 pOut = round(half * (syn_hop/ana_hop + 1))，去掉前置填充的零段
        const size_t p_out = static_cast<size_t>(std::round(static_cast<float>(bins - 1) * (static_cast<float>(syn_hop) / static_cast<float>(ana_hop) + 1.0f)));
        const size_t want = static_cast<size_t>(std::round(static_cast<float>(input.size()) * static_cast<float>(syn_hop) / static_cast<float>(ana_hop)));

        std::vector<float> result;
        result.reserve(want);
        for (size_t i = 0; i < want; ++i) {
            const size_t idx = p_out + i;
            if (idx >= out.size())
                break;
            const float n = norm[idx];
            result.push_back(n > 1e-4f ? out[idx] / n : out[idx]);
        }

        // ---- 音高移动：重采样 ----
        if (kp != 1.0f) {
            result = detail::LinearResample(result, kp);
        }

        return result;
    }

private:
    size_t frame_size_ = kDefaultFrameSize;
    size_t hop_size_ = kDefaultHopSize;
    float kt_ = 1.0f;
    float kp_ = 1.0f;

    std::vector<float> prev_;     // 上一帧分析相位
    std::vector<float> syn_prev_; // 上一帧合成相位
    bool first_ = true;
};

} // namespace qwqdsp_test
