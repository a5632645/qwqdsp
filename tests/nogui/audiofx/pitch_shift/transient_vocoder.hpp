#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <qwqdsp/spectral/real_fft.hpp>
#include <span>
#include <vector>

#include "phase_locked_vocoder.hpp"

namespace qwqdsp_test {

// ------------------------------------------------------------
// TransientVocoder
// ------------------------------------------------------------
/**
 * @brief 瞬态感知相位锁定声码器（transient-aware vocoder）。
 *
 * 移植自 ref_repos/stretch 的 stretch-transient（Röbel 2003）：
 * 在相位锁定声码器之上测量逐帧谱通量，检测到突变的起始（瞬态）时
 * 把合成相位重置为原始分析相位，保留鼓点/拨弦等攻击的锐度；
 * 非瞬态帧仍走相位传播 + 相位锁定。
 *
 * 瞬态判定（与 stretch 一致）：
 * - flux = Σ max(0, log1p(mag[k]) - log1p(prev_mag[k]))
 * - energy = Σ (0.5 + 0.5·k/half) · log1p(mag[k])
 * - norm_flux = flux / energy
 * - 触发：frames > 4 && cooldown == 0
 *   && norm_flux > mean + threshold·max(0.07, std) && norm_flux > mean·1.35
 * - 均值/方差用 EMA 更新（瞬态帧 α=0.3，其余 α=0.12）
 */
class TransientVocoder {
public:
    static constexpr size_t kDefaultFrameSize = 2048;
    static constexpr size_t kDefaultHopSize = 512;
    static constexpr float kDefaultThreshold = 1.5f;

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

    /** 设置瞬态检测阈值（越大越少触发） */
    void SetTransientThreshold(float t) noexcept {
        threshold_ = t;
    }

    /** 清理相位传播与瞬态统计状态 */
    void Reset() noexcept {
        prev_.clear();
        syn_prev_.clear();
        prev_mag_.clear();
        flux_mean_ = 0.0f;
        flux_var_ = 0.0f;
        frames_ = 0;
        cooldown_ = 0;
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

        // 前后各 La 零填充
        std::vector<float> padded(2 * La + input.size(), 0.0f);
        std::copy(input.begin(), input.end(), padded.begin() + static_cast<std::ptrdiff_t>(La));

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
            prev_mag_.assign(bins, 0.0f);
            flux_mean_ = 0.0f;
            flux_var_ = 0.0f;
            frames_ = 0;
            cooldown_ = 0;
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

            // ---- 瞬态检测 ----
            const bool is_transient = DetectTransient(mag);

            // ---- 相位传播：瞬态帧重置为分析相位，其余传播 + 锁定 ----
            if (first_ || is_transient) {
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
            prev_mag_ = mag;

            // ---- 综合（由幅度与相位重建） ----
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
    /** 逐帧谱通量瞬态检测（stretch transient.js 的 updateFluxStats + 双阈值） */
    bool DetectTransient(const std::vector<float>& mag) {
        const size_t half = mag.size() - 1;
        const float threshold = threshold_;

        bool is_transient = false;
        if (!first_) {
            float flux = 0.0f;
            float energy = 0.0f;
            for (size_t k = 0; k <= half; ++k) {
                const float weight = 0.5f + 0.5f * static_cast<float>(k) / std::max<size_t>(1, half);
                const float d = std::log1p(mag[k]) - std::log1p(prev_mag_[k]);
                if (d > 0.0f)
                    flux += d;
                energy += weight * std::log1p(mag[k]);
            }
            const float norm_flux = (energy > 1e-10f) ? flux / energy : 0.0f;
            const float mean = (frames_ == 0 && !has_mean_) ? norm_flux : flux_mean_;
            const float std = std::sqrt(flux_var_);
            // Std floor 0.07：稳定的和弦拍频 normFlux ≤ ~0.07，真实起始 ≥ ~0.19
            is_transient = frames_ > 4 && cooldown_ == 0 && norm_flux > mean + threshold * std::max(0.07f, std)
                           && norm_flux > mean * 1.35f;
            UpdateFluxStats(norm_flux, is_transient ? 0.3f : 0.12f);
            cooldown_ = is_transient ? 1 : (cooldown_ > 0 ? cooldown_ - 1 : 0);
        }

        ++frames_;
        return is_transient;
    }

    /** EMA 更新谱通量均值与方差 */
    void UpdateFluxStats(float value, float alpha) noexcept {
        if (!has_mean_) {
            flux_mean_ = value;
            flux_var_ = 0.0f;
            has_mean_ = true;
            return;
        }
        const float delta = value - flux_mean_;
        flux_mean_ += alpha * delta;
        flux_var_ = (1.0f - alpha) * (flux_var_ + alpha * delta * delta);
    }

    size_t frame_size_ = kDefaultFrameSize;
    size_t hop_size_ = kDefaultHopSize;
    float kt_ = 1.0f;
    float kp_ = 1.0f;
    float threshold_ = kDefaultThreshold;

    std::vector<float> prev_;     // 上一帧分析相位
    std::vector<float> syn_prev_; // 上一帧合成相位
    std::vector<float> prev_mag_; // 上一帧幅度
    float flux_mean_ = 0.0f;
    float flux_var_ = 0.0f;
    bool has_mean_ = false;
    size_t frames_ = 0;
    size_t cooldown_ = 0;
    bool first_ = true;
};

} // namespace qwqdsp_test
