#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <numbers>
#include <span>
#include <vector>

#include <qwqdsp/spectral/real_fft_adv.hpp>

#include "raylib.h"

/**
 * @brief 原始帧+窗 → 频率重分配 + 时间重分配 → magma 颜色输出
 *
 * 内部维护 subColumns 个子列缓存，每帧将能量按 (inst_freq, group_delay)
 * 2D 双线性分布到缓存，然后弹出最旧子列。
 *
 * X_h = FFT(x·w)
 * X_t = FFT(x[n-1]·w[n])           (频率重分配)
 * X_pf = roll(X_h,1), X_pf[0]=0    (时间重分配 → 群延迟)
 */
template <typename Colormap>
struct TfReassignmentFrame {
    void Init(int sampleRate, int fftSize, int hopSize, int zeroPad, int outputHeight, float freqMin, float freqMax,
              float dbFloor) noexcept {
        sampleRate_ = sampleRate;
        fftSize_ = fftSize;
        hopSize_ = hopSize;
        zeroPad_ = zeroPad;
        fftLen_ = fftSize * zeroPad;
        binSize_ = fftLen_ / 2 + 1;
        subColumns_ = fftSize / hopSize;
        outputHeight_ = outputHeight;
        freqMin_ = freqMin;
        freqMax_ = freqMax;
        dbFloor_ = dbFloor;
        logMin_ = std::log10(freqMin);
        logMax_ = std::log10(freqMax);
        subColScale_ = static_cast<float>(subColumns_);

        fft_.Init(fftLen_);

        fft_in_.resize(fftLen_, 0.0f);
        shift_in_.resize(fftSize_);
        X_h_.resize(binSize_);
        X_t_.resize(binSize_);
        X_pf_.resize(binSize_);

        col_buf_.resize(subColumns_ * outputHeight_, 0.0f);
        column_.resize(outputHeight_);
    }

    void Process(std::span<const float> raw_frame, std::span<const float> window,
                 std::span<const float> windowed_frame) noexcept {
        // ── 1. X_h = FFT(x * w) ──
        std::copy(windowed_frame.begin(), windowed_frame.end(), fft_in_.begin());
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_h_);

        // ── 2. X_t = FFT(x[n-1] * w[n]) ──
        shift_in_[0] = 0.0f;
        for (int i = 1; i < fftSize_; ++i)
            shift_in_[i] = raw_frame[i - 1] * window[i];
        std::copy(shift_in_.begin(), shift_in_.end(), fft_in_.begin());
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_t_);

        // ── 3. X_pf = roll(X_h, 1); X_pf[0] = 0 ──
        std::copy(X_h_.begin(), X_h_.end() - 1, X_pf_.begin() + 1);
        X_pf_.front() = {};

        // ── 4. 遍历每个 bin: 2D 重分配 ──
        const float two_pi = 2.0f * std::numbers::pi_v<float>;

        for (int k = 0; k < binSize_; ++k) {
            float mag_lin = std::abs(X_h_[k]);
            if (mag_lin < 1e-8f)
                continue;

            // ── 频率重分配: inst_freq_hz ──
            auto cross_t = X_h_[k] * std::conj(X_t_[k]);
            float inst_freq_norm = std::arg(cross_t) / two_pi;
            while (inst_freq_norm < 0.0f)
                inst_freq_norm += 1.0f;
            inst_freq_norm = std::fmod(inst_freq_norm, 1.0f);
            float inst_freq_hz = inst_freq_norm * sampleRate_;
            if (inst_freq_hz < freqMin_ || inst_freq_hz > freqMax_)
                continue;

            // ── 时间重分配: group_delay ∈ (-0.5, 0.5] ──
            auto cross_f = X_h_[k] * std::conj(X_pf_[k]);
            float arg_f = std::arg(cross_f) / two_pi;
            while (arg_f < 0.0f)
                arg_f += 1.0f;
            arg_f = std::fmod(arg_f, 1.0f);
            float group_delay = 0.5f - arg_f; // (-0.5, 0.5]

            // ── 映射到子列 + Y 像素 ──
            float logF = std::log10(inst_freq_hz);
            float norm = (logF - logMin_) / (logMax_ - logMin_);
            float y_pos = static_cast<float>(outputHeight_ - 1) * (1.0f - norm);

            float c_pos = (group_delay + 0.5f) * subColScale_;
            // ── 2D 双线性分布 ──
            y_pos = std::clamp(y_pos, 0.0f, static_cast<float>(outputHeight_ - 0.01f));
            int y_idx = static_cast<int>(std::floor(y_pos));
            bool y_last = (y_idx >= outputHeight_ - 1);

            c_pos = std::clamp(c_pos, 0.0f, static_cast<float>(subColumns_ - 0.01f));
            int c_idx = static_cast<int>(std::floor(c_pos));
            bool c_last = (c_idx >= subColumns_ - 1);

            // 处理边界特殊情况
            if (y_last && c_last) {
                col_buf_[c_idx * outputHeight_ + y_idx] += mag_lin;
            }
            else if (y_last) {
                float c_frac = c_pos - static_cast<float>(c_idx);
                col_buf_[c_idx * outputHeight_ + y_idx] += mag_lin * (1.0f - c_frac);
                col_buf_[(c_idx + 1) * outputHeight_ + y_idx] += mag_lin * c_frac;
            }
            else if (c_last) {
                float y_frac = y_pos - static_cast<float>(y_idx);
                col_buf_[c_idx * outputHeight_ + y_idx] += mag_lin * (1.0f - y_frac);
                col_buf_[c_idx * outputHeight_ + y_idx + 1] += mag_lin * y_frac;
            }
            else {
                float c_frac = c_pos - static_cast<float>(c_idx);
                float y_frac = y_pos - static_cast<float>(y_idx);
                float wgt = mag_lin;
                col_buf_[c_idx * outputHeight_ + y_idx] += wgt * (1.0f - c_frac) * (1.0f - y_frac);
                col_buf_[(c_idx + 1) * outputHeight_ + y_idx] += wgt * c_frac * (1.0f - y_frac);
                col_buf_[c_idx * outputHeight_ + y_idx + 1] += wgt * (1.0f - c_frac) * y_frac;
                col_buf_[(c_idx + 1) * outputHeight_ + y_idx + 1] += wgt * c_frac * y_frac;
            }
        }

        // ── 5. 弹出最旧子列 → dB → Color ──
        constexpr float kEps = 1e-12f;
        for (int y = 0; y < outputHeight_; ++y) {
            float dB = 20.0f * std::log10(col_buf_[y] + kEps);
            dB = std::clamp(dB, dbFloor_, 0.0f);
            int idx = static_cast<int>((dB - dbFloor_) / (-dbFloor_) * 255.0f);
            idx = std::clamp(idx, 0, 255);
            column_[y] = Colormap::kTable[idx];
        }

        // 左移: col_buf_[0..n-2] = col_buf_[1..n-1]
        std::memmove(col_buf_.data(), col_buf_.data() + outputHeight_,
                     (subColumns_ - 1) * outputHeight_ * sizeof(float));
        // 清零新列
        std::fill(col_buf_.begin() + (subColumns_ - 1) * outputHeight_, col_buf_.end(), 0.0f);
    }

    std::span<const Color> GetColumn() const noexcept {
        return {column_.data(), static_cast<size_t>(outputHeight_)};
    }

    int ColumnHeight() const noexcept {
        return outputHeight_;
    }
private:
    int sampleRate_{}, fftSize_{}, hopSize_{}, zeroPad_{}, fftLen_{}, binSize_{};
    int subColumns_{}, outputHeight_{};
    float freqMin_{}, freqMax_{}, logMin_{}, logMax_{}, dbFloor_{}, subColScale_{};

    qwqdsp_spectral::RealFftAdv fft_;
    std::vector<float> fft_in_;
    std::vector<float> shift_in_;
    std::vector<std::complex<float>> X_h_, X_t_, X_pf_;
    std::vector<float> col_buf_; // [subCol * outputHeight + y]
    std::vector<Color> column_;
};
