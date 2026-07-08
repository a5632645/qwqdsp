#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <numbers>
#include <span>
#include <vector>

#include <qwqdsp/spectral/real_fft.hpp>

#include "raylib.h"

/**
 * @brief 原始帧+窗 → 时间重分配（roll 法）→ 颜色输出
 *
 * 不做频率重分配，Y 轴固定在 bin 中心频率。
 * 群延迟通过 roll(X_h,1) 的相位差估计，能量按 2D 双线性分布到滚动缓存。
 *
 * X_h  = FFT(x·w)
 * X_pf = roll(X_h,1), X_pf[0]=0   (时间重分配 → 群延迟)
 */
template <typename Colormap>
struct TimeReassignmentFrame {
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
        X_h_.resize(binSize_);
        X_pf_.resize(binSize_);

        col_buf_.resize(subColumns_ * outputHeight_, 0.0f);
        column_.resize(outputHeight_);
    }

    void Process(std::span<const float> /*raw_frame*/, std::span<const float> /*window*/,
                 std::span<const float> windowed_frame) noexcept {
        // ── 1. X_h = FFT(x * w) ──
        std::copy(windowed_frame.begin(), windowed_frame.end(), fft_in_.begin());
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_h_);

        // ── 2. X_pf = roll(X_h, 1); X_pf[0] = 0 ──
        std::copy(X_h_.begin(), X_h_.end() - 1, X_pf_.begin() + 1);
        X_pf_.front() = {};

        // ── 3. 遍历每个 bin: 时间重分配 + 2D 分布 ──
        const float bin_hz = static_cast<float>(sampleRate_) / static_cast<float>(fftLen_);

        for (int k = 0; k < binSize_; ++k) {
            float mag_lin = std::abs(X_h_[k]);
            if (mag_lin < 1e-8f)
                continue;

            // ── 时间重分配: group_delay ∈ (-0.5, 0.5] ──
            auto cross_f = X_h_[k] * std::conj(X_pf_[k]);
            float arg_f = std::arg(cross_f) / (2.0f * std::numbers::pi_v<float>);
            while (arg_f < 0.0f)
                arg_f += 1.0f;
            arg_f = std::fmod(arg_f, 1.0f);
            float group_delay = 0.5f - arg_f; // (-0.5, 0.5]

            // ── Y 范围: bin k 覆盖 [k·bin_hz, (k+1)·bin_hz] ──
            float f_low = static_cast<float>(k) * bin_hz;
            float f_high = static_cast<float>(k + 1) * bin_hz;
            if (f_high < freqMin_ || f_low > freqMax_)
                continue;
            f_low = std::max(f_low, freqMin_);
            f_high = std::min(f_high, freqMax_);

            float n_low = (std::log10(f_low) - logMin_) / (logMax_ - logMin_);
            float n_high = (std::log10(f_high) - logMin_) / (logMax_ - logMin_);
            float y_top = static_cast<float>(outputHeight_ - 1) * (1.0f - n_high);
            float y_bot = static_cast<float>(outputHeight_ - 1) * (1.0f - n_low);
            y_top = std::clamp(y_top, 0.0f, static_cast<float>(outputHeight_ - 1));
            y_bot = std::clamp(y_bot, 0.0f, static_cast<float>(outputHeight_ - 1));

            int y0_idx = static_cast<int>(std::floor(y_top));
            int y1_idx = static_cast<int>(std::floor(y_bot));
            float y_range = y_bot - y_top;
            if (y_range < 1.0f)
                y_range = 1.0f;

            // ── 水平位置 ──
            float c_pos = (group_delay + 0.5f) * subColScale_;
            c_pos = std::clamp(c_pos, 0.0f, static_cast<float>(subColumns_ - 0.01f));
            int c_idx = static_cast<int>(std::floor(c_pos));
            float c_frac = c_pos - static_cast<float>(c_idx);
            bool c_last = (c_idx >= subColumns_ - 1);

            // ── 按频率区间比例分布到各 Y 行 ──
            for (int y = y0_idx; y <= y1_idx; ++y) {
                if (y < 0 || y >= outputHeight_)
                    continue;
                float row_y0 = std::max(static_cast<float>(y), y_top);
                float row_y1 = std::min(static_cast<float>(y + 1), y_bot);
                float y_weight = (row_y1 - row_y0) / y_range;

                float energy = mag_lin * y_weight;

                if (c_last) {
                    col_buf_[c_idx * outputHeight_ + y] += energy;
                }
                else {
                    col_buf_[c_idx * outputHeight_ + y] += energy * (1.0f - c_frac);
                    col_buf_[(c_idx + 1) * outputHeight_ + y] += energy * c_frac;
                }
            }
        }

        // ── 4. 弹出最旧子列 → dB → Color ──
        constexpr float kEps = 1e-12f;
        for (int y = 0; y < outputHeight_; ++y) {
            float dB = 20.0f * std::log10(col_buf_[y] + kEps);
            dB = std::clamp(dB, dbFloor_, 0.0f);
            int idx = static_cast<int>((dB - dbFloor_) / (-dbFloor_) * 255.0f);
            idx = std::clamp(idx, 0, 255);
            column_[y] = Colormap::kTable[idx];
        }

        std::memmove(col_buf_.data(), col_buf_.data() + outputHeight_,
                     (subColumns_ - 1) * outputHeight_ * sizeof(float));
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

    qwqdsp_spectral::RealFFT fft_;
    std::vector<float> fft_in_;
    std::vector<std::complex<float>> X_h_, X_pf_;
    std::vector<float> col_buf_; // [subCol * outputHeight + y]
    std::vector<Color> column_;
};
