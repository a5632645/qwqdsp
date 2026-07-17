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
 * @brief NC 方法 + 时间重分配（roll 法）
 *
 * 矩形窗 FFT，NC 值计算相邻 bin 负相关，群延迟取中间 bin (k + zeroPad_/2) 做时间重分配。
 *
 * X   = FFT(x)                          (矩形窗)
 * X_pf = roll(X, 1), X_pf[0] = 0        (时间重分配 → 群延迟)
 *
 * NC bin k 对应 bin 对 (k, k+zeroPad_)：
 *   ncSum = -(Re(X[k])·Re(X[k+zeroPad_]) + Im(X[k])·Im(X[k+zeroPad_]))
 *
 * Y 范围 [k·bin_hz, (k+zeroPad_)·bin_hz]，
 * 能量 (sqrt(ncSum)) 按群延迟水平分布到滚动缓存。
 *
 * @ref https://arxiv.org/html/2410.07982v3#bib.bib1
 */
template <typename Colormap>
struct NcTimeReassignmentFrame {
    void Init(int sampleRate, int fftSize, int hopSize, int zeroPad, int outputHeight, float freqMin, float freqMax,
              float dbFloor) noexcept {
        sampleRate_ = sampleRate;
        fftSize_ = fftSize;
        hopSize_ = hopSize;
        zeroPad_ = zeroPad;
        fftLen_ = fftSize * zeroPad;
        binSize_ = fftLen_ / 2 + 1;
        ncSize_ = binSize_ - zeroPad_;
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
        X_.resize(binSize_);
        X_pf_.resize(binSize_);

        col_buf_.resize(subColumns_ * outputHeight_, 0.0f);
        column_.resize(outputHeight_);
    }

    void Process(std::span<const float> raw_frame, std::span<const float> /*window*/,
                 std::span<const float> /*windowed_frame*/) noexcept {
        // ── 1. 矩形窗 FFT ──
        std::copy(raw_frame.begin(), raw_frame.end(), fft_in_.begin());
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_);

        // ── 2. 归一化增益 2/fftSize ──
        const float norm = 2.0f / static_cast<float>(fftSize_);
        for (int k = 0; k < binSize_; ++k) {
            X_[k] *= norm;
        }

        // ── 3. X_pf = roll(X, 1); X_pf[0] = 0 ──
        std::copy(X_.begin(), X_.end() - 1, X_pf_.begin() + 1);
        X_pf_.front() = {};

        // ── 4. 计算 NC + 时间重分配（只叠加，不做加权归一化）──
        constexpr float kEps = 1e-18f;
        const float bin_hz = static_cast<float>(sampleRate_) / static_cast<float>(fftLen_);

        for (int k = 0; k < ncSize_; ++k) {
            // ── NC gain（线性）──
            float ncSum = -(X_[k].real() * X_[k + zeroPad_].real() + X_[k].imag() * X_[k + zeroPad_].imag());
            if (ncSum < 0.0f)
                continue;
            float gain = std::sqrt(ncSum + kEps);

            // ── 时间重分配: 中间 bin (k + zeroPad_/2) 的群延迟 ──
            int mid = k + zeroPad_ / 2;
            auto cross_f = X_[mid] * std::conj(X_pf_[mid]);
            float arg_f = std::arg(cross_f) / (2.0f * std::numbers::pi_v<float>);
            while (arg_f < 0.0f)
                arg_f += 1.0f;
            arg_f = std::fmod(arg_f, 1.0f);
            float group_delay = 0.5f - arg_f; // (-0.5, 0.5]

            // ── Y 范围: NC bin k 覆盖 [k·bin_hz, (k+zeroPad_)·bin_hz] ──
            float f_low = static_cast<float>(k) * bin_hz;
            float f_high = static_cast<float>(k + zeroPad_) * bin_hz;
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

            // ── 按频率区间比例 + 水平双线性分布到缓存 ──
            for (int y = y0_idx; y <= y1_idx; ++y) {
                if (y < 0 || y >= outputHeight_)
                    continue;
                float row_y0 = std::max(static_cast<float>(y), y_top);
                float row_y1 = std::min(static_cast<float>(y + 1), y_bot);
                float y_weight = (row_y1 - row_y0) / y_range;

                float energy = gain * y_weight;

                if (c_last) {
                    col_buf_[c_idx * outputHeight_ + y] += energy;
                }
                else {
                    col_buf_[c_idx * outputHeight_ + y] += energy * (1.0f - c_frac);
                    col_buf_[(c_idx + 1) * outputHeight_ + y] += energy * c_frac;
                }
            }
        }

        // ── 5. 弹出最旧子列 → gain 转 dB → Color ──
        constexpr float kEpsDb = 1e-12f;
        for (int y = 0; y < outputHeight_; ++y) {
            float dB = 20.0f * std::log10(col_buf_[y] + kEpsDb);
            dB = std::clamp(dB, dbFloor_, 0.0f);
            int idx = static_cast<int>((dB - dbFloor_) / (-dbFloor_) * 255.0f);
            idx = std::clamp(idx, 0, 255);
            column_[y] = Colormap::kTable[idx];
        }

        std::move(col_buf_.begin() + outputHeight_, col_buf_.end(), col_buf_.begin());
        std::fill(col_buf_.begin() + (subColumns_ - 1) * outputHeight_, col_buf_.end(), 0.0f);
    }

    std::span<const Color> GetColumn() const noexcept {
        return {column_.data(), static_cast<size_t>(outputHeight_)};
    }

    int ColumnHeight() const noexcept {
        return outputHeight_;
    }
private:
    int sampleRate_{}, fftSize_{}, hopSize_{}, zeroPad_{}, fftLen_{}, binSize_{}, ncSize_{};
    int subColumns_{}, outputHeight_{};
    float freqMin_{}, freqMax_{}, logMin_{}, logMax_{}, dbFloor_{}, subColScale_{};

    qwqdsp_spectral::RealFftAdv fft_;
    std::vector<float> fft_in_;
    std::vector<std::complex<float>> X_, X_pf_;
    std::vector<float> col_buf_; // [subCol * outputHeight + y]
    std::vector<Color> column_;
};
