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

        // 缓冲在原始 FFT bin 网格上行数 = binSize_（每个 FFT bin 一行），宽度 = 子列数
        col_buf_.resize(subColumns_ * binSize_, 0.0f);
        colMax_.resize(outputHeight_);
        column_.resize(outputHeight_);
        rowToBinLo_.resize(outputHeight_);
        rowToBinHi_.resize(outputHeight_);

        // ── 预计算每个频率像素行 y 覆盖的 FFT bin 区间 [k_lo, k_hi] ──
        // 像素 y 覆盖 [f_bot, f_top)（y 越小频率越高）；FFT bin k 覆盖 [k·bin_hz, (k+1)·bin_hz]
        const float bin_hz_i = static_cast<float>(sampleRate_) / static_cast<float>(fftLen_);
        const float log_step_i = (logMax_ - logMin_) / static_cast<float>(outputHeight_);
        for (int y = 0; y < outputHeight_; ++y) {
            float f_top = std::pow(10.0f, logMax_ - static_cast<float>(y) * log_step_i);
            float f_bot = std::pow(10.0f, logMax_ - static_cast<float>(y + 1) * log_step_i);
            int k_top = static_cast<int>(std::round(f_top / bin_hz_i));
            int k_bot = static_cast<int>(std::round(f_bot / bin_hz_i));
            int k_lo = std::min(k_top, k_bot);
            int k_hi = std::max(k_top, k_bot);
            rowToBinLo_[y] = std::clamp(k_lo, 0, binSize_ - 1);
            rowToBinHi_[y] = std::clamp(k_hi, 0, binSize_ - 1);
            if (rowToBinLo_[y] > rowToBinHi_[y])
                std::swap(rowToBinLo_[y], rowToBinHi_[y]);
        }
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

        // ── 3. 遍历每个 bin: 时间重分配（水平子列），行 = 原始 FFT bin k ──
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

            // ── 水平位置（子列）——行 = 原始 FFT bin k，不做垂直 Y 分布 ──
            float c_pos = (group_delay + 0.5f) * subColScale_;
            c_pos = std::clamp(c_pos, 0.0f, static_cast<float>(subColumns_ - 0.01f));
            int c_idx = static_cast<int>(std::floor(c_pos));
            float c_frac = c_pos - static_cast<float>(c_idx);
            bool c_last = (c_idx >= subColumns_ - 1);

            if (c_last) {
                col_buf_[c_idx * binSize_ + k] += mag_lin;
            }
            else {
                col_buf_[c_idx * binSize_ + k] += mag_lin * (1.0f - c_frac);
                col_buf_[(c_idx + 1) * binSize_ + k] += mag_lin * c_frac;
            }
        }

        // ── 4. 弹出最旧子列（原始 FFT 网格）→ max-hold 映射频率网格 → Color ──
        constexpr float kEps = 1e-12f;
        const int subColStride = binSize_;

        std::fill(colMax_.begin(), colMax_.end(), dbFloor_);
        for (int y = 0; y < outputHeight_; ++y) {
            const int k_lo = rowToBinLo_[y];
            const int k_hi = rowToBinHi_[y];
            float best = dbFloor_;
            for (int k = k_lo; k <= k_hi; ++k) {
                const float v = col_buf_[k];   // 第 0 个子列、FFT bin k
                if (v > best)
                    best = v;
            }
            float dB = 20.0f * std::log10(best + kEps);
            dB = std::clamp(dB, dbFloor_, 0.0f);
            int idx = static_cast<int>((dB - dbFloor_) / (-dbFloor_) * 255.0f);
            idx = std::clamp(idx, 0, 255);
            column_[y] = Colormap::kTable[idx];
        }

        // 左移一个子列并清空最后一列
        std::move(col_buf_.begin() + subColStride, col_buf_.end(), col_buf_.begin());
        std::fill(col_buf_.end() - subColStride, col_buf_.end(), 0.0f);
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
    std::vector<std::complex<float>> X_h_, X_pf_;
    std::vector<float> col_buf_;   // [subCol * binSize_ + k]，原始 FFT bin 网格，行 = FFT bin
    std::vector<float> colMax_;    // [outputHeight_]，频率网格每行 max dB
    std::vector<Color> column_;
    std::vector<int> rowToBinLo_;  // [outputHeight_]，每行覆盖的 FFT bin 区间下界
    std::vector<int> rowToBinHi_;  // [outputHeight_]，每行覆盖的 FFT bin 区间上界
};
