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

        // 缓冲在原始 NC 网格上行数 = ncSize_（每个 NC bin 一行），宽度 = 子列数
        col_buf_.resize(subColumns_ * ncSize_, 0.0f);
        colMax_.resize(outputHeight_);
        column_.resize(outputHeight_);
        rowToBinLo_.resize(outputHeight_);
        rowToBinHi_.resize(outputHeight_);

        // ── 预计算每个频率像素行 y 覆盖的 NC bin 区间 [k_lo, k_hi] ──
        // 像素 y 覆盖 [f_bot, f_top)（y 越小频率越高）；NC bin k 覆盖 [k·bin_hz, (k+zeroPad_)·bin_hz]
        const float bin_hz_i = static_cast<float>(sampleRate_) / static_cast<float>(fftLen_);
        const float log_step_i = (logMax_ - logMin_) / static_cast<float>(outputHeight_);
        for (int y = 0; y < outputHeight_; ++y) {
            float f_top = std::pow(10.0f, logMax_ - static_cast<float>(y) * log_step_i);
            float f_bot = std::pow(10.0f, logMax_ - static_cast<float>(y + 1) * log_step_i);
            int k_top = static_cast<int>(std::round(f_top / bin_hz_i));
            int k_bot = static_cast<int>(std::round(f_bot / bin_hz_i));
            int k_lo = std::min(k_top, k_bot);
            int k_hi = std::max(k_top, k_bot);
            rowToBinLo_[y] = std::clamp(k_lo, 0, ncSize_ - 1);
            rowToBinHi_[y] = std::clamp(k_hi, 0, ncSize_ - 1);
            if (rowToBinLo_[y] > rowToBinHi_[y])
                std::swap(rowToBinLo_[y], rowToBinHi_[y]);
        }
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

            // ── 水平位置（子列）——时间重分配只做水平分布，行 = 原始 NC bin k ──
            float c_pos = (group_delay + 0.5f) * subColScale_;
            c_pos = std::clamp(c_pos, 0.0f, static_cast<float>(subColumns_ - 0.01f));
            int c_idx = static_cast<int>(std::floor(c_pos));
            float c_frac = c_pos - static_cast<float>(c_idx);
            bool c_last = (c_idx >= subColumns_ - 1);

            if (c_last) {
                col_buf_[c_idx * ncSize_ + k] += gain;
            }
            else {
                col_buf_[c_idx * ncSize_ + k] += gain * (1.0f - c_frac);
                col_buf_[(c_idx + 1) * ncSize_ + k] += gain * c_frac;
            }
        }

        // ── 5. 弹出最旧子列（原始 NC 网格）→ max-hold 映射到频率网格 → Color ──
        //    grid[k] = col_buf_[k]（col_buf_ 布局 [subCol * ncSize_ + k]，第 0 个子列即最旧列）
        constexpr float kEpsDb = 1e-12f;
        const int subColStride = ncSize_;

        std::fill(colMax_.begin(), colMax_.end(), dbFloor_);
        for (int y = 0; y < outputHeight_; ++y) {
            const int k_lo = rowToBinLo_[y];
            const int k_hi = rowToBinHi_[y];
            float best = dbFloor_;
            for (int k = k_lo; k <= k_hi; ++k) {
                const float v = col_buf_[k];   // 第 0 个子列、NC bin k
                if (v > best)
                    best = v;
            }
            float dB = 20.0f * std::log10(best + kEpsDb);
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
    int sampleRate_{}, fftSize_{}, hopSize_{}, zeroPad_{}, fftLen_{}, binSize_{}, ncSize_{};
    int subColumns_{}, outputHeight_{};
    float freqMin_{}, freqMax_{}, logMin_{}, logMax_{}, dbFloor_{}, subColScale_{};

    qwqdsp_spectral::RealFftAdv fft_;
    std::vector<float> fft_in_;
    std::vector<std::complex<float>> X_, X_pf_;
    std::vector<float> col_buf_;   // [subCol * ncSize_ + k]，原始 NC 网格，行 = NC bin
    std::vector<float> colMax_;    // [outputHeight_]，频率网格每行 max dB
    std::vector<Color> column_;
    std::vector<int> rowToBinLo_;  // [outputHeight_]，每行覆盖的 NC bin 区间下界
    std::vector<int> rowToBinHi_;  // [outputHeight_]，每行覆盖的 NC bin 区间上界
};
