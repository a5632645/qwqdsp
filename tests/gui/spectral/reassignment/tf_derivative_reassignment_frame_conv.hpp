#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <numbers>
#include <span>
#include <vector>

#include <qwqdsp/spectral/real_fft_adv.hpp>
#include <qwqdsp/window/blackman_harris.hpp>
#include <qwqdsp/window/blackman_harris_3term.hpp>
#include <qwqdsp/window/window.hpp>
#include <qwqdsp/window/helper.hpp>

#include "raylib.h"

/**
 * @brief 原始帧+窗 → 导数窗频率重分配 + 时间矩重分配 → magma 颜色输出
 *
 * 复制自 TfReassignmentFrame，但不使用 1-sample 时移和频率轴 roll。
 * 频率位置通过导数窗 x·dw 估计，时间位置通过时间加权窗 x·(n-center)w 估计。
 *
 * X_h  = FFT(x·w)
 * X_dh = FFT(x·dw)               (导数窗频率重分配)
 * X_th = FFT(x·(n-center)w[n])   (时间重分配)
 */
template <typename Colormap, float kMinWeight = 0.1f>
struct TfDerivativeReassignmentFrameConv {
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
        window_.resize(fftSize_);
        dwindow_.resize(fftSize_);
        twindow_.resize(fftSize_);
        X_h_.resize(binSize_);
        X_dh_.resize(binSize_);
        X_th_.resize(binSize_);
        X_tdh_.resize(binSize_);
        conv_.resize(binSize_);

        col_buf_.resize(subColumns_ * outputHeight_, 0.0f);
        weight_buf_.resize(subColumns_ * outputHeight_, 0.0f);
        column_.resize(outputHeight_);

        InitDerivativeWindow();
    }

    void Process(std::span<const float> raw_frame, std::span<const float> /*window*/,
                 std::span<const float> /*windowed_frame*/) noexcept {
        // ── 1. X_h = FFT(x * w) ──
        for (int i = 0; i < fftSize_; ++i)
            fft_in_[i] = raw_frame[i] * window_[i];
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_h_);

        // ── 2. X_dh = FFT(x * dw) ──
        for (int i = 0; i < fftSize_; ++i)
            fft_in_[i] = raw_frame[i] * dwindow_[i];
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_dh_);

        // ── 3. X_th = FFT(x * (n-center)w[n]) ──
        for (int i = 0; i < fftSize_; ++i)
            fft_in_[i] = raw_frame[i] * twindow_[i];
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_th_);

        // ── 3.5. X_tdh = FFT(x * tdw_l) ──
        for (int i = 0; i < fftSize_; ++i)
            fft_in_[i] = raw_frame[i] * tdwindow_[i];
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_tdh_);

        // ── 4. 收敛指示器 (conv, 不影响输出) ──
        {
            constexpr float kEpsConv = 1e-20f;
            const float scale_by = 2.0f * std::numbers::pi_v<float> / static_cast<float>(fftSize_);
            for (int k = 0; k < binSize_; ++k) {
                auto& xh = X_h_[k];
                auto& xth = X_th_[k];
                auto& xtdh = X_tdh_[k];

                float msq = std::max(std::norm(xh), kEpsConv);

                // term1 = Re(X_tdh · conj(X_h)) / |X_h|²
                float term1 = (xtdh.real() * xh.real() + xtdh.imag() * xh.imag()) / msq;

                // Xdh_l = X_dh / 2π  (dwindow_ 比 dw_l 大 2π 倍)
                // term2 = Re((X_th · Xdh_l) / X_h²)
                std::complex<float> xdh_l = X_dh_[k] / (2.0f * std::numbers::pi_v<float>);
                std::complex<float> num2 = xth * xdh_l;
                std::complex<float> den2 = xh * xh;
                den2.real(den2.real() + kEpsConv);
                float term2 = (num2 / den2).real();

                float c = std::abs(1.0f + scale_by * (term1 - term2));
                conv_[k] = std::clamp(c, 0.0f, 1.0f);
            }
        }

        // ── 5. 遍历每个 bin: 2D 重分配 ──
        constexpr float kEps = 1e-20f;
        const float bin_hz = static_cast<float>(sampleRate_) / static_cast<float>(fftLen_);

        for (int k = 0; k < binSize_; ++k) {
            const float mag_lin = std::abs(X_h_[k]);
            if (mag_lin < 1e-8f)
                continue;

            if (conv_[k] > 0.3f && conv_[k] < 0.7f) continue;

            const float mag_sq = std::max(std::norm(X_h_[k]), kEps);

            // ── 导数窗频率重分配: inst_freq_hz ──
            auto xdh = X_dh_[k] / (2.0f * std::numbers::pi_v<float>);
            auto const& xh = X_h_[k];
            float up = xdh.imag() * xh.real() - xdh.real() * xh.imag();
            float freq_c = -up / mag_sq * zeroPad_;
            float inst_freq_hz = (static_cast<float>(k) + freq_c) * bin_hz;
            if (inst_freq_hz < freqMin_ || inst_freq_hz > freqMax_)
                continue;

            // ── 时间矩重分配: group_delay ∈ [-0.5, 0.5] ──
            float time_num = xh.real() * X_th_[k].real() + xh.imag() * X_th_[k].imag();
            float time_offset = time_num / mag_sq;
            float group_delay = time_offset / static_cast<float>(fftSize_);
            group_delay = std::clamp(group_delay, -0.5f, 0.5f);

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

            if (y_last && c_last) {
                col_buf_[c_idx * outputHeight_ + y_idx] += mag_lin;
                weight_buf_[c_idx * outputHeight_ + y_idx] += 1.0f;
            }
            else if (y_last) {
                float c_frac = c_pos - static_cast<float>(c_idx);
                col_buf_[c_idx * outputHeight_ + y_idx] += mag_lin * (1.0f - c_frac);
                col_buf_[(c_idx + 1) * outputHeight_ + y_idx] += mag_lin * c_frac;
                weight_buf_[c_idx * outputHeight_ + y_idx] += 1.0f;
                weight_buf_[(c_idx + 1) * outputHeight_ + y_idx] += 1.0f;
            }
            else if (c_last) {
                float y_frac = y_pos - static_cast<float>(y_idx);
                col_buf_[c_idx * outputHeight_ + y_idx] += mag_lin * (1.0f - y_frac);
                col_buf_[c_idx * outputHeight_ + y_idx + 1] += mag_lin * y_frac;
                weight_buf_[c_idx * outputHeight_ + y_idx] += (1.0f - y_frac);
                weight_buf_[c_idx * outputHeight_ + y_idx + 1] += y_frac;
            }
            else {
                float c_frac = c_pos - static_cast<float>(c_idx);
                float y_frac = y_pos - static_cast<float>(y_idx);
                col_buf_[c_idx * outputHeight_ + y_idx] += mag_lin * (1.0f - c_frac) * (1.0f - y_frac);
                col_buf_[(c_idx + 1) * outputHeight_ + y_idx] += mag_lin * c_frac * (1.0f - y_frac);
                col_buf_[c_idx * outputHeight_ + y_idx + 1] += mag_lin * (1.0f - c_frac) * y_frac;
                col_buf_[(c_idx + 1) * outputHeight_ + y_idx + 1] += mag_lin * c_frac * y_frac;
                weight_buf_[c_idx * outputHeight_ + y_idx] += (1.0f - y_frac);
                weight_buf_[(c_idx + 1) * outputHeight_ + y_idx] += (1.0f - y_frac);
                weight_buf_[c_idx * outputHeight_ + y_idx + 1] += y_frac;
                weight_buf_[(c_idx + 1) * outputHeight_ + y_idx + 1] += y_frac;
            }
        }

        // ── 6. 弹出最旧子列 → dB → Color ──
        constexpr float kDbEps = 1e-12f;
        for (int y = 0; y < outputHeight_; ++y) {
            float dB;
            if (weight_buf_[y] < kMinWeight) {
                dB = dbFloor_;
            } else {
                float avg = col_buf_[y] / weight_buf_[y];
                dB = 20.0f * std::log10(avg + kDbEps);
            }
            dB = std::clamp(dB, dbFloor_, 0.0f);
            int idx = static_cast<int>((dB - dbFloor_) / (-dbFloor_) * 255.0f);
            idx = std::clamp(idx, 0, 255);
            column_[y] = Colormap::kTable[idx];
        }

        std::move(col_buf_.begin() + outputHeight_, col_buf_.end(), col_buf_.begin());
        std::move(weight_buf_.begin() + outputHeight_, weight_buf_.end(), weight_buf_.begin());
        std::fill(col_buf_.begin() + (subColumns_ - 1) * outputHeight_, col_buf_.end(), 0.0f);
        std::fill(weight_buf_.begin() + (subColumns_ - 1) * outputHeight_, weight_buf_.end(), 0.0f);
    }

    std::span<const Color> GetColumn() const noexcept {
        return {column_.data(), static_cast<size_t>(outputHeight_)};
    }

    int ColumnHeight() const noexcept {
        return outputHeight_;
    }
private:
    void InitDerivativeWindow() noexcept {
        qwqdsp_window::BlackmanHarrisThreeTerm::Window(window_, true);
        const float window_gain = qwqdsp_window::Helper::NormalizeGain(window_);
        for (float& v : window_)
            v *= window_gain;

        qwqdsp_window::BlackmanHarrisThreeTerm::DWindow(dwindow_);
        for (float& v : dwindow_)
            v *= window_gain;

        qwqdsp_window::Helper::TWindow(twindow_, window_);

        // tdw_l = dwindow_ · (n - (N-1)/2) / (2π)
        tdwindow_.resize(fftSize_);
        qwqdsp_window::Helper::TWindow(tdwindow_, dwindow_);
        float inv_two_pi = 1.0f / (2.0f * std::numbers::pi_v<float>);
        for (float& v : tdwindow_)
            v *= inv_two_pi;
    }

    int sampleRate_{}, fftSize_{}, hopSize_{}, zeroPad_{}, fftLen_{}, binSize_{};
    int subColumns_{}, outputHeight_{};
    float freqMin_{}, freqMax_{}, logMin_{}, logMax_{}, dbFloor_{}, subColScale_{};

    qwqdsp_spectral::RealFftAdv fft_;
    std::vector<float> fft_in_;
    std::vector<float> window_;
    std::vector<float> dwindow_;
    std::vector<float> twindow_;
    std::vector<std::complex<float>> X_h_, X_dh_, X_th_, X_tdh_;
    std::vector<float> conv_;
    std::vector<float> tdwindow_;
    std::vector<float> col_buf_; // [subCol * outputHeight + y]
    std::vector<float> weight_buf_;
    std::vector<Color> column_;
};
