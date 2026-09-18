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
#include <qwqdsp/window/helper.hpp>
#include <qwqdsp/window/window.hpp>

#include "raylib.h"

#include "log_reassign_grid.hpp"

/**
 * @brief Phase Vocoder 瞬时频率 + 收敛加权重分配 → magma 颜色输出
 *
 * 频率: Bernsee Phase Vocoder (相邻帧相位差)
 * 时间: 时间加权窗 x·(n-center)w
 * 收敛: Loris 混合偏导 (conv ∈ [0,1], 0=正弦, 1=脉冲);
 *       conv ∈ (0.3, 0.7) 的 bin 被丢弃 (过渡态不参与显示)
 * 显示: LogReassignGrid (log 行 → 线性子格 → 子格求和 / 行 max, 不加权)
 *
 * X_h   = FFT(x·w)                 (标准窗)
 * X_dh  = FFT(x·dw)                (导数窗, 仅用于 conv)
 * X_th  = FFT(x·(n-center)w[n])    (时间重分配)
 * X_tdh = FFT(x·(n-center)dw[n])   (仅用于 conv)
 */
template <typename Colormap>
struct TfPhaseVocoderReassignmentFrameConv {
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
        expct_ = std::numbers::pi_v<float> * 2.0f * static_cast<float>(hopSize) / static_cast<float>(fftLen_);
        osamp_ = static_cast<float>(fftLen_) / static_cast<float>(hopSize);

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
        lastPhase_.resize(binSize_, 0.0f);
        column_.resize(outputHeight_);

        grid_.Init(sampleRate, fftLen_, subColumns_, outputHeight, freqMin, freqMax, dbFloor);
        InitDerivativeWindow();
        grid_.SetWindow(window_, fft_);
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

        // ── 5. 遍历每个 bin: PV 频率 + 时间重分配 + 落到网格 ──
        constexpr float kEps = 1e-20f;
        const float bin_hz = static_cast<float>(sampleRate_) / static_cast<float>(fftLen_);
        const float two_pi = std::numbers::pi_v<float> * 2.0f;

        for (int k = 0; k < binSize_; ++k) {
            const float mag_lin = std::abs(X_h_[k]);
            if (mag_lin < 1e-8f)
                continue;

            if (conv_[k] > 0.3f && conv_[k] < 0.7f)
                continue;

            const float mag_sq = std::max(std::norm(X_h_[k]), kEps);

            // ── Phase Vocoder 瞬时频率 (Bernsee) ──
            auto const& xh = X_h_[k];
            float phase = std::arg(xh);
            float dphi = phase - lastPhase_[k];
            lastPhase_[k] = phase;
            dphi -= static_cast<float>(k) * expct_;
            float qpd = std::floor((dphi + std::numbers::pi_v<float>) / two_pi);
            dphi -= qpd * two_pi;
            float bin_dev = osamp_ * dphi / two_pi;
            float inst_freq_hz = (static_cast<float>(k) + bin_dev) * bin_hz;
            if (inst_freq_hz < freqMin_ || inst_freq_hz > freqMax_)
                continue;

            // ── 时间矩重分配: group_delay ∈ [-0.5, 0.5] ──
            float time_num = xh.real() * X_th_[k].real() + xh.imag() * X_th_[k].imag();
            float group_delay = std::clamp(time_num / (mag_sq * static_cast<float>(fftSize_)), -0.5f, 0.5f);

            grid_.Add(inst_freq_hz, group_delay, mag_lin);
        }

        grid_.Emit(Colormap::kTable, column_);
    }

    std::span<const Color> GetColumn() const noexcept {
        return {column_.data(), static_cast<size_t>(outputHeight_)};
    }

    int ColumnHeight() const noexcept {
        return outputHeight_;
    }

    /// @brief 窗相干增益 (自动标定)
    float Calibration() const noexcept {
        return grid_.Calibration();
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
    float freqMin_{}, freqMax_{}, expct_{}, osamp_{};

    qwqdsp_spectral::RealFftAdv fft_;
    std::vector<float> fft_in_;
    std::vector<float> window_;
    std::vector<float> dwindow_;
    std::vector<float> twindow_;
    std::vector<std::complex<float>> X_h_, X_dh_, X_th_, X_tdh_;
    std::vector<float> conv_;
    std::vector<float> lastPhase_;
    std::vector<float> tdwindow_;
    std::vector<Color> column_;
    LogReassignGrid grid_;
};
