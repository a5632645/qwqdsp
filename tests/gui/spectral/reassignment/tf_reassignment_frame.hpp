#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <vector>

#include <qwqdsp/spectral/real_fft_adv.hpp>

#include "raylib.h"

#include "log_reassign_grid.hpp"

/**
 * @brief 原始帧+窗 → 频率重分配 + 时间重分配 → magma 颜色输出
 *
 * 显示累加交给 LogReassignGrid: log 频率行 → 线性子格, 子格内求和 / 行内取 max,
 * 不加权, 幅度由窗函数自动标定。
 *
 * X_h = FFT(x·w)
 * X_t = FFT(x[n-1]·w[n])           (频率重分配)
 * X_pf = roll(X_h,1), X_pf[0]=0    (时间重分配 → 群延迟)
 */
template <typename Colormap, bool EnableFreqInterp = true>
struct TfReassignmentFrame {
    void Init(int sampleRate, int fftSize, int hopSize, int zeroPad, int outputHeight, float freqMin, float freqMax,
              float dbFloor) noexcept {
        sampleRate_ = sampleRate;
        fftSize_ = fftSize;
        subColumns_ = fftSize / hopSize;
        fftLen_ = fftSize * zeroPad;
        binSize_ = fftLen_ / 2 + 1;
        outputHeight_ = outputHeight;
        freqMin_ = freqMin;
        freqMax_ = freqMax;

        fft_.Init(fftLen_);

        fft_in_.resize(fftLen_, 0.0f);
        shift_in_.resize(fftSize_);
        X_h_.resize(binSize_);
        X_t_.resize(binSize_);
        X_pf_.resize(binSize_);
        column_.resize(outputHeight_);

        grid_.Init(sampleRate, fftLen_, subColumns_, outputHeight, freqMin, freqMax, dbFloor);
    }

    void Process(std::span<const float> raw_frame, std::span<const float> window,
                 std::span<const float> windowed_frame) noexcept {
        // ── 0. 首次调用: 由用户的窗算出幅度标定 ──
        if (!grid_.CalibrationReady())
            grid_.SetWindow(window, fft_);

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
            inst_freq_norm -= std::floor(inst_freq_norm);
            float inst_freq_hz = inst_freq_norm * sampleRate_;
            if (inst_freq_hz < freqMin_ || inst_freq_hz > freqMax_)
                continue;

            // ── 时间重分配: group_delay ∈ (-0.5, 0.5] ──
            auto cross_f = X_h_[k] * std::conj(X_pf_[k]);
            float arg_f = std::arg(cross_f) / two_pi;
            arg_f -= std::floor(arg_f);
            float group_delay = 0.5f - arg_f;

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

    /// @brief 窗相干增益 (单位纯音的重分配和 / 谱峰), 首次 Process 后有效
    float Calibration() const noexcept {
        return grid_.Calibration();
    }

private:
    int sampleRate_{}, fftSize_{}, subColumns_{}, zeroPad_{}, fftLen_{}, binSize_{}, outputHeight_{};
    float freqMin_{}, freqMax_{};

    qwqdsp_spectral::RealFftAdv fft_;
    std::vector<float> fft_in_;
    std::vector<float> shift_in_;
    std::vector<std::complex<float>> X_h_, X_t_, X_pf_;
    std::vector<Color> column_;
    LogReassignGrid<EnableFreqInterp> grid_;
};
