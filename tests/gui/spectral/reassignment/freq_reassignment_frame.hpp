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
 * @brief 对原始帧+窗做 FFT → 频率重分配 → magma 颜色输出
 *
 * 接收 (raw_frame, window, windowed_frame)。
 * X_h = FFT(x·w)，X_t = FFT(x[n-1]·w[n])。
 * 互谱相位差 → 瞬时频率 → 能量按 LogReassignGrid 的
 * "log 行 → 线性子格 → 子格求和 / 行 max (不加权)" 分布。
 *
 * 无时间重分配 (subColumns = 1), 因此不接收 hopSize。
 */
template <typename Colormap>
struct FreqReassignmentFrame {
    void Init(int sampleRate, int fftSize, int zeroPad, int outputHeight, float freqMin, float freqMax,
              float dbFloor) noexcept {
        sampleRate_ = sampleRate;
        fftSize_ = fftSize;
        zeroPad_ = zeroPad;
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
        column_.resize(outputHeight_);

        grid_.Init(sampleRate, fftLen_, 1, outputHeight, freqMin, freqMax, dbFloor);
    }

    void Process(std::span<const float> raw_frame, std::span<const float> window,
                 std::span<const float> windowed_frame) noexcept {
        if (!grid_.CalibrationReady())
            grid_.SetWindow(window, fft_);

        // ── 1. 零填充后 FFT: X_h = FFT(x * w, n=fftLen) ──
        std::copy(windowed_frame.begin(), windowed_frame.end(), fft_in_.begin());
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_h_);

        // ── 2. 时移 1 样本再乘窗 → FFT: X_t = FFT(x[n-1] * w[n], n=fftLen) ──
        //     注意: 不能移动已加窗信号 (x*w)[n-1], 必须移动原始帧 x[n-1] 再乘 w[n]
        shift_in_[0] = 0.0f;
        for (int i = 1; i < fftSize_; ++i)
            shift_in_[i] = raw_frame[i - 1] * window[i];
        std::copy(shift_in_.begin(), shift_in_.end(), fft_in_.begin());
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_t_);

        // ── 3. 遍历每个 bin: 互谱 → 瞬时频率 → 落到网格 ──
        const float two_pi = 2.0f * std::numbers::pi_v<float>;

        for (int k = 0; k < binSize_; ++k) {
            float mag_lin = std::abs(X_h_[k]);
            if (mag_lin < 1e-8f)
                continue;

            auto cross = X_h_[k] * std::conj(X_t_[k]);
            float inst_freq_norm = std::arg(cross) / two_pi;
            inst_freq_norm -= std::floor(inst_freq_norm);
            float inst_freq_hz = inst_freq_norm * sampleRate_;
            if (inst_freq_hz < freqMin_ || inst_freq_hz > freqMax_)
                continue;

            grid_.Add(inst_freq_hz, 0.0f, mag_lin);
        }

        grid_.Emit(Colormap::kTable, column_);
    }

    std::span<const Color> GetColumn() const noexcept {
        return {column_.data(), static_cast<size_t>(outputHeight_)};
    }

    int ColumnHeight() const noexcept {
        return outputHeight_;
    }

    /// @brief 窗相干增益 (自动标定), 首次 Process 后有效
    float Calibration() const noexcept {
        return grid_.Calibration();
    }

private:
    int sampleRate_{}, fftSize_{}, zeroPad_{}, fftLen_{}, binSize_{}, outputHeight_{};
    float freqMin_{}, freqMax_{};

    qwqdsp_spectral::RealFftAdv fft_;
    std::vector<float> fft_in_;
    std::vector<float> shift_in_;
    std::vector<std::complex<float>> X_h_, X_t_;
    std::vector<Color> column_;
    LogReassignGrid grid_;
};
