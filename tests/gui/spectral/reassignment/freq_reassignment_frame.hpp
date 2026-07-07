#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <vector>

#include <qwqdsp/spectral/real_fft.hpp>

#include "raylib.h"

/**
 * @brief 对原始帧+窗做 FFT → 频率重分配 → magma 颜色输出
 *
 * 接收 (raw_frame, window, windowed_frame)。
 * X_h = FFT(x·w)，X_t = FFT(x[n-1]·w[n])。
 * 互谱相位差 → 瞬时频率 → 能量线性分布到输出网格。
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
        dbFloor_ = dbFloor;
        logMin_ = std::log10(freqMin);
        logMax_ = std::log10(freqMax);

        fft_.Init(fftLen_);

        fft_in_.resize(fftLen_, 0.0f);
        shift_in_.resize(fftSize_);
        X_h_.resize(binSize_);
        X_t_.resize(binSize_);
        col_lin_.resize(outputHeight_, 0.0f);
        column_.resize(outputHeight_);
    }

    void Process(std::span<const float> raw_frame, std::span<const float> window,
                 std::span<const float> windowed_frame) noexcept {
        // ── 1. 零填充后 FFT: X_h = FFT(x * w, n=fftLen) ──
        std::copy(windowed_frame.begin(), windowed_frame.end(), fft_in_.begin());
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_h_);

        // ── 2. 时移 1 样本再乘窗 → 零填充后 FFT: X_t = FFT(x[n-1] * w[n], n=fftLen) ──
        //     注意: 不能移动已加窗信号 (x*w)[n-1], 必须移动原始帧 x[n-1] 再乘 w[n]
        shift_in_[0] = 0.0f;
        for (int i = 1; i < fftSize_; ++i)
            shift_in_[i] = raw_frame[i - 1] * window[i];
        std::copy(shift_in_.begin(), shift_in_.end(), fft_in_.begin());
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_t_);

        // ── 3. 遍历每个 bin: 互谱 → 瞬时频率 → 能量线性分布 ──
        std::fill(col_lin_.begin(), col_lin_.end(), 0.0f);
        const float two_pi = 2.0f * std::numbers::pi_v<float>;

        for (int k = 0; k < binSize_; ++k) {
            // cross[k] = X_h[k] · conj(X_t[k])
            auto cross = X_h_[k] * std::conj(X_t_[k]);

            // mag_lin = |X_h|  (线性幅度归一化)
            float mag_lin = std::abs(X_h_[k]);
            if (mag_lin < 1e-8f)
                continue;

            // inst_freq_norm = mod(angle(cross) / 2π, 1)  → [0, 1)
            float inst_freq_norm = std::arg(cross) / two_pi;
            while (inst_freq_norm < 0.0f)
                inst_freq_norm += 1.0f;
            inst_freq_norm = std::fmod(inst_freq_norm, 1.0f);

            // inst_freq_hz = inst_freq_norm · sr
            float inst_freq_hz = inst_freq_norm * sampleRate_;
            if (inst_freq_hz < freqMin_ || inst_freq_hz > freqMax_)
                continue;

            // 瞬时频率 → Y 像素浮点位置 (对数坐标)
            float logF = std::log10(inst_freq_hz);
            float norm = (logF - logMin_) / (logMax_ - logMin_);
            float y_pos = static_cast<float>(outputHeight_ - 1) * (1.0f - norm);

            // 线性分配给相邻两行 (先 clamp y_pos 再算 y_idx/y_frac)
            y_pos = std::clamp(y_pos, 0.0f, static_cast<float>(outputHeight_ - 1));
            int y_idx = static_cast<int>(std::floor(y_pos));
            if (y_idx >= outputHeight_ - 1) {
                col_lin_[outputHeight_ - 1] += mag_lin;
                continue;
            }
            float y_frac = y_pos - static_cast<float>(y_idx);

            col_lin_[y_idx] += mag_lin * (1.0f - y_frac);
            col_lin_[y_idx + 1] += mag_lin * y_frac;
        }

        // ── 4. 线性幅度 → dB → magma 颜色 ──
        constexpr float kEps = 1e-12f;
        for (int y = 0; y < outputHeight_; ++y) {
            float dB = 20.0f * std::log10(col_lin_[y] + kEps);
            dB = std::clamp(dB, dbFloor_, 0.0f);
            int idx = static_cast<int>((dB - dbFloor_) / (-dbFloor_) * 255.0f);
            idx = std::clamp(idx, 0, 255);
            column_[y] = Colormap::kTable[idx];
        }
    }

    std::span<const Color> GetColumn() const noexcept {
        return {column_.data(), static_cast<size_t>(outputHeight_)};
    }

    int ColumnHeight() const noexcept {
        return outputHeight_;
    }
private:
    int sampleRate_{}, fftSize_{}, zeroPad_{}, fftLen_{}, binSize_{};
    int outputHeight_{};
    float freqMin_{}, freqMax_{}, logMin_{}, logMax_{}, dbFloor_{};

    qwqdsp_spectral::RealFFT fft_;
    std::vector<float> fft_in_;
    std::vector<float> shift_in_;
    std::vector<std::complex<float>> X_h_, X_t_;
    std::vector<float> col_lin_;
    std::vector<Color> column_;
};
