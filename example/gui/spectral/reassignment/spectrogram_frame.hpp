#pragma once

#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

#include <qwqdsp/spectral/real_fft_adv.hpp>

#include "raylib.h"

/**
 * @brief 对一帧已加窗的时域数据做 FFT → 频率映射 → magma 颜色输出
 *
 * 核心处理方法，不管理音频缓冲和环回。
 */
template <typename Colormap>
struct SpectrogramFrame {
    void Init(int sampleRate, int fftSize, int zeroPad, int height, float freqMin, float freqMax,
              float dbFloor) noexcept {
        sample_rate_ = sampleRate;
        fft_size_ = fftSize;
        zero_pad_ = zeroPad;
        fft_len_ = fftSize * zeroPad;
        bin_size_ = fft_len_ / 2 + 1;
        freq_min_ = freqMin;
        freq_max_ = freqMax;
        db_floor_ = dbFloor;
        height_ = height;

        fft_.Init(fft_len_);

        log_min_ = std::log10(freq_min_);
        log_max_ = std::log10(freq_max_);

        fft_in_.resize(fft_len_, 0.0f);
        gain_.resize(bin_size_);
        column_.resize(height_);
        colMax_.resize(height_);
        row_to_bin_.resize(height_);

        // ── 预计算每行 → 最近单 bin ──
        const float bin_spacing = static_cast<float>(sample_rate_) / static_cast<float>(fft_len_);
        const float log_step = (log_max_ - log_min_) / static_cast<float>(height_);

        for (int y = 0; y < height_; ++y) {
            float f_center = std::pow(10.0f, log_max_ - (y + 0.5f) * log_step);
            int bin = static_cast<int>(std::round(f_center / bin_spacing));
            row_to_bin_[y] = std::clamp(bin, 0, bin_size_ - 1);
        }

        // ── 自适应 crossover ──
        const float pixelSlope = std::log(10.0f) * (log_max_ - log_min_) / static_cast<float>(height_ - 1);
        const float crossoverHz = bin_spacing / pixelSlope;
        crossoverY_ = static_cast<int>((1.0f - (std::log10(crossoverHz) - log_min_) / (log_max_ - log_min_))
                                       * static_cast<float>(height_ - 1));
        crossoverY_ = std::clamp(crossoverY_, 0, height_ - 1);
    }

    void Process(std::span<const float> /*raw_frame*/, std::span<const float> /*window*/,
                 std::span<const float> windowed_frame) noexcept {
        // 零填充后 FFT
        std::copy(windowed_frame.begin(), windowed_frame.end(), fft_in_.begin());
        std::fill(fft_in_.begin() + fft_size_, fft_in_.end(), 0.0f);
        fft_.FFTGainPhase(fft_in_, gain_);

        const float bin_spacing = static_cast<float>(sample_rate_) / static_cast<float>(fft_len_);
        const float log_step = (log_max_ - log_min_) / static_cast<float>(height_);

        // ── Pass 1 — 低频 nearest-bin ──
        for (int y = crossoverY_; y < height_; ++y) {
            int bin = row_to_bin_[y];
            colMax_[y] = 20.0f * std::log10(gain_[bin] + 1e-12f);
        }

        // ── Pass 2 — 高频 max-hold ──
        for (int y = 0; y < crossoverY_; ++y) {
            // 像素 y 覆盖 [f_bot, f_top)（y 越小频率越高，故 f_top > f_bot，对应 bin_bot < bin_top）
            float f_top = std::pow(10.0f, log_max_ - static_cast<float>(y) * log_step);
            float f_bot = std::pow(10.0f, log_max_ - static_cast<float>(y + 1) * log_step);
            int bin_top = static_cast<int>(std::round(f_top / bin_spacing));
            int bin_bot = static_cast<int>(std::round(f_bot / bin_spacing));
            bin_top = std::clamp(bin_top, 0, bin_size_ - 1);
            bin_bot = std::clamp(bin_bot, 0, bin_top);

            float maxDb = -1e12f;
            for (int b = bin_bot; b <= bin_top; ++b) {
                float dB = 20.0f * std::log10(gain_[b] + 1e-12f);
                if (dB > maxDb)
                    maxDb = dB;
            }
            colMax_[y] = maxDb;
        }

        // ── colMax → 颜色 ──
        for (int y = 0; y < height_; ++y) {
            float dB = std::clamp(colMax_[y], db_floor_, 0.0f);
            int idx = static_cast<int>((dB - db_floor_) / (-db_floor_) * 255.0f);
            idx = std::clamp(idx, 0, 255);
            column_[y] = Colormap::kTable[idx];
        }
    }

    std::span<const Color> GetColumn() const noexcept {
        return {column_.data(), static_cast<size_t>(height_)};
    }

    int ColumnHeight() const noexcept {
        return height_;
    }
private:
    qwqdsp_spectral::RealFftAdv fft_;
    std::vector<float> fft_in_;
    std::vector<float> gain_;
    std::vector<float> colMax_;
    std::vector<Color> column_;
    std::vector<int> row_to_bin_;
    float freq_min_{}, freq_max_{}, log_min_{}, log_max_{}, db_floor_{};
    int bin_size_{}, height_{}, sample_rate_{}, fft_size_{}, zero_pad_{}, fft_len_{}, crossoverY_{};
};
