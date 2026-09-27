#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <vector>

#include <qwqdsp/spectral/real_fft_adv.hpp>

#include "raylib.h"

/**
 * @brief NC (Negative Correlation) 方法谱帧
 *
 * 对原始帧做矩形窗 FFT，计算相邻 bin 之间的负相关系数：
 *   ncSum = -(Re(X[k])·Re(X[k+1]) + Im(X[k])·Im(X[k+1]))
 *
 * ncSum < 0 时输出 dbFloor，否则输出 20·log10(sqrt(ncSum))。
 * 每行按频率映射到最近 NC bin 对。
 *
 * @ref https://arxiv.org/html/2410.07982v3#bib.bib1
 *
 */
template <typename Colormap>
struct NcReassignmentFrame {
    void Init(int sampleRate, int fftSize, int zeroPad, int outputHeight, float freqMin, float freqMax,
              float dbFloor) noexcept {
        sampleRate_ = sampleRate;
        fftSize_ = fftSize;
        zeroPad_ = zeroPad;
        fftLen_ = fftSize * zeroPad;
        binSize_ = fftLen_ / 2 + 1;
        ncSize_ = binSize_ - zeroPad_; // NC 值数量 = 零填充后有效 bin 对数
        outputHeight_ = outputHeight;
        freqMin_ = freqMin;
        freqMax_ = freqMax;
        dbFloor_ = dbFloor;
        logMin_ = std::log10(freqMin);
        logMax_ = std::log10(freqMax);

        fft_.Init(fftLen_);

        fft_in_.resize(fftLen_, 0.0f);
        X_.resize(binSize_);
        nc_.resize(ncSize_);
        colMax_.resize(outputHeight_);
        column_.resize(outputHeight_);
        rowToBin_.resize(outputHeight_);

        // ── 预计算每行 → 最近 NC bin ──
        const float bin_hz = static_cast<float>(sampleRate_) / static_cast<float>(fftLen_);
        const float logStep = (logMax_ - logMin_) / static_cast<float>(outputHeight_);

        for (int y = 0; y < outputHeight_; ++y) {
            float f_center = std::pow(10.0f, logMax_ - (y + 0.5f) * logStep);
            float centerRatio = f_center / bin_hz - 0.5f * static_cast<float>(zeroPad_);
            rowToBin_[y] = static_cast<int>(std::round(centerRatio));
            rowToBin_[y] = std::clamp(rowToBin_[y], 0, ncSize_ - 1);
        }

        // ── 自适应 crossover: 像素宽度 == NC bin 宽度的频率 ──
        const float pixelSlope = std::log(10.0f) * (logMax_ - logMin_) / static_cast<float>(outputHeight_ - 1);
        const float binWidth = static_cast<float>(zeroPad_) * bin_hz;
        const float crossoverHz = binWidth / pixelSlope;
        crossoverBin_ = static_cast<int>(std::round(crossoverHz / bin_hz - 0.5f * static_cast<float>(zeroPad_)));
        crossoverBin_ = std::clamp(crossoverBin_, 0, ncSize_ - 1);
        float crNorm = (std::log10(crossoverHz) - logMin_) / (logMax_ - logMin_);
        crossoverY_ = static_cast<int>((1.0f - crNorm) * static_cast<float>(outputHeight_ - 1));
        crossoverY_ = std::clamp(crossoverY_, 0, outputHeight_ - 1);
    }

    void Process(std::span<const float> raw_frame, std::span<const float> /*window*/,
                 std::span<const float> /*windowed_frame*/) noexcept {
        // ── 1. 矩形窗 FFT（原始帧直接零填充）──
        std::copy(raw_frame.begin(), raw_frame.end(), fft_in_.begin());
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_);

        // ── 2. 归一化增益 2/fftSize ──
        const float norm = 2.0f / static_cast<float>(fftSize_);
        for (int k = 0; k < binSize_; ++k) {
            X_[k] *= norm;
        }

        // ── 3. 预计算所有 NC bin 的 dB 值 ──
        constexpr float kEps = 1e-18f;
        for (int k = 0; k < ncSize_; ++k) {
            float ncSum = -(X_[k].real() * X_[k + zeroPad_].real() + X_[k].imag() * X_[k + zeroPad_].imag());
            if (ncSum < 0.0f) {
                nc_[k] = dbFloor_;
            }
            else {
                float gain = std::sqrt(ncSum + kEps);
                nc_[k] = 20.0f * std::log10(gain);
            }
            nc_[k] = std::clamp(nc_[k], dbFloor_, 0.0f);
        }

        // ── 4. Pass 1 — 低频像素驱动: y=crossoverY_..bottom 查预计算表 ──
        const float bin_hz = static_cast<float>(sampleRate_) / static_cast<float>(fftLen_);

        std::fill(colMax_.begin(), colMax_.end(), dbFloor_);

        for (int y = crossoverY_; y < outputHeight_; ++y) {
            colMax_[y] = nc_[rowToBin_[y]];
        }

        // ── 5. Pass 2 — 高频 bin 驱动: k=crossoverBin_..ncSize_-1 max-hold ──
        for (int k = crossoverBin_; k < ncSize_; ++k) {
            float f_low = static_cast<float>(k) * bin_hz;
            float f_high = static_cast<float>(k + zeroPad_) * bin_hz;
            if (f_high < freqMin_ || f_low > freqMax_)
                continue;
            f_low = std::max(f_low, freqMin_);
            f_high = std::min(f_high, freqMax_);

            float n_low = (std::log10(f_low) - logMin_) / (logMax_ - logMin_);
            float n_high = (std::log10(f_high) - logMin_) / (logMax_ - logMin_);
            int y_start = static_cast<int>((1.0f - n_high) * static_cast<float>(outputHeight_ - 1));
            int y_end = static_cast<int>((1.0f - n_low) * static_cast<float>(outputHeight_ - 1));
            y_start = std::clamp(y_start, 0, outputHeight_ - 1);
            y_end = std::clamp(y_end, 0, outputHeight_ - 1);

            float ncDb = nc_[k];
            for (int y = y_start; y <= y_end; ++y) {
                if (ncDb > colMax_[y])
                    colMax_[y] = ncDb;
            }
        }

        // ── 4. colMax → 颜色 ──
        for (int y = 0; y < outputHeight_; ++y) {
            int idx = static_cast<int>((colMax_[y] - dbFloor_) / (-dbFloor_) * 255.0f);
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
    int sampleRate_{}, fftSize_{}, zeroPad_{}, fftLen_{}, binSize_{}, ncSize_{}, crossoverBin_{}, crossoverY_{};
    int outputHeight_{};
    float freqMin_{}, freqMax_{}, logMin_{}, logMax_{}, dbFloor_{};

    qwqdsp_spectral::RealFftAdv fft_;
    std::vector<float> fft_in_;
    std::vector<std::complex<float>> X_;
    std::vector<float> nc_;
    std::vector<float> colMax_;
    std::vector<int> rowToBin_;
    std::vector<Color> column_;
};
