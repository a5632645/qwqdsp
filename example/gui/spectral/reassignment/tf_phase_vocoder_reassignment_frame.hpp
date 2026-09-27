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

#include "raylib.h"

#include "log_reassign_grid.hpp"

/**
 * @brief Phase Vocoder 瞬时频率 + 时间矩重分配 → 颜色输出
 *
 * EnablePeakFilter=true 时启用 Loris 式重分配最小值旁瓣剔除。
 * 显示累加交给 LogReassignGrid (log 行 → 线性子格 → 子格求和 / 行 max, 不加权)。
 *
 * X_h  = FFT(x·w)                           (标准窗)
 * X_th = FFT(x·(n-center)w[n])              (时间重分配)
 *
 * 频率: Bernsee Phase Vocoder
 *   Δφ = angle(X_h[m]) - lastPhase
 *   Δφ -= k · 2π · hop / N          (减去 bin 预期相位增量)
 *   Δφ_wrapped → [-π, π]
 *   bin_dev = N/hop · Δφ_wrapped / (2π)
 *   f_inst = (k + bin_dev) · sr / N
 *
 * 时间: 时间加权窗法
 *   group_delay = Re(conj(X_h)·X_th) / (|X_h|² · N)
 */
template <typename Colormap, bool EnablePeakFilter = false, bool EnableFreqInterp = true>
struct TfPhaseVocoderReassignmentFrame {
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
        twindow_.resize(fftSize_);
        X_h_.resize(binSize_);
        X_th_.resize(binSize_);
        lastPhase_.resize(binSize_, 0.0f);
        column_.resize(outputHeight_);

        if constexpr (EnablePeakFilter) {
            freq_c_arr_.resize(binSize_);
            sidelobe_mask_.resize(binSize_);
        }

        grid_.Init(sampleRate, fftLen_, subColumns_, outputHeight, freqMin, freqMax, dbFloor);
        InitWindowImpl();
        grid_.SetWindow(window_, fft_);
    }

    void Process(std::span<const float> raw_frame, std::span<const float> /*window*/,
                 std::span<const float> /*windowed_frame*/) noexcept {
        constexpr float kEps = 1e-20f;
        const float bin_hz = static_cast<float>(sampleRate_) / static_cast<float>(fftLen_);
        const float two_pi = std::numbers::pi_v<float> * 2.0f;

        // ── 1. X_h = FFT(x * w) ──
        for (int i = 0; i < fftSize_; ++i)
            fft_in_[i] = raw_frame[i] * window_[i];
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_h_);

        // ── 2. X_th = FFT(x * (n-center)w[n]) ──
        for (int i = 0; i < fftSize_; ++i)
            fft_in_[i] = raw_frame[i] * twindow_[i];
        std::fill(fft_in_.begin() + fftSize_, fft_in_.end(), 0.0f);
        fft_.FFT(fft_in_, X_th_);

        // ── 3a. 第一趟: 计算 bin_dev (仅 EnablePeakFilter 需提前算) ──
        if constexpr (EnablePeakFilter) {
            for (int k = 0; k < binSize_; ++k) {
                freq_c_arr_[k] = 0.0f;
                if (std::abs(X_h_[k]) < 1e-8f)
                    continue;
                float phase = std::arg(X_h_[k]);
                float dphi = phase - lastPhase_[k];
                lastPhase_[k] = phase;
                dphi -= static_cast<float>(k) * expct_;
                float qpd = std::floor((dphi + std::numbers::pi_v<float>) / two_pi);
                dphi -= qpd * two_pi;
                freq_c_arr_[k] = osamp_ * dphi / two_pi;
            }
        }

        // ── 3b. 旁瓣检测 (仅 EnablePeakFilter) ──
        if constexpr (EnablePeakFilter) {
            std::fill(sidelobe_mask_.begin(), sidelobe_mask_.end(), 0);
            for (int j = 1; j < binSize_ - 2; ++j) {
                float f_rs_j = static_cast<float>(j) + freq_c_arr_[j];
                float f_rs_j1 = static_cast<float>(j + 1) + freq_c_arr_[j + 1];
                if (f_rs_j > static_cast<float>(j) && f_rs_j1 < static_cast<float>(j + 1)) {
                    if ((f_rs_j - static_cast<float>(j)) > (static_cast<float>(j + 1) - f_rs_j1))
                        sidelobe_mask_[j] = 1;
                    else
                        sidelobe_mask_[j + 1] = 1;
                }
            }
        }

        // ── 3c. 遍历每个 bin → 时间重分配 + 落到网格 ──
        for (int k = 0; k < binSize_; ++k) {
            if constexpr (EnablePeakFilter) {
                if (sidelobe_mask_[k])
                    continue;
            }

            const float mag_lin = std::abs(X_h_[k]);
            if (mag_lin < 1e-8f)
                continue;

            const float mag_sq = std::max(std::norm(X_h_[k]), kEps);

            // ── Phase Vocoder 瞬时频率 (Bernsee) ──
            float bin_dev;
            if constexpr (EnablePeakFilter) {
                bin_dev = freq_c_arr_[k];
            }
            else {
                float phase = std::arg(X_h_[k]);
                float dphi = phase - lastPhase_[k];
                lastPhase_[k] = phase;
                dphi -= static_cast<float>(k) * expct_;
                float qpd = std::floor((dphi + std::numbers::pi_v<float>) / two_pi);
                dphi -= qpd * two_pi;
                bin_dev = osamp_ * dphi / two_pi;
            }
            float inst_freq_hz = (static_cast<float>(k) + bin_dev) * bin_hz;

            // ── 时间矩重分配: group_delay ∈ [-0.5, 0.5] ──
            auto const& xh = X_h_[k];
            float time_num = xh.real() * X_th_[k].real() + xh.imag() * X_th_[k].imag();
            float time_offset = time_num / mag_sq;
            float group_delay = std::clamp(time_offset / static_cast<float>(fftSize_), -0.5f, 0.5f);

            if (inst_freq_hz < freqMin_ || inst_freq_hz > freqMax_)
                continue;

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
    void InitWindowImpl() noexcept {
        qwqdsp_window::BlackmanHarrisThreeTerm::Window(window_, true);
        const float window_gain = qwqdsp_window::Helper::NormalizeGain(window_);
        for (float& v : window_)
            v *= window_gain;

        qwqdsp_window::Helper::TWindow(twindow_, window_);
    }

    int sampleRate_{}, fftSize_{}, hopSize_{}, zeroPad_{}, fftLen_{}, binSize_{};
    int subColumns_{}, outputHeight_{};
    float freqMin_{}, freqMax_{}, expct_{}, osamp_{};

    qwqdsp_spectral::RealFftAdv fft_;
    std::vector<float> fft_in_;
    std::vector<float> window_;
    std::vector<float> twindow_;
    std::vector<std::complex<float>> X_h_, X_th_;
    std::vector<float> lastPhase_;
    std::vector<float> freq_c_arr_;
    std::vector<uint8_t> sidelobe_mask_;
    std::vector<Color> column_;
    LogReassignGrid<EnableFreqInterp> grid_;
};
