#pragma once
#include "../../spectral/real_fft_adv.hpp"
#include "../pitch.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

namespace qwqdsp_pitch {

// ------------------------------------------------------------
// PyinCore
// ------------------------------------------------------------
// 概率 YIN (pYIN) 核心算法。
// 对每帧音频计算差函数 → CMNDF → 多阈值概率分布 → 返回多个基频候选。
// 参考: Mauch & Dixon 的 pYIN 算法。

struct PyinCandidate {
    float pitch_hz;    // 基频 (Hz)
    float probability; // 概率
};

class PyinCore {
public:
    void Init(float sample_rate, int block_size) {
        fs_ = sample_rate;
        block_size_ = block_size;
        yin_buffer_size_ = block_size / 2;
        yin_buffer_.resize(yin_buffer_size_);
        powers_.resize(yin_buffer_size_);
        temp_.resize(block_size);
        peak_prob_.resize(yin_buffer_size_);
        fft_.Init(block_size);
        fft1_.resize(fft_.NumBins());
        fft2_.resize(fft_.NumBins());
        SetMinPitch(min_pitch_);
        SetMaxPitch(max_pitch_);
    }

    void SetMinPitch(float hz) noexcept {
        min_pitch_ = hz;
        max_tau_ = std::min(static_cast<int>(std::round(fs_ / hz)), yin_buffer_size_ - 1);
    }

    void SetMaxPitch(float hz) noexcept {
        max_pitch_ = hz;
        min_tau_ = std::max(static_cast<int>(std::round(fs_ / hz)), 2);
    }

    // --------------------------------------------------------
    // Process
    // --------------------------------------------------------
    // 处理一帧音频，返回候选基频列表（已按概率降序排列）。
    // 候选不超过 max_candidates 个，概率低于 min_prob 的被过滤。
    std::vector<PyinCandidate> Process(std::span<const float> block, int max_candidates = 5, float min_prob = 0.01f) {
        FastDifference(block);
        CcumulativeMeanNormalizedDifference();
        YinProb(static_cast<size_t>(min_tau_), static_cast<size_t>(max_tau_));

        std::vector<PyinCandidate> candidates;
        for (size_t i = 0; i < yin_buffer_size_; ++i) {
            if (peak_prob_[i] > min_prob && i > 0) {
                float tau = ParabolicInterpolationFromBuffer(static_cast<int>(i));
                float pitch = fs_ / tau;
                if (pitch >= 20.0f && pitch <= 2000.0f) {
                    candidates.push_back({pitch, static_cast<float>(peak_prob_[i])});
                }
            }
        }

        // 按概率降序排列
        std::sort(candidates.begin(), candidates.end(),
                  [](const PyinCandidate& a, const PyinCandidate& b) { return a.probability > b.probability; });

        if (candidates.size() > static_cast<size_t>(max_candidates)) {
            candidates.resize(max_candidates);
        }

        return candidates;
    }

    // 直接访问归一化后的差函数数组（用于调试或自定义处理）
    std::span<const float> YinBuffer() const noexcept {
        return yin_buffer_;
    }

    std::span<const double> PeakProbability() const noexcept {
        return peak_prob_;
    }
private:
    // --------------------------------------------------------
    // fastDifference — FFT 加速差函数
    // --------------------------------------------------------
    void FastDifference(std::span<const float> block) noexcept {
        int const max_tal = yin_buffer_size_;

        // ----- power terms -----
        std::fill(powers_.begin(), powers_.end(), 0.0f);
        for (int i = 0; i < max_tal; ++i) {
            powers_[0] += block[i] * block[i];
        }
        for (int tau = 1; tau < max_tal; ++tau) {
            powers_[tau] =
                powers_[tau - 1] - block[tau - 1] * block[tau - 1] + block[tau + max_tal] * block[tau + max_tal];
        }

        // ----- FFT of full block -----
        fft_.FFT(block, fft1_);

        // ----- kernel: reversed x[0:N/2], zero-padded to N -----
        for (int i = 0; i < max_tal; ++i) {
            temp_[i] = block[max_tal - 1 - i];
            temp_[i + max_tal] = 0.0f;
        }
        fft_.FFT(temp_, fft2_);

        // ----- complex multiply (convolution in freq domain) -----
        size_t const num_bins = fft_.NumBins();
        for (size_t i = 0; i < num_bins; ++i) {
            float re = fft1_[i].real() * fft2_[i].real() - fft1_[i].imag() * fft2_[i].imag();
            float im = fft1_[i].imag() * fft2_[i].real() + fft1_[i].real() * fft2_[i].imag();
            fft1_[i] = {re, im};
        }
        fft_.IFFT(temp_, fft1_);

        // ----- difference function d(tau) -----
        for (int i = 0; i < max_tal; ++i) {
            yin_buffer_[i] = powers_[0] + powers_[i] - 2.0f * temp_[max_tal - 1 + i];
            if (yin_buffer_[i] < 0.0f) {
                yin_buffer_[i] = 0.0f; // 数值误差保护
            }
        }
    }

    // --------------------------------------------------------
    // cumulativeMeanNormalizedDifference — CMNDF
    // --------------------------------------------------------
    void CcumulativeMeanNormalizedDifference() noexcept {
        yin_buffer_[0] = 1.0f;
        float running_sum = 0.0f;
        for (int tau = 1; tau < yin_buffer_size_; ++tau) {
            running_sum += yin_buffer_[tau];
            if (running_sum != 0.0f) {
                yin_buffer_[tau] *= static_cast<float>(tau) / running_sum;
            }
            else {
                yin_buffer_[tau] = 1.0f;
            }
        }
    }

    // --------------------------------------------------------
    // yinProb — 多阈值概率估计
    // --------------------------------------------------------
    // 使用 100 个 threshold 对差函数局部最小值累积概率权重。
    void YinProb(size_t min_tau, size_t max_tau) noexcept {
        std::fill(peak_prob_.begin(), peak_prob_.end(), 0.0);

        size_t min_ind = 0;
        float min_val = std::numeric_limits<float>::max();
        float sum_prob = 0.0f;

        size_t tau = min_tau;
        while (tau + 1 < max_tau) {
            if (yin_buffer_[tau] < kThresholds[99] && yin_buffer_[tau + 1] < yin_buffer_[tau]) {
                // walk downhill to local minimum
                while (tau + 1 < max_tau && yin_buffer_[tau + 1] < yin_buffer_[tau]) {
                    ++tau;
                }
                // tau is now at a local minimum
                if (yin_buffer_[tau] < min_val && tau > 2) {
                    min_val = yin_buffer_[tau];
                    min_ind = tau;
                }
                // accumulate probability from all thresholds below this minimum's value
                int thresh_idx = 99;
                while (thresh_idx >= 0 && kThresholds[thresh_idx] > yin_buffer_[tau]) {
                    peak_prob_[tau] += kBetaDist2[thresh_idx];
                    --thresh_idx;
                }
                sum_prob += static_cast<float>(peak_prob_[tau]);
                ++tau;
            }
            else {
                ++tau;
            }
        }

        // 归一化概率
        if (sum_prob > 0.0f && min_ind > 0 && min_ind < max_tau) {
            double peak_min = peak_prob_[min_ind];
            for (size_t i = min_tau; i < max_tau; ++i) {
                peak_prob_[i] = peak_prob_[i] / sum_prob * peak_min;
            }
            // 将剩余概率分配给全局最小值（保证总和为 1）
            double non_peak_prob = 1.0;
            for (size_t i = min_tau; i < max_tau; ++i) {
                non_peak_prob -= peak_prob_[i];
            }
            if (non_peak_prob > 0.0) {
                peak_prob_[min_ind] += non_peak_prob * kMinWeight;
            }
        }
    }

    // --------------------------------------------------------
    // parabolicInterpolation — 抛物线插值
    // --------------------------------------------------------
    float ParabolicInterpolationFromBuffer(int tau_estimate) const noexcept {
        if (tau_estimate <= 0 || tau_estimate >= static_cast<int>(yin_buffer_size_) - 1) {
            return static_cast<float>(tau_estimate);
        }
        float s0 = yin_buffer_[tau_estimate - 1];
        float s1 = yin_buffer_[tau_estimate];
        float s2 = yin_buffer_[tau_estimate + 1];
        double adjustment = (s2 - s0) / (2.0 * (2.0 * s1 - s2 - s0));
        if (std::abs(adjustment) > 1.0) {
            adjustment = 0.0;
        }
        return static_cast<float>(tau_estimate) + static_cast<float>(adjustment);
    }

    // --------------------------------------------------------
    // 静态数据
    // --------------------------------------------------------

    // 100 个等间距 threshold: 0.01, 0.02, ..., 1.00
    static constexpr float kThresholds[100] = {
        0.01f, 0.02f, 0.03f, 0.04f, 0.05f, 0.06f, 0.07f, 0.08f, 0.09f, 0.10f, 0.11f, 0.12f, 0.13f, 0.14f, 0.15f,
        0.16f, 0.17f, 0.18f, 0.19f, 0.20f, 0.21f, 0.22f, 0.23f, 0.24f, 0.25f, 0.26f, 0.27f, 0.28f, 0.29f, 0.30f,
        0.31f, 0.32f, 0.33f, 0.34f, 0.35f, 0.36f, 0.37f, 0.38f, 0.39f, 0.40f, 0.41f, 0.42f, 0.43f, 0.44f, 0.45f,
        0.46f, 0.47f, 0.48f, 0.49f, 0.50f, 0.51f, 0.52f, 0.53f, 0.54f, 0.55f, 0.56f, 0.57f, 0.58f, 0.59f, 0.60f,
        0.61f, 0.62f, 0.63f, 0.64f, 0.65f, 0.66f, 0.67f, 0.68f, 0.69f, 0.70f, 0.71f, 0.72f, 0.73f, 0.74f, 0.75f,
        0.76f, 0.77f, 0.78f, 0.79f, 0.80f, 0.81f, 0.82f, 0.83f, 0.84f, 0.85f, 0.86f, 0.87f, 0.88f, 0.89f, 0.90f,
        0.91f, 0.92f, 0.93f, 0.94f, 0.95f, 0.96f, 0.97f, 0.98f, 0.99f, 1.00f};

    static constexpr float kBetaDist2[100] = {
        0.012614f, 0.022715f, 0.030646f, 0.036712f, 0.041184f, 0.044301f, 0.046277f, 0.047298f, 0.047528f, 0.047110f,
        0.046171f, 0.044817f, 0.043144f, 0.041231f, 0.039147f, 0.036950f, 0.034690f, 0.032406f, 0.030133f, 0.027898f,
        0.025722f, 0.023624f, 0.021614f, 0.019704f, 0.017900f, 0.016205f, 0.014621f, 0.013148f, 0.011785f, 0.010530f,
        0.009377f, 0.008324f, 0.007366f, 0.006497f, 0.005712f, 0.005005f, 0.004372f, 0.003806f, 0.003302f, 0.002855f,
        0.002460f, 0.002112f, 0.001806f, 0.001539f, 0.001307f, 0.001105f, 0.000931f, 0.000781f, 0.000652f, 0.000542f,
        0.000449f, 0.000370f, 0.000303f, 0.000247f, 0.000201f, 0.000162f, 0.000130f, 0.000104f, 0.000082f, 0.000065f,
        0.000051f, 0.000039f, 0.000030f, 0.000023f, 0.000018f, 0.000013f, 0.000010f, 0.000007f, 0.000005f, 0.000004f,
        0.000003f, 0.000002f, 0.000001f, 0.000001f, 0.000001f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f,
        0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f,
        0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f};

    static constexpr float kMinWeight = 0.01f;

    // --------------------------------------------------------
    // 成员变量
    // --------------------------------------------------------
    float fs_{};
    int block_size_{};
    int yin_buffer_size_{};
    float min_pitch_{60.0f};
    float max_pitch_{900.0f};
    int min_tau_{2};
    int max_tau_{};

    std::vector<float> yin_buffer_;
    std::vector<float> powers_;
    std::vector<float> temp_;
    std::vector<double> peak_prob_;

    qwqdsp_spectral::RealFftAdv fft_;
    std::vector<std::complex<float>> fft1_;
    std::vector<std::complex<float>> fft2_;
};

} // namespace qwqdsp_pitch
