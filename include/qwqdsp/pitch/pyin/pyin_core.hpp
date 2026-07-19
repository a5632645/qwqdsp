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
        fastDifference(block);
        cumulativeMeanNormalizedDifference();
        yinProb(static_cast<size_t>(min_tau_), static_cast<size_t>(max_tau_));

        std::vector<PyinCandidate> candidates;
        for (size_t i = 0; i < yin_buffer_size_; ++i) {
            if (peak_prob_[i] > min_prob && i > 0) {
                float tau = parabolicInterpolationFromBuffer(static_cast<int>(i));
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

    void SetThresholdDistribution(int dist) noexcept {
        thresh_distr_ = std::clamp(dist, 0, 7);
    }
private:
    // --------------------------------------------------------
    // fastDifference — FFT 加速差函数
    // --------------------------------------------------------
    // 与 Java FastYin / pYIN YinUtil::fastDifference 算法一致。
    void fastDifference(std::span<const float> block) noexcept {
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
    void cumulativeMeanNormalizedDifference() noexcept {
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
    void yinProb(size_t min_tau, size_t max_tau) noexcept {
        std::fill(peak_prob_.begin(), peak_prob_.end(), 0.0);

        float const* dist = kThresholdDistributions[thresh_distr_];

        if (max_tau > static_cast<size_t>(yin_buffer_size_)) {
            max_tau = static_cast<size_t>(yin_buffer_size_);
        }
        if (min_tau < 2)
            min_tau = 2;
        if (min_tau >= max_tau)
            return;

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
                    peak_prob_[tau] += dist[thresh_idx];
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
    float parabolicInterpolationFromBuffer(int tau_estimate) const noexcept {
        if (tau_estimate <= 0 || tau_estimate >= static_cast<int>(yin_buffer_size_) - 1) {
            return static_cast<float>(tau_estimate);
        }
        float s0 = yin_buffer_[tau_estimate - 1];
        float s1 = yin_buffer_[tau_estimate];
        float s2 = yin_buffer_[tau_estimate + 1];
        double adjustment = (s2 - s0) / (2.0 * (2.0 * s1 - s2 - s0));
        if (std::abs(adjustment) > 1.0)
            adjustment = 0.0;
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

    // 8 种 threshold 分布 (100 bins each)
    // 0: Uniform    1: Beta(μ=0.10)  2: Beta(μ=0.15)  3: Beta(μ=0.20)
    // 4: Beta(μ=0.30)  5: Single@0.10  6: Single@0.15  7: Single@0.20
    static constexpr float kUniformDist[100] = {
        0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f,
        0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f,
        0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f,
        0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f,
        0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f,
        0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f,
        0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f,
        0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f,
        0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f,
        0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f,
        0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f, 0.0100000f,
        0.0100000f};
    static constexpr float kBetaDist1[100] = {
        0.028911f, 0.048656f, 0.061306f, 0.068539f, 0.071703f, 0.071877f, 0.069915f, 0.066489f, 0.062117f, 0.057199f,
        0.052034f, 0.046844f, 0.041786f, 0.036971f, 0.032470f, 0.028323f, 0.024549f, 0.021153f, 0.018124f, 0.015446f,
        0.013096f, 0.011048f, 0.009275f, 0.007750f, 0.006445f, 0.005336f, 0.004397f, 0.003606f, 0.002945f, 0.002394f,
        0.001937f, 0.001560f, 0.001250f, 0.000998f, 0.000792f, 0.000626f, 0.000492f, 0.000385f, 0.000300f, 0.000232f,
        0.000179f, 0.000137f, 0.000104f, 0.000079f, 0.000060f, 0.000045f, 0.000033f, 0.000024f, 0.000018f, 0.000013f,
        0.000009f, 0.000007f, 0.000005f, 0.000003f, 0.000002f, 0.000002f, 0.000001f, 0.000001f, 0.000001f, 0.000000f,
        0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f,
        0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f,
        0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f,
        0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f};
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
    static constexpr float kBetaDist3[100] = {
        0.006715f, 0.012509f, 0.017463f, 0.021655f, 0.025155f, 0.028031f, 0.030344f, 0.032151f, 0.033506f, 0.034458f,
        0.035052f, 0.035331f, 0.035332f, 0.035092f, 0.034643f, 0.034015f, 0.033234f, 0.032327f, 0.031314f, 0.030217f,
        0.029054f, 0.027841f, 0.026592f, 0.025322f, 0.024042f, 0.022761f, 0.021489f, 0.020234f, 0.019002f, 0.017799f,
        0.016630f, 0.015499f, 0.014409f, 0.013362f, 0.012361f, 0.011407f, 0.010500f, 0.009641f, 0.008830f, 0.008067f,
        0.007351f, 0.006681f, 0.006056f, 0.005475f, 0.004936f, 0.004437f, 0.003978f, 0.003555f, 0.003168f, 0.002814f,
        0.002492f, 0.002199f, 0.001934f, 0.001695f, 0.001481f, 0.001288f, 0.001116f, 0.000963f, 0.000828f, 0.000708f,
        0.000603f, 0.000511f, 0.000431f, 0.000361f, 0.000301f, 0.000250f, 0.000206f, 0.000168f, 0.000137f, 0.000110f,
        0.000088f, 0.000070f, 0.000055f, 0.000043f, 0.000033f, 0.000025f, 0.000019f, 0.000014f, 0.000010f, 0.000007f,
        0.000005f, 0.000004f, 0.000002f, 0.000002f, 0.000001f, 0.000001f, 0.000000f, 0.000000f, 0.000000f, 0.000000f,
        0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f};
    static constexpr float kBetaDist4[100] = {
        0.003996f, 0.007596f, 0.010824f, 0.013703f, 0.016255f, 0.018501f, 0.020460f, 0.022153f, 0.023597f, 0.024809f,
        0.025807f, 0.026607f, 0.027223f, 0.027671f, 0.027963f, 0.028114f, 0.028135f, 0.028038f, 0.027834f, 0.027535f,
        0.027149f, 0.026687f, 0.026157f, 0.025567f, 0.024926f, 0.024240f, 0.023517f, 0.022763f, 0.021983f, 0.021184f,
        0.020371f, 0.019548f, 0.018719f, 0.017890f, 0.017062f, 0.016241f, 0.015428f, 0.014627f, 0.013839f, 0.013068f,
        0.012315f, 0.011582f, 0.010870f, 0.010181f, 0.009515f, 0.008874f, 0.008258f, 0.007668f, 0.007103f, 0.006565f,
        0.006053f, 0.005567f, 0.005107f, 0.004673f, 0.004264f, 0.003880f, 0.003521f, 0.003185f, 0.002872f, 0.002581f,
        0.002312f, 0.002064f, 0.001835f, 0.001626f, 0.001434f, 0.001260f, 0.001102f, 0.000959f, 0.000830f, 0.000715f,
        0.000612f, 0.000521f, 0.000440f, 0.000369f, 0.000308f, 0.000254f, 0.000208f, 0.000169f, 0.000136f, 0.000108f,
        0.000084f, 0.000065f, 0.000050f, 0.000037f, 0.000027f, 0.000019f, 0.000014f, 0.000009f, 0.000006f, 0.000004f,
        0.000002f, 0.000001f, 0.000001f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f};
    static constexpr float kSingle10[100] = {
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 1.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f};
    static constexpr float kSingle15[100] = {
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 1.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f};
    static constexpr float kSingle20[100] = {
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 1.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f,
        0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f, 0.00000f};

    static constexpr float const* kThresholdDistributions[8] = {kUniformDist, kBetaDist1, kBetaDist2, kBetaDist3,
                                                                kBetaDist4,   kSingle10,  kSingle15,  kSingle20};

    static constexpr float kMinWeight = 0.01f;

    // --------------------------------------------------------
    // 成员变量
    // --------------------------------------------------------
    float fs_{};
    int block_size_{};
    int yin_buffer_size_{};
    int thresh_distr_{2}; // 默认 Beta(μ=0.15)
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
