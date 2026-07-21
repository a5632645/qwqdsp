#pragma once
#include "../../spectral/real_fft.hpp"
#include "../pitch.hpp"
#include "pyin_core.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <span>
#include <vector>

namespace qwqdsp_pitch {

// ------------------------------------------------------------
// PmpmCore
// ------------------------------------------------------------
// 概率 MPM (pMPM) 核心算法。
// 对每帧计算 NSDF 自相关 → 峰值拾取 → 多 cutoff 扫描 → 概率累积。
// 参考: Mauch & Dixon pYIN 思路 + McLeod MPM 方法

class PmpmCore {
public:
    void Init(float sample_rate, int block_size) {
        fs_ = sample_rate;
        block_size_ = block_size;
        fft_.Init(block_size * 2);
        fft_in_.resize(block_size * 2);
        fft_out_.resize(block_size * 2 + 2);
        autocorr_.resize(block_size);
        SetMinPitch(min_pitch_);
        SetMaxPitch(max_pitch_);
    }

    void SetMinPitch(float hz) noexcept {
        min_pitch_ = hz;
        max_tau_ = std::min(static_cast<int>(std::round(fs_ / hz)), block_size_);
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
        ComputeAutocorrelation(block);
        auto peaks = PeakPicking();

        // ---- 抛物线插值所有有效峰值 ----
        estimates_.clear();
        float highest_amplitude = -std::numeric_limits<float>::max();
        for (int idx : peaks) {
            float val = autocorr_[static_cast<size_t>(idx)];
            highest_amplitude = std::max(highest_amplitude, val);
            if (val > kSmallCutoff) {
                auto est = ParabolicInterpolation(autocorr_, static_cast<size_t>(idx));
                estimates_.push_back(est);
                highest_amplitude = std::max(highest_amplitude, est.second);
            }
        }

        if (estimates_.empty())
            return {};

        // ---- 多 cutoff 概率累积 ----
        std::map<int, float> tau_prob; // tau → 累积概率

        for (int n = 0; n < kNumCutoffs; ++n) {
            float cutoff = kCutoffBegin + static_cast<float>(n) * kCutoffStep;
            float actual_cutoff = cutoff * highest_amplitude;

            int period_tau = -1;
            for (auto const& est : estimates_) {
                if (est.second >= actual_cutoff) {
                    period_tau = static_cast<int>(std::round(est.first));
                    break;
                }
            }

            float weight = (period_tau > 0) ? 1.0f : kPa;
            tau_prob[period_tau] += weight * kProbDist;
        }

        // ---- 转换为基频候选 ----
        std::vector<PyinCandidate> candidates;
        for (auto const& [tau, prob] : tau_prob) {
            if (tau <= 0 || prob < min_prob)
                continue;
            float tau_interp = ParabolicInterpolationFromBuffer(static_cast<size_t>(tau));
            float pitch = fs_ / tau_interp;
            if (pitch >= 20.0f && pitch <= 2000.0f) {
                candidates.push_back({pitch, prob});
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
private:
    // --------------------------------------------------------
    // computeAutocorrelation — FFT 加速自相关
    // --------------------------------------------------------
    void ComputeAutocorrelation(std::span<const float> block) noexcept {
        // 零填充到 2*block_size
        std::copy_n(block.begin(), block.size(), fft_in_.begin());
        std::fill_n(fft_in_.begin() + block.size(), block.size(), 0.0f);

        // FFT → 功率谱
        fft_.FFT(fft_in_.data(), fft_out_.data());
        size_t const num_bins = fft_.GetFFTSize() / 2 + 1;
        for (size_t i = 0; i < num_bins; ++i) {
            float re = fft_out_[2 * i];
            float im = fft_out_[2 * i + 1];
            fft_out_[2 * i] = re * re + im * im;
            fft_out_[2 * i + 1] = 0.0f;
        }

        // IFFT → 自相关
        fft_.IFFT(fft_out_.data(), fft_in_.data());

        // 提取前 block_size 个有效滞后
        for (int i = 0; i < block_size_; ++i) {
            autocorr_[static_cast<size_t>(i)] = fft_in_[static_cast<size_t>(i)];
        }
    }

    // --------------------------------------------------------
    // peakPicking — 在自相关函数中寻找正峰值
    // --------------------------------------------------------
    std::vector<int> const& PeakPicking() noexcept {
        peaks_.clear();

        int const size = max_tau_;
        int pos = min_tau_;

        // 跳过初始正区
        while (pos < (size - 1) / 3 && autocorr_[static_cast<size_t>(pos)] > 0) {
            ++pos;
        }
        while (pos < size - 1 && autocorr_[static_cast<size_t>(pos)] <= 0.0f) {
            ++pos;
        }

        if (pos == 0) {
            pos = 1;
        }

        int cur_max_pos = 0;
        while (pos < size - 1) {
            if (autocorr_[static_cast<size_t>(pos)] > autocorr_[static_cast<size_t>(pos) - 1]
                && autocorr_[static_cast<size_t>(pos)] >= autocorr_[static_cast<size_t>(pos) + 1]
                && (cur_max_pos == 0
                    || autocorr_[static_cast<size_t>(pos)] > autocorr_[static_cast<size_t>(cur_max_pos)])) {
                cur_max_pos = pos;
            }
            ++pos;
            if (pos < size - 1 && autocorr_[static_cast<size_t>(pos)] <= 0.0f) {
                if (cur_max_pos > 0) {
                    peaks_.push_back(cur_max_pos);
                    cur_max_pos = 0;
                }
                while (pos < size - 1 && autocorr_[static_cast<size_t>(pos)] <= 0.0f) {
                    ++pos;
                }
            }
        }
        if (cur_max_pos > 0) {
            peaks_.push_back(cur_max_pos);
        }

        return peaks_;
    }

    // --------------------------------------------------------
    // parabolicInterpolation
    // --------------------------------------------------------
    static std::pair<float, float> ParabolicInterpolation(std::vector<float> const& array, size_t x) noexcept {
        if (x < 1) {
            size_t x_adj = (array[x] <= array[x + 1]) ? x : x + 1;
            return {static_cast<float>(x_adj), array[x_adj]};
        }
        if (x >= array.size() - 1) {
            size_t x_adj = (array[x] <= array[x - 1]) ? x : x - 1;
            return {static_cast<float>(x_adj), array[x_adj]};
        }

        float den = array[x + 1] + array[x - 1] - 2.0f * array[x];
        if (den == 0.0f) {
            return {static_cast<float>(x), array[x]};
        }
        float delta = array[x - 1] - array[x + 1];
        float tau = static_cast<float>(x) + delta / (2.0f * den);
        float amp = array[x] - delta * delta / (8.0f * den);
        return {tau, amp};
    }

    float ParabolicInterpolationFromBuffer(size_t tau_estimate) const noexcept {
        if (tau_estimate <= 0 || tau_estimate >= autocorr_.size() - 1) {
            return static_cast<float>(tau_estimate);
        }
        float s0 = autocorr_[tau_estimate - 1];
        float s1 = autocorr_[tau_estimate];
        float s2 = autocorr_[tau_estimate + 1];
        double adj = static_cast<double>(s2 - s0) / (2.0 * (2.0 * s1 - s2 - s0));
        if (std::abs(adj) > 1.0) {
            adj = 0.0;
        }
        return static_cast<float>(static_cast<double>(tau_estimate) + adj);
    }

    static constexpr float kCutoffBegin = 0.80f;
    static constexpr float kCutoffStep = 0.01f;
    static constexpr int kNumCutoffs = 20;
    static constexpr float kProbDist = 1.0f / static_cast<float>(kNumCutoffs); // 0.05
    static constexpr float kPa = 0.01f;
    static constexpr float kSmallCutoff = 0.5f;
    static constexpr float kMinPitchDefault = 50.0f;
    static constexpr float kMaxPitchDefault = 500.0f;

    float fs_{};
    int block_size_{};
    int min_tau_{2};
    int max_tau_{};
    float min_pitch_{kMinPitchDefault};
    float max_pitch_{kMaxPitchDefault};

    qwqdsp_spectral::RealFFT fft_;
    std::vector<float> fft_in_;
    std::vector<float> fft_out_;
    std::vector<float> autocorr_;
    std::vector<int> peaks_;
    std::vector<std::pair<float, float>> estimates_;
};

} // namespace qwqdsp_pitch
