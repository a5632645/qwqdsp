#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>

#include "../filter/biquad.hpp"
#include "../filter/rbj.hpp"

namespace qwqdsp_pitch {

// ------------------------------------------------------------
//  SimpleVoicingDetector
// ------------------------------------------------------------
//  基于高低频能量比的实时清/浊音检测器。
//
//  算法：clamp(1 − hp_scale × hp_msq / (lp_msq + 1e-5), 0, 1)
//  使用均值平方（Mean Square）而非 RMS，省去开方。
//
//  固定参数:
//    Threshold  = 0.5      LP Freq  = 500 Hz
//    HP Freq    = 5000 Hz  RMS Tau  = 5 ms
//  运行时参数:
//    Sensitivity（实质 HP Scale，默认 1.0）
//    RMS Tau
// ------------------------------------------------------------

struct SimpleVoicingDetector {
    static constexpr float kThreshold = 0.5f;

    static constexpr float kLpFreq = 500.0f;
    static constexpr float kHpFreq = 5000.0f;
    static constexpr float kRmsTau = 0.005f; // 5 ms

    SimpleVoicingDetector() noexcept = default;

    void Init(float sample_rate) noexcept {
        sample_rate_ = sample_rate;
        UpdateRmsAlpha();

        float const kQ = std::numbers::sqrt2_v<float> / 2.0f;
        float const w_lp = 2.0f * std::numbers::pi_v<float> * kLpFreq / sample_rate;
        float const w_hp = 2.0f * std::numbers::pi_v<float> * kHpFreq / sample_rate;

        qwqdsp_filter::RBJ coeff;
        coeff.Lowpass(w_lp, kQ);
        lpf_.Set(coeff.ToBiquadCoeff());
        coeff.Highpass(w_hp, kQ);
        hpf_.Set(coeff.ToBiquadCoeff());
    }

    void Reset() noexcept {
        lpf_.Reset();
        hpf_.Reset();
        full_mean_sq_ = 0.0f;
        lp_mean_sq_ = 0.0f;
        hp_mean_sq_ = 0.0f;
        last_prob_ = 0.0f;
    }

    void ProcessSample(float sample) noexcept {
        full_mean_sq_ = rms_alpha_ * sample * sample + (1.0f - rms_alpha_) * full_mean_sq_;

        float lp = lpf_.Tick(sample);
        lp_mean_sq_ = rms_alpha_ * lp * lp + (1.0f - rms_alpha_) * lp_mean_sq_;

        float hp = hpf_.Tick(sample);
        hp_mean_sq_ = rms_alpha_ * hp * hp + (1.0f - rms_alpha_) * hp_mean_sq_;
    }

    float FrameResult() noexcept {
        float prob = std::clamp(1.0f - (hp_scale_ * hp_mean_sq_) / (lp_mean_sq_ + 1e-5f), 0.0f, 1.0f);
        last_prob_ = prob;
        return prob;
    }

    float LastProbability() const noexcept {
        return last_prob_;
    }

    // 调节敏感度（实质为 HP Scale，默认 1.0，越大越易判为清音）
    void SetSensitivity(float v) noexcept {
        hp_scale_ = v;
    }

    // 设置 RMS 时间常数（秒），默认 0.005 (5 ms)
    void SetRmsTau(float tau_sec) noexcept {
        rms_tau_ = tau_sec;
        UpdateRmsAlpha();
    }
private:
    void UpdateRmsAlpha() noexcept {
        rms_alpha_ = 1.0f - std::exp(-1.0f / (rms_tau_ * sample_rate_));
    }

    qwqdsp_filter::Biquad lpf_;
    qwqdsp_filter::Biquad hpf_;

    float full_mean_sq_ = 0.0f;
    float lp_mean_sq_ = 0.0f;
    float hp_mean_sq_ = 0.0f;

    float rms_alpha_ = 1.0f;
    float sample_rate_ = 0.0f;
    float rms_tau_ = 0.005f;
    float hp_scale_ = 1.0f;
    float last_prob_ = 0.0f;
};

} // namespace qwqdsp_pitch
