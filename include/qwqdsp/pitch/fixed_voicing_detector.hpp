#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>

#include "../filter/biquad.hpp"
#include "../filter/rbj.hpp"

namespace qwqdsp_pitch {

// ------------------------------------------------------------
//  FixedVoicingDetector
// ------------------------------------------------------------
//  基于高低频能量比的实时清/浊音检测器（全编译期参数）。
//
//  所有参数 constexpr 固定，无运行时可调接口。
//
//  固定参数:
//    Threshold  = 0.5      LP Scale  = 1.0       HP Scale  = 1.0
//    LP Freq    = 500 Hz   HP Freq   = 5000 Hz   Delta     = 0.0001
//    Center     = 0.0      Width     = 0.1       RMS Tau   = 5 ms
// ------------------------------------------------------------

struct FixedVoicingDetector {
    static constexpr float kThreshold = 0.5f;

    static constexpr float kLpFreq = 500.0f;
    static constexpr float kHpFreq = 5000.0f;
    static constexpr float kLpScale = 1.0f;

    static constexpr float kDelta = 0.0001f;
    static constexpr float kCenter = 0.0f;
    static constexpr float kWidth = 0.1f;
    static constexpr float kRmsTau = 0.005f; // 5 ms

    static constexpr float kSilenceRms = 0.002f;

    FixedVoicingDetector() noexcept = default;

    void Init(float sample_rate) noexcept {
        float alpha = 1.0f - std::exp(-1.0f / (kRmsTau * sample_rate));
        rms_alpha_ = alpha;

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
        full_mean_sq_ = 0.0;
        lp_mean_sq_ = 0.0;
        hp_mean_sq_ = 0.0;
        last_prob_ = 0.0f;
    }

    void ProcessSample(float sample) noexcept {
        double s = static_cast<double>(sample);
        full_mean_sq_ = rms_alpha_ * s * s + (1.0 - rms_alpha_) * full_mean_sq_;

        float lp = lpf_.Tick(sample);
        double ld = static_cast<double>(lp);
        lp_mean_sq_ = rms_alpha_ * ld * ld + (1.0 - rms_alpha_) * lp_mean_sq_;

        float hp = hpf_.Tick(sample);
        double hd = static_cast<double>(hp);
        hp_mean_sq_ = rms_alpha_ * hd * hd + (1.0 - rms_alpha_) * hp_mean_sq_;
    }

    float FrameResult(float& out_energy_ratio, float& out_rms) noexcept {
        out_rms = static_cast<float>(std::sqrt(full_mean_sq_));

        float lp_val = static_cast<float>(std::sqrt(lp_mean_sq_));
        float hp_val = static_cast<float>(std::sqrt(hp_mean_sq_));
        out_energy_ratio = (hp_scale_ * hp_val) / (kLpScale * lp_val + kDelta);

        float prob = ProbabilityFromRatio(out_rms, out_energy_ratio);
        last_prob_ = prob;
        return prob;
    }

    float LastProbability() const noexcept {
        return last_prob_;
    }

    /// 调节敏感度（实质为 HP Scale，默认 1.0，越大越易判为清音）
    void SetSensitivity(float v) noexcept {
        hp_scale_ = v;
    }
private:
    static float ProbabilityFromRatio(float rms, float ratio) noexcept {
        if (rms < kSilenceRms)
            return 0.0f;
        float log_r = std::log10(ratio + 1e-10f);
        float prob = 1.0f - std::clamp((log_r + kCenter) / kWidth, 0.0f, 1.0f);
        return std::clamp(prob, 0.0f, 1.0f);
    }

    qwqdsp_filter::Biquad lpf_;
    qwqdsp_filter::Biquad hpf_;

    double full_mean_sq_ = 0.0;
    double lp_mean_sq_ = 0.0;
    double hp_mean_sq_ = 0.0;

    double rms_alpha_ = 1.0;
    float hp_scale_ = 1.0f;
    float last_prob_ = 0.0f;
};

} // namespace qwqdsp_pitch
