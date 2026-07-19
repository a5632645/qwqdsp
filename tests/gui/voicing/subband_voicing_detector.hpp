#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>

#include <qwqdsp/filter/biquad.hpp>
#include <qwqdsp/filter/rbj.hpp>
#include <qwqdsp/pitch/fixed_voicing_detector.hpp>

#include "rms_tracker.hpp"

// ------------------------------------------------------------
//  RatioToProbability
// ------------------------------------------------------------
//  将能量比映射到浊音概率（不含静音检测）。
//  prob = 1 − clamp((log10(ratio) + center) / width, 0, 1)
// ------------------------------------------------------------
inline float RatioToProbability(float energy_ratio, float center, float width) noexcept {
    float log_r = std::log10(energy_ratio + 1e-10f);
    float prob = 1.0f - std::clamp((log_r + center) / width, 0.0f, 1.0f);
    return std::clamp(prob, 0.0f, 1.0f);
}

// ------------------------------------------------------------
//  SubbandVoicingDetector
// ------------------------------------------------------------
//  基于高低频能量比的实时清/浊音检测器。
//
//  算法：
//    时域并联 Biquad LPF@800Hz + HPF@3kHz (RBJ, Q=√2/2)
//    短时 RMS → ratio = (hp_scale * RMS_high) / (lp_scale * RMS_low + delta)
//    log10(ratio) 使用 center/width 映射到 [0,1] 浊音概率
//
//  可调参数 (mutex 保护):
//    setLpScale(v)    — 低频权重 (默认 1.0)
//    setHpScale(v)    — 高频权重 (默认 1.0)
//    setDelta(v)      — 分母底噪 (默认 1e-10)
//    setCenter(v)     — log10(ratio) 中点偏移 (默认 1.0)
//    setWidth(v)      — V/UV 过渡区宽度 (默认 2.5)
//
//  用法：
//    SubbandVoicingDetector detector(sample_rate);
//    for (each sample)     detector.processSample(sample);
//    if (frame full)       auto result = detector.frameResult();
// ------------------------------------------------------------

struct SubbandVoicingDetector {
    /// @param sample_rate 采样率 (Hz)
    SubbandVoicingDetector(float sample_rate) noexcept
        : sample_rate_(sample_rate) {
        setLpFreq(800.0f);
        setHpFreq(3000.0f);
        setRmsTau(0.03f);
    }

    /// 重置累积和滤波器状态
    void reset() noexcept {
        lpf_.Reset();
        hpf_.Reset();
        full_rms_.reset();
        lp_rms_.reset();
        hp_rms_.reset();
        last_prob_ = 0.0f;
    }

    /// 每采样处理：Biquad 滤波 + 短时 RMS 更新
    void processSample(float sample) noexcept {
        full_rms_.addSample(sample);

        float lp = lpf_.Tick(sample);
        lp_rms_.addSample(lp);

        float hp = hpf_.Tick(sample);
        hp_rms_.addSample(hp);
    }

    // --------------------------------------------------------
    //  可调参数 (mutex 保护)
    // --------------------------------------------------------
    void setLpScale(float v) noexcept {
        lp_scale_ = v;
    }
    void setHpScale(float v) noexcept {
        hp_scale_ = v;
    }
    void setDelta(float v) noexcept {
        delta_ = v;
    }
    void setCenter(float v) noexcept {
        center_ = v;
    }
    void setWidth(float v) noexcept {
        width_ = v;
    }
    void setRmsTau(float tau_sec) noexcept {
        full_rms_.setTimeConstant(tau_sec, sample_rate_);
        lp_rms_.setTimeConstant(tau_sec, sample_rate_);
        hp_rms_.setTimeConstant(tau_sec, sample_rate_);
    }
    void setLpFreq(float hz) noexcept {
        float w = 2.0f * std::numbers::pi_v<float> * hz / sample_rate_;
        qwqdsp_filter::RBJ coeff;
        coeff.Lowpass(w, kQ);
        lpf_.Set(coeff.ToBiquadCoeff());
        lpf_.Reset();
    }
    void setHpFreq(float hz) noexcept {
        float w = 2.0f * std::numbers::pi_v<float> * hz / sample_rate_;
        qwqdsp_filter::RBJ coeff;
        coeff.Highpass(w, kQ);
        hpf_.Set(coeff.ToBiquadCoeff());
        hpf_.Reset();
    }

    float center() const noexcept {
        return center_;
    }
    float width() const noexcept {
        return width_;
    }

    // --------------------------------------------------------
    //  frameResult
    // --------------------------------------------------------
    //  读取当前短时 RMS 值计算能量比和浊音概率。
    //  不重置 RMS 状态（连续运行）。
    //
    //  @param rms        输出：全频段短时线性 RMS
    //  @return            浊音概率 [0, 1]
    // --------------------------------------------------------
    float frameResult(float& out_energy_ratio, float& out_rms) noexcept {
        out_rms = full_rms_.value();

        float lp_val = lp_rms_.value();
        float hp_val = hp_rms_.value();

        out_energy_ratio = (hp_scale_ * hp_val) / (lp_scale_ * lp_val + delta_);
        float prob = probabilityFromRatio(out_rms, out_energy_ratio, center_, width_);
        last_prob_ = prob;
        return prob;
    }

    /// 获取最近一次 frameResult 的概率
    float lastProbability() const noexcept {
        return last_prob_;
    }
private:
    // --------------------------------------------------------
    //  probabilityFromRatio
    // --------------------------------------------------------
    //  将能量比映射到浊音概率。
    //  低 ratio (<1) → 浊音 → 高概率
    //  高 ratio (>1) → 清音 → 低概率
    //  静音 (rms < 0.002) → 概率 0
    // --------------------------------------------------------
    static float probabilityFromRatio(float rms, float energy_ratio, float center, float width) noexcept {
        if (rms < 0.002f)
            return 0.0f;
        return RatioToProbability(energy_ratio, center, width);
    }

    static constexpr float kQ = std::numbers::sqrt2_v<float> / 2.0f; // Butterworth

    qwqdsp_filter::Biquad lpf_;
    qwqdsp_filter::Biquad hpf_;

    // 短时 RMS（连续运行，不重置）
    RmsTracker full_rms_;
    RmsTracker lp_rms_;
    RmsTracker hp_rms_;

    float sample_rate_;

    // 可调参数
    float lp_scale_ = 1.0f;
    float hp_scale_ = 1.0f;
    float delta_ = 1e-10f;
    float center_ = 1.0f;
    float width_ = 2.5f;

    // 最近一次概率
    float last_prob_ = 0.0f;
};
