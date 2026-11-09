#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>

#include "maximizer_dsp.hpp"
#include "qwqdsp/convert.hpp"
#include "qwqdsp/misc/smoother.hpp"

namespace hybrid_maximizer_demo {

// ------------------------------------------------------------
//  FeedbackLimiter — 检测输出电平的反馈式限制器
//
//  检测的是上一个输出样本的峰值，因此没有 look-ahead，
//  峰值到来时增益只能事后下降，会产生 overshoot。
//  这正是后面还要接 feedforward brickwall limiter 的原因。
// ------------------------------------------------------------
class FeedbackLimiter {
public:
    void Init(float sample_rate) noexcept {
        sample_rate_ = sample_rate;
        Reset();
    }

    void Reset() noexcept {
        gain_reduction_db_ = 0.0f;
        last_out_l_ = 0.0f;
        last_out_r_ = 0.0f;
    }

    void SetParameters(float threshold_db, float ratio, float attack_ms,
                       float release_ms) noexcept {
        threshold_db_ = threshold_db;
        ratio_ = std::max(ratio, 1.0f);
        attack_coeff_ = qwqdsp_misc::ExpSmoother::ComputeSmoothFactor(attack_ms, sample_rate_);
        release_coeff_ = qwqdsp_misc::ExpSmoother::ComputeSmoothFactor(release_ms, sample_rate_);
    }

    std::pair<float, float> ProcessSample(float l, float r) noexcept {
        // 反馈检测：用上一个输出样本判定当前增益
        float const level = std::max(std::abs(last_out_l_), std::abs(last_out_r_));
        float const level_db = 20.0f * std::log10(level + kEpsilon);
        float const over_db = level_db - threshold_db_;
        float const target_reduction_db = over_db > 0.0f ? over_db * (1.0f - 1.0f / ratio_) : 0.0f;

        float const coeff = target_reduction_db > gain_reduction_db_ ? attack_coeff_ : release_coeff_;
        gain_reduction_db_ = target_reduction_db + coeff * (gain_reduction_db_ - target_reduction_db);

        float const gain = qwqdsp::convert::Db2Gain(-gain_reduction_db_);
        last_out_l_ = l * gain;
        last_out_r_ = r * gain;
        return {last_out_l_, last_out_r_};
    }

    /// @brief 当前增益（dB，<=0）
    float GetGainReductionDb() const noexcept {
        return -gain_reduction_db_;
    }

private:
    static constexpr float kEpsilon = 1.0e-10f;

    float sample_rate_{};
    float threshold_db_{-6.0f};
    float ratio_{20.0f};
    float attack_coeff_{};
    float release_coeff_{};
    float gain_reduction_db_{};
    float last_out_l_{};
    float last_out_r_{};
};

// ------------------------------------------------------------
//  HybridMaximizer — feedback limiter + feedforward brickwall
//
//  input gain -> feedback limiter -> feedforward look-ahead limiter
//  反馈级负责平滑的模拟式动作，前馈级负责保证最终 ceiling。
// ------------------------------------------------------------
class HybridMaximizer {
public:
    void Init(float sample_rate) noexcept {
        feedback_.Init(sample_rate);
        final_.Init(sample_rate);
        input_gain_smoother_.SetSmoothTime(20.0f, sample_rate);
        input_gain_smoother_.SetTargetImmediately(1.0f);
        Reset();
    }

    void Reset() noexcept {
        feedback_.Reset();
        final_.Reset();
        input_gain_smoother_.Reset();
        feedback_gain_reduction_db_ = 0.0f;
        final_gain_reduction_db_ = 0.0f;
    }

    void SetParameters(float input_gain_db, float threshold_db, float ratio,
                       float attack_ms, float release_ms, float ceiling_db,
                       float lookahead_ms) noexcept {
        input_gain_smoother_.SetTarget(qwqdsp::convert::Db2Gain(input_gain_db));
        feedback_.SetParameters(threshold_db, ratio, attack_ms, release_ms);
        final_.SetParameters(0.0f, ceiling_db, release_ms, lookahead_ms);
    }

    void Process(float* buffer, size_t frame_count) noexcept {
        for (size_t i = 0; i < frame_count; ++i) {
            float const input_gain = input_gain_smoother_.Tick();
            auto const [l, r] = feedback_.ProcessSample(
                buffer[2 * i] * input_gain,
                buffer[2 * i + 1] * input_gain);
            buffer[2 * i] = l;
            buffer[2 * i + 1] = r;
        }

        final_.Process(buffer, frame_count);

        feedback_gain_reduction_db_ = feedback_.GetGainReductionDb();
        final_gain_reduction_db_ = final_.GetGainReductionDb();
    }

    float GetFeedbackGainReductionDb() const noexcept {
        return feedback_gain_reduction_db_;
    }

    float GetFinalGainReductionDb() const noexcept {
        return final_gain_reduction_db_;
    }

private:
    FeedbackLimiter feedback_;
    maximizer_demo::StereoMaximizer final_;
    qwqdsp_misc::ExpSmoother input_gain_smoother_;
    float feedback_gain_reduction_db_{};
    float final_gain_reduction_db_{};
};

} // namespace hybrid_maximizer_demo
