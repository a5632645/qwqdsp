#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "qwqdsp/convert.hpp"
#include "qwqdsp/misc/smoother.hpp"

namespace maximizer_demo {

// ------------------------------------------------------------
//  StereoMaximizer — input gain + 立体声 look-ahead 峰值限制器
//
//  链路：
//    input gain -> 峰值检测 -> look-ahead 最小增益 -> attack/release 平滑
//                -> 延迟后的信号 * gain -> ceiling
//
//  look-ahead 窗口内取“所需增益”的最小值，保证峰值到达输出时
//  增益已经提前压下；attack 时间取 look-ahead / 8，使增益在
//  峰值到达前基本收敛。检测为 sample peak，不做 true-peak 过采样。
// ------------------------------------------------------------
class StereoMaximizer {
public:
    static constexpr size_t kMaxLookaheadSamples = 1024;
    static constexpr float kMaxLookaheadMs = 10.0f;

    void Init(float sample_rate) noexcept {
        sample_rate_ = sample_rate;
        input_gain_smoother_.SetSmoothTime(20.0f, sample_rate);
        input_gain_smoother_.SetTargetImmediately(1.0f);
        Reset();
    }

    void Reset() noexcept {
        input_gain_smoother_.Reset();
        ResetDelay();
    }

    void SetParameters(float input_gain_db, float ceiling_db, float release_ms,
                       float lookahead_ms) noexcept {
        input_gain_smoother_.SetTarget(qwqdsp::convert::Db2Gain(input_gain_db));
        ceiling_ = qwqdsp::convert::Db2Gain(ceiling_db);
        release_coeff_ = qwqdsp_misc::ExpSmoother::ComputeSmoothFactor(release_ms, sample_rate_);

        float const lookahead = std::clamp(lookahead_ms, 0.1f, kMaxLookaheadMs);
        size_t lookahead_samples =
            static_cast<size_t>(std::lround(lookahead * sample_rate_ / 1000.0f));
        lookahead_samples = std::clamp<size_t>(lookahead_samples, 1, kMaxLookaheadSamples);
        if (lookahead_samples != lookahead_samples_) {
            lookahead_samples_ = lookahead_samples;
            ResetDelay();
        }

        float const attack_ms = std::max(0.05f, lookahead / 8.0f);
        attack_coeff_ = qwqdsp_misc::ExpSmoother::ComputeSmoothFactor(attack_ms, sample_rate_);
    }

    /// @brief 原地处理交错双声道样本
    void Process(float* buffer, size_t frame_count) noexcept {
        for (size_t i = 0; i < frame_count; ++i) {
            float const input_gain = input_gain_smoother_.Tick();
            float const l = buffer[2 * i] * input_gain;
            float const r = buffer[2 * i + 1] * input_gain;

            float const peak = std::max(std::abs(l), std::abs(r));
            float required_gain = 1.0f;
            if (peak > ceiling_) {
                required_gain = ceiling_ / peak;
            }
            PushRequiredGain(required_gain);

            float const target_gain = deque_[deque_head_].value;
            float const coeff = target_gain < gain_ ? attack_coeff_ : release_coeff_;
            gain_ = target_gain + coeff * (gain_ - target_gain);

            float out_l = l;
            float out_r = r;
            if (lookahead_samples_ > 0) {
                out_l = delay_l_[delay_pos_];
                out_r = delay_r_[delay_pos_];
                delay_l_[delay_pos_] = l;
                delay_r_[delay_pos_] = r;
                delay_pos_ = (delay_pos_ + 1) % lookahead_samples_;
            }

            buffer[2 * i] = out_l * gain_;
            buffer[2 * i + 1] = out_r * gain_;
            gain_reduction_db_ = 20.0f * std::log10(std::max(gain_, 1.0e-6f));
        }
    }

    /// @brief 当前增益（dB，<=0）
    float GetGainReductionDb() const noexcept {
        return gain_reduction_db_;
    }

private:
    struct GainEntry {
        uint64_t index;
        float value;
    };

    static constexpr size_t kDequeCapacity = kMaxLookaheadSamples + 2;

    void ResetDelay() noexcept {
        delay_l_.fill(0.0f);
        delay_r_.fill(0.0f);
        delay_pos_ = 0;
        deque_head_ = 0;
        deque_size_ = 0;
        sample_index_ = 0;
        gain_ = 1.0f;
        gain_reduction_db_ = 0.0f;
    }

    void PushRequiredGain(float value) noexcept {
        while (deque_size_ > 0) {
            size_t const back = (deque_head_ + deque_size_ - 1) % kDequeCapacity;
            if (deque_[back].value < value) {
                break;
            }
            --deque_size_;
        }

        size_t const pos = (deque_head_ + deque_size_) % kDequeCapacity;
        deque_[pos] = GainEntry{sample_index_, value};
        ++deque_size_;

        uint64_t const window = static_cast<uint64_t>(lookahead_samples_) + 1;
        while (deque_size_ > 0 && deque_[deque_head_].index + window <= sample_index_) {
            deque_head_ = (deque_head_ + 1) % kDequeCapacity;
            --deque_size_;
        }

        ++sample_index_;
    }

    float sample_rate_{};
    float ceiling_{1.0f};
    float attack_coeff_{};
    float release_coeff_{};
    float gain_{1.0f};
    float gain_reduction_db_{};
    size_t lookahead_samples_{};
    size_t delay_pos_{};

    std::array<float, kMaxLookaheadSamples> delay_l_{};
    std::array<float, kMaxLookaheadSamples> delay_r_{};

    std::array<GainEntry, kDequeCapacity> deque_{};
    size_t deque_head_{};
    size_t deque_size_{};
    uint64_t sample_index_{};

    qwqdsp_misc::ExpSmoother input_gain_smoother_;
};

} // namespace maximizer_demo
