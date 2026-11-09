#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "qwqdsp/convert.hpp"
#include "qwqdsp/misc/smoother.hpp"

namespace maximizer_demo {

enum class MaximizerMode {
    Iir,  // sliding minimum + one-pole attack/release
    Fir,  // sliding minimum + release + cascaded-box FIR attack smoothing
};

template <MaximizerMode Mode>
class StereoMaximizerT;

using StereoMaximizer = StereoMaximizerT<MaximizerMode::Iir>;
using FirMaximizer = StereoMaximizerT<MaximizerMode::Fir>;

// ------------------------------------------------------------
//  BoxFilter — FIR 移动平均，用于 FIR 模式的 attack smoothing
// ------------------------------------------------------------
class BoxFilter {
public:
    static constexpr size_t kMaxLength = 1024;

    void SetLength(size_t length) noexcept {
        length_ = std::clamp<size_t>(length, 1, kMaxLength);
        buffer_.fill(1.0f);
        sum_ = static_cast<float>(length_);
        pos_ = 0;
    }

    void Reset(float value) noexcept {
        length_ = std::max<size_t>(length_, 1);
        buffer_.fill(value);
        sum_ = value * static_cast<float>(length_);
        pos_ = 0;
    }

    float Tick(float x) noexcept {
        sum_ += x - buffer_[pos_];
        buffer_[pos_] = x;
        pos_ = (pos_ + 1) % length_;
        return sum_ / static_cast<float>(length_);
    }

private:
    std::array<float, kMaxLength> buffer_{};
    size_t length_{1};
    size_t pos_{};
    float sum_{1.0f};
};

// ------------------------------------------------------------
//  GainSmoother — 编译期选择的增益平滑策略
// ------------------------------------------------------------
template <MaximizerMode Mode>
class GainSmoother;

// 原来的实现：一阶 attack/release
template <>
class GainSmoother<MaximizerMode::Iir> {
public:
    void Init(float sample_rate) noexcept {
        sample_rate_ = sample_rate;
        Reset();
    }

    void Reset() noexcept {
        gain_ = 1.0f;
    }

    void SetParameters(size_t attack_samples, size_t /*hold_samples*/,
                       float release_ms, float sample_rate) noexcept {
        sample_rate_ = sample_rate;
        float const attack_ms = 1000.0f * static_cast<float>(attack_samples) / sample_rate;
        float const attack = std::max(0.05f, attack_ms / 12.0f);
        attack_coeff_ = qwqdsp_misc::ExpSmoother::ComputeSmoothFactor(attack, sample_rate);
        release_coeff_ = qwqdsp_misc::ExpSmoother::ComputeSmoothFactor(release_ms, sample_rate);
    }

    float Tick(float target_gain) noexcept {
        float const coeff = target_gain < gain_ ? attack_coeff_ : release_coeff_;
        gain_ = target_gain + coeff * (gain_ - target_gain);
        return gain_;
    }

private:
    float sample_rate_{};
    float attack_coeff_{};
    float release_coeff_{};
    float gain_{1.0f};
};

// Signalsmith 风格：release 在前，FIR attack smoothing 在后
template <>
class GainSmoother<MaximizerMode::Fir> {
public:
    void Init(float sample_rate) noexcept {
        sample_rate_ = sample_rate;
        Reset();
    }

    void Reset() noexcept {
        released_ = 1.0f;
        box1_.Reset(1.0f);
        box2_.Reset(1.0f);
    }

    void SetParameters(size_t attack_samples, size_t /*hold_samples*/,
                       float release_ms, float sample_rate) noexcept {
        sample_rate_ = sample_rate;
        release_coeff_ = qwqdsp_misc::ExpSmoother::ComputeSmoothFactor(
            std::max(release_ms, 0.1f), sample_rate);

        // 两个 box 串联，总 FIR 长度约为 attack_samples
        size_t const first = std::max<size_t>(1, attack_samples / 2);
        size_t const second = std::max<size_t>(1, attack_samples - first);
        if (first != first_length_ || second != second_length_) {
            first_length_ = first;
            second_length_ = second;
            box1_.SetLength(first);
            box2_.SetLength(second);
        }
    }

    float Tick(float target_gain) noexcept {
        // 立即下降，指数恢复
        released_ += (target_gain - released_) * release_coeff_;
        released_ = std::min(released_, target_gain);

        // FIR attack smoothing
        return box2_.Tick(box1_.Tick(released_));
    }

private:
    float sample_rate_{};
    float release_coeff_{};
    float released_{1.0f};
    size_t first_length_{};
    size_t second_length_{};
    BoxFilter box1_;
    BoxFilter box2_;
};

// ------------------------------------------------------------
//  StereoMaximizerT — input gain + look-ahead 峰值限制器
//
//  公共部分：
//    peak -> required gain -> moving minimum -> [Mode 平滑] -> 延迟信号 * gain
//
//  IIR 模式：moving minimum + 一阶 attack/release
//  FIR 模式：moving minimum + release + cascaded box FIR
// ------------------------------------------------------------
template <MaximizerMode Mode>
class StereoMaximizerT {
public:
    static constexpr size_t kMaxLookaheadSamples = 1024;
    static constexpr size_t kMaxHoldSamples = 4096;
    static constexpr float kMaxLookaheadMs = 10.0f;
    static constexpr float kMaxHoldMs = 80.0f;

    void Init(float sample_rate) noexcept {
        sample_rate_ = sample_rate;
        input_gain_smoother_.SetSmoothTime(20.0f, sample_rate);
        input_gain_smoother_.SetTargetImmediately(1.0f);
        smoother_.Init(sample_rate);
        Reset();
    }

    void Reset() noexcept {
        input_gain_smoother_.Reset();
        ResetDelay();
        smoother_.Reset();
    }

    void SetParameters(float input_gain_db, float ceiling_db, float release_ms,
                       float lookahead_ms, float hold_ms = 0.0f) noexcept {
        input_gain_smoother_.SetTarget(qwqdsp::convert::Db2Gain(input_gain_db));
        ceiling_ = qwqdsp::convert::Db2Gain(ceiling_db);

        float const lookahead = std::clamp(lookahead_ms, 0.1f, kMaxLookaheadMs);
        float const hold = std::clamp(hold_ms, 0.0f, kMaxHoldMs);
        size_t const attack_samples = std::clamp<size_t>(
            static_cast<size_t>(std::lround(lookahead * sample_rate_ / 1000.0f)),
            1, kMaxLookaheadSamples);
        size_t const hold_samples = std::min<size_t>(
            static_cast<size_t>(std::lround(hold * sample_rate_ / 1000.0f)),
            kMaxHoldSamples);

        if (attack_samples != attack_samples_ || hold_samples != hold_samples_) {
            attack_samples_ = attack_samples;
            hold_samples_ = hold_samples;
            if constexpr (Mode == MaximizerMode::Iir) {
                window_samples_ = attack_samples_ + 1;
            }
            else {
                window_samples_ = attack_samples_ + hold_samples_ + 1;
            }
            ResetDelay();
            smoother_.Reset();
        }

        smoother_.SetParameters(attack_samples_, hold_samples_, release_ms, sample_rate_);
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
            float const gain = smoother_.Tick(target_gain);

            float out_l = l;
            float out_r = r;
            if (attack_samples_ > 0) {
                out_l = delay_l_[delay_pos_];
                out_r = delay_r_[delay_pos_];
                delay_l_[delay_pos_] = l;
                delay_r_[delay_pos_] = r;
                delay_pos_ = (delay_pos_ + 1) % attack_samples_;
            }

            buffer[2 * i] = out_l * gain;
            buffer[2 * i + 1] = out_r * gain;
            gain_reduction_db_ = 20.0f * std::log10(std::max(gain, 1.0e-6f));
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

    static constexpr size_t kDequeCapacity = kMaxLookaheadSamples + kMaxHoldSamples + 2;

    void ResetDelay() noexcept {
        delay_l_.fill(0.0f);
        delay_r_.fill(0.0f);
        delay_pos_ = 0;
        deque_head_ = 0;
        deque_size_ = 0;
        sample_index_ = 0;
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

        uint64_t const window = static_cast<uint64_t>(window_samples_);
        while (deque_size_ > 0 && deque_[deque_head_].index + window <= sample_index_) {
            deque_head_ = (deque_head_ + 1) % kDequeCapacity;
            --deque_size_;
        }

        ++sample_index_;
    }

    float sample_rate_{};
    float ceiling_{1.0f};
    size_t attack_samples_{};
    size_t hold_samples_{};
    size_t window_samples_{1};
    size_t delay_pos_{};

    std::array<float, kMaxLookaheadSamples> delay_l_{};
    std::array<float, kMaxLookaheadSamples> delay_r_{};

    std::array<GainEntry, kDequeCapacity> deque_{};
    size_t deque_head_{};
    size_t deque_size_{};
    uint64_t sample_index_{};

    qwqdsp_misc::ExpSmoother input_gain_smoother_;
    GainSmoother<Mode> smoother_;
    float gain_reduction_db_{};
};

} // namespace maximizer_demo
