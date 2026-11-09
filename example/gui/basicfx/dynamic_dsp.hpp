#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "qwqdsp/convert.hpp"
#include "qwqdsp/misc/smoother.hpp"

namespace dynamic_demo {

enum class DynamicsMode : int {
    Compressor = 0,
    Expander,
    Gate,
    kNumModes
};

// ------------------------------------------------------------
//  DynamicsProcessor — block 峰值检测的 feed-forward 动态处理器
//
//  每个 block：
//    1. 取 L/R 联动峰值
//    2. 按当前模式（compressor / expander / gate）计算目标增益
//    3. 逐样本用攻击/释放一阶平滑把增益推向目标
//
//  compressor：超过 threshold 后按 ratio 压缩
//  expander  ：低于 threshold 后按 ratio 向下扩展，最多衰减 range
//  gate      ：低于 threshold 后衰减 range，带 3 dB 回差
// ------------------------------------------------------------
class DynamicsProcessor {
public:
    void Init(float sample_rate) noexcept {
        sample_rate_ = sample_rate;
        makeup_smoother_.SetSmoothTime(20.0f, sample_rate);
        makeup_smoother_.SetTargetImmediately(1.0f);
        Reset();
    }

    void Reset() noexcept {
        current_gain_ = 1.0f;
        gate_open_ = true;
        block_gain_reduction_db_ = 0.0f;
        makeup_smoother_.Reset();
    }

    void SetParameters(DynamicsMode mode, float threshold_db, float ratio, float range_db,
                       float attack_ms, float release_ms, float makeup_db) noexcept {
        mode_ = mode;
        threshold_db_ = threshold_db;
        ratio_ = std::max(ratio, 1.0f);
        range_db_ = std::max(range_db, 0.0f);
        attack_coeff_ = qwqdsp_misc::ExpSmoother::ComputeSmoothFactor(attack_ms, sample_rate_);
        release_coeff_ = qwqdsp_misc::ExpSmoother::ComputeSmoothFactor(release_ms, sample_rate_);

        // gate 不叠加 makeup
        makeup_smoother_.SetTarget(mode == DynamicsMode::Gate
                                       ? 1.0f
                                       : qwqdsp::convert::Db2Gain(makeup_db));
    }

    /// @brief 原地处理交错双声道样本
    void Process(float* buffer, size_t frame_count) noexcept {
        // 1. block 内双声道联合峰值
        float peak = 0.0f;
        size_t const num_samples = frame_count * 2;
        for (size_t i = 0; i < num_samples; ++i) {
            peak = std::max(peak, std::abs(buffer[i]));
        }

        // 2. gain computer
        float const level_db = 20.0f * std::log10(std::max(peak, 1.0e-6f));
        float const target_gain = qwqdsp::convert::Db2Gain(ComputeTargetGainDb(level_db));

        // compressor 的 attack 是增益下降；gate/expander 的 attack 是增益上升
        bool const attack = (mode_ == DynamicsMode::Compressor)
                                ? target_gain < current_gain_
                                : target_gain > current_gain_;
        float const coeff = attack ? attack_coeff_ : release_coeff_;

        // 3. 逐样本平滑增益并施加
        for (size_t i = 0; i < frame_count; ++i) {
            current_gain_ = target_gain + coeff * (current_gain_ - target_gain);

            float const gain = current_gain_ * makeup_smoother_.Tick();
            buffer[2 * i] *= gain;
            buffer[2 * i + 1] *= gain;
        }

        block_gain_reduction_db_ = 20.0f * std::log10(std::max(current_gain_, 1.0e-6f));
    }

    float GetGainReductionDb() const noexcept {
        return block_gain_reduction_db_;
    }

private:
    static constexpr float kGateHysteresisDb = 3.0f;

    float ComputeTargetGainDb(float level_db) noexcept {
        switch (mode_) {
        case DynamicsMode::Compressor: {
            float const over_db = level_db - threshold_db_;
            if (over_db > 0.0f) {
                // ratio:1 -> 超出 threshold 的部分按 1 - 1/ratio 衰减
                return over_db * (1.0f / ratio_ - 1.0f);
            }
            return 0.0f;
        }
        case DynamicsMode::Expander: {
            float const under_db = threshold_db_ - level_db;
            if (under_db > 0.0f) {
                // 低于 threshold 后按 ratio:1 向下扩展
                return -std::min(under_db * (ratio_ - 1.0f), range_db_);
            }
            return 0.0f;
        }
        case DynamicsMode::Gate: {
            if (gate_open_) {
                if (level_db < threshold_db_ - kGateHysteresisDb) {
                    gate_open_ = false;
                }
            }
            else {
                if (level_db > threshold_db_) {
                    gate_open_ = true;
                }
            }
            return gate_open_ ? 0.0f : -range_db_;
        }
        default:
            return 0.0f;
        }
    }

    float sample_rate_{};
    DynamicsMode mode_{DynamicsMode::Compressor};
    float threshold_db_{-24.0f};
    float ratio_{4.0f};
    float range_db_{60.0f};
    float attack_coeff_{};
    float release_coeff_{};
    float current_gain_{1.0f};
    float block_gain_reduction_db_{};
    bool gate_open_{true};
    qwqdsp_misc::ExpSmoother makeup_smoother_;
};

} // namespace dynamic_demo
