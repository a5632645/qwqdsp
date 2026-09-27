#pragma once

#include <cmath>

// ------------------------------------------------------------
// RmsTracker
// ------------------------------------------------------------
// 短时 RMS（指数滑动平均）。
//   y[n] = α * x²[n] + (1-α) * y[n-1]
//   value = sqrt(y[n])
// ------------------------------------------------------------

struct RmsTracker {
    /// 配置时间常数
    /// @param tau_sec    积分时间 (秒)，如 0.03 = 30ms
    /// @param sample_rate 采样率 (Hz)
    void setTimeConstant(float tau_sec, float sample_rate) noexcept {
        alpha_ = 1.0f - std::exp(-1.0f / (tau_sec * sample_rate));
    }

    /// 重置
    void reset() noexcept {
        mean_sq_ = 0.0;
    }

    /// 添加一个采样（运行短时 RMS）
    void addSample(float sample) noexcept {
        double sq = static_cast<double>(sample) * static_cast<double>(sample);
        mean_sq_ = static_cast<double>(alpha_) * sq + (1.0 - static_cast<double>(alpha_)) * mean_sq_;
    }

    /// 返回当前短时线性 RMS
    float value() const noexcept {
        return static_cast<float>(std::sqrt(mean_sq_));
    }

    float MeanSq() const noexcept {
        return mean_sq_;
    }
private:
    double mean_sq_ = 0.0;
    double alpha_ = 1.0;
};
