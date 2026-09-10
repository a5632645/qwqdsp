#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <span>
#include <vector>

#include <qwqdsp/filter/window_fir.hpp>
#include <qwqdsp/fx/polyphase_resample_fir.hpp>
#include <qwqdsp/window/kaiser.hpp>

namespace swift_f0_rt {

/**
 * @brief 整数倍抽取器（48kHz → 16kHz）
 *
 * 抗混叠 FIR 用 Kaiser 窗加权的 sinc 设计，直流增益归一化。
 * 内部通过 PolyphaseDownsamplerFir 维护跨块滤波器状态，可逐块流式调用。
 */
class Decimator {
public:
    /**
     * @param factor         抽取倍数
     * @param atten_db       阻带衰减 (dB)，决定 Kaiser beta
     * @param taps_per_phase 每个相位的抽头数（总长度 = factor * taps_per_phase）
     */
    void Init(int factor, float atten_db, int taps_per_phase) noexcept {
        factor_ = factor;
        group_.fill(0.0f);
        pending_ = 0;

        int const num_taps = factor * taps_per_phase;
        delay_samples_ = (num_taps - 1) / 2;
        std::vector<float> coeff(static_cast<size_t>(num_taps));
        // 截止略低于 target Nyquist，为过渡带留余量
        float const cutoff = std::numbers::pi_v<float> / static_cast<float>(factor) * 0.92f;
        qwqdsp_filter::WindowFIR::Lowpass(coeff, cutoff);
        qwqdsp_window::Kaiser::ApplyWindow(coeff, qwqdsp_window::Kaiser::Beta(atten_db), false);
        qwqdsp_filter::WindowFIR::Normalize(coeff);
        fir_.Init(coeff, factor);
    }

    void Reset() noexcept {
        fir_.Reset();
        group_.fill(0.0f);
        pending_ = 0;
    }

    /// @brief FIR 群延迟（输入侧样本数），用于时间对齐
    int64_t DelaySamples() const noexcept {
        return delay_samples_;
    }

    /**
     * @brief 流式抽取，输出样本数 = floor((pending + in.size()) / factor)
     * @note out 会被清空后追加，容量保留，不在内部做动态分配
     */
    void Process(std::span<const float> in, std::vector<float>& out) {
        out.clear();
        size_t i = 0;
        while (i < in.size()) {
            size_t const need = static_cast<size_t>(factor_) - pending_;
            size_t const take = std::min(need, in.size() - i);
            std::copy_n(in.begin() + static_cast<std::ptrdiff_t>(i), take,
                        group_.begin() + static_cast<std::ptrdiff_t>(pending_));
            pending_ += take;
            i += take;

            if (pending_ == static_cast<size_t>(factor_)) {
                out.push_back(fir_.Tick(std::span<float>{group_.data(), static_cast<size_t>(factor_)}));
                pending_ = 0;
            }
        }
    }

private:
    qwqdsp_fx::PolyphaseDownsamplerFir fir_;
    std::array<float, 8> group_{};
    size_t pending_{};
    int factor_{};
    int64_t delay_samples_{};
};

} // namespace swift_f0_rt
