#pragma once
#include <cmath>

namespace qwqdsp_filter {
class OnePoleTPT {
public:
    void Reset() noexcept {
        s_ = 0;
    }

    void MakeLowpass(float w) noexcept {
        auto k = std::tan(w / 2);
        g_ = k / (1 + k);
    }

    void MakePass() noexcept {
        g_ = 1;
    }

    void CopyFrom(OnePoleTPT const& other) noexcept {
        g_ = other.g_;
    }

    /**
     * @note 多模式输出
     *       hp = x - lp
     *       ap = lp - hp or ap = 2*lp - x
     */
    float TickLowpass(float x) noexcept {
        float const delta = g_ * (x - s_);
        s_ += delta;
        float const y = s_;
        s_ += delta;
        return y;
    }

    float TickHighpass(float x) noexcept {
        return x - TickLowpass(x);
    }

    float TickAllpass(float x) noexcept {
        return 2 * TickLowpass(x) - x;
    }
private:
    float g_{};
    float s_{};
};
} // namespace qwqdsp_filter
