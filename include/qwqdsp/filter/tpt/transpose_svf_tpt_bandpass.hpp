#pragma once
#include <cmath>
#include <complex>

namespace qwqdsp_filter {
/// @brief 转置 TPT LP2BP 带通：积分器替换为带通积分器后的转置 SVF 结构
/// @note 低通原型为 a*s^2 + b*s + c，极点共轭成对，零点必须位于虚轴上
class TransposeSvfTPTBandpass {
public:
    /// @param[in] pole 左半平面极点，与共轭成对
    /// @param[in] zero 虚轴零点，与共轭成对，必须为纯虚数
    /// @param[in] k 原型增益
    TransposeSvfTPTBandpass(std::complex<double> pole, std::complex<double> zero, double k) noexcept {
        float const wc = static_cast<float>(std::abs(pole));
        prototype_.wc = wc;
        prototype_.r2 = static_cast<float>(-2.0 * pole.real() / wc);
        float const wz = static_cast<float>(std::abs(zero));
        if (wz > 0.0f) {
            prototype_.hp_mix = static_cast<float>(k);
            prototype_.bp_mix = 0.0f;
            prototype_.lp_mix = static_cast<float>(k) * (wz * wz) / (wc * wc);
        }
        else {
            prototype_.hp_mix = 0.0f;
            prototype_.bp_mix = static_cast<float>(k) / wc;
            prototype_.lp_mix = 0.0f;
        }
    }

    void Reset() noexcept {
        sb1_ = 0;
        sl1_ = 0;
        sb2_ = 0;
        sl2_ = 0;
    }

    /// @brief 设置数字角频率边沿
    /// @param[in] w1 下边沿，范围 (0, w2)，单位 rad/sample
    /// @param[in] w2 上边沿，范围 (w1, pi)，单位 rad/sample
    /// @note 内部按 tan(w/2) 预扭曲，再解 LP->BP 积分器替换后的 ZDF 反馈
    void Set(float w1, float w2) noexcept {
        float const aw1 = std::tan(w1 / 2);
        float const aw2 = std::tan(w2 / 2);
        float const bw = aw2 - aw1;
        float const w0 = std::sqrt(aw1 * aw2);

        g_ = w0;
        s_ = prototype_.wc * bw / w0;
        r2_ = prototype_.r2;
        r2s_ = r2_ * s_;
        // 2x2 ZDF 系统: [at bt; -bt dt] [b1; b2] = [u; v]
        at_ = 1 + g_ * g_;
        bt_ = s_ * g_;
        dt_ = 1 + g_ * g_ + r2s_ * g_;
        inv_d_ = 1 / (at_ * dt_ + bt_ * bt_);
    }

    /// @brief 处理一个采样，三路输入共用同一混合信号
    float Tick(float x) noexcept {
        return Tick(prototype_.hp_mix * x, prototype_.bp_mix * x, prototype_.lp_mix * x);
    }

    /// @brief 处理一个采样，三路输入独立馈入
    float Tick(float x_hp, float x_bp, float x_lp) noexcept {
        float const u = g_ * (x_lp - x_hp) - g_ * sl1_ + sb1_;
        float const v = g_ * (x_bp - r2_ * x_hp) - g_ * sl2_ + sb2_;
        float const b1 = (dt_ * u - bt_ * v) * inv_d_;
        float const b2 = (at_ * v + bt_ * u) * inv_d_;
        float const l1 = g_ * b1 + sl1_;
        float const l2 = g_ * b2 + sl2_;

        sb1_ = 2 * b1 - sb1_;
        sl1_ = 2 * l1 - sl1_;
        sb2_ = 2 * b2 - sb2_;
        sl2_ = 2 * l2 - sl2_;

        return x_hp + s_ * b2;
    }
private:
    struct Prototype {
        float wc;
        float r2;
        float hp_mix;
        float bp_mix;
        float lp_mix;
    };

    Prototype prototype_{};
    float g_{};
    float s_{};
    float r2_{};
    float r2s_{};
    float at_{};
    float bt_{};
    float dt_{};
    float inv_d_{};
    // 两个带通积分器各用零阻尼 SVF 实现：b 为 BP 输出，l 为 LP 输出
    float sb1_{};
    float sl1_{};
    float sb2_{};
    float sl2_{};
};
} // namespace qwqdsp_filter
