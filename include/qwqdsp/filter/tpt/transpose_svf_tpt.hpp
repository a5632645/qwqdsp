#pragma once
#include <cmath>

namespace qwqdsp_filter {
/// @brief 转置 TPT 状态变量滤波器，三路输入可独立混合，默认低通
class TransposeSvfTPT {
public:
    void Reset() noexcept {
        s1_ = 0;
        s2_ = 0;
    }

    /// @brief 设置归一化频率与阻尼
    /// @param[in] w 归一化角频率，单位为 rad/sample，范围 [0, pi]
    /// @param[in] r2 2R = 1/Q，大于等于 2 时极点为实数
    void SetCoeffSVF(float w, float r2) noexcept {
        SetAnalogCoeff(std::tan(w / 2), r2);
    }

    /// @brief 直接从模拟域设置积分器增益与阻尼
    /// @param[in] g 模拟角频率，等于 tan(w/2)
    /// @param[in] r2 2R = 1/Q
    void SetAnalogCoeff(float g, float r2) noexcept {
        g_ = g;
        r2_ = r2;
        g2_ = g_ * g_;
        gr2g2_ = g_ * r2_ + g2_;
        d_ = 1 / (1 + r2_ * g_ + g2_);
        UpdateMix();
    }

    /// @brief 设置归一化频率与品质因数
    /// @param[in] w 归一化角频率，单位为 rad/sample，范围 [0, pi]
    /// @param[in] Q 品质因数，必须大于 0
    void SetCoeffQ(float w, float Q) noexcept {
        SetCoeffSVF(w, 1 / Q);
    }

    /// @brief 设置三路输入混合系数
    /// @param[in] hp 高通输入系数
    /// @param[in] bp 带通输入系数
    /// @param[in] lp 低通输入系数
    /// @note 与普通 SVF 的 hp/bp/lp 输出混合系数同义，转置结构会在内部与积分器增益合并
    void SetInputMix(float hp, float bp, float lp) noexcept {
        hp_mix_ = hp;
        bp_mix_ = bp;
        lp_mix_ = lp;
        UpdateMix();
    }

    /// @brief 处理一个采样，三路输入共用同一混合信号
    /// @param[in] x 输入采样
    /// @return 混合输出
    float Tick(float x) noexcept {
        float const y = (ff_ * x + s1_ + s2_) * d_;
        float const t1 = s1_;
        s1_ = t1 + 2 * (g2_lp_ * x - g2_ * y);
        s2_ += 2 * (g_bp_ * x + g2_lp_ * x - gr2g2_ * y + t1);
        return y;
    }

    /// @brief 处理一个采样，三路输入独立馈入
    /// @param[in] x_hp 高通输入
    /// @param[in] x_bp 带通输入
    /// @param[in] x_lp 低通输入
    /// @return 混合输出
    float Tick(float x_hp, float x_bp, float x_lp) noexcept {
        float const y = (x_hp + g_ * x_bp + g2_ * x_lp + s1_ + s2_) * d_;
        float const t1 = s1_;
        s1_ = t1 + 2 * (g2_ * x_lp - g2_ * y);
        s2_ += 2 * (g_ * x_bp + g2_ * x_lp - gr2g2_ * y + t1);
        return y;
    }
private:
    void UpdateMix() noexcept {
        ff_ = hp_mix_ + g_ * bp_mix_ + g2_ * lp_mix_;
        g_bp_ = g_ * bp_mix_;
        g2_lp_ = g2_ * lp_mix_;
    }

    float g_{};
    float g2_{};
    float r2_{};
    float gr2g2_{};
    float d_{};
    float hp_mix_{};
    float bp_mix_{};
    float lp_mix_{1};
    float ff_{1};
    float g_bp_{};
    float g2_lp_{};
    float s1_{};
    float s2_{};
};
} // namespace qwqdsp_filter
