#pragma once
#include "burg_lp.hpp"

#include <algorithm>
#include <cassert>
#include <complex>
#include <numbers>
#include <span>
#include <vector>

namespace qwqdsp_adaptive {
/**
 * @brief 全通扭曲的 Burg LPC（warped Burg）
 *
 * 与 @ref BurgLP 的唯一区别：每级递推里先把**后向误差**过一遍一阶全通
 *
 *     A(z) = (a + z^-1) / (1 + a * z^-1)
 *
 * 再做相关与更新，等价于在**扭曲后的频率轴**上做 Burg。
 * `a > 0` 把频率轴朝低频压缩（低频段被拉长、分辨率提高），`a < 0` 反之；
 * `a = 0` 时全通退化成 1 拍延迟，等价于普通 Burg（边界处理略有差别，见 @note）。
 *
 * 得到的全极点模型仍是 `1/A(z)`（A(z) 由反射系数经 `Lattice2Tf_KeepK` 得到），
 * 但**算频响时要把 z^-1 换成 A(z)**：
 *
 *     H(ω) = 1 / A(A(e^{-jw}))
 *
 * 否则画出来的包络不是这条扭曲轴上的模型。
 *
 * @note 与 @ref BurgLP 的边界处理不同：本类每级都对整帧递推（和上游块处理实现一致），
 *       `BurgLP` 每级会把长度缩短 1 个采样。a=0 时低阶系数一致，高阶会有 ~1e-2 量级差异。
 * @ref 移植自 green_vocoder 的 `BlockBurgLPC` 的**分析部分**（那块代码还含声码器用的
 *      平滑/增益/格型合成，不属于本算法）
 */
class WarpedBurgLP {
public:
    void Init(size_t block_len) {
        eb_.resize(block_len);
    }

    /**
     * @brief 设置扭曲用的一阶全通系数
     * @param allpass_coeff a ∈ (-1, 1)；内部会夹到 ±0.999，避免极点贴单位圆
     */
    void SetWarp(float allpass_coeff) noexcept {
        warp_ = std::clamp(allpass_coeff, -0.999f, 0.999f);
    }

    /**
     * @brief 设置岭回归常数
     * @param ridge k = -2*up / (down + ridge)，防止 down≈0 时病态；0 表示关闭
     * @note 上游块处理实现用 1e-4
     */
    void SetRidge(float ridge) noexcept {
        ridge_ = std::max(ridge, 0.0f);
    }

    /**
     * @brief 逐级扭曲 Burg 递推，求出各阶反射系数
     * @param x 输入信号。**会被就地修改**（递推的 upgoing 结果直接写在 x 上），
     *          需要保留原信号时调用方先自行拷贝
     * @param latticek 输出的反射系数，长度 = 想要的最大阶数
     * @note 需要 `x.size() >= latticek.size()`；内部缓冲由 `Init(block_len)` 分配，
     *       应保证 `block_len >= x.size()`
     */
    void Process(std::span<float> x, std::span<float> latticek) noexcept {
        assert(eb_.size() >= x.size());
        assert(x.size() >= latticek.size());

        std::copy(x.begin(), x.end(), eb_.begin());
        for (auto& k : latticek) {
            float up{};
            float down{};
            float allpass_state{};
            for (size_t i = 0; i < x.size(); ++i) {
                // 一阶全通 A(z) = (a + z^-1) / (1 + a*z^-1)
                float const y = warp_ * eb_[i] + allpass_state;
                allpass_state = eb_[i] - warp_ * y;
                eb_[i] = y;

                up += x[i] * y;
                down += x[i] * x[i];
                down += y * y;
            }
            k = -2.0f * up / (down + ridge_);

            for (size_t i = 0; i < x.size(); ++i) {
                float const upgo = x[i] + eb_[i] * k;
                float const downgo = eb_[i] + x[i] * k;
                x[i] = upgo;
                eb_[i] = downgo;
            }
        }
    }

    /**
     * @brief 扭曲轴上某频点的全通延迟算子 A(e^{-jw})
     * @param freq_hz 频率
     * @param sample_rate 采样率
     * @return A(z) = (a + z^-1) / (1 + a*z^-1) 在 z = e^{jw} 处的值
     * @note 由它把 A(z) 的 z^-1 替换掉即可得到扭曲后的频响
     */
    [[nodiscard]] std::complex<float> WarpDelay(float freq_hz, float sample_rate) const noexcept {
        auto const zinv = std::polar(1.0f, -2.0f * std::numbers::pi_v<float> * freq_hz / sample_rate);
        return (warp_ + zinv) / (1.0f + warp_ * zinv);
    }

private:
    std::vector<float> eb_;
    float warp_{0.0f};
    float ridge_{0.0f};
};
} // namespace qwqdsp_adaptive
