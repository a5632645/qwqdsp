#pragma once
#include "biquad_coeff.hpp"
#include <cmath>
#include <numbers>

namespace qwqdsp_filter {
/**
 * @ref https://www.w3.org/TR/audio-eq-cookbook/
 */
struct RBJ {
    float b0;
    float b1;
    float b2;
    float a1;
    float a2;

    [[nodiscard]]
    BiquadCoeff ToBiquadCoeff() const noexcept {
        return BiquadCoeff{b0, b1, b2, a1, a2};
    }

    /**
     * @note 使用数字倍频程表示带宽具有一阶prewarp Q
     */
    static float DigitalOctave2AnalogQ(float w, float octave) noexcept {
        auto a = std::numbers::ln2_v<float> * 0.5f * octave * w / std::sin(w);
        return 0.5f / std::sinh(a);
    }

    static float DigitalBW2AnalogQ(float w, float bw) noexcept {
        auto f0 = w - bw * 0.5f;
        auto f1 = w + bw * 0.5f;
        auto octave = f1 / f0;
        return DigitalOctave2AnalogQ(w, octave);
    }

    void Lowpass(float w, float Q) noexcept {
        auto a = std::sin(w) / (2 * Q);
        auto cosw = std::cos(w);
        b0 = (1 - cosw) / 2.0f;
        b1 = 1 - cosw;
        b2 = b0;
        a1 = -2 * cosw;
        a2 = 1 - a;
        float inva0 = 1.0f / (1 + a);
        b0 *= inva0;
        b1 *= inva0;
        b2 *= inva0;
        a1 *= inva0;
        a2 *= inva0;
    }

    [[nodiscard]]
    BiquadCoeff Dicimate(size_t dicimate) noexcept {
        if (dicimate == 1) {
            return kBiquadPassthrough;
        }
        else {
            Lowpass(std::numbers::pi_v<float> / static_cast<float>(dicimate), std::numbers::sqrt2_v<float> / 2);
            return ToBiquadCoeff();
        }
    }

    void Highpass(float w, float Q) noexcept {
        auto a = std::sin(w) / (2 * Q);
        auto cosw = std::cos(w);
        b0 = (1 + cosw) / 2.0f;
        b1 = -(1 + cosw);
        b2 = b0;
        a1 = -2 * cosw;
        a2 = 1 - a;
        float inva0 = 1.0f / (1 + a);
        b0 *= inva0;
        b1 *= inva0;
        b2 *= inva0;
        a1 *= inva0;
        a2 *= inva0;
    }

    /**
     * |H(z=exp(jw))| = Q
     */
    void Bandpass(float w, float Q) noexcept {
        auto a = std::sin(w) / (2 * Q);
        auto cosw = std::cos(w);
        b0 = Q * a;
        b1 = 0;
        b2 = -Q * a;
        a1 = -2 * cosw;
        a2 = 1 - a;
        float inva0 = 1.0f / (1 + a);
        b0 *= inva0;
        b1 *= inva0;
        b2 *= inva0;
        a1 *= inva0;
        a2 *= inva0;
    }

    /**
     * |H(z=exp(jw))| = 1
     */
    void BandpassKeep0(float w, float Q) noexcept {
        auto a = std::sin(w) / (2 * Q);
        auto cosw = std::cos(w);
        b0 = a;
        b1 = 0;
        b2 = -a;
        a1 = -2 * cosw;
        a2 = 1 - a;
        float inva0 = 1.0f / (1 + a);
        b0 *= inva0;
        b1 *= inva0;
        b2 *= inva0;
        a1 *= inva0;
        a2 *= inva0;
    }

    void BandpassKeep0Precise(float w1, float w2) noexcept {
        float f1 = std::tan(w1 / 2);
        float f2 = std::tan(w2 / 2);
        float f0 = std::sqrt(f1 * f2);
        float Q = f0 / std::abs(f1 - f2);
        float w = 2 * std::atan(f0);
        BandpassKeep0(w, Q);
    }

    void Peak(float w, float Q, float g) noexcept {
        auto a = std::sin(w) / (2 * Q);
        auto cosw = std::cos(w);
        auto A = std::pow(10.0f, g / 40.0f);
        b0 = 1 + a * A;
        b1 = -2 * cosw;
        b2 = 1 - a * A;
        a1 = -2 * cosw;
        a2 = 1 - a / A;
        float inva0 = 1.0f / (1 + a / A);
        b0 *= inva0;
        b1 *= inva0;
        b2 *= inva0;
        a1 *= inva0;
        a2 *= inva0;
    }

    void Lowshelf(float w, float Q, float g) noexcept {
        auto a = std::sin(w) / (2 * Q);
        auto cosw = std::cos(w);
        auto A = std::pow(10.0f, g / 40.0f);
        auto sqrtA = std::pow(10.0f, g / 80.0f);
        b0 = A * ((A + 1) - (A - 1) * cosw + 2 * sqrtA * a);
        b1 = 2 * A * ((A - 1) - (A + 1) * cosw);
        b2 = A * ((A + 1) - (A - 1) * cosw - 2 * sqrtA * a);
        a1 = -2 * ((A - 1) + (A + 1) * cosw);
        a2 = (A + 1) + (A - 1) * cosw - 2 * sqrtA * a;
        float inva0 = 1.0f / ((A + 1) + (A - 1) * cosw + 2 * sqrtA * a);
        b0 *= inva0;
        b1 *= inva0;
        b2 *= inva0;
        a1 *= inva0;
        a2 *= inva0;
    }

    void HighShelf(float w, float Q, float g) noexcept {
        auto a = std::sin(w) / (2 * Q);
        auto cosw = std::cos(w);
        auto A = std::pow(10.0f, g / 40.0f);
        auto sqrtA = std::pow(10.0f, g / 80.0f);
        b0 = A * ((A + 1) + (A - 1) * cosw + 2 * sqrtA * a);
        b1 = -2 * A * ((A - 1) + (A + 1) * cosw);
        b2 = A * ((A + 1) + (A - 1) * cosw - 2 * sqrtA * a);
        a1 = 2 * ((A - 1) - (A + 1) * cosw);
        a2 = (A + 1) - (A - 1) * cosw - 2 * sqrtA * a;
        float inva0 = 1.0f / ((A + 1) - (A - 1) * cosw + 2 * sqrtA * a);
        b0 *= inva0;
        b1 *= inva0;
        b2 *= inva0;
        a1 *= inva0;
        a2 *= inva0;
    }

    void Notch(float w, float Q) noexcept {
        auto a = std::sin(w) / (2 * Q);
        auto cosw = std::cos(w);
        b0 = 1;
        b1 = -2 * cosw;
        b2 = 1;
        a1 = -2 * cosw;
        a2 = 1 - a;
        float inva0 = 1.0f / (1 + a);
        b0 *= inva0;
        b1 *= inva0;
        b2 *= inva0;
        a1 *= inva0;
        a2 *= inva0;
    }

    void Allpass(float w, float Q) noexcept {
        auto a = std::sin(w) / (2 * Q);
        auto cosw = std::cos(w);
        b0 = 1 - a;
        b1 = -2 * cosw;
        b2 = 1 + a;
        a1 = -2 * cosw;
        a2 = 1 - a;
        float inva0 = 1.0f / (1 + a);
        b0 *= inva0;
        b1 *= inva0;
        b2 *= inva0;
        a1 *= inva0;
        a2 *= inva0;
    }

    /**
     * @brief 倾斜架: 直流 -g/2 dB, 奈奎斯特 +g/2 dB, 中心频率处 0dB
     * @param w 数字中心角频率 (rad/sample)
     * @param Q 品质因子
     * @param g 两端到中心的增益落差 (dB), 正值表示抬高频压低频
     * @note 模拟原型 Tiltshelf(s) 与 Highshelf(s) 只差一个常数因子:
     *       Highshelf(s) = sqrt_A^2 * Tiltshelf(s), sqrt_A = 10^(g/80)
     *       双线性变换是线性的, 所以数字域同样有 H_ts(z) = H_hs(z) / sqrt_A^2,
     *       即只要把 HighShelf 的**分子**除以 sqrt_A^2, 分母不变。
     *       (已用符号推导 + scipy.signal.bilinear 交叉验证, 详见推导注释)
     */
    void Tiltshelf(float w, float Q, float g) noexcept {
        HighShelf(w, Q, g);
        float const A = std::pow(10.0f, g / 40.0f);
        float const inv_A = 1.0f / A;
        b0 *= inv_A;
        b1 *= inv_A;
        b2 *= inv_A;
    }

    // ------------------------------------------------------------
    // 单极点 (一阶): Q 对一阶没有意义, 所以不接受 Q 参数
    // ------------------------------------------------------------
    // 推导: 令 K = tan(w/2) 为预畸变量, 代入双线性变换
    //       s = (1 - z^-1) / (K * (1 + z^-1))
    // 到 analog_responce.hpp 的模拟原型里, 化简即得下面各式。
    // 增益参数 g_db 一律按 10^(g_db/40) 转成 sqrt_A, 与模拟原型的约定一致
    // (HighshelfOnepole/LowshelfOnepole/TiltshelfOnepole 的 @param sqrt_A 即此值)。
    // 等价关系: 直流与奈奎斯特之间的落差是 g_db, 中心频率处增益为 g_db/2。

    /**
     * @brief 一阶低通: 直流 0dB, 奈奎斯特 -inf
     * @param w 数字截止角频率 (rad/sample)
     * @note 一3dB 点落在 w 处
     */
    void LowpassOnepole(float w) noexcept {
        float const k = std::tan(w * 0.5f);
        float const inv = 1.0f / (k + 1);
        float const kn = k * inv;
        b0 = kn;
        b1 = kn;
        b2 = 0;
        a1 = (k - 1) * inv;
        a2 = 0;
    }

    /**
     * @brief 一阶高通: 奈奎斯特 0dB, 直流 -inf
     * @param w 数字截止角频率 (rad/sample)
     * @note 一3dB 点落在 w 处
     */
    void HighpassOnepole(float w) noexcept {
        float const k = std::tan(w * 0.5f);
        float const inv = 1.0f / (k + 1);
        b0 = inv;
        b1 = -inv;
        b2 = 0;
        a1 = (k - 1) * inv;
        a2 = 0;
    }

    /**
     * @brief 一阶全通: 幅度恒为 1, 相位在 w 处为 -90 度
     * @param w 数字转折角频率 (rad/sample)
     */
    void AllpassOnepole(float w) noexcept {
        float const k = std::tan(w * 0.5f);
        float const inv = 1.0f / (k + 1);
        float const a = (k - 1) * inv;
        b0 = a;
        b1 = 1;
        b2 = 0;
        a1 = a;
        a2 = 0;
    }

    /**
     * @brief 一阶低架: 直流增益 = g, 奈奎斯特 0dB, 中心处 g/2
     * @param w 数字中心角频率 (rad/sample)
     * @param g 直流处相对奈奎斯特的增益 (dB)
     */
    void LowshelfOnepole(float w, float g) noexcept {
        float const k = std::tan(w * 0.5f);
        float const sqrtA = std::pow(10.0f, g / 40.0f);
        float const inv = 1.0f / (k + sqrtA);
        b0 = sqrtA * (k * sqrtA + 1) * inv;
        b1 = sqrtA * (k * sqrtA - 1) * inv;
        b2 = 0;
        a1 = (k - sqrtA) * inv;
        a2 = 0;
    }

    /**
     * @brief 一阶高架: 奈奎斯特处增益 = g, 直流 0dB, 中心处 g/2
     * @param w 数字中心角频率 (rad/sample)
     * @param g 奈奎斯特处相对直流的增益 (dB)
     */
    void HighshelfOnepole(float w, float g) noexcept {
        float const k = std::tan(w * 0.5f);
        float const sqrtA = std::pow(10.0f, g / 40.0f);
        float const inv = 1.0f / (k * sqrtA + 1);
        b0 = sqrtA * (k + sqrtA) * inv;
        b1 = sqrtA * (k - sqrtA) * inv;
        b2 = 0;
        a1 = (k * sqrtA - 1) * inv;
        a2 = 0;
    }

    /**
     * @brief 一阶倾斜架: 直流 -g/2 dB, 奈奎斯特 +g/2 dB, 中心 0dB
     * @param w 数字中心角频率 (rad/sample)
     * @param g 两端到中心的增益落差 (dB), 正值表示抬高频压低频
     */
    void TiltshelfOnepole(float w, float g) noexcept {
        float const k = std::tan(w * 0.5f);
        float const sqrtA = std::pow(10.0f, g / 40.0f);
        float const inv = 1.0f / (k * sqrtA + 1);
        b0 = (k + sqrtA) * inv;
        b1 = (k - sqrtA) * inv;
        b2 = 0;
        a1 = (k * sqrtA - 1) * inv;
        a2 = 0;
    }
};
} // namespace qwqdsp_filter
