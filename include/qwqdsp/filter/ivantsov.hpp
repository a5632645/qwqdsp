#pragma once
#include "biquad_coeff.hpp"
#include <cmath>
#include <numbers>

namespace qwqdsp_filter {
/**
 * @brief Ivantsov 的 "理想双线性" (decramped) 双二阶设计器
 *
 * @ref Yuriy Ivantsov, "On the ideal bilinear and biquadratic digital filter"
 *      https://ivantsovy.com/research/paper1.pdf
 *
 * 与 RBJ(双线性变换, 走预畸变轴)的区别: 本设计器贴合**直接轴**上的模拟原型,
 * 也就是让 |H_d(e^jw)| 尽量等于 |H_a(j*w)|, 所以高频段比 RBJ 更接近原型。
 * 代价是它并非精确的代数映射, 而是在若干点上匹配原型的拟合:
 *   - fc 处精确匹配(预畸变后);
 *   - 越靠近奈奎斯特, 等效的匹配频率越会偏离直接轴 —— 大约 f=20kHz(fs=48k)
 *     处等效频率约为 1.1*f。这是算法固有的性质, 不是实现误差。
 *
 * 所有函数的 Q 与 notebook 里的 zeta 是同一个阻尼参数: zeta = 1/(2Q)。
 *
 * sigma 是算法的形状参数(无量纲), 控制去畸变的"力度": 它不改变 fc 处的匹配,
 * 但会改变从直流到奈奎斯特之间的频率轴弯曲程度。论文与 notebook 里常用
 * 2.0(notebook 默认)、2.5、pi*sqrt(2/3)(约 2.5651); 不同 fc 下最优值不同,
 * 所以做成参数而不是写死常量。
 *
 * @note 参考的 notebook: qwqdsp/notebooks/ivantsovy.ipynb
 */
class Ivantsov {
public:
    /// sigma 的默认值, 取自 notebook 的默认参数
    static constexpr float kDefaultSigma = 2.0f;

    // ------------------------------------------------------------
    // 一阶
    // ------------------------------------------------------------

    /**
     * @brief 一阶低通
     * @param wc 数字截止角频率 (rad/sample)
     * @param sigma 算法形状参数
     */
    BiquadCoeff LowpassOnepole(float wc, float sigma = kDefaultSigma) noexcept {
        double const omega = PrewarpedOmega(wc, sigma);
        double const alpha = Phi(0.0, sigma);
        double const beta = Phi(omega, sigma);
        double const gain = 1.0 / (1.0 + alpha);
        return AssembleOnepole(gain, alpha, beta);
    }

    /**
     * @brief 一阶高通
     * @param wc 数字截止角频率 (rad/sample)
     * @param sigma 算法形状参数
     */
    BiquadCoeff HighpassOnepole(float wc, float sigma = kDefaultSigma) noexcept {
        double const omega = PrewarpedOmega(wc, sigma);
        double const alpha = kPhiAtInfinity;
        double const beta = Phi(omega, sigma);
        double const gain = omega / (2.0 * kPi);
        return AssembleOnepole(gain, alpha, beta);
    }

    /**
     * @brief 一阶全通
     * @param wc 数字转折角频率 (rad/sample)
     * @param sigma 算法形状参数
     */
    BiquadCoeff AllpassOnepole(float wc, float sigma = kDefaultSigma) noexcept {
        double const omega = PrewarpedOmega(wc, sigma);
        double const beta = Phi(omega, sigma);
        double const alpha = (beta != 0.0) ? 1.0 / beta : 0.0;
        double const gain = 1.0 / (1.0 + alpha);
        return AssembleOnepole(gain, alpha, beta);
    }

    /**
     * @brief 一阶低架
     * @param wc 数字转折角频率 (rad/sample)
     * @param db 直流处相对奈奎斯特的增益 (dB)
     * @param sigma 算法形状参数
     */
    BiquadCoeff LowshelfOnepole(float wc, float db, float sigma = kDefaultSigma) noexcept {
        double const omega = PrewarpedOmega(wc, sigma);
        double const g = std::pow(10.0, static_cast<double>(db) / 20.0);
        double const sqrt_g = std::sqrt(g);
        double const alpha = Phi(omega / sqrt_g, sigma);
        double const beta = Phi(omega * sqrt_g, sigma);
        double const gain = g / (1.0 + alpha);
        return AssembleOnepole(gain, alpha, beta);
    }

    /**
     * @brief 一阶高架
     * @param wc 数字转折角频率 (rad/sample)
     * @param db 奈奎斯特处相对直流的增益 (dB)
     * @param sigma 算法形状参数
     */
    BiquadCoeff HighshelfOnepole(float wc, float db, float sigma = kDefaultSigma) noexcept {
        double const omega = PrewarpedOmega(wc, sigma);
        double const g = std::pow(10.0, static_cast<double>(db) / 20.0);
        double const sqrt_g = std::sqrt(g);
        double const alpha = Phi(omega * sqrt_g, sigma);
        double const beta = Phi(omega / sqrt_g, sigma);
        double const gain = 1.0 / (1.0 + alpha);
        return AssembleOnepole(gain, alpha, beta);
    }

    /**
     * @brief 一阶倾斜架: 直流 -db/2, 奈奎斯特 +db/2, 转折处 0dB
     * @param wc 数字转折角频率 (rad/sample)
     * @param db 两端到中心的落差 (dB)
     * @param sigma 算法形状参数
     * @note 模拟原型上 HighshelfOnepole = sqrt_A * TiltshelfOnepole
     *       (sqrt_A = 10^(db/40), 见 analog_responce.hpp);
     *       滤波器的分子乘常数等于整个 |H| 乘该常数, 所以直接把分子除以 sqrt_A。
     *       缩放后与模拟原型的偏差与 HighshelfOnepole 逐点相同。
     */
    BiquadCoeff TiltshelfOnepole(float wc, float db, float sigma = kDefaultSigma) noexcept {
        double const omega = PrewarpedOmega(wc, sigma);
        double const g = std::pow(10.0, static_cast<double>(db) / 20.0);
        double const sqrt_g = std::sqrt(g);
        double const alpha = Phi(omega * sqrt_g, sigma);
        double const beta = Phi(omega / sqrt_g, sigma);
        double const gain = 1.0 / ((1.0 + alpha) * sqrt_g);
        return AssembleOnepole(gain, alpha, beta);
    }

    // ------------------------------------------------------------
    // 二阶
    // ------------------------------------------------------------

    /**
     * @brief 二阶低通
     * @param wc 数字截止角频率 (rad/sample)
     * @param Q 品质因子 (notebook 里的 zeta = 1/(2Q))
     * @param sigma 算法形状参数
     */
    BiquadCoeff Lowpass(float wc, float Q, float sigma = kDefaultSigma) noexcept {
        double const omega = PrewarpedOmega(wc, sigma);
        double const alpha1 = 2.0 * Phi(0.0, sigma);
        double const alpha2 = Phi(0.0, sigma) * Phi(0.0, sigma);
        double const gain = 1.0 / (1.0 + alpha1 + alpha2);
        return AssembleSecondOrder(gain, alpha1, alpha2, wc, Q, sigma);
    }

    /**
     * @brief 二阶高通
     * @param wc 数字截止角频率 (rad/sample)
     * @param Q 品质因子 (notebook 里的 zeta = 1/(2Q))
     * @param sigma 算法形状参数
     */
    BiquadCoeff Highpass(float wc, float Q, float sigma = kDefaultSigma) noexcept {
        double const omega = PrewarpedOmega(wc, sigma);
        double const alpha1 = 2.0 * kPhiAtInfinity;
        double const alpha2 = kPhiAtInfinity * kPhiAtInfinity;
        double const gain = omega * omega / (4.0 * kPi * kPi);
        return AssembleSecondOrder(gain, alpha1, alpha2, wc, Q, sigma);
    }

    /**
     * @brief 二阶带通
     * @param wc 数字中心角频率 (rad/sample)
     * @param Q 品质因子 (notebook 里的 zeta = 1/(2Q))
     * @param sigma 算法形状参数
     * @note 峰值增益恒为 1 (0dB), 与 Q 无关 —— 对应 RBJ::BandpassKeep0 / MatchBiquad::NormBandpass
     *       的归一化约定, 而不是 RBJ::Bandpass(峰值增益 = Q)那一套
     */
    BiquadCoeff Bandpass(float wc, float Q, float sigma = kDefaultSigma) noexcept {
        double const omega = PrewarpedOmega(wc, sigma);
        double const zeta = ZetaOf(Q);
        double const phi0 = Phi(0.0, sigma);
        double const alpha1 = phi0 + kPhiAtInfinity;
        double const alpha2 = phi0 * kPhiAtInfinity;
        double const gain = (omega / (2.0 * kPi)) * (2.0 * zeta / (2.0 + alpha1));
        return AssembleSecondOrder(gain, alpha1, alpha2, wc, Q, sigma);
    }

    /**
     * @brief 二阶带通, 峰值增益 = Q (与 RBJ::Bandpass / MatchBiquad::Bandpass 同约定)
     * @param wc 数字中心角频率 (rad/sample)
     * @param Q 品质因子 (notebook 里的 zeta = 1/(2Q))
     * @param sigma 算法形状参数
     * @note 模拟原型上 Bandpass = Q * NormBandpass (见 analog_responce.hpp),
     *       所以把归一化带通的分子乘以 Q 即可; 缩放后与模拟原型的偏差逐点不变。
     */
    BiquadCoeff BandpassQ(float wc, float Q, float sigma = kDefaultSigma) noexcept {
        double const omega = PrewarpedOmega(wc, sigma);
        double const zeta = ZetaOf(Q);
        double const phi0 = Phi(0.0, sigma);
        double const alpha1 = phi0 + kPhiAtInfinity;
        double const alpha2 = phi0 * kPhiAtInfinity;
        double const gain = (omega / (2.0 * kPi)) * (2.0 * zeta / (2.0 + alpha1)) * static_cast<double>(Q);
        return AssembleSecondOrder(gain, alpha1, alpha2, wc, Q, sigma);
    }

    /**
     * @brief 二阶陷波
     * @param wc 数字陷波中心角频率 (rad/sample)
     * @param Q 品质因子 (notebook 里的 zeta = 1/(2Q))
     * @param sigma 算法形状参数
     * @note 用论文式(3.9)的简化分子: 零点精确落在单位圆上
     */
    BiquadCoeff Notch(float wc, float Q, float sigma = kDefaultSigma) noexcept {
        double const omega = PrewarpedOmega(wc, sigma);
        double const num = omega * omega - kPi * kPi - SigmaSquare(sigma);
        double const den = omega * omega + kPi * kPi - SigmaSquare(sigma);
        double const alpha1 = -2.0 * (num / den);
        double const alpha2 = 1.0;
        double const gain = 1.0 / (1.0 + alpha1 + alpha2);
        return AssembleSecondOrder(gain, alpha1, alpha2, wc, Q, sigma);
    }

    /**
     * @brief 二阶全通 (分子恰好是分母的反序, 幅度恒为 1)
     * @param wc 数字转折角频率 (rad/sample)
     * @param Q 品质因子 (notebook 里的 zeta = 1/(2Q))
     * @param sigma 算法形状参数
     */
    BiquadCoeff Allpass(float wc, float Q, float sigma = kDefaultSigma) noexcept {
        double const zeta = ZetaOf(Q);
        double const omega = PrewarpedOmega(wc, sigma);
        double const beta1 = Phi1(omega, zeta, sigma);
        double const beta2 = Phi2(omega, zeta, sigma);
        double const alpha1 = beta1 / beta2;
        double const alpha2 = 1.0 / beta2;
        double const gain = 1.0 / (1.0 + alpha1 + alpha2);
        return AssembleSecondOrder(gain, alpha1, alpha2, wc, Q, sigma);
    }

    /**
     * @brief 二阶峰值 (band shelf): 中心频率处增益 = db, 直流与奈奎斯特处 ≈ 0dB
     * @param wc 数字中心角频率 (rad/sample)
     * @param Q 品质因子 (notebook 里的 zeta = 1/(2Q))
     * @param db 中心频率处增益 (dB)
     * @param sigma 算法形状参数
     */
    BiquadCoeff Peaking(float wc, float Q, float db, float sigma = kDefaultSigma) noexcept {
        double const zeta = ZetaOf(Q);
        double const omega = PrewarpedOmega(wc, sigma);
        double const g = std::pow(10.0, static_cast<double>(db) / 20.0);
        double const sqrt_g = std::sqrt(g);
        // 注意: peaking 是把 zeta 按 sqrt(g) 缩放(增益作用在阻尼上),
        // 而 shelf 是把 omega 按 g^0.25 缩放 —— 两者不是同一套缩放
        double const alpha1 = Phi1(omega, zeta * sqrt_g, sigma);
        double const alpha2 = Phi2(omega, zeta * sqrt_g, sigma);
        double const gain = 1.0 / (1.0 + alpha1 + alpha2);
        double const beta1 = Phi1(omega, zeta / sqrt_g, sigma);
        double const beta2 = Phi2(omega, zeta / sqrt_g, sigma);
        return Assemble(gain, alpha1, alpha2, beta1, beta2);
    }

    /**
     * @brief 二阶低架: 直流处增益 = db, 奈奎斯特处 ≈ 0dB
     * @param wc 数字转折角频率 (rad/sample)
     * @param Q 品质因子 (notebook 里的 zeta = 1/(2Q))
     * @param db 直流处增益 (dB)
     * @param sigma 算法形状参数
     */
    BiquadCoeff Lowshelf(float wc, float Q, float db, float sigma = kDefaultSigma) noexcept {
        double const zeta = ZetaOf(Q);
        double const omega = PrewarpedOmega(wc, sigma);
        double const g = std::pow(10.0, static_cast<double>(db) / 20.0);
        double const sqrt_sqrt_g = std::pow(g, 0.25);
        double const alpha1 = Phi1(omega / sqrt_sqrt_g, zeta, sigma);
        double const alpha2 = Phi2(omega / sqrt_sqrt_g, zeta, sigma);
        double const beta1 = Phi1(omega * sqrt_sqrt_g, zeta, sigma);
        double const beta2 = Phi2(omega * sqrt_sqrt_g, zeta, sigma);
        double const gain = g / (1.0 + alpha1 + alpha2);
        return Assemble(gain, alpha1, alpha2, beta1, beta2);
    }

    /**
     * @brief 二阶高架: 奈奎斯特处增益 = db, 直流处 ≈ 0dB
     * @param wc 数字转折角频率 (rad/sample)
     * @param Q 品质因子 (notebook 里的 zeta = 1/(2Q))
     * @param db 奈奎斯特处增益 (dB)
     * @param sigma 算法形状参数
     */
    BiquadCoeff Highshelf(float wc, float Q, float db, float sigma = kDefaultSigma) noexcept {
        double const zeta = ZetaOf(Q);
        double const omega = PrewarpedOmega(wc, sigma);
        double const g = std::pow(10.0, static_cast<double>(db) / 20.0);
        double const sqrt_sqrt_g = std::pow(g, 0.25);
        double const alpha1 = Phi1(omega * sqrt_sqrt_g, zeta, sigma);
        double const alpha2 = Phi2(omega * sqrt_sqrt_g, zeta, sigma);
        double const gain = 1.0 / (1.0 + alpha1 + alpha2);
        double const beta1 = Phi1(omega / sqrt_sqrt_g, zeta, sigma);
        double const beta2 = Phi2(omega / sqrt_sqrt_g, zeta, sigma);
        return Assemble(gain, alpha1, alpha2, beta1, beta2);
    }

    /**
     * @brief 二阶倾斜架: 直流 -db/2, 奈奎斯特 +db/2, 中心处 0dB
     * @param wc 数字中心角频率 (rad/sample)
     * @param Q 品质因子 (notebook 里的 zeta = 1/(2Q))
     * @param db 两端到中心的落差 (dB), 正值表示抬高频压低频
     * @param sigma 算法形状参数
     * @note 模拟原型上 Highshelf = A * Tiltshelf (A = 10^(db/40), 见 analog_responce.hpp),
     *       且两者极点零点完全相同, 只差这个常数因子; 滤波器的分子乘常数等于整个 |H|
     *       乘该常数, 所以直接把分子除以 A。缩放后与模拟原型的偏差与 Highshelf 逐点相同。
     */
    BiquadCoeff Tiltshelf(float wc, float Q, float db, float sigma = kDefaultSigma) noexcept {
        double const zeta = ZetaOf(Q);
        double const omega = PrewarpedOmega(wc, sigma);
        double const g = std::pow(10.0, static_cast<double>(db) / 20.0);
        double const sqrt_sqrt_g = std::pow(g, 0.25);
        double const alpha1 = Phi1(omega * sqrt_sqrt_g, zeta, sigma);
        double const alpha2 = Phi2(omega * sqrt_sqrt_g, zeta, sigma);
        double const gain = 1.0 / ((1.0 + alpha1 + alpha2) * std::sqrt(g));
        double const beta1 = Phi1(omega / sqrt_sqrt_g, zeta, sigma);
        double const beta2 = Phi2(omega / sqrt_sqrt_g, zeta, sigma);
        return Assemble(gain, alpha1, alpha2, beta1, beta2);
    }
private:
    static constexpr double kPi = std::numbers::pi_v<double>;
    /// 一阶原型在 x -> +inf 时的 phi 值, 即 phi(inf, sigma) 的极限
    static constexpr double kPhiAtInfinity = -1.0;

    /// Q 到 notebook 里的阻尼系数 zeta 的换算
    static inline double ZetaOf(float Q) noexcept {
        return 1.0 / (2.0 * static_cast<double>(Q));
    }

    static inline double SigmaSquare(float sigma) noexcept {
        double const s = sigma;
        return s * s;
    }

    /**
     * @brief 把设计频率预畸变到算法的 omega 域
     * @param wc 数字角频率 (rad/sample)
     * @param sigma 算法形状参数
     * @return omega (无量纲)
     * @note 直接用 wc 表达: 论文里的 pi/tan(pi*fc/fs) 就是 pi/tan(wc/2)
     * @note wc > pi (fc 超过奈奎斯特) 时退化为 2*pi/wc, 与 notebook 保持一致
     */
    static inline double PrewarpedOmega(float wc, float sigma) noexcept {
        double const w = wc;
        if (w > kPi) {
            return 2.0 * kPi / w;
        }
        double const s = sigma;
        double const cot = kPi / std::tan(w / 2.0);
        return std::sqrt(s * s + cot * cot);
    }

    /**
     * @brief 一阶原型的 phi 函数
     * @param x 频率 (omega 域)
     * @param sigma 算法形状参数
     * @return phi(x, sigma) ∈ (-1, 1)
     */
    static inline double Phi(double x, float sigma) noexcept {
        double const s = sigma;
        double const value = std::sqrt(x * x + s * s);
        return (kPi - value) / (kPi + value);
    }

    /// nu(x, y, sigma): 二阶原型的辅助量
    static inline double Nu(double x, double y, float sigma) noexcept {
        double const s = sigma;
        double const x2 = x * x;
        double const s2 = s * s;
        return std::sqrt(x2 * x2 + 2.0 * x2 * s2 * (2.0 * y * y - 1.0) + s2 * s2);
    }

    /// kappa(x, y, sigma): 二阶原型的辅助量
    static inline double Kappa(double x, double y, float sigma) noexcept {
        double const s = sigma;
        return x * x * (2.0 * y * y - 1.0) + s * s;
    }

    /// 二阶原型的 phi1
    static inline double Phi1(double x, double y, float sigma) noexcept {
        double const v = Nu(x, y, sigma);
        double const k = Kappa(x, y, sigma);
        double const root = std::sqrt(v + k);
        double const den = kPi * kPi + kPi * std::numbers::sqrt2_v<double> * root + v;
        return (2.0 * kPi * kPi - 2.0 * v) / den;
    }

    /// 二阶原型的 phi2
    static inline double Phi2(double x, double y, float sigma) noexcept {
        double const v = Nu(x, y, sigma);
        double const k = Kappa(x, y, sigma);
        double const root = std::sqrt(v + k);
        double const den = kPi * kPi + kPi * std::numbers::sqrt2_v<double> * root + v;
        return (kPi * kPi - kPi * std::numbers::sqrt2_v<double> * root + v) / den;
    }

    /**
     * @brief 由分母多项式与增益装配二阶系数
     * @param gain 增益因子 G
     * @param num1 分子一次项 (论文里的 alpha1)
     * @param num2 分子二次项 (论文里的 alpha2)
     * @param wc 数字角频率 (rad/sample)
     * @param Q 品质因子
     * @param sigma 算法形状参数
     * @return b0..a2
     * @note 分母对所有二阶类型都相同, 只由 (omega, zeta, sigma) 决定
     */
    static inline BiquadCoeff AssembleSecondOrder(double gain, double num1, double num2, float wc, float Q,
                                                  float sigma) noexcept {
        double const zeta = ZetaOf(Q);
        double const omega = PrewarpedOmega(wc, sigma);
        double const den1 = Phi1(omega, zeta, sigma);
        double const den2 = Phi2(omega, zeta, sigma);
        return Assemble(gain, num1, num2, den1, den2);
    }

    /**
     * @brief 装配二阶系数: b = G*(1+den1+den2)*[1, num1, num2], a = [1, den1, den2]
     * @param gain 增益因子 G
     * @param num1 分子一次项
     * @param num2 分子二次项
     * @param den1 分母一次项
     * @param den2 分母二次项
     * @return b0..a2
     */
    static inline BiquadCoeff Assemble(double gain, double num1, double num2, double den1, double den2) noexcept {
        double const scale = gain * (1.0 + den1 + den2);
        return DoubleBiquadCoeff{scale, scale * num1, scale * num2, den1, den2}.ToFloat();
    }

    /**
     * @brief 装配一阶系数(存进二阶结构, 高次项为 0)
     * @param gain 增益因子 G
     * @param alpha 分子一次项
     * @param beta 分母一次项
     * @return b0,b1,a1, 其余为 0
     */
    static inline BiquadCoeff AssembleOnepole(double gain, double alpha, double beta) noexcept {
        double const scale = gain * (1.0 + beta);
        return DoubleBiquadCoeff{scale, scale * alpha, 0.0, beta, 0.0}.ToFloat();
    }
};
} // namespace qwqdsp_filter
