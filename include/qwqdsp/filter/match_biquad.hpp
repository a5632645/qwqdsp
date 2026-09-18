#pragma once
#include "analog_responce.hpp"
#include "biquad_coeff.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <numeric>
/**
 * @ref 全通 https://apulsoft.ch/blog/matched-allpass/
 * @ref 复极点shelf https://vicanek.de/articles/2poleShelvingFits.pdf
 * @ref 单极点shelf https://vicanek.de/articles/ShelvingFits.pdf
 * @ref 经典 https://vicanek.de/articles/BiquadFits.pdf
 */
namespace qwqdsp_filter {
class MatchBiquad {
public:
    // -------------------- onepole --------------------
    BiquadCoeff HighshelfOnepole(float wc, float db) noexcept {
        auto G = std::pow(10.0, db / 20.0);
        auto fc = wc / kPi;
        auto alpha = 2.0 / (kPi * kPi) * (1.0 + 1.0 / (G * fc * fc)) - 0.5;
        auto beta = 2.0 / (kPi * kPi) * (1.0 + G / (fc * fc)) - 0.5;
        auto a1 = -alpha / (1.0 + alpha + ClampSqrt(1.0 + 2.0 * alpha));
        auto b = -beta / (1.0 + beta + ClampSqrt(1.0 + 2.0 * beta));
        auto b0 = (1.0 + a1) / (1.0 + b);
        auto b1 = b * b0;
        return DoubleBiquadCoeff{b0, b1, 0, a1, 0}.ToFloat();
    }

    BiquadCoeff TiltshelfOnepole(float wc, float db) noexcept {
        auto sqrt_G = std::pow(10.0f, -db / 40.0f);
        auto r = HighshelfOnepole(wc, db);
        r.b0 *= sqrt_G;
        r.b1 *= sqrt_G;
        r.b2 *= sqrt_G;
        return r;
    }

    BiquadCoeff LowshelfOnepole(float wc, float db) noexcept {
        auto G = std::pow(10.0f, db / 20.0f);
        auto r = HighshelfOnepole(wc, -db);
        r.b0 *= G;
        r.b1 *= G;
        r.b2 *= G;
        return r;
    }

    /**
     * @note match pi/2和双线性变换一样
     */
    template <double kMatchPhase = std::numbers::pi / 4>
    BiquadCoeff AllpassOnepole(float wc) noexcept {
        static auto const cos_a = std::cos(kMatchPhase);
        static auto const f_a = std::sqrt((1 - cos_a) / (1 + cos_a));
        auto w_a = f_a * wc;
        auto c_0 = (std::cos(w_a) - cos_a) / (std::cos(kMatchPhase + w_a) - 1);
        return DoubleBiquadCoeff{c_0, 1, 0, c_0, 0}.ToFloat();
    }

    /**
     * @brief 单极点低通: 直流 0dB, 奈奎斯特处与模拟原型相等
     * @param wc 数字截止角频率 (rad/sample)
     *
     * 与 HighshelfOnepole 用同一套匹配框架 (见类头引的 Vicanek 单极点论文):
     *   |H|^2 = (1 + beta*phi) / (1 + alpha*phi),  phi = 1 - cos(pi*f),  f 以奈奎斯特为单位
     * 三个条件:
     *   (i)   直流增益 = 1              -> b0 + b1 = 1 + a1
     *   (ii)  奈奎斯特处匹配模拟原型      -> (1+2beta)/(1+2alpha) = |H_a(1)|^2
     *   (iii) 低频处匹配到 2 阶          -> beta - alpha = -(2/pi^2) * (d|H_a|^2/d phi)|_0
     * 对 |H_a|^2 = 1/(1+(f/fc)^2) 解得
     *   alpha = 2/pi^2 * (1 + 1/fc^2) - 1/2,   beta = 2/pi^2 - 1/2  (beta 与 fc 无关)
     * 其中 fc = wc/pi 是归一化到奈奎斯特的截止频率, 恰好就是预畸变后的值:
     *   2/pi^2 * (1 + 1/fc^2) = 2/(pi^2 * fc^2) + 2/pi^2 = 2/wc^2 + 2/pi^2
     *
     * @note 与两极点 Lowpass 一样, 这是"尽量贴合"而非精确映射: 误差集中在最高一个
     *       倍频程(如 fc=1kHz 时峰值 0.96dB 在 18kHz 附近), 奈奎斯特处为 0,
     *       1kHz 以下不超过 0.003dB。要压低全频段峰值可把匹配点从奈奎斯特移开
     *       (见论文式 12 对 shelf 的处理), 本实现取奈奎斯特以与 HighshelfOnepole 一致。
     */
    BiquadCoeff LowpassOnepole(float wc) noexcept {
        double const fc = wc / kPi;
        double const alpha = 2.0 / (kPi * kPi) * (1.0 + 1.0 / (fc * fc)) - 0.5;
        double const beta = 2.0 / (kPi * kPi) - 0.5;
        auto a1 = -alpha / (1.0 + alpha + ClampSqrt(1.0 + 2.0 * alpha));
        auto b = -beta / (1.0 + beta + ClampSqrt(1.0 + 2.0 * beta));
        auto b0 = (1.0 + a1) / (1.0 + b);
        auto b1 = b * b0;
        return DoubleBiquadCoeff{b0, b1, 0, a1, 0}.ToFloat();
    }

    /**
     * @brief 单极点高通: 直流 0, 奈奎斯特处与模拟原型相等
     * @param wc 数字截止角频率 (rad/sample)
     *
     * 极点与 LowpassOnepole 相同(转折频率一样), 分子零点落在直流(z=1)。
     * 模拟原型 |H_a|^2 = (f/fc)^2/(1+(f/fc)^2) 与低通是**功率互补**的
     * (|H_lp|^2 + |H_hp|^2 = 1), 于是取
     *   |H_hp|^2 = 1 - |H_lp|^2 = (alpha - beta)*phi / (1 + alpha*phi)
     * 数字实现里零点在直流对应 b1 = -b0, 其分子为 2*phi, 故
     *   2*b0^2/(1+a1)^2 = alpha - beta = 2/wc^2   ->   b0 = (1 + a1)/wc
     * 这样直流精确为 0, 且由于互补关系, 全频段偏差比低通小得多。
     */
    BiquadCoeff HighpassOnepole(float wc) noexcept {
        double const fc = wc / kPi;
        double const alpha = 2.0 / (kPi * kPi) * (1.0 + 1.0 / (fc * fc)) - 0.5;
        auto a1 = -alpha / (1.0 + alpha + ClampSqrt(1.0 + 2.0 * alpha));
        auto b0 = (1.0 + a1) / wc;
        return DoubleBiquadCoeff{b0, -b0, 0, a1, 0}.ToFloat();
    }

    // -------------------- twopole --------------------
    BiquadCoeff Lowpass(float wc, float Q) noexcept {
        auto r = ImpluseInvarant(wc, Q);
        auto phi = GetPhi(wc);
        auto A = GetA(r);
        auto R = Q * Q * std::inner_product(phi.begin(), phi.end(), A.begin(), double{});
        auto B0 = A[0];
        auto B1 = (R - B0 * phi[0]) / phi[1];
        B1 = std::max(B1, 0.0);
        auto sqrt_B0 = 1 + r.a1 + r.a2;
        auto sqrt_B1 = std::sqrt(B1);
        r.b0 = (sqrt_B0 + sqrt_B1) / 2;
        r.b1 = sqrt_B0 - r.b0;
        r.b2 = 0;
        return r.ToFloat();
    }

    BiquadCoeff Highpass(float wc, float Q) noexcept {
        auto r = ImpluseInvarant(wc, Q);
        auto phi = GetPhi(wc);
        auto A = GetA(r);
        auto R = std::inner_product(phi.begin(), phi.end(), A.begin(), double{});
        R = std::max(R, 0.0);
        r.b0 = std::sqrt(R) * Q / (4 * phi[1]);
        r.b1 = -2 * r.b0;
        r.b2 = r.b0;
        return r.ToFloat();
    }

    BiquadCoeff NormBandpass(float wc, float Q) noexcept {
        auto r = ImpluseInvarant(wc, Q);
        auto phi = GetPhi(wc);
        auto A = GetA(r);
        auto R1 = std::inner_product(phi.begin(), phi.end(), A.begin(), double{});
        auto R2 = -A[0] + A[1] + 4 * (phi[0] - phi[1]) * A[2];
        auto B2 = (R1 - R2 * phi[1]) / (4 * phi[1] * phi[1]);
        auto B1 = R2 + 4 * (phi[1] - phi[0]) * B2;
        r.b1 = -std::sqrt(std::max(B1, 0.0)) / 2;
        r.b0 = (std::sqrt(std::max(0.0, B2 + B1 / 4)) - r.b1) / 2;
        r.b2 = -r.b0 - r.b1;
        return r.ToFloat();
    }

    BiquadCoeff Bandpass(float wc, float Q) noexcept {
        auto r = ImpluseInvarant(wc, Q);
        auto phi = GetPhi(wc);
        auto A = GetA(r);
        auto R1 = Q * Q * std::inner_product(phi.begin(), phi.end(), A.begin(), double{});
        auto R2 = Q * Q * (-A[0] + A[1] + 4 * (phi[0] - phi[1]) * A[2]);
        auto B2 = (R1 - R2 * phi[1]) / (4 * phi[1] * phi[1]);
        auto B1 = R2 + 4 * (phi[1] - phi[0]) * B2;
        r.b1 = -std::sqrt(std::max(B1, 0.0)) / 2;
        r.b0 = (std::sqrt(std::max(0.0, B2 + B1 / 4)) - r.b1) / 2;
        r.b2 = -r.b0 - r.b1;
        return r.ToFloat();
    }

    BiquadCoeff Notch(float wc, float Q) noexcept {
        auto r = ImpluseInvarant(wc, Q);
        auto phi = GetPhi(wc);
        auto A = GetA(r);
        auto B0 = A[0];
        auto B2 = (-B0) / (4 * phi[1] * phi[1]);
        auto B1 = B0 + 4 * (phi[1] - phi[0]) * B2;
        Solveb(r, B0, B1, B2);
        return r.ToFloat();
    }

    BiquadCoeff Peaking(float wc, float Q, float db) noexcept {
        auto G = std::pow(10.0f, db / 20.0f);
        // 原型的极点是 wn = wc, Qp = A*Q (A = 10^(db/40) = sqrt(G));
        // 注意与 shelf 族的区别: shelf 原型的增益作用在 wn 上(wn = wc*sqrt(A)), Q 不变,
        // 而 peaking 原型的增益作用在 Q 上, wn 不变。
        auto r = ImpluseInvarant(wc, Q * std::sqrt(G));
        auto phi = GetPhi(wc);
        auto A = GetA(r);
        auto R1 = G * G * std::inner_product(phi.begin(), phi.end(), A.begin(), double{});
        auto R2 = G * G * (-A[0] + A[1] + 4 * (phi[0] - phi[1]) * A[2]);
        auto B0 = A[0];
        auto B2 = (R1 - R2 * phi[1] - B0) / (4 * phi[1] * phi[1]);
        auto B1 = R2 + B0 + 4 * (phi[1] - phi[0]) * B2;
        Solveb(r, B0, B1, B2);
        return r.ToFloat();
    }

    /**
     * @brief 高架滤波器 (奈奎斯特处增益 = db, 直流 0dB, wc 处 db/2)
     * @param wc 数字中心角频率 (rad/sample)
     * @param Q 品质因子
     * @param db 奈奎斯特处相对直流的增益 (dB)
     */
    BiquadCoeff Highshelf(float wc, float Q, float db) noexcept {
        AnalogResponce analog_;
        auto sqrt_G = std::pow(10.0f, db / 80.0f);
        auto r = ImpluseInvarant(wc * sqrt_G, Q);
        return FitNumerator(r, wc, [&](double w) {
            return static_cast<double>(std::norm(analog_.Highshelf(static_cast<float>(w), wc, Q, sqrt_G)));
        });
    }

    /**
     * @brief 低架滤波器 (直流处增益 = db, 奈奎斯特 0dB, wc 处 db/2)
     * @param wc 数字中心角频率 (rad/sample)
     * @param Q 品质因子
     * @param db 直流处相对奈奎斯特的增益 (dB)
     */
    BiquadCoeff Lowshelf(float wc, float Q, float db) noexcept {
        AnalogResponce analog_;
        auto sqrt_G = std::pow(10.0f, db / 80.0f);
        auto r = ImpluseInvarant(wc / sqrt_G, Q);
        return FitNumerator(r, wc, [&](double w) {
            return static_cast<double>(std::norm(analog_.Lowshelf(static_cast<float>(w), wc, Q, sqrt_G)));
        });
    }

    /**
     * @brief 倾斜架 (直流 -db/2, 奈奎斯特 +db/2, wc 处 0dB)
     * @param wc 数字中心角频率 (rad/sample)
     * @param Q 品质因子
     * @param db 两端到中心的增益落差 (dB)
     */
    BiquadCoeff Tiltshelf(float wc, float Q, float db) noexcept {
        AnalogResponce analog_;
        auto sqrt_G = std::pow(10.0f, db / 80.0f);
        auto r = ImpluseInvarant(wc * sqrt_G, Q);
        return FitNumerator(r, wc, [&](double w) {
            return static_cast<double>(std::norm(analog_.Tiltshelf(static_cast<float>(w), wc, Q, sqrt_G)));
        });
    }

    BiquadCoeff Allpass(float wc, float Q) noexcept {
        auto f_c = wc / (kPi * 2);
        auto w_c = wc;
        // zeta damping factor
        auto zeta = 1 / (2 * Q);
        auto zetasq = zeta * zeta;

        // match phases at f_c and the point where the phase is -a.
        auto a = std::clamp(std::sqrt(zeta / (2 * f_c)), 0.01, 1.5);
        auto cos_a = std::cos(a);
        auto A = cos_a + 1;
        auto sin_a_h = std::sqrt(0.5 - 0.5 * cos_a);
        auto cos_a_h = std::sqrt(0.5 + 0.5 * cos_a);

        auto R = (A * (zetasq - 1) + 2) * A;
        auto B = -2 * zetasq * A + 2 * zeta * std::sqrt(R) + A - 2;
        auto f_a = std::sqrt(B / (cos_a - 1));

        // create digital allpass through two points
        auto w_a = 2 * kPi * std::min(0.499, f_c * f_a);

        auto cos_w_c = std::cos(w_c);
        auto sin_w_a = std::sin(w_a), cos_w_a = std::cos(w_a);

        auto C = sin_a_h * (cos_w_a - cos_w_c);
        auto D = sin_w_a * cos_a_h;

        auto bot = -1 / (C + D);
        auto c_0 = (C - D) * bot;
        auto c_1 = 2 * cos_w_c * D * bot;
        return DoubleBiquadCoeff{c_0, c_1, 1, c_1, c_0}.ToFloat();
    }
private:
    static constexpr auto kPi = std::numbers::pi_v<double>;
    /// 三频点拟合的最大重试次数 (每次把最高匹配点向奈奎斯特移动)
    static constexpr int kMaxFitTrial = 40;

    /**
     * @brief 脉冲响应不变的极点映射
     * @param wc 数字截止角频率 (rad/sample)
     * @param Q 极点品质因子
     * @return 只有 a1/a2 非零的二阶节
     * @note 全程用 double 计算: 返回类型本来就是 DoubleBiquadCoeff, 而下游要算
     *       1+a1+a2 = |D(1)|, 当 wc 很小时它只有 ~wc^2 量级, 是灾难性相消; 若这里
     *       用 float 算 a1/a2, 该量只剩两三位有效数字, 低频段(如 fc=20Hz)的
     *       零极点会偏到把响应算错几十 dB。参数仍是 float 精度, 进来先提升。
     */
    static inline DoubleBiquadCoeff ImpluseInvarant(float wc, float Q) noexcept {
        double const wc_d = wc;
        double const zeta = 1.0 / (2.0 * static_cast<double>(Q));
        double const exp_qwc = std::exp(-zeta * wc_d);
        double a1 = 0.0;
        if (zeta <= 1.0) {
            a1 = -2.0 * exp_qwc * std::cos(std::sqrt(1.0 - zeta * zeta) * wc_d);
        }
        else {
            a1 = -2.0 * exp_qwc * std::cosh(std::sqrt(zeta * zeta - 1.0) * wc_d);
        }
        double const a2 = exp_qwc * exp_qwc;
        return DoubleBiquadCoeff{0, 0, 0, a1, a2};
    }

    /**
     * @brief 由 |N|^2 的三个约束解出分子系数 b0/b1/b2
     * @param c 输入输出的系数(输入用 a1/a2, 输出写入 b0/b1/b2)
     * @param B0 |N(1)|^2 (直流处分子模平方)
     * @param B1 |N(-1)|^2 (奈奎斯特处分子模平方)
     * @param B2 由 wc 处目标幅度反推的量, 满足 b0*b2 = -B2/4
     * @note 参数用 double 而不是 float: b0 由 W^2+B2 解出, 低频时二者几乎相等,
     *       残差很小, 用 float 传参会把残差的有效位数吃掉, b0/b2 就失真了
     */
    static inline void Solveb(DoubleBiquadCoeff& c, double B0, double B1, double B2) noexcept {
        auto sqrt_B0 = ClampSqrt(B0);
        auto sqrt_B1 = ClampSqrt(B1);
        auto W = (sqrt_B0 + sqrt_B1) / 2;
        c.b0 = (W + ClampSqrt(W * W + B2)) / 2;
        c.b1 = (sqrt_B0 - sqrt_B1) / 2;
        c.b2 = -B2 / (4 * c.b0);
    }

    static inline std::array<double, 3> GetA(DoubleBiquadCoeff const& c) noexcept {
        auto a1 = c.a1;
        auto a2 = c.a2;
        return {X2(1 + a1 + a2), X2(1 - a1 + a2), -4 * a2};
    }

    /**
     * @brief 判断 {B0,B1,B2} 是否对应一个真实可实现的分子
     * @param B 分子模平方的三个系数
     * @return 可行则 true
     * @note 约束来自 b0/b1/b2 必须有实数解: 要求 B0,B1 >= 0 且 W^2+B2 >= 0,
     *       其中 W = (sqrt(B0)+sqrt(B1))/2
     */
    static inline bool FeasibleNumerator(std::array<double, 3> const& B) noexcept {
        if (!(B[0] > 0.0) || !(B[1] > 0.0)) {
            return false;
        }
        auto const W = (std::sqrt(B[0]) + std::sqrt(B[1])) / 2.0;
        return W * W + B[2] > 0.0;
    }

    /**
     * @brief 解 3x3 线性方程组 (按行主元选主元)
     * @param m 系数矩阵的三行
     * @param v 右端项
     * @param out 解
     * @return 矩阵非奇异则 true
     */
    static inline bool Solve3x3(std::array<std::array<double, 3>, 3> const& m, std::array<double, 3> const& v,
                                std::array<double, 3>& out) noexcept {
        auto const det3 = [](std::array<std::array<double, 3>, 3> const& a) {
            return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                 + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
        };
        double const d = det3(m);
        // 用相对判据: 三个匹配点的 phi 行范数量级接近, 除以 1 即可
        if (!(std::abs(d) > 1e-14)) {
            return false;
        }
        for (size_t col = 0; col < 3; ++col) {
            auto mc = m;
            for (size_t row = 0; row < 3; ++row) {
                mc[row][col] = v[row];
            }
            out[col] = det3(mc) / d;
        }
        return true;
    }

    /**
     * @brief 解出高/低/倾斜架的分子, 极点保持不变
     * @param r 已确定的极点 (a1/a2 有效)
     * @param wc 数字中心角频率 (rad/sample)
     * @param target_mag2 目标幅度平方 (输入数字角频率, 返回 |H_a|^2)
     * @return 完整系数
     *
     * 数学依据 (滤波器的 |H|^2 是 phi 的二次型):
     *   |H(e^jw)|^2 = [B0*p0 + B1*p1 + B2*p2] / [A0*p0 + A1*p1 + A2*p2]
     *   其中 (p0,p1,p2) = GetPhi(w) = (1-s, s, 4s(1-s)), s = sin^2(w/2)
     * 极点确定后 A 与分母固定, 于是"在 3 个频点上让 |H|^2 等于模拟原型"就是关于
     * {B0,B1,B2} 的 3x3 线性方程组:
     *   B . phi(w_k) = |H_a(w_k)|^2 * A . phi(w_k),  k = 1,2,3
     *
     * 关于可行性: 方程组的解不一定对应真实的 b0/b1/b2。由 Solveb 的形式
     *   b0 = (W + sqrt(W^2 + B2)) / 2,  W = (sqrt(B0) + sqrt(B1)) / 2
     * 可知必须有 B0 > 0, B1 > 0, W^2 + B2 > 0, 否则 b0 不是实数(原实现就是在这里
     * 被 ClampSqrt 静默夹断, 把直流与奈奎斯特的增益约束破坏掉)。
     *
     * 不可行为何会发生: 极点频率 = wc * sqrt_G (highshelf/tiltshelf) 或 wc / sqrt_G
     * (lowshelf), sqrt_G = 10^(db/80)。增益越大、fc 越高 sqrt_G 越大, 极点频率就越
     * 容易超过奈奎斯特: highshelf/tiltshelf 在 fc > fs/(2*sqrt_G) 时越界(+18dB 时约
     * 14.3kHz)。此时可选分子在整个频段上能提供的 |N|^2 形状已无法同时命中直流、
     * 奈奎斯特与 wc 三个目标, 三条件互相矛盾, 解出的 B 落到不可行区。
     *
     * 处理办法: 先照常试原来的三点闭式解; 只有它不可行时, 才把第三个约束点从 wc
     * 起逐步向奈奎斯特方向重新分布并重解, 直到得到可行的 B。这样既保住了原有的
     * 精确匹配(绝大多数参数下与修复前逐位一致), 又不会静默破坏约束。
     *
     * @ref 相关的二极点 shelf 拟合背景: https://vicanek.de/articles/2poleShelvingFits.pdf
     */
    template <typename TargetMag2>
    static inline BiquadCoeff FitNumerator(DoubleBiquadCoeff r, double wc, TargetMag2&& target_mag2) noexcept {
        auto const A = GetA(r);
        auto const phi = GetPhi(wc);

        // ----- 首选: 三点闭式解 (直流/奈奎斯特/wc 精确) -----
        double const B0 = A[0] * target_mag2(0.0);
        double const B1 = A[1] * target_mag2(kPi);
        double const R1 = target_mag2(wc) * std::inner_product(phi.begin(), phi.end(), A.begin(), double{});
        double const B2 = (R1 - B0 * phi[0] - B1 * phi[1]) / phi[2];

        if (FeasibleNumerator({B0, B1, B2})) {
            Solveb(r, B0, B1, B2);
            return r.ToFloat();
        }

        // ----- 退化: 重分布约束点后重解 3x3 -----
        // 直流点保留(它决定整体的增益基准), 另两点向奈奎斯特方向重新分布,
        // 相当于放弃"wc 处精确", 换取一个可行的解。
        double const w_low = 0.0;
        double w_mid = 0.5 * wc;
        double w_high = wc;
        std::array<double, 3> B{};
        bool solved = false;

        for (int trial = 0; trial < kMaxFitTrial; ++trial) {
            auto const p_low = GetPhi(w_low);
            auto const p_mid = GetPhi(w_mid);
            auto const p_high = GetPhi(w_high);
            std::array<std::array<double, 3>, 3> const m{p_low, p_mid, p_high};
            std::array<double, 3> const rhs{
                target_mag2(w_low) * std::inner_product(p_low.begin(), p_low.end(), A.begin(), double{}),
                target_mag2(w_mid) * std::inner_product(p_mid.begin(), p_mid.end(), A.begin(), double{}),
                target_mag2(w_high) * std::inner_product(p_high.begin(), p_high.end(), A.begin(), double{}),
            };
            std::array<double, 3> cand{};
            if (Solve3x3(m, rhs, cand) && FeasibleNumerator(cand)) {
                B = cand;
                solved = true;
                break;
            }
            w_mid = 0.5 * (w_mid + w_high);
            w_high = 0.5 * (w_high + kPi);
        }

        if (solved) {
            Solveb(r, B[0], B[1], B[2]);
        }
        else {
            // 实在解不出: 用原来的夹断结果, 至少给出一个能用的滤波器
            Solveb(r, B0, B1, B2);
        }
        return r.ToFloat();
    }

    static inline std::array<double, 3> GetPhi(double w) noexcept {
        auto sin2 = X2(std::sin(w / 2));
        auto phi0 = 1 - sin2;
        auto phi1 = sin2;
        return {phi0, phi1, 4 * phi0 * phi1};
    }

    static inline constexpr double X2(double x) noexcept {
        return x * x;
    }

    static inline double ClampSqrt(double x) noexcept {
        return std::sqrt(std::max(x, 0.0));
    }
};
} // namespace qwqdsp_filter
