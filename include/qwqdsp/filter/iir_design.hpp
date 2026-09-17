#pragma once
#include "biquad_coeff.hpp"
#include <array>
#include <cassert>
#include <complex>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <vector>

namespace qwqdsp_filter {
struct IIRDesign {
    static constexpr auto pi = std::numbers::pi;

    struct ZPK {
        std::optional<std::complex<double>> z; // 如果Null则在无穷远处
        std::complex<double> p;
        double k;
    };

    /**
     * @brief 复数按实数等比缩放
     * @param a 被缩放的复数
     * @param b 缩放系数
     * @return a * b
     */
    static std::complex<double> ScaleComplex(const std::complex<double>& a, double b) {
        return {a.real() * b, a.imag() * b};
    }

    // --------------------------------------------------------------------------------
    // 原型滤波器
    //
    // 所有原型都是模拟低通原型, 每个 ZPK 对应一个二阶节(一对共轭极点),
    // 因此阶数 = 2 * num_filter(总是偶数)。原型的 "1rad/sec" 是模拟角频率,
    // 要得到数字滤波器还需要经过 ProtyleTo* 映射与 Bilinear 离散化。
    // --------------------------------------------------------------------------------

    /**
     * @brief 巴特沃斯原型, 通带最平坦, (1)rad/sec 处 -3.01dB, 直流增益 0dB
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @note 零点全在无穷远处, 调用前请保证 ret[i].z 为空
     * @note 逐节 k = 1
     */
    static void Butterworth(std::span<ZPK> ret, size_t num_filter) {
        assert(ret.size() >= num_filter);

        size_t n = 2 * num_filter;
        size_t i = 0;
        for (size_t k = 1; k <= num_filter; ++k) {
            double phi = (2.0 * static_cast<double>(k) - 1.0) * pi / (2.0 * static_cast<double>(n));
            ret[i].p = std::complex{-std::sin(phi), std::cos(phi)};
            ret[i].k = 1;
            ++i;
        }
    }

    /**
     * @brief 切比雪夫 I 型原型, 通带等波纹, (1)rad/sec 处 -ripple dB
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param ripple 通带纹波(dB, >0), 通带在 0dB 与 -ripple dB 之间等波纹
     * @param even_pole_modify 偶数阶修正(见 @ref 链接): 把通带参考电平从 -ripple dB
     *        抬到 0dB(直流增益 1), 通带变成 [0, -ripple]dB; (1)rad/sec 处仍是 -ripple dB
     * @note 零点全在无穷远处, 调用前请保证 ret[i].z 为空
     * @note 逐节 k = |p|^2; 不修正时 section 0 再除以 sqrt(1 + eps^2)
     * @ref https://en.wikipedia.org/wiki/Chebyshev_filter
     */
    static void Chebyshev1(std::span<ZPK> ret, size_t num_filter, double ripple, bool even_pole_modify) {
        assert(ret.size() >= num_filter);

        size_t n = 2 * num_filter;
        size_t i = 0;
        double eps = std::sqrt(std::pow(10.0, ripple / 10.0) - 1.0);
        double A = 1.0 / static_cast<double>(n) * std::asinh(1.0 / eps);
        double k_re = std::sinh(A);
        double k_im = std::cosh(A);
        double first_pole = std::cos(pi * (static_cast<double>(n) - 1.0) / (2.0 * static_cast<double>(n)));
        first_pole = first_pole * first_pole;
        for (size_t k = 1; k <= num_filter; ++k) {
            double phi = (2.0 * static_cast<double>(k) - 1.0) * pi / (2.0 * static_cast<double>(n));
            if (even_pole_modify) {
                auto pole = std::complex{-std::sin(phi) * k_re, std::cos(phi) * k_im};
                ret[i].p = std::sqrt((pole * pole + first_pole) / (1.0 - first_pole));
            }
            else {
                ret[i].p = std::complex{-std::sin(phi) * k_re, std::cos(phi) * k_im};
            }
            ret[i].k = std::norm(ret[i].p);
            ++i;
        }
        if (!even_pole_modify) {
            ret[0].k /= std::sqrt(1.0f + eps * eps);
        }
    }

    /**
     * @brief 切比雪夫 II 型(逆切比雪夫)原型, 通带最平坦, 阻带等波纹
     *
     * (1)rad/sec 处 -3.01dB, 直流增益 0dB, 阻带等波纹深度为 -ripple dB。
     * 形状与 scipy/MATLAB 的 cheby2 相同, 区别只是归一化位置: 这里把 -3.01dB 点
     * 放在 (1)rad/sec, 而 scipy 的 Wn 是阻带边沿(第一瓣峰值处)。
     *
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param ripple 阻带衰减(dB, <0, 例如 -40), 注意与 Chebyshev1 的正的通带纹波相反
     * @param even_order_modify 偶数阶修正: 让最后一对极点不再有有限零点
     * @note 零点在虚轴上(有限频率)
     * @note 阻带第一个等波纹峰值在模拟角频率 cosh(acosh(sqrt(S))/n) rad/sec 处,
     *       其中 S = 10^(-ripple/10) - 1, n = 2 * num_filter
     * @ref https://en.wikipedia.org/wiki/Chebyshev_filter
     * @ref https://en.wikipedia.org/wiki/Chebyshev_nodes#Even_order_modified_Chebyshev_nodes
     */
    static void Chebyshev2(std::span<ZPK> ret, size_t num_filter, double ripple, bool even_order_modify) {
        assert(ret.size() >= num_filter);

        size_t n = 2 * num_filter;
        size_t i = 0;
        double eps = 1.0 / std::sqrt(std::pow(10.0, -ripple / 10.0) - 1.0);
        double A = 1.0 / static_cast<double>(n) * std::asinh(1.0 / eps);
        double scale =
            1.0 / std::cosh(std::acosh(std::sqrt(std::pow(10.0, -ripple / 10.0) - 1.0)) / static_cast<double>(n));
        double k_re = std::sinh(A) * scale;
        double k_im = std::cosh(A) * scale;

        double first_pole = std::cos(pi * (static_cast<double>(n) - 1.0) / (2.0 * static_cast<double>(n)));
        first_pole = first_pole * first_pole;
        // 最接近0的零点
        double const first_zero =
            std::cos((static_cast<double>(n) / 2.0 - 1.0 + 0.5) * std::numbers::pi_v<double> / static_cast<double>(n));
        for (size_t k = 1; k <= num_filter; ++k) {
            double phi = (2.0 * static_cast<double>(k) - 1.0) * pi / (2.0 * static_cast<double>(n));
            if (!even_order_modify) {
                ret[i].z = 1.0 / std::complex{0.0, std::cos(phi) * scale};
                ret[i].p = 1.0 / std::complex{-std::sin(phi) * k_re, std::cos(phi) * k_im};
            }
            else {
                auto pole = std::complex{-std::sin(phi) * k_re, std::cos(phi) * k_im};
                ret[i].p = 1.0 / std::sqrt((pole * pole + first_pole) / (1.0 - first_pole));
                if (k != num_filter) {
                    // 最靠近0的切比雪夫多项式的零点被映射到0，所以零点在无穷远处不赋值
                    double const zero = std::cos(phi);
                    double const tt = std::sqrt(
                        std::max(0.0, (zero * zero - first_zero * first_zero) / (1.0 - first_zero * first_zero)));
                    ret[i].z = 1.0 / std::complex{0.0, tt * scale};
                }
            }
            if (ret[i].z) {
                ret[i].k = std::norm(ret[i].p) / std::norm(*ret[i].z);
            }
            else {
                ret[i].k = std::norm(ret[i].p);
            }
            ++i;
        }
    }

    struct EllipticHelper {
        static constexpr double kError = std::numeric_limits<double>::epsilon();
        static constexpr bool kUseStd = false;

        EllipticHelper(double k0) {
            while (k0 > kError) {
                k_.push_back(k0);
                double kdot = std::sqrt(1.0 - k0 * k0);
                kdot_.push_back(kdot);
                k0 = (1.0 - kdot) / (1.0 + kdot);
            }
            k_.push_back(k0);
            kdot_.push_back(Kdot(k0));
        }

        double CompleteIntegral() {
#ifndef __APPLE__
            if constexpr (kUseStd) {
                return std::comp_ellint_1(k_.front());
            }
            else {
#endif
                const double a = std::sqrt(1.0 - 1.0e-3);
                const double k0 = k_.front();
                if (k0 <= a) {
                    double km = pi / 2.0;
                    for (size_t i = 1; i < k_.size(); ++i) {
                        double const k = k_[i];
                        km = km * (1 + k);
                    }
                    return km;
                }
                else {
                    const double kdot = kdot_.front();
                    double L = -std::log(kdot / 4.0);
                    return L + (L - 1.0) * (kdot * kdot) / 2.0;
                }
#ifndef __APPLE__
            }
#endif
        }

        std::complex<double> Cd(std::complex<double> u) {
            auto cdm = std::cos(u * pi / 2.0);
            for (size_t i = 0; i < k_.size() - 1; ++i) {
                double const k = k_[k_.size() - i - 1];
                cdm = (1.0 + k) / (1.0 / cdm + k * cdm);
            }
            return cdm;
        }

        double Cd(double u) {
            auto cdm = std::cos(u * pi / 2.0);
            for (size_t i = 0; i < k_.size() - 1; ++i) {
                double const k = k_[k_.size() - i - 1];
                cdm = (1.0 + k) / (1.0 / cdm + k * cdm);
            }
            return cdm;
        }

        std::complex<double> Sn(std::complex<double> u) {
            auto snm = std::sin(u * (pi / 2.0));
            for (size_t i = 0; i < k_.size() - 1; ++i) {
                double const k = k_[k_.size() - 1 - i];
                snm = (1.0 + k) / (1.0 / snm + k * snm);
            }
            return snm;
        }

        double Sn(double u) {
            auto snm = std::sin(u * (pi / 2.0));
            for (size_t i = 0; i < k_.size() - 1; ++i) {
                double const k = k_[k_.size() - 1 - i];
                snm = (1.0 + k) / (1.0 / snm + k * snm);
            }
            return snm;
        }

        std::complex<double> ArcSn(std::complex<double> sn0) {
            for (size_t i = 1; i < k_.size(); ++i) {
                sn0 = 2.0 * sn0 / ((1.0 + k_[i]) * (1.0 + std::sqrt(1.0 - k_[i - 1] * k_[i - 1] * sn0 * sn0)));
            }
            return std::asin(sn0) * (2.0 / pi);
        }

        double ArcSn(double sn0) {
            for (size_t i = 1; i < k_.size(); ++i) {
                sn0 = 2.0 * sn0 / ((1.0 + k_[i]) * (1.0 + std::sqrt(1.0 - k_[i - 1] * k_[i - 1] * sn0 * sn0)));
            }
            return std::asin(sn0) * (2.0 / pi);
        }

        static double Kdot(double k) {
            return std::sqrt(1.0 - k * k);
        }
    private:
        std::vector<double> k_;
        std::vector<double> kdot_;
    };

    /**
     * @brief cephes 椭圆函数工具(移植自 ref_repos/cephes-master/ellf)
     *
     * 与 EllipticHelper 的区别:
     * - 完全积分 K 用 minimax 多项式逼近(cephes ellpk), 不做 Landen 下降;
     * - 由 nome 的 theta 级数反解模数(cephes cay), 在 k 逼近 1 时收敛最快,
     *   不会像 Landen 上升/下降那样退化或死循环;
     * - sn/cn/dn 用 AGM + 反向递推(cephes ellpj), 参数是弧度;
     * - 不完全积分 F(phi|m) 用 AGM(cephes ellik)。
     *
     * @note cephes 的 ellpk(x) 里 x = 1 - m, 即 ellpk(x) = K(m)
     */
    struct EllipticHelperCephes {
        /// cephes polevl: 降幂 Horner 求值
        static double Polevl(double x, std::span<double const> coef) {
            double ans = coef[0];
            for (size_t i = 1; i < coef.size(); ++i) {
                ans = ans * x + coef[i];
            }
            return ans;
        }

        /**
         * @brief 完全椭圆积分 K(m), 参数 m = 1 - m1
         * @param m1 互补参数, 即 K 的奇点在 m1 = 0
         * @return K(m)
         */
        static double Ellpk(double m1) {
            static constexpr std::array<double, 11> kP{
                1.37982864606273237150E-4, 2.28025724005875567385E-3, 7.97404013220415179367E-3,
                9.85821379021226008714E-3, 6.87489687449949877925E-3, 6.18901033637687613229E-3,
                8.79078273952743772254E-3, 1.49380448916805252718E-2, 3.08851465246711995998E-2,
                9.65735902811690126535E-2, 1.38629436111989062502E0};
            static constexpr std::array<double, 11> kQ{
                2.94078955048598507511E-5, 9.14184723865917226571E-4, 5.94058303753167793257E-3,
                1.54850516649762399335E-2, 2.39089602715924892727E-2, 3.01204715227604046988E-2,
                3.73774314173823228969E-2, 4.88280347570998239232E-2, 7.03124996963957469739E-2,
                1.24999999999870820058E-1, 4.99999999999999999821E-1};
            constexpr double kMachep = 2.220446049250313e-16;
            if (m1 > kMachep) {
                return Polevl(m1, kP) - std::log(m1) * Polevl(m1, kQ);
            }
            // log(4) - log(sqrt(m1))
            return 1.3862943611198906188E0 - 0.5 * std::log(m1);
        }

        /**
         * @brief 由 nome 反解模数 k = sqrt(m)
         * @param q nome, 0 < q < 1
         * @return 模数 k; q 太大时会返回 >= 1 的值, 由调用者判无效
         * @note cephes cay(): (2K/pi)^(1/2) m^(1/4) 用 theta 级数展开后反解
         */
        static double Cay(double q) {
            constexpr double kMachep = 2.220446049250313e-16;
            double a = 1.0;
            double b = 1.0;
            double r = 1.0;
            double p = q;
            for (int iter = 0; iter < 10000; ++iter) {
                r *= p;
                a += 2.0 * r;
                double t1 = std::abs(r / a);
                r *= p;
                b += r;
                p *= q;
                double const t2 = std::abs(r / b);
                if (t2 > t1) {
                    t1 = t2;
                }
                if (!(t1 > kMachep)) {
                    break;
                }
            }
            a = b / a;
            return 4.0 * std::sqrt(q) * a * a;
        }

        /**
         * @brief 不完全椭圆积分 F(phi|m)
         * @param phi 幅角(弧度)
         * @param m 参数 m = k^2
         * @return F(phi|m)
         */
        static double Ellik(double phi, double m) {
            constexpr double kMachep = 2.220446049250313e-16;
            if (m == 0.0) {
                return phi;
            }
            double const a0 = 1.0 - m;
            if (a0 == 0.0) {
                if (std::abs(phi) >= pi / 2.0) {
                    return std::numeric_limits<double>::max();
                }
                return std::log(std::tan((pi / 2.0 + phi) / 2.0));
            }
            int npio2 = static_cast<int>(std::floor(phi / (pi / 2.0)));
            if (npio2 & 1) {
                npio2 += 1;
            }
            double K = 0.0;
            if (npio2) {
                K = Ellpk(a0);
                phi = phi - static_cast<double>(npio2) * (pi / 2.0);
            }
            int sign = 0;
            if (phi < 0.0) {
                phi = -phi;
                sign = -1;
            }
            double b = std::sqrt(a0);
            double t = std::tan(phi);
            bool transformed = false;
            double result = 0.0;
            if (std::abs(t) > 10.0) {
                double const e = 1.0 / (b * t);
                if (std::abs(e) < 10.0) {
                    if (npio2 == 0) {
                        K = Ellpk(a0);
                    }
                    result = K - Ellik(std::atan(e), m);
                    transformed = true;
                }
            }
            if (!transformed) {
                double a = 1.0;
                double c = std::sqrt(m);
                int d = 1;
                int mod = 0;
                while (std::abs(c / a) > kMachep) {
                    double const ratio = b / a;
                    phi = phi + std::atan(t * ratio) + static_cast<double>(mod) * pi;
                    mod = static_cast<int>((phi + pi / 2.0) / pi);
                    t = t * (1.0 + ratio) / (1.0 - ratio * t * t);
                    c = (a - b) / 2.0;
                    double const next_b = std::sqrt(a * b);
                    a = (a + b) / 2.0;
                    b = next_b;
                    d += d;
                }
                result = (std::atan(t) + static_cast<double>(mod) * pi) / (static_cast<double>(d) * a);
            }
            if (sign < 0) {
                result = -result;
            }
            return result + static_cast<double>(npio2) * K;
        }

        /**
         * @brief 实参数下的 Jacobi 椭圆函数 sn/cn/dn
         * @param u 自变量(弧度)
         * @param m 参数 m = k^2, 需要 0 <= m <= 1
         * @param sn 输出 sn(u, k)
         * @param cn 输出 cn(u, k)
         * @param dn 输出 dn(u, k)
         */
        static void Ellpj(double u, double m, double& sn, double& cn, double& dn) {
            constexpr double kMachep = 2.220446049250313e-16;
            if (m < 1.0e-9) {
                double const t = std::sin(u);
                double const b = std::cos(u);
                double const ai = 0.25 * m * (u - t * b);
                sn = t - ai * b;
                cn = b + ai * t;
                dn = 1.0 - 0.5 * m * t * t;
                return;
            }
            if (m >= 0.9999999999) {
                double const ai = 0.25 * (1.0 - m);
                double const b = std::cosh(u);
                double const t = std::tanh(u);
                double const ph = 1.0 / b;
                double const twon = b * std::sinh(u);
                sn = t + ai * (twon - u) / (b * b);
                double const ai2 = ai * t * ph;
                cn = ph - ai2 * (twon - u);
                dn = ph + ai2 * (twon + u);
                return;
            }
            std::array<double, 9> a{};
            std::array<double, 9> c{};
            a[0] = 1.0;
            double b = std::sqrt(1.0 - m);
            double twon = 1.0;
            c[0] = std::sqrt(m);
            int i = 0;
            while (std::abs(c[i] / a[i]) > kMachep && i <= 7) {
                double const ai = a[i];
                ++i;
                c[i] = (ai - b) / 2.0;
                double const next_b = std::sqrt(ai * b);
                a[i] = (ai + b) / 2.0;
                b = next_b;
                twon *= 2.0;
            }
            // 反向递推
            double phi = twon * a[i] * u;
            do {
                double const t = c[i] * std::sin(phi) / a[i];
                b = phi;
                phi = (std::asin(t) + phi) / 2.0;
            } while (--i > 0);
            sn = std::sin(phi);
            double const t = std::cos(phi);
            cn = t;
            dn = t / std::cos(phi - b);
        }
    };

    // qwqfixme 偶数极点零点修改
    /**
     * @brief 椭圆(考尔)原型, 通带与阻带都等波纹
     *
     * 通带边沿在 (1)rad/sec 且增益为 -db_passband dB, 直流增益也是 -db_passband dB;
     * 阻带边沿在 (1)rad/sec 的 1/k 倍处(k 是椭圆模数, 由阶数与两个纹波决定),
     * 从那里起等波纹深度为 -db_stopband dB。相同阶数下过渡带最陡。
     *
     * 模数 k 用 cephes 的做法(由 nome 的 theta 级数反解)计算, 不使用 Landen 下降
     * 序列, 因此在很陡的规格下也不会陷入死循环; 旧的 Orfanidis 式实现保留在
     * EllipticLanden 里。
     *
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param db_passband 通带纹波(dB, >0), 通带边沿电平
     * @param db_stopband 阻带衰减(dB, >0), 需要大于 db_passband
     * @return 成功返回 true
     * @retval false db_stopband <= db_passband, 或该规格要求的模数在 double 下就是 1
     *         (规格过陡, 阻带边沿与通带边沿重合; 例如 order=16, 通带 6dB, 阻带 10dB)
     * @note 零点在虚轴上(有限频率)
     * @ref Gray & Markel, "A Computer Program for Designing Digital Elliptic Filters",
     *      IEEE Trans. ASSP, Dec. 1976; 以及 cephes 的 ellf.c
     */
    [[nodiscard]] static bool Elliptic(std::span<ZPK> ret, size_t num_filter, double db_passband, double db_stopband) {
        assert(ret.size() >= num_filter);

        double const eps_passband = std::sqrt(std::pow(10.0, db_passband / 10.0) - 1.0);
        double const eps_stopband = std::sqrt(std::pow(10.0, db_stopband / 10.0) - 1.0);
        if (!(eps_stopband > eps_passband)) {
            return false;
        }
        size_t const N = 2 * num_filter;
        double const m1 = eps_passband / eps_stopband; // 选择性(第二模数)
        double const m1_sq = m1 * m1;
        // 设计模数: q = exp(-pi*K(k1')/(N*K(k1))), k = cay(q)
        double const Kk1 = EllipticHelperCephes::Ellpk(1.0 - m1_sq);
        double const Kpk1 = EllipticHelperCephes::Ellpk(m1_sq);
        double const q = std::exp(-pi * Kpk1 / (static_cast<double>(N) * Kk1));
        double const k = EllipticHelperCephes::Cay(q);
        if (!(k > 0.0 && k < 1.0)) {
            return false;
        }
        double const m = k * k;
        double const Kk = EllipticHelperCephes::Ellpk(1.0 - m);
        // u = F(atan(1/eps_p) | k1'^2) * K(k) / (N*K(k1))
        double const u = EllipticHelperCephes::Ellik(std::atan(1.0 / eps_passband), 1.0 - m1_sq) * Kk
                       / (static_cast<double>(N) * Kk1);
        double sn1 = 0.0;
        double cn1 = 0.0;
        double dn1 = 0.0;
        EllipticHelperCephes::Ellpj(u, 1.0 - m, sn1, cn1, dn1);

        for (size_t i = 0; i < num_filter; ++i) {
            auto& s = ret[i];
            double const arg = static_cast<double>(N - 1 - 2 * static_cast<int>(i)) * Kk / static_cast<double>(N);
            double sn = 0.0;
            double cn = 0.0;
            double dn = 0.0;
            EllipticHelperCephes::Ellpj(arg, m, sn, cn, dn);
            // 零点在虚轴上
            s.z = std::complex{0.0, 1.0 / (k * sn)};
            // 极点: Gray & Markel 的实值公式
            double const r = k * sn * sn1;
            double const den = cn1 * cn1 + r * r;
            s.p = std::complex{-cn * dn * sn1 * cn1 / den, sn * dn1 / den};
            s.k = std::norm(s.p) / std::norm(*s.z);
        }
        ret[0].k /= std::sqrt(1.0 + eps_passband * eps_passband);
        return true;
    }

    /**
     * @brief 椭圆(考尔)原型, Orfanidis 式实现(Landen 序列求模数)
     *
     * 与 Elliptic 的数学定义相同, 但模数用 k' = (k1')^N * prod sn(u_i)^4 的
     * Landen 路线求, 需要构造多级 Landen 序列。
     *
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param db_passband 通带纹波(dB, >0), 通带边沿电平
     * @param db_stopband 阻带衰减(dB, >0), 需要大于 db_passband
     * @return 成功返回 true
     * @retval false 参数非法, 或求出的互补模数不在 (0, 1)(此时继续做 Landen 下降
     *         会死循环, 所以直接返回失败)
     * @warning 深阻带 + 小通带纹波时(k1 很小)精度不如 Elliptic(nome 反解),
     *          例如 order=12, 通带 0.05dB, 阻带 140dB 时模数误差可达 1e-2
     * @note 供对照/兼容使用, 新代码建议直接用 Elliptic
     */
    [[nodiscard]] static bool EllipticLanden(std::span<ZPK> ret, size_t num_filter, double db_passband,
                                             double db_stopband) {
        assert(ret.size() >= num_filter);

        double const eps_passband = std::sqrt(std::pow(10.0, db_passband / 10.0) - 1.0);
        double const eps_stopband = std::sqrt(std::pow(10.0, db_stopband / 10.0) - 1.0);
        double const k1 = eps_passband / eps_stopband;
        if (!(k1 > 0.0 && k1 < 1.0)) {
            return false;
        }
        size_t const N = 2 * num_filter;
        auto const k1dot = EllipticHelper::Kdot(k1);
        if (!(k1dot > 0.0 && k1dot < 1.0)) {
            return false;
        }
        // ellipdeg: k' = (k1')^N * prod sn(u_i)^4
        double kdot = 0.0;
        {
            EllipticHelper helper{k1dot};
            double const f1 = std::pow(k1dot, N);
            std::complex<double> back{1.0, 0.0};
            for (size_t i = 1; i <= num_filter; ++i) {
                auto ui = (2.0 * static_cast<double>(i) - 1.0) / static_cast<double>(N);
                back *= std::pow(helper.Sn(ui), 4.0);
            }
            kdot = f1 * std::real(back);
        }
        if (!(kdot > 0.0 && kdot < 1.0)) {
            return false;
        }
        double const k = EllipticHelper::Kdot(kdot);
        // 注意: 这里必须检验 k 本身。kdot 只要是"非零的很小值"(例如 2e-9),
        // sqrt(1 - kdot^2) 就会在 double 下舍入成精确的 1.0, 于是下面的
        // EllipticHelper{k} 会卡在 k0 恒为 1 的 Landen 下降里死循环。
        if (!(k > 0.0 && k < 1.0)) {
            return false;
        }

        EllipticHelper helper{k};
        EllipticHelper helper1{k1};
        // ArcSn 返回的就是以 K(k1) 归一化的自变量, 所以这里只除以 N;
        // 若再除一次 K(k1) 会让极点整体偏移, 通带纹波变成向上凸
        auto const v0 =
            std::complex{0.0, -1.0} * helper1.ArcSn(std::complex{0.0, 1.0} / eps_passband) / static_cast<double>(N);
        for (size_t i = 0; i < num_filter; ++i) {
            auto& s = ret[i];
            auto ui = (2.0 * static_cast<double>(i + 1) - 1.0) / static_cast<double>(N);
            // zero
            auto epsi = helper.Cd(ui);
            s.z = std::complex<double>{0.0, 1.0} / (k * epsi);
            // pole
            s.p = std::complex{0.0, 1.0} * helper.Cd(ui - v0 * std::complex{0.0, 1.0});
            s.k = std::norm(s.p) / std::norm(*s.z);
        }
        ret[0].k /= std::sqrt(1.0 + eps_passband * eps_passband);
        return true;
    }

    // --------------------------------------------------------------------------------
    // 滤波器映射
    // --------------------------------------------------------------------------------
    /**
     * @brief 原型 -> 低通: 极点与零点整体乘以 omega
     * @param analog 被就地修改的原型节
     * @param num_filter 参与的极点对数, 只处理前 num_filter 节
     * @param omega 目标通带边沿的模拟角频率(rad/sec)
     * @note 原型 (1)rad/sec 处的边沿被搬到 omega
     * @note 无有限零点(零点在无穷远)的节 k 乘 omega^2, 有零点的节 k 不变,
     *       这样直流增益与原型相同
     */
    static void ProtyleToLowpass(std::span<ZPK> analog, size_t num_filter, double omega) {
        assert(analog.size() >= num_filter);

        for (size_t i = 0; i < num_filter; ++i) {
            const auto& s = analog[i];
            ZPK lps = analog[i];
            double gain = omega * omega;
            lps.p = ScaleComplex(s.p, omega);
            if (s.z) {
                lps.z = ScaleComplex(*s.z, omega);
                gain = 1.0;
            }
            analog[i] = lps;
            analog[i].k *= gain;
        }
    }

    /**
     * @brief 原型 -> 高通: 做 s -> omega/s 的频率反转
     * @param protyle 被就地修改的原型节
     * @param num_filter 参与的极点对数, 只处理前 num_filter 节
     * @param omega 目标通带边沿的模拟角频率(rad/sec)
     * @note 极点为 omega/p; 原型有有限零点时零点变成 omega/z, 否则在原点补一个零点
     * @note k 按 1/|p|^2(有零点时再乘 |z|^2) 归一化, 使高频增益与原型直流增益相同
     */
    static void ProtyleToHighpass(std::span<ZPK> protyle, size_t num_filter, double omega) {
        assert(protyle.size() >= num_filter);

        for (size_t i = 0; i < num_filter; ++i) {
            const auto& s = protyle[i];
            ZPK hps = protyle[i];
            double gain = 1.0 / std::norm(s.p);
            hps.p = omega / s.p;
            if (s.z) {
                gain *= std::norm(*s.z);
                hps.z = omega / *s.z;
            }
            else {
                hps.z = 0;
            }
            protyle[i] = hps;
            protyle[i].k *= gain;
        }
    }

    /**
     * @brief 原型 -> 带通, 用中心频率与 Q 指定带宽
     *
     * 做 u = (s^2 + wo^2) / (s * bw) 的映射(bw = wo / Q), 原型 (1)rad/sec 处的边沿
     * 变成 s^2 - bw*s + wo^2 = 0 的两个根, 两个边沿的几何平均是 wo、差是 bw。
     *
     * @param protyle 被就地修改的原型节
     * @param num_filter 原型极点对数, 处理前 num_filter 节
     * @param wo 带通中心角频率(rad/sec), 两个边沿的几何平均
     * @param Q 中心频率 / 带宽, 需要 > 0
     * @note 节数会翻倍: 结果写入 [0, num_filter) 与 [num_filter, 2*num_filter) 两半,
     *       因此 protyle 至少要有 2 * num_filter 个元素
     * @note 无零点的原型节会拆成"原点零点"与"无穷远零点"两半
     */
    static void ProtyleToBandpass(std::span<ZPK> protyle, size_t num_filter, double wo, double Q) {
        assert(protyle.size() >= num_filter * 2);

        double bw = wo / Q;
        double wo2 = wo * wo;

        for (size_t i = 0; i < num_filter; ++i) {
            ZPK s = protyle[i];
            ZPK bp1;
            ZPK bp2;
            bp1.k = s.k;
            bp2.k = 1;

            if (s.z) {
                auto z_delta = std::sqrt((*s.z * *s.z * bw * bw) - 4.0 * wo2);
                bp1.z = ScaleComplex((*s.z * bw) + z_delta, 0.5);
                bp2.z = ScaleComplex((*s.z * bw) - z_delta, 0.5);
            }
            else {
                bp1.z = 0.0;
                bp2.z = std::nullopt;
                bp1.k *= bw;
                bp2.k *= bw;
            }

            auto p_delta = std::sqrt((s.p * s.p * bw * bw) - 4.0 * wo2);
            bp1.p = ScaleComplex((s.p * bw) + p_delta, 0.5);
            bp2.p = ScaleComplex((s.p * bw) - p_delta, 0.5);

            protyle[i] = bp1;
            protyle[i + num_filter] = bp2;
        }
    }

    /**
     * @brief 原型 -> 带通, 用两个边沿频率指定带宽
     * @param protyle 被就地修改的原型节
     * @param num_filter 原型极点对数, 处理前 num_filter 节
     * @param w1 低边沿的模拟角频率(rad/sec), w1 < w2
     * @param w2 高边沿的模拟角频率(rad/sec)
     * @note 带宽取 bw = w2 - w1; 节数翻倍, protyle 至少要有 2 * num_filter 个元素
     * @note 与 ProtyleToBandpass(wo, Q) 的区别是这里直接给两个边沿,
     *       且两半节的 k 都取 sqrt(原型节 k)
     */
    static void ProtyleToBandpass2(std::span<ZPK> protyle, size_t num_filter, double w1, double w2) {
        assert(protyle.size() >= num_filter * 2);

        double bw = w2 - w1;
        for (size_t i = 0; i < num_filter; ++i) {
            ZPK bp1;
            ZPK bp2;
            ZPK s = protyle[i];
            bp2.k = bp1.k = std::sqrt(s.k);
            if (s.z) {
                auto p_delta = std::sqrt(s.p * s.p * bw * bw - 4.0 * w1 * w2);
                auto z_delta = std::sqrt(*s.z * *s.z * bw * bw - 4.0 * w1 * w2);
                bp1.p = ScaleComplex(s.p * bw + p_delta, 0.5);
                bp2.p = ScaleComplex(s.p * bw - p_delta, 0.5);
                bp1.z = ScaleComplex(*s.z * bw + z_delta, 0.5);
                bp2.z = ScaleComplex(*s.z * bw - z_delta, 0.5);
            }
            else {
                auto delta = std::sqrt(s.p * s.p * bw * bw - 4.0 * w1 * w2);
                bp1.p = ScaleComplex(s.p * bw + delta, 0.5);
                bp2.p = ScaleComplex(s.p * bw - delta, 0.5);
                bp1.z = 0;
                bp1.k *= bw;
                bp2.k *= bw;
            }
            protyle[i] = bp1;
            protyle[i + num_filter] = bp2;
        }
    }

    /**
     * @brief 原型 -> 带阻, 用中心频率与 Q 指定带宽
     *
     * 先做原型 -> 高通(bw = wo / Q), 再做高通 -> 带阻; 原型 (1)rad/sec 处的边沿
     * 变成 s^2 - bw*s + wo^2 = 0 的两个根, 两个边沿的几何平均是 wo、差是 bw。
     *
     * @param protyle 被就地修改的原型节
     * @param num_filter 原型极点对数, 处理前 num_filter 节
     * @param wo 带阻中心角频率(rad/sec), 也是陷波零点所在频率
     * @param Q 中心频率 / 带宽, 需要 > 0
     * @note 节数会翻倍: 结果写入 [0, num_filter) 与 [num_filter, 2*num_filter) 两半,
     *       因此 protyle 至少要有 2 * num_filter 个元素
     */
    static void ProtyleToBandstop(std::span<ZPK> protyle, size_t num_filter, double wo, double Q) {
        assert(protyle.size() >= num_filter * 2);

        double bw = wo / Q;
        // 只遍历原型节, 结果写入 i 与 i + num_filter 两半
        for (size_t i = 0; i < num_filter; ++i) {
            // prototype -> highpass at bw
            ZPK s = protyle[i];
            {
                auto const& ss = protyle[i];
                s.k /= std::norm(ss.p);
                s.p = bw / ss.p;
                if (ss.z) {
                    s.k *= std::norm(*ss.z);
                    s.z = bw / *ss.z;
                }
                else {
                    s.z = 0;
                }
            }
            // highpass -> bandstop
            ZPK bp1;
            ZPK bp2;
            bp1.k = bp2.k = std::sqrt(s.k);
            if (s.z) {
                auto p_delta = std::sqrt(s.p * s.p - 4.0 * wo * wo);
                auto z_delta = std::sqrt(*s.z * *s.z - 4.0 * wo * wo);
                bp1.p = ScaleComplex(s.p + p_delta, 0.5);
                bp2.p = ScaleComplex(s.p - p_delta, 0.5);
                bp1.z = ScaleComplex(*s.z + z_delta, 0.5);
                bp2.z = ScaleComplex(*s.z - z_delta, 0.5);
            }
            else {
                auto delta = std::sqrt(s.p * s.p - 4.0 * wo * wo);
                bp1.p = ScaleComplex(s.p + delta, 0.5);
                bp2.p = ScaleComplex(s.p - delta, 0.5);
                bp1.z = 0;
            }
            protyle[i] = bp1;
            protyle[i + num_filter] = bp2;
        }
    }

    /**
     * @brief 原型 -> 带阻, 用两个边沿频率指定带宽
     * @param protyle 被就地修改的原型节
     * @param num_filter 原型极点对数, 处理前 num_filter 节
     * @param w1 低边沿的模拟角频率(rad/sec), w1 < w2
     * @param w2 高边沿的模拟角频率(rad/sec)
     * @note 带宽取 bw = w2 - w1; 节数翻倍, protyle 至少要有 2 * num_filter 个元素
     * @note 与 ProtyleToBandstop(wo, Q) 的区别是这里直接给两个边沿
     */
    static void ProtyleToBandstop2(std::span<ZPK> protyle, size_t num_filter, double w1, double w2) {
        assert(protyle.size() >= 2 * num_filter);

        double bw = w2 - w1;
        // 只遍历原型节, 结果写入 i 与 i + num_filter 两半
        for (size_t i = 0; i < num_filter; ++i) {
            ZPK s = protyle[i];
            {
                auto const& ss = protyle[i];
                s.k /= std::norm(ss.p);
                s.p = bw / ss.p;
                if (ss.z) {
                    s.k *= std::norm(*ss.z);
                    s.z = bw / *ss.z;
                }
                else {
                    s.z = 0;
                }
            }

            ZPK bp1;
            ZPK bp2;
            bp1.k = bp2.k = std::sqrt(s.k);
            if (s.z) {
                auto p_delta = std::sqrt(s.p * s.p - 4.0 * w1 * w2);
                auto z_delta = std::sqrt(*s.z * *s.z - 4.0 * w1 * w2);
                bp1.p = ScaleComplex(s.p + p_delta, 0.5);
                bp2.p = ScaleComplex(s.p - p_delta, 0.5);
                bp1.z = ScaleComplex(*s.z + z_delta, 0.5);
                bp2.z = ScaleComplex(*s.z - z_delta, 0.5);
            }
            else {
                auto delta = std::sqrt(s.p * s.p - 4.0 * w1 * w2);
                bp1.p = ScaleComplex(s.p + delta, 0.5);
                bp2.p = ScaleComplex(s.p - delta, 0.5);
                bp1.z = 0;
            }
            protyle[i] = bp1;
            protyle[i + num_filter] = bp2;
        }
    }

    // --------------------------------------------------------------------------------
    // 离散化
    // --------------------------------------------------------------------------------
    /**
     * @brief 双线性变换, 模拟 ZPK -> 数字 ZPK
     *
     * 做 s = 2*fs*(1 - z^-1)/(1 + z^-1): 极点为 (k+p)/(k-p), k = 2*fs;
     * 有限零点为 (k+z)/(k-z), 无穷远零点变成 z = -1(Nyquist)。
     *
     * @param analog 被就地修改的模拟域零极点节
     * @param fs 采样率(Hz)
     * @note 每节的 k 按 |k-z|^2/|k-p|^2(无零点时 1/|k-p|^2)缩放, 保持频响幅度不变
     * @note 变换后每节一定带有有限零点, 可以直接交给 TfToBiquad
     */
    static void Bilinear(std::span<ZPK> analog, double fs) {
        std::complex k = 2.0 * fs;
        for (size_t i = 0; i < analog.size(); ++i) {
            const ZPK& s = analog[i];
            ZPK z;
            if (s.z) {
                z.p = (k + s.p) / (k - s.p);
                z.z = (k + *s.z) / (k - *s.z);
                z.k = s.k * std::real((k - *s.z) * (k - std::conj(*s.z)) / (k - s.p) / (k - std::conj(s.p)));
            }
            else {
                z.p = (k + s.p) / (k - s.p);
                z.z = -1;
                z.k = s.k / std::real((k - s.p) * (k - std::conj(s.p)));
            }
            analog[i] = z;
        }
    }

    /**
     * @brief 数字域 ZPK 转成双二阶系数
     *
     * 每节映射为 b0 = k, b1 = -2*k*Re(z), b2 = k*|z|^2, a1 = -2*Re(p), a2 = |p|^2,
     * 也就是 H(z) = k(1 - z*z^-1)(1 - z'*z^-1) / ((1 - p*z^-1)(1 - p'*z^-1))。
     *
     * @param digital 数字域零极点节
     * @param biquad 输出的双二阶系数, 至少 digital.size() 个
     * @pre 每节都必须有有限零点(通常先调用 Bilinear, 无穷远零点会被写成 z = -1),
     *      否则解引用空的 optional 是未定义行为
     * @warning 系数按 float 存放: 高选择性/窄陷波设计的实测频响会与设计值有偏差,
     *          偏差量随极点靠近单位圆而变大
     */
    static void TfToBiquad(std::span<ZPK> digital, std::span<BiquadCoeff> biquad) {
        assert(biquad.size() >= digital.size());

        size_t num_filter = digital.size();
        for (size_t i = 0; i < num_filter; ++i) {
            const auto& z = digital[i];
            float k = static_cast<float>(z.k);
            float b0 = static_cast<float>(k);
            float b1 = static_cast<float>(-k * 2.0 * std::real(*z.z));
            float b2 = static_cast<float>(k * std::norm(*z.z));
            float a1 = static_cast<float>(-2.0 * std::real(z.p));
            float a2 = static_cast<float>(std::norm(z.p));

            auto& coeff = biquad[i];
            coeff.a1 = a1;
            coeff.a2 = a2;
            coeff.b0 = b0;
            coeff.b1 = b1;
            coeff.b2 = b2;
        }
    }

    /**
     * @brief 数字频率的预畸变, 得到双线性变换前应使用的模拟角频率
     * @param freq 数字频率(Hz)
     * @param fs 采样率(Hz)
     * @return 模拟角频率(rad/sec), 即 2*fs*tan(pi*freq/fs)
     * @note 逆变换为 freq = fs/pi * atan(omega/(2*fs))
     */
    static double Digital2AnalogW(double freq, double fs) {
        return 2 * fs * std::tan(freq * pi / fs);
    }
};
} // namespace qwqdsp_filter
