#pragma once

#include <cmath>
#include <limits>

namespace qwqdsp_cephes {

namespace detail {

/**
 * @brief 多项式求值（首项系数在数组中）
 */
static inline double polevl(double x, const double* coef, int N) noexcept {
    double ans;
    const double* p = coef;
    ans = *p++;
    int i = N;
    do
        ans = ans * x + *p++;
    while (--i);
    return ans;
}

/**
 * @brief 多项式求值（首项系数为 1.0，省略）
 */
static inline double p1evl(double x, const double* coef, int N) noexcept {
    double ans;
    const double* p = coef;
    ans = x + *p++;
    int i = N - 1;
    do
        ans = ans * x + *p++;
    while (--i);
    return ans;
}

} // namespace detail

/**
 * @brief 椭圆积分与雅可比椭圆函数（双精度）
 */
struct Elliptic {
    /**
     * @brief K(m) — 第一类完全椭圆积分
     * @details 参数 m1 = 1 - m，将 m=1 处的对数奇点移至原点以保证精度。
     *          相对误差: IEEE 域 [0,1], 峰值 2.5e-16, RMS 6.8e-17
     */
    static double ellpk(double m1) noexcept {
        static constexpr double kP[] = {
            1.37982864606273237150E-4, 2.28025724005875567385E-3, 7.97404013220415179367E-3, 9.85821379021226008714E-3,
            6.87489687449949877925E-3, 6.18901033637687613229E-3, 8.79078273952743772254E-3, 1.49380448916805252718E-2,
            3.08851465246711995998E-2, 9.65735902811690126535E-2, 1.38629436111989062502E0,
        };
        static constexpr double kQ[] = {
            2.94078955048598507511E-5, 9.14184723865917226571E-4, 5.94058303753167793257E-3, 1.54850516649762399335E-2,
            2.39089602715924892727E-2, 3.01204715227604046988E-2, 3.73774314173823228969E-2, 4.88280347570998239232E-2,
            7.03124996963957469739E-2, 1.24999999999870820058E-1, 4.99999999999999999821E-1,
        };
        static constexpr double kC1 = 1.3862943611198906188E0; // log(4)

        if (m1 < 0.0 || m1 > 1.0)
            return 0.0;
        if (m1 > 2.2204460492503131E-16) // MACHEP
            return detail::polevl(m1, kP, 10) - std::log(m1) * detail::polevl(m1, kQ, 10);
        if (m1 == 0.0)
            return std::numeric_limits<double>::infinity(); // MAXNUM
        return kC1 - 0.5 * std::log(m1);
    }

    /**
     * @brief E(m) — 第二类完全椭圆积分
     * @details 参数 m1 = 1 - m。使用多项式逼近 P(x) - x·log(x)·Q(x)。
     *          相对误差: IEEE 域 [0,1], 峰值 2.1e-16, RMS 7.3e-17
     */
    static double ellpe(double m1) noexcept {
        static constexpr double kP[] = {
            1.53552577301013293365E-4, 2.50888492163602060990E-3, 8.68786816565889628429E-3, 1.07350949056076193403E-2,
            7.77395492516787092951E-3, 7.58395289413514708519E-3, 1.15688436810574127319E-2, 2.18317996015557253103E-2,
            5.68051945617860553470E-2, 4.43147180560990850618E-1, 1.00000000000000000299E0,
        };
        static constexpr double kQ[] = {
            3.27954898576485872656E-5, 1.00962792679356715133E-3, 6.50609489976927491433E-3, 1.68862163993311317300E-2,
            2.61769742454493659583E-2, 3.34833904888224918614E-2, 4.27180926518931511717E-2, 5.85936634471101055642E-2,
            9.37499997197644278445E-2, 2.49999999999888314361E-1,
        };

        if (m1 <= 0.0 || m1 > 1.0) {
            if (m1 == 0.0)
                return 1.0;
            return 0.0;
        }
        return detail::polevl(m1, kP, 10) - std::log(m1) * (m1 * detail::polevl(m1, kQ, 9));
    }

    /**
     * @brief F(φ\\m) — 第一类不完全椭圆积分
     * @details 使用算术-几何平均（AGM）算法。
     *          相对误差: IEEE 域 φ∈[-10,10], m∈[0,1], 峰值 7.4e-16, RMS 1.0e-16
     */
    static double ellik(double phi, double m) noexcept {
        static constexpr double kPIO2 = 1.57079632679489661923;
        static constexpr double kPI = 3.14159265358979323846;
        static constexpr double kMACHEP = 2.2204460492503131E-16;

        if (m == 0.0)
            return phi;
        double a = 1.0 - m;
        if (a == 0.0) {
            if (std::abs(phi) >= kPIO2)
                return std::numeric_limits<double>::infinity(); // MAXNUM
            return std::log(std::tan((kPIO2 + phi) / 2.0));
        }

        int npio2 = (int)std::floor(phi / kPIO2);
        if (npio2 & 1)
            npio2 += 1;
        double K = 0.0;
        if (npio2)
            K = ellpk(a);
        phi = phi - npio2 * kPIO2;

        int sign = 0;
        if (phi < 0.0) {
            phi = -phi;
            sign = -1;
        }

        double b = std::sqrt(a);
        double t = std::tan(phi);

        if (std::abs(t) > 10.0) {
            double e = 1.0 / (b * t);
            if (std::abs(e) < 10.0) {
                e = std::atan(e);
                if (npio2 == 0)
                    K = ellpk(a);
                double temp = K - ellik(e, m);
                if (sign < 0)
                    temp = -temp;
                temp += npio2 * K;
                return temp;
            }
        }

        double c = std::sqrt(m);
        a = 1.0;
        int d = 1;
        int mod = 0;

        while (std::abs(c / a) > kMACHEP) {
            double temp = b / a;
            phi = phi + std::atan(t * temp) + mod * kPI;
            mod = (int)((phi + kPIO2) / kPI);
            t = t * (1.0 + temp) / (1.0 - temp * t * t);
            c = (a - b) / 2.0;
            temp = std::sqrt(a * b);
            a = (a + b) / 2.0;
            b = temp;
            d += d;
        }

        double temp = (std::atan(t) + mod * kPI) / (d * a);
        if (sign < 0)
            temp = -temp;
        temp += npio2 * K;
        return temp;
    }

    /**
     * @brief E(φ\\m) — 第二类不完全椭圆积分
     * @details 使用算术-几何平均（AGM）算法。
     *          相对误差: IEEE 域 φ∈[-10,10], m∈[0,1], 峰值 3.3e-15, RMS 1.4e-16
     */
    static double ellie(double phi, double m) noexcept {
        static constexpr double kPIO2 = 1.57079632679489661923;
        static constexpr double kPI = 3.14159265358979323846;
        static constexpr double kMACHEP = 2.2204460492503131E-16;

        if (m == 0.0)
            return phi;

        double lphi = phi;
        int npio2 = (int)std::floor(lphi / kPIO2);
        if (npio2 & 1)
            npio2 += 1;
        lphi = lphi - npio2 * kPIO2;

        int sign = 1;
        if (lphi < 0.0) {
            lphi = -lphi;
            sign = -1;
        }

        double a = 1.0 - m;
        double E = ellpe(a);
        if (a == 0.0) {
            double temp = std::sin(lphi);
            if (sign < 0)
                temp = -temp;
            temp += npio2 * E;
            return temp;
        }

        double t = std::tan(lphi);
        double b = std::sqrt(a);

        if (std::abs(t) > 10.0) {
            double e = 1.0 / (b * t);
            if (std::abs(e) < 10.0) {
                e = std::atan(e);
                double temp = E + m * std::sin(lphi) * std::sin(e) - ellie(e, m);
                if (sign < 0)
                    temp = -temp;
                temp += npio2 * E;
                return temp;
            }
        }

        double c = std::sqrt(m);
        a = 1.0;
        int d = 1;
        double e = 0.0;
        int mod = 0;

        while (std::abs(c / a) > kMACHEP) {
            double temp = b / a;
            lphi = lphi + std::atan(t * temp) + mod * kPI;
            mod = (int)((lphi + kPIO2) / kPI);
            t = t * (1.0 + temp) / (1.0 - temp * t * t);
            c = (a - b) / 2.0;
            temp = std::sqrt(a * b);
            a = (a + b) / 2.0;
            b = temp;
            d += d;
            e += c * std::sin(lphi);
        }

        double temp = E / ellpk(1.0 - m);
        temp *= (std::atan(t) + mod * kPI) / (d * a);
        temp += e;

        if (sign < 0)
            temp = -temp;
        temp += npio2 * E;
        return temp;
    }

    /**
     * @brief Jacobi 椭圆函数 sn(u|m), cn(u|m), dn(u|m)
     * @details 使用算术-几何平均（AGM）算法计算。当 m 接近 0 或 1 时有特例处理。
     *          sn 绝对误差: IEEE, 峰值 4.1e-15, RMS 4.6e-16
     *          cn 绝对误差: IEEE, 峰值 3.6e-15, RMS 4.4e-16
     *          dn 绝对误差: IEEE, 峰值 1.3e-12, RMS 1.8e-14
     *          phi 相对误差: IEEE, 峰值 9.2e-16, RMS 1.4e-16
     * @param [in] u  实参数
     * @param [in] m  参数（0 ≤ m ≤ 1）
     * @param [out] sn  sn(u|m)
     * @param [out] cn  cn(u|m)
     * @param [out] dn  dn(u|m)
     * @param [out] ph  幅角 φ
     * @return true 成功, false 域错误
     */
    static bool ellpj(double u, double m, double* sn, double* cn, double* dn, double* ph) noexcept {
        static constexpr double kPIO2 = 1.57079632679489661923;
        static constexpr double kMACHEP = 2.2204460492503131E-16;

        if (m < 0.0 || m > 1.0) {
            *sn = 0.0;
            *cn = 0.0;
            *dn = 0.0;
            *ph = 0.0;
            return false;
        }

        if (m < 1.0e-9) {
            double t = std::sin(u);
            double b = std::cos(u);
            double ai = 0.25 * m * (u - t * b);
            *sn = t - ai * b;
            *cn = b + ai * t;
            *ph = u - ai;
            *dn = 1.0 - 0.5 * m * t * t;
            return true;
        }

        if (m >= 0.9999999999) {
            double ai = 0.25 * (1.0 - m);
            double b = std::cosh(u);
            double t = std::tanh(u);
            double phi = 1.0 / b;
            double twon = b * std::sinh(u);
            *sn = t + ai * (twon - u) / (b * b);
            *ph = 2.0 * std::atan(std::exp(u)) - kPIO2 + ai * (twon - u) / b;
            ai *= t * phi;
            *cn = phi - ai * (twon - u);
            *dn = phi + ai * (twon + u);
            return true;
        }

        /* AGM scale */
        double a[9], c[9];
        a[0] = 1.0;
        double b = std::sqrt(1.0 - m);
        c[0] = std::sqrt(m);
        double twon = 1.0;
        int i = 0;

        while (std::abs(c[i] / a[i]) > kMACHEP) {
            if (i > 7) {
                goto done;
            }
            double ai = a[i];
            ++i;
            c[i] = (ai - b) / 2.0;
            double t = std::sqrt(ai * b);
            a[i] = (ai + b) / 2.0;
            b = t;
            twon *= 2.0;
        }

done:
        /* backward recurrence */
        double phi = twon * a[i] * u;
        double phi_prev = phi;
        do {
            double t = c[i] * std::sin(phi) / a[i];
            phi_prev = phi;
            phi = (std::asin(t) + phi) / 2.0;
        } while (--i);

        *sn = std::sin(phi);
        double t = std::cos(phi);
        *cn = t;
        *dn = t / std::cos(phi - phi_prev);
        *ph = phi;
        return true;
    }
};

} // namespace qwqdsp_cephes
