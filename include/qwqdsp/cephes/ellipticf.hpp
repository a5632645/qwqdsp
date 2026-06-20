#pragma once

#include <cmath>

namespace qwqdsp_cephes {

namespace detail {

/**
 * @brief 多项式求值（单精度，首项系数在数组中）
 */
static inline float polevlf(float x, const float* coef, int N) noexcept {
    float ans;
    const float* p = coef;
    ans = *p++;
    int i = N;
    do
        ans = ans * x + *p++;
    while (--i);
    return ans;
}

/**
 * @brief 多项式求值（单精度，首项系数为 1.0，省略）
 */
static inline float p1evlf(float x, const float* coef, int N) noexcept {
    float ans;
    const float* p = coef;
    ans = x + *p++;
    int i = N - 1;
    do
        ans = ans * x + *p++;
    while (--i);
    return ans;
}

} // namespace detail

/**
 * @brief 椭圆积分与雅可比椭圆函数（单精度）
 */
struct Ellipticf {
    /**
     * @brief K(m) — 第一类完全椭圆积分（单精度）
     * @details 参数 m1 = 1 - m。
     *          相对误差: IEEE 域 [0,1], 峰值 1.3e-7, RMS 3.4e-8
     */
    static float ellpk(float m1) noexcept {
        static constexpr float kP[] = {
            1.37982864606273237150E-4f, 2.28025724005875567385E-3f, 7.97404013220415179367E-3f,
            9.85821379021226008714E-3f, 6.87489687449949877925E-3f, 6.18901033637687613229E-3f,
            8.79078273952743772254E-3f, 1.49380448916805252718E-2f, 3.08851465246711995998E-2f,
            9.65735902811690126535E-2f, 1.38629436111989062502E0f,
        };
        static constexpr float kQ[] = {
            2.94078955048598507511E-5f, 9.14184723865917226571E-4f, 5.94058303753167793257E-3f,
            1.54850516649762399335E-2f, 2.39089602715924892727E-2f, 3.01204715227604046988E-2f,
            3.73774314173823228969E-2f, 4.88280347570998239232E-2f, 7.03124996963957469739E-2f,
            1.24999999999870820058E-1f, 4.99999999999999999821E-1f,
        };
        static constexpr float kC1 = 1.3862943611198906188E0f; // log(4)

        if (m1 < 0.0f || m1 > 1.0f)
            return 0.0f;
        if (m1 > 5.9604644775390625E-8f) // MACHEPF
            return detail::polevlf(m1, kP, 10) - std::log(m1) * detail::polevlf(m1, kQ, 10);
        if (m1 == 0.0f)
            return 3.4028234663852885981e38f; // MAXNUMF
        return kC1 - 0.5f * std::log(m1);
    }

    /**
     * @brief E(m) — 第二类完全椭圆积分（单精度）
     * @details 参数 m1 = 1 - m。使用多项式逼近 P(x) - x·log(x)·Q(x)。
     *          相对误差: IEEE 域 [0,1], 峰值 1.1e-7, RMS 3.9e-8
     */
    static float ellpe(float m1) noexcept {
        static constexpr float kP[] = {
            1.53552577301013293365E-4f, 2.50888492163602060990E-3f, 8.68786816565889628429E-3f,
            1.07350949056076193403E-2f, 7.77395492516787092951E-3f, 7.58395289413514708519E-3f,
            1.15688436810574127319E-2f, 2.18317996015557253103E-2f, 5.68051945617860553470E-2f,
            4.43147180560990850618E-1f, 1.00000000000000000299E0f,
        };
        static constexpr float kQ[] = {
            3.27954898576485872656E-5f, 1.00962792679356715133E-3f, 6.50609489976927491433E-3f,
            1.68862163993311317300E-2f, 2.61769742454493659583E-2f, 3.34833904888224918614E-2f,
            4.27180926518931511717E-2f, 5.85936634471101055642E-2f, 9.37499997197644278445E-2f,
            2.49999999999888314361E-1f,
        };

        if (m1 <= 0.0f || m1 > 1.0f) {
            if (m1 == 0.0f)
                return 1.0f;
            return 0.0f;
        }
        return detail::polevlf(m1, kP, 10) - std::log(m1) * (m1 * detail::polevlf(m1, kQ, 9));
    }

    /**
     * @brief F(φ\\m) — 第一类不完全椭圆积分（单精度）
     * @details 使用算术-几何平均（AGM）算法。
     *          相对误差: IEEE 域 φ∈[0,2], m∈[0,1], 峰值 2.9e-7, RMS 5.8e-8
     */
    static float ellik(float phi, float m) noexcept {
        static constexpr float kPIF = 3.141592653589793238f;
        static constexpr float kPIO2F = 1.5707963267948966192f;
        static constexpr float kMACHEPF = 5.9604644775390625E-8f;

        if (m == 0.0f)
            return phi;

        int sign = 0;
        if (phi < 0.0f) {
            phi = -phi;
            sign = -1;
        }

        float a = 1.0f;
        float b = 1.0f - m;
        if (b == 0.0f)
            return std::log(std::tan(0.5f * (kPIO2F + phi)));
        b = std::sqrt(b);
        float c = std::sqrt(m);
        int d = 1;
        float t = std::tan(phi);
        int mod = (int)((phi + kPIO2F) / kPIF);

        while (std::abs(c / a) > kMACHEPF) {
            float temp = b / a;
            phi = phi + std::atan(t * temp) + mod * kPIF;
            mod = (int)((phi + kPIO2F) / kPIF);
            t = t * (1.0f + temp) / (1.0f - temp * t * t);
            c = (a - b) / 2.0f;
            temp = std::sqrt(a * b);
            a = (a + b) / 2.0f;
            b = temp;
            d += d;
        }

        float result = (std::atan(t) + mod * kPIF) / (d * a);
        if (sign < 0)
            result = -result;
        return result;
    }

    /**
     * @brief E(φ\\m) — 第二类不完全椭圆积分（单精度）
     * @details 使用算术-几何平均（AGM）算法。
     *          相对误差: IEEE 域 φ∈[0,2], m∈[0,1], 峰值 4.5e-7, RMS 7.4e-8
     */
    static float ellie(float phi, float m) noexcept {
        static constexpr float kPIF = 3.141592653589793238f;
        static constexpr float kPIO2F = 1.5707963267948966192f;
        static constexpr float kMACHEPF = 5.9604644775390625E-8f;

        if (m == 0.0f)
            return phi;
        if (m == 1.0f)
            return std::sin(phi);

        float lphi = phi;
        if (lphi < 0.0f)
            lphi = -lphi;

        float a = 1.0f;
        float b = 1.0f - m;
        b = std::sqrt(b);
        float c = std::sqrt(m);
        int d = 1;
        float e = 0.0f;
        float t = std::tan(lphi);
        int mod = (int)((lphi + kPIO2F) / kPIF);

        while (std::abs(c / a) > kMACHEPF) {
            float temp = b / a;
            lphi = lphi + std::atan(t * temp) + mod * kPIF;
            mod = (int)((lphi + kPIO2F) / kPIF);
            t = t * (1.0f + temp) / (1.0f - temp * t * t);
            c = (a - b) / 2.0f;
            temp = std::sqrt(a * b);
            a = (a + b) / 2.0f;
            b = temp;
            d += d;
            e += c * std::sin(lphi);
        }

        b = 1.0f - m; // restore b = complementary modulus
        float result = ellpe(b) / ellpk(b);
        result *= (std::atan(t) + mod * kPIF) / (d * a);
        result += e;

        if (phi < 0.0f)
            result = -result;
        return result;
    }

    /**
     * @brief Jacobi 椭圆函数 sn(u|m), cn(u|m), dn(u|m)（单精度）
     * @details 使用算术-几何平均（AGM）算法。当 m 接近 0 或 1 时有特例处理。
     *          sn 绝对误差: IEEE, 峰值 1.7e-6, RMS 2.2e-7
     *          cn 绝对误差: IEEE, 峰值 1.6e-6, RMS 2.2e-7
     *          dn 绝对误差: IEEE, 峰值 1.4e-3, RMS 1.9e-5
     *          phi 相对误差: IEEE, 峰值 3.9e-7, RMS 6.7e-8
     * @param [in] u  实参数
     * @param [in] m  参数（0 ≤ m ≤ 1）
     * @param [out] sn  sn(u|m)
     * @param [out] cn  cn(u|m)
     * @param [out] dn  dn(u|m)
     * @param [out] ph  幅角 φ
     * @return true 成功, false 域错误
     */
    static bool ellpjf(float u, float m, float* sn, float* cn, float* dn, float* ph) noexcept {
        static constexpr float kPIO2F = 1.5707963267948966192f;
        static constexpr float kMACHEPF = 5.9604644775390625E-8f;

        if (m < 0.0f || m > 1.0f) {
            *sn = 0.0f;
            *cn = 0.0f;
            *dn = 0.0f;
            *ph = 0.0f;
            return false;
        }

        if (m < 1.0e-5f) {
            float t = std::sin(u);
            float b = std::cos(u);
            float ai = 0.25f * m * (u - t * b);
            *sn = t - ai * b;
            *cn = b + ai * t;
            *ph = u - ai;
            *dn = 1.0f - 0.5f * m * t * t;
            return true;
        }

        if (m >= 0.99999f) {
            float ai = 0.25f * (1.0f - m);
            float b = std::cosh(u);
            float t = std::tanh(u);
            float phi = 1.0f / b;
            float twon = b * std::sinh(u);
            *sn = t + ai * (twon - u) / (b * b);
            *ph = 2.0f * std::atan(std::exp(u)) - kPIO2F + ai * (twon - u) / b;
            ai *= t * phi;
            *cn = phi - ai * (twon - u);
            *dn = phi + ai * (twon + u);
            return true;
        }

        /* AGM scale */
        float a[10], c[10];
        a[0] = 1.0f;
        float b = std::sqrt(1.0f - m);
        c[0] = std::sqrt(m);
        float twon = 1.0f;
        int i = 0;

        while (std::abs(c[i] / a[i]) > kMACHEPF) {
            if (i > 8)
                break;
            float ai = a[i];
            ++i;
            c[i] = (ai - b) / 2.0f;
            float t = std::sqrt(ai * b);
            a[i] = (ai + b) / 2.0f;
            b = t;
            twon += twon;
        }

        /* backward recurrence */
        float phi = twon * a[i] * u;
        float phi_prev = phi;
        do {
            float t = c[i] * std::sin(phi) / a[i];
            phi_prev = phi;
            phi = (std::asin(t) + phi) / 2.0f;
        } while (--i);

        *sn = std::sin(phi);
        float t = std::cos(phi);
        *cn = t;
        *dn = t / std::cos(phi - phi_prev);
        *ph = phi;
        return true;
    }
};

} // namespace qwqdsp_cephes
