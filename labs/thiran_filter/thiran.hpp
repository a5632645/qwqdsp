#pragma once
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <span>
#include <vector>

/// @brief Thiran 分数延迟滤波器两种形式的系数公式与频域评估
/// @note 全通形式 H(z) = z^-n A(1/z)/A(z)，全极点形式 1/A(z)，两者共用同一个分母向量
namespace qwqdsp_thiran_lab {

inline constexpr double kPi = std::numbers::pi;

/// @brief 二项式系数，用连乘避免阶乘溢出
/// @param[in] n 上标
/// @param[in] k 下标
/// @return C(n,k)，k > n 时为 0
inline double Binomial(std::size_t n, std::size_t k) noexcept {
    if (k > n) {
        return 0.0;
    }
    double result = 1.0;
    for (std::size_t i = 1; i <= k; ++i) {
        result = result * static_cast<double>(n - i + 1) / static_cast<double>(i);
    }
    return result;
}

/// @brief 计算 Thiran 全通形式的分母系数 a[0..n]
/// @param[in] n 阶数
/// @param[in] delta 相对整数延迟的偏移，直流群延迟为 n+delta
/// @return 系数向量，a[0] 恒为 1
/// @note 公式在 delta = -1 处奇异，此时极零点落在单位圆上
inline std::vector<double> AllpassCoef(std::size_t n, double delta) {
    std::vector<double> a(n + 1, 0.0);
    a[0] = 1.0;
    for (std::size_t k = 1; k <= n; ++k) {
        double prod = 1.0;
        for (std::size_t i = 0; i <= n; ++i) {
            prod *= (delta + static_cast<double>(i)) / (delta + static_cast<double>(k) + static_cast<double>(i));
        }
        a[k] = (k % 2 == 0 ? 1.0 : -1.0) * Binomial(n, k) * prod;
    }
    return a;
}

/// @brief 计算 Thiran 全极点形式 1/A(z) 的分母系数 a[0..n]
/// @param[in] n 阶数
/// @param[in] delta 直流群延迟为 delta，公式中的乘积项使用 2*delta
/// @return 系数向量，a[0] 恒为 1
/// @note 公式在 delta = -0.5 处奇异，此时极点落在单位圆上
inline std::vector<double> AllpoleCoef(std::size_t n, double delta) {
    std::vector<double> a(n + 1, 0.0);
    a[0] = 1.0;
    for (std::size_t k = 1; k <= n; ++k) {
        double prod = 1.0;
        for (std::size_t i = 0; i <= n; ++i) {
            double const scaled = 2.0 * delta;
            prod *= (scaled + static_cast<double>(i)) / (scaled + static_cast<double>(k) + static_cast<double>(i));
        }
        a[k] = (k % 2 == 0 ? 1.0 : -1.0) * Binomial(n, k) * prod;
    }
    return a;
}

/// @brief 取全通形式的分子系数（分母的反序）
/// @param[in] a 分母系数
/// @return 分子系数 b[k] = a[n-k]
inline std::vector<double> AllpassNumerator(std::span<const double> a) {
    return std::vector<double>(a.rbegin(), a.rend());
}

/// @brief 计算 H(e^{jw}) = Σ b_k e^{-jkw} / Σ a_k e^{-jkw}
inline std::complex<double> Response(std::span<const double> b, std::span<const double> a, double w) {
    std::complex<double> num = 0.0;
    std::complex<double> den = 0.0;
    for (std::size_t k = 0; k < b.size(); ++k) {
        num += b[k] * std::polar(1.0, -static_cast<double>(k) * w);
    }
    for (std::size_t k = 0; k < a.size(); ++k) {
        den += a[k] * std::polar(1.0, -static_cast<double>(k) * w);
    }
    return num / den;
}

namespace detail {
/// @brief 多项式 Σ c_k e^{-jkw} 对群延迟的贡献
inline double GroupDelayOf(std::span<const double> c, double w) {
    std::complex<double> sum = 0.0;
    std::complex<double> weighted = 0.0;
    for (std::size_t k = 0; k < c.size(); ++k) {
        std::complex<double> const term = c[k] * std::polar(1.0, -static_cast<double>(k) * w);
        sum += term;
        weighted += static_cast<double>(k) * term;
    }
    return (weighted / sum).real();
}
} // namespace detail

/// @brief 计算群延迟 tau_g(w) = Im[-(H'/H)]，可直接用于全通或全极点
inline double GroupDelay(std::span<const double> b, std::span<const double> a, double w) {
    return detail::GroupDelayOf(b, w) - detail::GroupDelayOf(a, w);
}

/// @brief 计算相延迟 tau_p(w) = -arg H(w)/w
inline double PhaseDelay(std::span<const double> b, std::span<const double> a, double w) {
    return -std::arg(Response(b, a, w)) / w;
}

/// @brief 求 A(z) = Σ a_k z^-k 的极点，即 Σ a_k z^{n-k} 的根（Durand-Kerner）
/// @param[in] a 系数向量，a[0] = 1
/// @return 极点；A(z) 退化为常数时返回空向量
inline std::vector<std::complex<double>> Poles(std::span<const double> a) {
    std::size_t const n = a.size() - 1;
    std::vector<std::complex<double>> roots;
    if (n == 0) {
        return roots;
    }
    bool trivial = true;
    for (std::size_t k = 1; k <= n; ++k) {
        if (a[k] != 0.0) {
            trivial = false;
            break;
        }
    }
    if (trivial) {
        return roots;
    }

    double bound = 1.0;
    for (std::size_t k = 1; k <= n; ++k) {
        bound = std::max(bound, 1.0 + std::abs(a[k]));
    }
    roots.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        double const angle = 2.0 * kPi * (static_cast<double>(i) + 0.4) / static_cast<double>(n);
        roots[i] = std::polar(bound, angle);
    }

    auto const eval = [a, n](std::complex<double> z) {
        std::complex<double> value = a[0];
        for (std::size_t k = 1; k <= n; ++k) {
            value = value * z + a[k];
        }
        return value;
    };

    for (int iter = 0; iter < 200; ++iter) {
        double max_step = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            std::complex<double> denom = 1.0;
            for (std::size_t j = 0; j < n; ++j) {
                if (j != i) {
                    denom *= roots[i] - roots[j];
                }
            }
            std::complex<double> const step = eval(roots[i]) / denom;
            roots[i] -= step;
            max_step = std::max(max_step, std::abs(step));
        }
        if (max_step < 1e-14) {
            break;
        }
    }
    return roots;
}

/// @brief 返回最大的极点模长，小于 1 表示稳定
inline double MaxPoleRadius(std::span<const double> a) {
    double radius = 0.0;
    for (auto const& root : Poles(a)) {
        radius = std::max(radius, std::abs(root));
    }
    return radius;
}

/// @brief 计算两组正数在双对数坐标下的最小二乘斜率
/// @param[in] x 自变量
/// @param[in] y 因变量，与 x 等长
/// @return 斜率；点数少于 2 时返回 0
inline double LogLogSlope(std::span<const double> x, std::span<const double> y) {
    std::size_t const n = std::min(x.size(), y.size());
    if (n < 2) {
        return 0.0;
    }
    double sx = 0.0;
    double sy = 0.0;
    double sxy = 0.0;
    double sxx = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double const lx = std::log(x[i]);
        double const ly = std::log(y[i]);
        sx += lx;
        sy += ly;
        sxy += lx * ly;
        sxx += lx * lx;
    }
    double const count = static_cast<double>(n);
    return (count * sxy - sx * sy) / (count * sxx - sx * sx);
}

/// @brief 在 [0, wmax] 上按梯形法积分
/// @param[in] f 被积函数
/// @param[in] steps 区间数
inline double Trapezoid(auto const& f, double wmax, std::size_t steps = 2000) {
    double const dw = wmax / static_cast<double>(steps);
    double sum = 0.5 * (f(0.0) + f(wmax));
    for (std::size_t i = 1; i < steps; ++i) {
        sum += f(dw * static_cast<double>(i));
    }
    return sum * dw;
}

} // namespace qwqdsp_thiran_lab
