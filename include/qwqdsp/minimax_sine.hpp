#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>

namespace qwqdsp::minimax_sine {

// ------------------------------------------------------------
// 内部工具
// ------------------------------------------------------------
namespace detail {
/**
 * @brief 把任意相位分数折到 [-1/4, 1/4]，对应 sin(2πx) 的一个单调四分周期
 *
 * 文档折叠 f(x) = max(min(x, 1/2-x), -1/2-x)，先折到 [-1/2, 1/2]（减去最近的
 * 整数）再折到一个四分周期。这里采用无分支的通用形式：仅一次 floor + 两对
 * min/max。相比按输入限定 [0,1) 的分支版，本函数接受任意相位（可为负或 ≥1），
 * 且在编译器自动向量化时因无分支而吞吐与分支版持平（实测相等），故保留通用性。
 *
 * @param x 任意相位分数（仅依赖周期性，无需预先归一化）
 * @return [-1/4, 1/4] 区间内的折纸结果
 */
template <class T>
constexpr T FoldFraction(T x) noexcept {
    T const xr = x - std::floor(x + T(0.5));
    return std::max(std::min(xr, T(0.5) - xr), T(-0.5) - xr);
}

/**
 * @brief Horner 求值 sin(2πx) 多项式（只用奇次幂，系数从常数项到高次）
 *
 * 文档写法：x1*(c0 + x2*(c1 + x2*(c2 + ...)))，x1 = x，x2 = x*x。
 *
 * @tparam T 求值类型（float：尾数约 24 位；double：接近文档声明的 36 位精度）
 * @param c 系数数组（c[0] 为常数项，渐次到最高次）
 * @param x 已折到 [0,1/4] 的自变量（x1 = x）
 */
template <class T, std::size_t N>
T Horner(const std::array<double, N>& c, T x) noexcept {
    T const x2 = x * x;
    T r = static_cast<T>(c[N - 1]);
    for (std::size_t i = N - 1; i-- > 0;)
        r = r * x2 + static_cast<T>(c[i]);
    return x * r;
}
} // namespace detail

// ------------------------------------------------------------
// 系数表（sin(2πx)，仅奇次幂，最小化 [0,1/4] 上的最大相对误差）
// 来源：Lasse Schlör, "Fast MiniMax Polynomial Approximations of Sine and Cosine"
// @ref https://publik-void.github.io/sin-cos-approximations/
// 相对误差声明：degree 7 ≈ 9.39e-7，degree 9 ≈ 5.31e-9，degree 13 ≈ 6.18e-14。
// ------------------------------------------------------------

/** @brief degree 7 系数（|误差| ≤ 9.4e-7） */
inline constexpr std::array kSin7Coeffs{
    6.28317940663383263695539184084148414,
    -41.3389424984438475211600845004648977,
    81.3953586300215790478948871165643136,
    -71.4746942820704835483055527883527435,
};

/** @brief degree 9 系数（|误差| ≤ 5.4e-9） */
inline constexpr std::array kSin9Coeffs{
    6.28318527379078585274731929079414949,
    -41.3416774783915252855640244027643612,
    81.6022312427274226421465134076212909,
    -76.5749921819992128192000934020817094,
    39.7109181438058471453004860893416233,
};

/** @brief degree 13 系数（|误差| ≤ 6.3e-14） */
inline constexpr std::array kSin13Coeffs{
    6.28318530717919415440631052356951227,
    -41.3417022398184912491504586563309009,
    81.6052491334177909789178729942153114,
    -76.7058464941280158505651312164454235,
    42.0581028415940046209613080938107769,
    -15.0810317173017800774891418165142071,
    3.66346472229432872653352143098494556,
};

// ------------------------------------------------------------
// 对外模板函数：输入相位分数 x ∈ [0,1)，输出 sin(2πx) ∈ [-1,1]
// ------------------------------------------------------------

/**
 * @brief sin(2πx)，degree 7，最大误差 ≈ 9.4e-7
 * @tparam T 类型（float 或 double）
 * @param x 相位分数（任意实数，内部按周期折纸；[0,1) 为一整周期）
 * @return [-1, 1]
 */
template <class T>
constexpr T Sin7(T x) noexcept {
    T const f = detail::FoldFraction(static_cast<T>(x));
    return detail::Horner(kSin7Coeffs, f);
}

/**
 * @brief sin(2πx)，degree 9，最大误差 ≈ 5.3e-9（T=double 时）
 * @tparam T 类型（float：~1e-7；double：~5e-9）
 * @param x 相位分数（任意实数，内部按周期折纸；[0,1) 为一整周期）
 * @return [-1, 1]
 */
template <class T>
constexpr T Sin9(T x) noexcept {
    T const f = detail::FoldFraction(static_cast<T>(x));
    return detail::Horner(kSin9Coeffs, f);
}

/**
 * @brief sin(2πx)，degree 13，最大误差 ≈ 6.2e-14（T=double 时）
 * @tparam T 类型（float：~1e-7；double：~6e-14）
 * @param x 相位分数（任意实数，内部按周期折纸；[0,1) 为一整周期）
 * @return [-1, 1]
 */
template <class T>
constexpr T Sin13(T x) noexcept {
    T const f = detail::FoldFraction(static_cast<T>(x));
    return detail::Horner(kSin13Coeffs, f);
}

} // namespace qwqdsp::minimax_sine
