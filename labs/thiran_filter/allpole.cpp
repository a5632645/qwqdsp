/// @brief 研究 Thiran 全极点分数延迟滤波器的性质及其与全通形式的关系
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <iostream>
#include <numeric>
#include <span>
#include <utility>
#include <vector>

#include "thiran.hpp"

namespace {

using qwqdsp_thiran_lab::AllpassCoef;
using qwqdsp_thiran_lab::AllpassNumerator;
using qwqdsp_thiran_lab::AllpoleCoef;
using qwqdsp_thiran_lab::GroupDelay;
using qwqdsp_thiran_lab::kPi;
using qwqdsp_thiran_lab::MaxPoleRadius;
using qwqdsp_thiran_lab::Response;
using qwqdsp_thiran_lab::Trapezoid;

constexpr std::size_t kMaxOrder = 6;
constexpr std::array<double, 7> kDeltas{0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0};
constexpr std::array<double, 1> kUnit{1.0};
constexpr std::array<std::pair<std::size_t, double>, 6> kMagnitudeCases{{{2, 0.5}, {2, 2.0}, {4, 0.5}, {4, 2.0}, {6, 0.5}, {6, 5.0}}};

int g_failures = 0;

void Check(bool ok, std::string_view what) {
    if (!ok) {
        ++g_failures;
        std::cout << std::format("FAIL: {}\n", what);
    }
}

/// @brief 输出系数与直流群延迟，供 thiran_check.py 逐点比对
void DumpCoefs() {
    std::cout << "--- 系数 a_k = (-1)^k C(n,k) prod (2d+i)/(2d+k+i)，被延迟的是 1/A(z) ---\n";
    for (std::size_t n = 1; n <= kMaxOrder; ++n) {
        for (double d : kDeltas) {
            auto const a = AllpoleCoef(n, d);
            for (std::size_t k = 1; k <= n; ++k) {
                std::cout << std::format("COEF P N={} D={:.17g} K={} VALUE={:.17g}\n", n, d, k, a[k]);
            }
            std::cout << std::format("TAU0 P N={} D={:.17g} TAU={:.17g}\n", n, d, GroupDelay(kUnit, a, 0.0));
        }
    }
}

/// @brief 直流群延迟、直流增益、极点模长、带内群延迟符号
void CheckProperties() {
    std::cout << "\n--- 性质表（1/A(z)，幅度不受控）---\n";
    std::cout << std::format("{:>3} {:>6} {:>11} {:>11} {:>10} {:>9} {:>8} {:>10} {:>10}\n",
                             "N", "delta", "tau_g(0)", "gain(0)", "max|p|", "stable", "tau_g(pi)", "|H(pi)|", "tau_p(0)");
    for (std::size_t n = 1; n <= kMaxOrder; ++n) {
        for (double d : kDeltas) {
            auto const a = AllpoleCoef(n, d);
            double const tau0 = GroupDelay(kUnit, a, 0.0);
            double const gain_dc = std::abs(Response(kUnit, a, 0.0));
            double const radius = MaxPoleRadius(a);
            double const tau_pi = GroupDelay(kUnit, a, kPi);
            double const gain_pi = std::abs(Response(kUnit, a, kPi));
            double const tau_p0 = qwqdsp_thiran_lab::PhaseDelay(kUnit, a, 1e-6);

            std::cout << std::format("{:>3} {:>6.2f} {:>11.6f} {:>11.6f} {:>10.4f} {:>9} {:>10.4f} {:>10.4f} {:>10.6f}\n",
                                     n, d, tau0, gain_dc, radius, radius < 1.0 ? "yes" : "no", tau_pi, gain_pi, tau_p0);

            Check(std::abs(tau0 - d) < 1e-9, std::format("n={} d={:+.2f}: 直流群延迟 != delta", n, d));
            Check(std::abs(tau_p0 - d) < 1e-4, std::format("n={} d={:+.2f}: 直流相延迟 != delta", n, d));
            Check(tau_pi < 0.0, std::format("n={} d={:+.2f}: 奈奎斯特处群延迟应为负", n, d));
            Check(gain_dc >= 1.0 - 1e-12, std::format("n={} d={:+.2f}: 直流增益应不小于 1", n, d));
            Check(std::abs(gain_dc - 1.0 / std::accumulate(a.begin(), a.end(), 0.0)) < 1e-12,
                  std::format("n={} d={:+.2f}: 直流增益 != 1/A(1)", n, d));
        }
    }
}

/// @brief 群延迟误差在直流附近的幂次（最平坦阶数）
void CheckFlatness() {
    std::cout << "\n--- 最平坦阶数：tau_g(w)-delta 随 w 的幂次 ---\n";
    std::vector<double> ws;
    std::vector<double> errs;
    for (std::size_t n = 1; n <= 4; ++n) {
        double const d = 0.5;
        auto const a = AllpoleCoef(n, d);
        ws.clear();
        errs.clear();
        for (std::size_t i = 1; i <= 32; ++i) {
            double const w = std::pow(10.0, -4.0 + 3.0 * static_cast<double>(i) / 32.0);
            double const err = std::abs(GroupDelay(kUnit, a, w) - d);
            if (err > 1e-13) {
                ws.push_back(w);
                errs.push_back(err);
            }
        }
        double const slope = qwqdsp_thiran_lab::LogLogSlope(ws, errs);
        std::cout << std::format("N={} delta=+0.50  测量幂次 {:.3f}   预期 2N={}\n", n, slope, 2 * n);
        Check(ws.size() >= 4 && std::abs(slope - 2.0 * static_cast<double>(n)) < 0.3,
              std::format("n={}: 最平坦幂次 != {}", n, 2 * n));
    }
}

/// @brief 同一分母向量下全极点与全通两种形式的定量关系
void CheckRelation() {
    std::cout << "\n--- 关系：全极点参数 delta 等价于全通参数 2*delta ---\n";
    std::cout << std::format("{:>3} {:>8} {:>14} {:>16} {:>16} {:>16}\n",
                             "N", "delta", "coef diff", "tau_pole(0)", "tau_allpass(0)", "N+2*delta");
    for (std::size_t n = 1; n <= kMaxOrder; ++n) {
        for (double d : kDeltas) {
            auto const a_pole = AllpoleCoef(n, d);
            auto const a_pass = AllpassCoef(n, 2.0 * d);
            double coef_diff = 0.0;
            for (std::size_t k = 0; k <= n; ++k) {
                coef_diff = std::max(coef_diff, std::abs(a_pole[k] - a_pass[k]));
            }
            auto const b_pass = AllpassNumerator(a_pass);
            double const tau_pole = GroupDelay(kUnit, a_pole, 0.0);
            double const tau_pass = GroupDelay(b_pass, a_pass, 0.0);

            std::cout << std::format("{:>3} {:>8.2f} {:>14.3e} {:>16.6f} {:>16.6f} {:>16.6f}\n",
                                     n, d, coef_diff, tau_pole, tau_pass, static_cast<double>(n) + 2.0 * d);
            Check(coef_diff < 1e-15, std::format("n={} d={:+.2f}: 两种公式的系数不相等", n, d));
            Check(std::abs(tau_pass - (static_cast<double>(n) + 2.0 * d)) < 1e-9,
                  std::format("n={} d={:+.2f}: 同系数的全通直流群延迟 != n+2*delta", n, d));
            Check(std::abs(tau_pass - (static_cast<double>(n) + 2.0 * tau_pole)) < 1e-9,
                  std::format("n={} d={:+.2f}: 全通与全极点的直流群延迟关系不成立", n, d));
        }
    }
}

/// @brief 幅度响应：低通、增益随 delta 与阶数上升、不受控
void CheckMagnitude() {
    std::cout << "\n--- 幅度：|1/A(e^{jw})|（低通，增益无控制手段）---\n";
    std::cout << std::format("{:>3} {:>6} {:>10} {:>10} {:>10} {:>10} {:>10}\n",
                             "N", "delta", "|H(0)|", "0.25pi", "0.5pi", "0.75pi", "|H(pi)|");
    for (auto const& [n, d] : kMagnitudeCases) {
        auto const a = AllpoleCoef(n, d);
        double const g0 = std::abs(Response(kUnit, a, 0.0));
        double const g1 = std::abs(Response(kUnit, a, 0.25 * kPi));
        double const g2 = std::abs(Response(kUnit, a, 0.5 * kPi));
        double const g3 = std::abs(Response(kUnit, a, 0.75 * kPi));
        double const g4 = std::abs(Response(kUnit, a, kPi));
        std::cout << std::format("{:>3} {:>6.2f} {:>10.4f} {:>10.4f} {:>10.4f} {:>10.4f} {:>10.4f}\n",
                                 n, d, g0, g1, g2, g3, g4);
        Check(g0 > g2 && g2 > g3 && g3 > g4, std::format("n={} d={:+.2f}: 幅度不是单调低通", n, d));
    }
}

/// @brief 极点模长随 delta 的变化，标出稳定边界
void CheckStability() {
    std::cout << "\n--- 稳定域：最大极点模长 vs delta（公式在 delta=-0.5 处奇异）---\n";
    for (std::size_t n = 1; n <= kMaxOrder; ++n) {
        std::cout << std::format("N={}:", n);
        for (int i = -9; i <= 10; ++i) {
            double const d = static_cast<double>(i) / 10.0;
            if (i == -5) {
                std::cout << std::format(" {:.1f}:奇异", d);
                continue;
            }
            auto const a = AllpoleCoef(n, d);
            double const radius = MaxPoleRadius(a);
            std::cout << std::format(" {:.1f}:{:.4f}", d, radius);
            Check((radius < 1.0) == (d > -0.5), std::format("n={} d={:+.1f}: 稳定域边界 != delta>-0.5", n, d));
        }
        std::cout << "\n";
    }
}

} // namespace

int main() {
    DumpCoefs();
    CheckProperties();
    CheckFlatness();
    CheckRelation();
    CheckMagnitude();
    CheckStability();
    std::cout << std::format("\nfailures: {}\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
