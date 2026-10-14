/// @brief 研究 Thiran 全通分数延迟滤波器的性质
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <iostream>
#include <vector>

#include <qwqdsp/filter/thiran_filter.hpp>

#include "thiran.hpp"

namespace {

using qwqdsp_thiran_lab::AllpassCoef;
using qwqdsp_thiran_lab::AllpassNumerator;
using qwqdsp_thiran_lab::GroupDelay;
using qwqdsp_thiran_lab::kPi;
using qwqdsp_thiran_lab::MaxPoleRadius;
using qwqdsp_thiran_lab::Response;
using qwqdsp_thiran_lab::Trapezoid;

constexpr std::size_t kMaxOrder = 6;
constexpr std::array<double, 8> kDeltas{-0.9, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0};

int g_failures = 0;

void Check(bool ok, std::string_view what) {
    if (!ok) {
        ++g_failures;
        std::cout << std::format("FAIL: {}\n", what);
    }
}

std::vector<double> const kScanW = [] {
    std::vector<double> w(1025);
    for (std::size_t i = 0; i < w.size(); ++i) {
        w[i] = kPi * static_cast<double>(i) / static_cast<double>(w.size() - 1);
    }
    return w;
}();

/// @brief 输出系数与直流群延迟，供 thiran_check.py 逐点比对
void DumpCoefs() {
    std::cout << "--- 系数 a_k = (-1)^k C(n,k) prod (d+i)/(d+k+i)，分子取分母反序 ---\n";
    for (std::size_t n = 1; n <= kMaxOrder; ++n) {
        for (double d : kDeltas) {
            auto const a = AllpassCoef(n, d);
            auto const b = AllpassNumerator(a);
            for (std::size_t k = 1; k <= n; ++k) {
                std::cout << std::format("COEF AP N={} D={:.17g} K={} VALUE={:.17g}\n", n, d, k, a[k]);
            }
            std::cout << std::format("TAU0 AP N={} D={:.17g} TAU={:.17g}\n", n, d, GroupDelay(b, a, 0.0));
        }
    }
}

/// @brief 直流群延迟、幅度平坦度、平均群延迟、极点模长
void CheckProperties() {
    std::cout << "\n--- 性质表（分子为分母反序，故幅度恒为 1）---\n";
    std::cout << std::format("{:>3} {:>7} {:>11} {:>11} {:>11} {:>10} {:>9} {:>7}\n",
                             "N", "delta", "tau_g(0)", "N+delta", "|H|-1max", "mean_tau", "max|p|", "stable");
    for (std::size_t n = 1; n <= kMaxOrder; ++n) {
        for (double d : kDeltas) {
            auto const a = AllpassCoef(n, d);
            auto const b = AllpassNumerator(a);
            double const tau0 = GroupDelay(b, a, 0.0);
            double const delay = static_cast<double>(n) + d;

            double max_dev = 0.0;
            for (double w : kScanW) {
                max_dev = std::max(max_dev, std::abs(std::abs(Response(b, a, w)) - 1.0));
            }
            double const mean_tau = Trapezoid([&](double w) { return GroupDelay(b, a, w); }, kPi, 4096) / kPi;
            double const radius = MaxPoleRadius(a);

            std::cout << std::format("{:>3} {:>+7.2f} {:>11.6f} {:>11.6f} {:>11.3e} {:>10.6f} {:>9.4f} {:>7}\n",
                                     n, d, tau0, delay, max_dev, mean_tau, radius, radius < 1.0 ? "yes" : "no");

            Check(std::abs(tau0 - delay) < 1e-9, std::format("n={} d={:+.2f}: 直流群延迟 != n+delta", n, d));
            Check(max_dev < 1e-13, std::format("n={} d={:+.2f}: |H| 偏离 1 超过 1e-13", n, d));
            Check(std::abs(mean_tau - static_cast<double>(n)) < 1e-9,
                  std::format("n={} d={:+.2f}: 平均群延迟 != n", n, d));
            Check((radius < 1.0) == (d > -1.0), std::format("n={} d={:+.2f}: 稳定域边界 != delta>-1", n, d));
        }
    }
}

/// @brief 群延迟误差在直流附近的幂次（最平坦阶数）
void CheckFlatness() {
    std::cout << "\n--- 最平坦阶数：tau_g(w)-(N+d) 随 w 的幂次 ---\n";
    std::vector<double> ws;
    std::vector<double> errs;
    for (std::size_t n = 1; n <= 3; ++n) {
        double const d = 0.5;
        auto const a = AllpassCoef(n, d);
        auto const b = AllpassNumerator(a);
        ws.clear();
        errs.clear();
        for (std::size_t i = 1; i <= 32; ++i) {
            double const w = std::pow(10.0, -4.0 + 3.0 * static_cast<double>(i) / 32.0);
            double const err = std::abs(GroupDelay(b, a, w) - (static_cast<double>(n) + d));
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

/// @brief 群延迟偏差：符号恒定、带内均值、正负 delta 的高频形态差异
void CheckDelayDeviation() {
    std::cout << "\n--- 群延迟偏差 tau_g(w)-D（D=N+delta）---\n";
    std::cout << std::format("{:>3} {:>7} {:>12} {:>12} {:>12} {:>13} {:>9}\n",
                             "N", "delta", "dev(0.5pi)", "dev(0.75pi)", "dev(pi)", "mean|dev|", "sign");
    for (std::size_t n = 1; n <= kMaxOrder; ++n) {
        for (double d : std::array<double, 4>{-0.5, -0.25, 0.25, 0.5}) {
            auto const a = AllpassCoef(n, d);
            auto const b = AllpassNumerator(a);
            double const delay = static_cast<double>(n) + d;
            double const dev_a = GroupDelay(b, a, 0.5 * kPi) - delay;
            double const dev_b = GroupDelay(b, a, 0.75 * kPi) - delay;
            double const dev_c = GroupDelay(b, a, kPi) - delay;
            double const mean_abs = Trapezoid([&](double w) { return std::abs(GroupDelay(b, a, w) - delay); }, kPi, 4096) / kPi;

            double lo = 0.0;
            double hi = 0.0;
            for (double w : kScanW) {
                double const dev = GroupDelay(b, a, w) - delay;
                lo = std::min(lo, dev);
                hi = std::max(hi, dev);
            }
            std::cout << std::format("{:>3} {:>+7.2f} {:>12.5f} {:>12.5f} {:>12.5f} {:>13.6f} {:>9}\n",
                                     n, d, dev_a, dev_b, dev_c, mean_abs, lo > -1e-9 ? ">=" : "<=");
            Check(lo > -1e-9 || hi < 1e-9, std::format("n={} d={:+.2f}: 偏差在带内变号", n, d));
            Check(std::abs(mean_abs - std::abs(d)) < 1e-6,
                  std::format("n={} d={:+.2f}: 平均绝对偏差 != |delta|", n, d));
        }
    }
}

/// @brief 极点模长随 delta 的变化，标出稳定边界
void CheckStability() {
    std::cout << "\n--- 稳定域：最大极点模长 vs delta（公式在 delta=-1 处奇异）---\n";
    for (std::size_t n = 1; n <= kMaxOrder; ++n) {
        std::cout << std::format("N={}:", n);
        for (int i = -15; i <= 5; ++i) {
            double const d = static_cast<double>(i) / 10.0;
            if (i == -10) {
                std::cout << std::format(" {:.1f}:奇异", d);
                continue;
            }
            auto const a = AllpassCoef(n, d);
            double const radius = MaxPoleRadius(a);
            std::cout << std::format(" {:.1f}:{:.4f}", d, radius);
            Check((radius < 1.0) == (d > -1.0), std::format("n={} d={:+.1f}: 稳定域边界 != delta>-1", n, d));
        }
        std::cout << "\n";
    }
}

/// @brief 交付的 ThiranFilter 与本文件的双精度全通模型逐点比对（黑盒：只看脉冲响应）
template <std::size_t N>
void CheckLibrary(float delay) {
    qwqdsp::ThiranFilter<N> filter;
    filter.Reset();
    std::size_t const int_delay = filter.Make(delay);

    double const shifted = static_cast<double>(delay) + 0.5;
    double const frac = shifted - std::floor(shifted);
    double const delta = frac - 0.5;
    int expected_int = static_cast<int>(std::lround(shifted - 0.5 - delta - static_cast<double>(N)));
    if (expected_int < 0) {
        expected_int = 0;
    }

    auto const a = AllpassCoef(N, delta);
    auto const b = AllpassNumerator(a);

    constexpr std::size_t kLen = 256;
    std::vector<double> x_history(N + 1, 0.0);
    std::vector<double> y_history(N + 1, 0.0);
    double max_diff = 0.0;
    for (std::size_t t = 0; t < kLen; ++t) {
        double const x = t == 0 ? 1.0 : 0.0;
        for (std::size_t k = N; k >= 1; --k) {
            x_history[k] = x_history[k - 1];
            y_history[k] = y_history[k - 1];
        }
        x_history[0] = x;
        double y = 0.0;
        for (std::size_t k = 0; k <= N; ++k) {
            y += b[k] * x_history[k];
        }
        for (std::size_t k = 1; k <= N; ++k) {
            y -= a[k] * y_history[k];
        }
        y_history[0] = y;
        double const got = static_cast<double>(filter.Tick(static_cast<float>(x)));
        max_diff = std::max(max_diff, std::abs(got - y));
    }

    double const tau_total = static_cast<double>(int_delay) + GroupDelay(b, a, 0.0);
    std::cout << std::format("LIB N={} delay={:>5.2f} int={} delta={:+.6f} total={:.5f} maxdiff={:.3e}\n",
                             N, delay, int_delay, delta, tau_total, max_diff);

    Check(int_delay == static_cast<std::size_t>(expected_int),
          std::format("N={} delay={}: Make 返回的整数延迟 != 拆分规则", N, delay));
    Check(delta >= -0.5 && delta < 0.5, std::format("N={} delay={}: delta 未归一到 [-0.5,0.5)", N, delay));
    if (static_cast<double>(delay) >= static_cast<double>(N) - 0.5) {
        Check(std::abs(tau_total - static_cast<double>(delay)) < 1e-4,
              std::format("N={} delay={}: 总延迟 != 请求延迟", N, delay));
    }
    else {
        Check(int_delay == 0, std::format("N={} delay={}: 低于 N-0.5 时整数延迟应被夹到 0", N, delay));
    }
    Check(max_diff < 1e-5, std::format("N={} delay={}: 脉冲响应与双精度参考差 {} 超过 1e-5", N, delay, max_diff));
}

/// @brief 覆盖请求延迟的不同小数位置（含落在 0.5 两侧与低于 N-0.5 的情形）
void CheckLibraryAll() {
    std::cout << "\n--- 与本库 ThiranFilter 的黑盒比对 ---\n";
    CheckLibrary<1>(0.25f);
    CheckLibrary<1>(0.80f);
    CheckLibrary<1>(3.40f);
    CheckLibrary<2>(0.60f);
    CheckLibrary<2>(1.50f);
    CheckLibrary<2>(4.25f);
    CheckLibrary<3>(1.40f);
    CheckLibrary<3>(2.55f);
    CheckLibrary<3>(6.10f);
    CheckLibrary<4>(0.30f);
    CheckLibrary<4>(2.50f);
    CheckLibrary<4>(3.60f);
    CheckLibrary<4>(5.75f);
    CheckLibrary<4>(9.40f);
}

} // namespace

int main() {
    DumpCoefs();
    CheckProperties();
    CheckFlatness();
    CheckDelayDeviation();
    CheckStability();
    CheckLibraryAll();
    std::cout << std::format("\nfailures: {}\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
