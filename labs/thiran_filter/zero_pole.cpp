/// @brief 验证线性相位零点与 Thiran 全极点组合低通的延迟和稳定性
#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <iostream>
#include <numeric>
#include <string_view>
#include <utility>

#include "thiran.hpp"

namespace {

using qwqdsp_thiran_lab::AllpoleCoef;
using qwqdsp_thiran_lab::GroupDelay;
using qwqdsp_thiran_lab::kPi;
using qwqdsp_thiran_lab::MaxPoleRadius;
using qwqdsp_thiran_lab::Response;
using qwqdsp_thiran_lab::ZeroPoleLowpassCoef;

struct Case {
    std::size_t zeros;
    std::size_t poles;
    double tau0;
};

constexpr std::array<Case, 6> kCases{{
    {0, 4, 1.0},
    {1, 4, 1.0},
    {2, 4, 2.0},
    {4, 4, 3.0},
    {2, 6, 2.0},
    {3, 6, 3.0},
}};

int g_failures = 0;

void Check(bool ok, std::string_view what) {
    if (!ok) {
        ++g_failures;
        std::cout << std::format("FAIL: {}\n", what);
    }
}

void CheckCascade(Case const& c) {
    auto [b, a] = ZeroPoleLowpassCoef(c.zeros, c.poles, c.tau0);
    double const delta = c.tau0 - 0.5 * static_cast<double>(c.zeros);

    Check(std::abs(GroupDelay(b, a, 0.0) - c.tau0) < 1e-9,
          std::format("zeros={} poles={} tau={}: 直流群延迟不等于目标", c.zeros, c.poles, c.tau0));
    Check(std::abs(Response(b, a, 0.0) - 1.0) < 1e-12,
          std::format("zeros={} poles={} tau={}: 直流增益未归一", c.zeros, c.poles, c.tau0));
    if (c.zeros > 0) {
        Check(std::abs(Response(b, a, kPi)) < 1e-9,
              std::format("zeros={} poles={} tau={}: z=-1 处零点不为零", c.zeros, c.poles, c.tau0));
    }
    if (delta > 0.0) {
        Check(MaxPoleRadius(a) < 1.0,
              std::format("zeros={} poles={} tau={}: 极点不稳定", c.zeros, c.poles, c.tau0));
    }
    if (c.zeros == 0) {
        auto const a_ref = AllpoleCoef(c.poles, delta);
        double const sum_a_ref = std::accumulate(a_ref.begin(), a_ref.end(), 0.0);
        Check(std::abs(b[0] - sum_a_ref) < 1e-12,
              "zeros=0: 分子常量未归一为 A(1)");
    }
}

void CheckFlatness() {
    constexpr std::size_t zeros = 2;
    constexpr std::size_t poles = 4;
    constexpr double tau0 = 2.0;
    auto [b, a] = ZeroPoleLowpassCoef(zeros, poles, tau0);

    constexpr double w1 = 0.05;
    constexpr double w2 = 0.10;
    double const e1 = std::abs(GroupDelay(b, a, w1) - tau0);
    double const e2 = std::abs(GroupDelay(b, a, w2) - tau0);
    double const slope = std::log(e2 / e1) / std::log(w2 / w1);
    Check(slope > 7.5 && slope < 8.5,
          std::format("低频群延迟误差幂次 {:.3f}，期望 2*poles=8", slope));

    auto const a_ref = AllpoleCoef(poles, tau0 - 0.5 * zeros);
    std::array<double, 1> const unit{1.0};
    double const e_ref = std::abs(GroupDelay(unit, a_ref, w1) - (tau0 - 0.5 * zeros));
    Check(std::abs(e1 - e_ref) < 1e-12,
          std::format("线性相位分子改变了分母误差：{:.3e} vs {:.3e}", e1, e_ref));
}

void PrintTable() {
    std::cout << "--- z=-1 零点 + Thiran 全极点：直流延迟、极点模长、带内幅度 ---\n";
    std::cout << std::format("{:>5} {:>5} {:>6} {:>7} {:>9} {:>12} {:>12}\n",
                             "zeros", "poles", "tau0", "max|p|", "err@0.1pi", "|H|@0.5pi", "|H|@0.9pi");
    for (Case const& c : kCases) {
        auto [b, a] = ZeroPoleLowpassCoef(c.zeros, c.poles, c.tau0);
        double const err = std::abs(GroupDelay(b, a, 0.1 * kPi) - c.tau0);
        double const mag_mid = std::abs(Response(b, a, 0.5 * kPi));
        double const mag_high = std::abs(Response(b, a, 0.9 * kPi));
        std::cout << std::format("{:>5} {:>5} {:>6.1f} {:>7.4f} {:>12.3e} {:>12.3e} {:>12.3e}\n",
                                 c.zeros, c.poles, c.tau0, MaxPoleRadius(a), err, mag_mid, mag_high);
    }
}

} // namespace

int main() {
    PrintTable();
    CheckFlatness();
    for (Case const& c : kCases) {
        CheckCascade(c);
    }
    std::cout << std::format("\nfailures: {}\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
