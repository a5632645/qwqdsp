// ------------------------------------------------------------
// minimax_sine 多项式正弦/余弦近似精度校验
// ------------------------------------------------------------
#include <cmath>
#include <cstdio>
#include <numbers>
#include <vector>

#include <qwqdsp/minimax_sine.hpp>

namespace mm = qwqdsp::minimax_sine;

// 各阶数声明的最大绝对误差上限（double 求值，文档值）
static constexpr double kSin7MaxErr = 1.5e-6;
static constexpr double kSin9MaxErr = 1.0e-8;
static constexpr double kSin13MaxErr = 6.0e-13;

template <class F>
static double MaxAbsErr(F approx, int n, double lo, double hi) noexcept {
    double worst = 0.0;
    for (int i = 0; i <= n; ++i) {
        double const x = lo + (hi - lo) * static_cast<double>(i) / static_cast<double>(n);
        double const true_v = std::sin(std::numbers::pi_v<double> * 2.0 * x);
        double const err = std::fabs(static_cast<double>(approx(static_cast<double>(x))) - true_v);
        worst = std::max(worst, err);
    }
    return worst;
}

int main() {
    int failures = 0;
    constexpr int kSteps = 1 << 16; // 65k 采样点足以暴露多项式误差

    auto check = [&](const char* name, double err, double tol) {
        bool const ok = err <= tol;
        std::printf("%-16s max_abs_err = %.3e  (tol %.1e)  %s\n", name, err, tol, ok ? "OK" : "FAIL");
        if (!ok)
            ++failures;
    };

    // ── 精度测试：跨 [−0.5, 1.5]，验证任意相位下的折纸 + 多项式 ──
    constexpr double kLo = -0.5, kHi = 1.5; // 覆盖折纸边界（0.5、1.0）与出界相位
    check("Sin7<double>", MaxAbsErr([](double v) { return mm::Sin7(v); }, kSteps, kLo, kHi), kSin7MaxErr);
    check("Sin9<double>", MaxAbsErr([](double v) { return mm::Sin9(v); }, kSteps, kLo, kHi), kSin9MaxErr);
    check("Sin13<double>", MaxAbsErr([](double v) { return mm::Sin13(v); }, kSteps, kLo, kHi), kSin13MaxErr);

    // ── float 求值：物理极限（float 尾数约 24 位，绝对误差 ~1e-7 ~ 1e-6）──
    check("Sin7<float>", MaxAbsErr([](float v) { return mm::Sin7(v); }, kSteps, kLo, kHi), 1.5e-6);
    check("Sin9<float>", MaxAbsErr([](float v) { return mm::Sin9(v); }, kSteps, kLo, kHi), 1.5e-6);
    check("Sin13<float>", MaxAbsErr([](float v) { return mm::Sin13(v); }, kSteps, kLo, kHi), 1.5e-6);

    // ── 关键点校验（含折纸边界与出界相位）──
    struct Pt { double x; double want; };
    const Pt pts[] = {{0.0, 0.0}, {0.25, 1.0}, {0.5, 0.0}, {0.75, -1.0},
                      {1.0, 0.0}, {1.25, 1.0}, {-0.25, -1.0}, {2.5, 0.0}};
    for (auto [x, want] : pts) {
        double const got = mm::Sin13<double>(x);
        bool const ok = std::fabs(got - want) < 1e-12;
        std::printf("keypoint: sin(2pi*%5.2f) = %15.12f (want %5.1f)  %s\n", x, got, want, ok ? "OK" : "FAIL");
        if (!ok)
            ++failures;
    }

    // ── 周期性与连续规约：sin 应在 [0,1) 上连续，无跳变 ──
    double prev = mm::Sin13<double>(0.0);
    double worst_step = 0.0;
    constexpr int kN = 1 << 16;
    for (int i = 1; i <= kN; ++i) {
        double const x = static_cast<double>(i) / static_cast<double>(kN);
        double const y = mm::Sin13<double>(x);
        worst_step = std::max(worst_step, std::fabs(y - prev));
        prev = y;
    }
    bool const cont_ok = worst_step < 2e-4;
    std::printf("continuity: max_sample_step = %.3e  %s\n", worst_step, cont_ok ? "OK" : "FAIL");
    if (!cont_ok)
        ++failures;

    if (failures == 0)
        std::printf("ALL OK\n");
    else
        std::printf("%d FAILURES\n", failures);
    return failures == 0 ? 0 : 1;
}
