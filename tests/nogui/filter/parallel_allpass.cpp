// 并行全通（双路全通）低通设计的数值探针 / 测试。
//
// 对 butter / cheby1 / cheby2 / ellip 四族、多个奇数阶与截止频率：
//   1. 两条全通链的 |A_i(e^jw)| = 1（全通）；
//   2. |H|^2 + |Hc|^2 = 1（功率互补），H = (A0+A1)/2 低通，Hc = (A0-A1)/2 高通；
//   3. |H| 与各族**解析幅频定义式**逐点一致；
//   4. ellip 额外校验：通带纹波触底、阻带峰值触顶、阻带边沿由模数 k 给出、
//      阻带内传输零点个数 = (order-1)/2；
//   5. 退化规格（stop_db <= pass_db、规格过陡）返回 false 而不是静默给出错误结果。
//
// 失败返回非零退出码。设计公式与分解规则的推导/交叉验证见 tests/hybrid/tf2ca。
#include <qwqdsp/filter/iir_design.hpp>
#include <qwqdsp/filter/parallel_allpass.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <format>
#include <iostream>
#include <numbers>
#include <string_view>
#include <vector>

namespace {
int g_failures = 0;

// 系数按 float 存，高 Q 极点(|p| 接近 1)处的全通偏差约 eps/(1-|p|)，故留 1e-4；
// 真正的算法错误（极点分组错）会给出 O(1) 的偏差，这个门限足以区分。
constexpr double kAllpassTol = 1e-4;
constexpr double kMagnitudeTol = 1e-3;
constexpr double kSpecTol = 5e-3;

void Check(bool ok, std::string_view what, double got, double limit) {
    if (!ok) {
        ++g_failures;
        std::cout << std::format("  FAIL {}: got {:.3e}, limit {:.3e}\n", what, got, limit);
    }
}

/// 切比雪夫多项式 T_n(x)，|x| > 1 时用 cosh 分支
double ChebyshevT(size_t n, double x) {
    if (std::abs(x) <= 1.0) {
        return std::cos(static_cast<double>(n) * std::acos(x));
    }
    return std::cosh(static_cast<double>(n) * std::acosh(std::abs(x)));
}

struct AllpassResponse {
    std::complex<double> up;
    std::complex<double> down;
};

AllpassResponse Evaluate(qwqdsp_filter::ParallelAllpass& filter, double w) {
    auto const z = std::complex<float>(static_cast<float>(std::cos(w)), static_cast<float>(std::sin(w)));
    auto const [up, down] = filter.GetResponce(z);
    return {up, down};
}

double Magnitude(qwqdsp_filter::ParallelAllpass& filter, double w) {
    auto const [up, down] = Evaluate(filter, w);
    return std::abs(0.5 * (up + down));
}

/// 扫频检查全通性与功率互补，返回 {全通偏差, 互补偏差}
std::pair<double, double> CheckAllpassShape(qwqdsp_filter::ParallelAllpass& filter) {
    constexpr size_t kNumPoints = 2000;
    double flat = 0.0;
    double complementary = 0.0;
    for (size_t i = 0; i <= kNumPoints; ++i) {
        double const w = std::numbers::pi_v<double> * static_cast<double>(i) / static_cast<double>(kNumPoints);
        auto const [up, down] = Evaluate(filter, w);
        flat = std::max({flat, std::abs(std::abs(up) - 1.0), std::abs(std::abs(down) - 1.0)});
        auto const lp = 0.5 * (up + down);
        auto const hp = 0.5 * (up - down);
        complementary = std::max(complementary, std::abs(std::norm(lp) + std::norm(hp) - 1.0));
    }
    return {flat, complementary};
}

/// 扫频比较 |H| 与解析参考
double CompareMagnitude(qwqdsp_filter::ParallelAllpass& filter, auto&& reference) {
    constexpr size_t kNumPoints = 2000;
    double worst = 0.0;
    for (size_t i = 0; i <= kNumPoints; ++i) {
        double const w = std::numbers::pi_v<double> * static_cast<double>(i) / static_cast<double>(kNumPoints);
        worst = std::max(worst, std::abs(Magnitude(filter, w) - reference(w)));
    }
    return worst;
}

struct Case {
    std::string_view name;
    size_t order;
    double cutoff;
    double pass_db;
    double stop_db;
};

void ReportCase(Case const& c, double err) {
    Check(err < kMagnitudeTol, std::format("{} magnitude", c.name), err, kMagnitudeTol);
    std::cout << std::format("  {:<28s} max|dH| = {:.2e}\n", c.name, err);
}

void CheckShape(Case const& c, qwqdsp_filter::ParallelAllpass& filter) {
    auto const [flat, comp] = CheckAllpassShape(filter);
    Check(flat < kAllpassTol, std::format("{} allpass", c.name), flat, kAllpassTol);
    Check(comp < kAllpassTol, std::format("{} complementary", c.name), comp, kAllpassTol);
}

void RunButterworth(Case const& c) {
    qwqdsp_filter::ParallelAllpass filter;
    filter.BuildButterworth(c.order, static_cast<float>(c.cutoff));
    CheckShape(c, filter);
    double const n = static_cast<double>(c.order);
    double const tp = std::tan(c.cutoff / 2.0);
    ReportCase(c, CompareMagnitude(filter, [=](double w) {
        return 1.0 / std::sqrt(1.0 + std::pow(std::tan(w / 2.0) / tp, 2.0 * n));
    }));
}

void RunChebyshev1(Case const& c) {
    qwqdsp_filter::ParallelAllpass filter;
    filter.BuildChebyshev1(c.order, static_cast<float>(c.cutoff), static_cast<float>(c.pass_db));
    CheckShape(c, filter);
    double const eps = std::sqrt(std::pow(10.0, c.pass_db / 10.0) - 1.0);
    double const tp = std::tan(c.cutoff / 2.0);
    ReportCase(c, CompareMagnitude(filter, [=](double w) {
        double const t = ChebyshevT(c.order, std::tan(w / 2.0) / tp);
        return 1.0 / std::sqrt(1.0 + eps * eps * t * t);
    }));
}

void RunChebyshev2(Case const& c) {
    qwqdsp_filter::ParallelAllpass filter;
    filter.BuildChebyshev2(c.order, static_cast<float>(c.cutoff), static_cast<float>(c.stop_db));
    CheckShape(c, filter);
    double const eps = 1.0 / std::sqrt(std::pow(10.0, c.stop_db / 10.0) - 1.0);
    double const ts = std::tan(c.cutoff / 2.0);
    ReportCase(c, CompareMagnitude(filter, [=](double w) {
        double const t = ChebyshevT(c.order, ts / std::tan(w / 2.0));
        return 1.0 / std::sqrt(1.0 + 1.0 / (eps * eps * t * t));
    }));
}

/// 由阶数与两个纹波求出椭圆模数 k（构造时的定义式：q = exp(-pi K(k')/(n K(k))), k = cay(q)）
double EllipticModulus(size_t order, double pass_db, double stop_db) {
    using Helper = qwqdsp_filter::IIRDesign::EllipticHelperCephes;
    double const eps_p = std::sqrt(std::pow(10.0, pass_db / 10.0) - 1.0);
    double const eps_s = std::sqrt(std::pow(10.0, stop_db / 10.0) - 1.0);
    double const m1_sq = (eps_p / eps_s) * (eps_p / eps_s);
    double const Kk1 = Helper::Ellpk(1.0 - m1_sq);
    double const Kpk1 = Helper::Ellpk(m1_sq);
    double const q = std::exp(-std::numbers::pi_v<double> * Kpk1 / (static_cast<double>(order) * Kk1));
    return Helper::Cay(q);
}

void RunElliptic(Case const& c) {
    qwqdsp_filter::ParallelAllpass filter;
    bool const built = filter.BuildElliptic(c.order, static_cast<float>(c.cutoff), static_cast<float>(c.pass_db),
                                            static_cast<float>(c.stop_db));
    Check(built, std::format("{} build", c.name), built ? 1.0 : 0.0, 1.0);
    if (!built) {
        return;
    }
    CheckShape(c, filter);

    // 阻带边沿：模拟原型阻带边沿在 1/k，数字域 w_s = 2 atan(tan(w/2)/k)
    double const k = EllipticModulus(c.order, c.pass_db, c.stop_db);
    double const omega_stop = 2.0 * std::atan(std::tan(c.cutoff / 2.0) / k);
    double const pass_floor = std::pow(10.0, -c.pass_db / 20.0);
    double const stop_ceil = std::pow(10.0, -c.stop_db / 20.0);

    // 通带：幅度落在 [pass_floor, 1]，且真的触到纹波下限（等波纹定义）
    constexpr size_t kNumPass = 800;
    double pass_min = 1.0;
    double pass_max = 0.0;
    for (size_t i = 0; i <= kNumPass; ++i) {
        double const w = c.cutoff * static_cast<double>(i) / static_cast<double>(kNumPass);
        double const mag = Magnitude(filter, w);
        pass_min = std::min(pass_min, mag);
        pass_max = std::max(pass_max, mag);
    }
    Check(pass_max <= 1.0 + kAllpassTol, std::format("{} passband max", c.name), pass_max, 1.0);
    Check(pass_min >= pass_floor - kSpecTol, std::format("{} passband min", c.name), pass_min, pass_floor);
    Check(std::abs(pass_min - pass_floor) < kSpecTol, std::format("{} passband ripple floor", c.name),
          std::abs(pass_min - pass_floor), kSpecTol);

    // 阻带：峰值既不超 stop_ceil 也要触到 stop_ceil
    constexpr size_t kNumStop = 20000;
    std::vector<double> stop(kNumStop + 1);
    double stop_peak = 0.0;
    for (size_t i = 0; i <= kNumStop; ++i) {
        double const w = omega_stop
                       + (std::numbers::pi_v<double> - omega_stop) * static_cast<double>(i)
                             / static_cast<double>(kNumStop);
        stop[i] = Magnitude(filter, w);
        stop_peak = std::max(stop_peak, stop[i]);
    }
    Check(stop_peak <= stop_ceil + kAllpassTol, std::format("{} stopband peak", c.name), stop_peak, stop_ceil);
    Check(std::abs(stop_peak - stop_ceil) < kSpecTol, std::format("{} stopband ripple ceiling", c.name),
          std::abs(stop_peak - stop_ceil), kSpecTol);
    Check(std::abs(Magnitude(filter, omega_stop) - stop_ceil) < kSpecTol,
          std::format("{} stopband edge level", c.name), std::abs(Magnitude(filter, omega_stop) - stop_ceil),
          kSpecTol);

    // 传输零点：窗口内取到极小、且明显低于阻带峰值。
    // 在整段 (0, pi) 上扫：通带幅度远高于门限不会被计入，而第一个零点可能紧贴
    // 阻带边沿，若从 omega_stop 起扫会被窗口边界漏掉。近 pi 的那个零点(分子在
    // z=-1)落在尾部保护窗口内，不计。
    constexpr size_t kWindow = 8;
    size_t nulls = 0;
    double const floor_threshold = 0.5 * stop_peak;
    std::vector<double> all(kNumStop + 1);
    for (size_t i = 0; i <= kNumStop; ++i) {
        all[i] = Magnitude(filter, std::numbers::pi_v<double> * static_cast<double>(i) / static_cast<double>(kNumStop));
    }
    for (size_t i = kWindow; i + kWindow <= kNumStop; ++i) {
        if (all[i] >= floor_threshold) {
            continue;
        }
        bool is_min = true;
        for (size_t j = 1; j <= kWindow; ++j) {
            if (all[i] >= all[i - j] || all[i] >= all[i + j]) {
                is_min = false;
                break;
            }
        }
        if (is_min) {
            ++nulls;
        }
    }
    size_t const expected_nulls = (c.order - 1) / 2;
    Check(nulls == expected_nulls, std::format("{} transmission zeros", c.name), static_cast<double>(nulls),
          static_cast<double>(expected_nulls));
    std::cout << std::format("  {:<28s} ripple=[{:.4f},{:.4f}]dB stop={:.2f}dB nulls={}\n", c.name,
                             20.0 * std::log10(pass_min), 20.0 * std::log10(pass_max),
                             20.0 * std::log10(stop_peak), nulls);
}

void RunFailurePaths() {
    qwqdsp_filter::ParallelAllpass filter;
    // stop_db <= pass_db：定义上非法，必须拒绝
    Check(!filter.BuildElliptic(5, 0.3f, 40.0f, 20.0f), "elliptic rejects stop<=pass", 0.0, 0.0);
    // 规格过陡（阻带边沿贴住通带边沿，模数在 double 下就是 1）：必须拒绝而不是死循环/静默给错
    Check(!filter.BuildElliptic(25, 0.3f, 6.0f, 6.6f), "elliptic rejects too-steep", 0.0, 0.0);
    Check(!filter.BuildElliptic(25, 0.3f, 3.0f, 3.3f), "elliptic rejects too-steep (narrow)", 0.0, 0.0);
}
} // namespace

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::vector<Case> const cases{
        {"butter  n=3 w=0.3", 3, 0.3, 0.0, 0.0},
        {"butter  n=7 w=0.5", 7, 0.5, 0.0, 0.0},
        {"butter  n=9 w=0.2", 9, 0.2, 0.0, 0.0},
        {"cheby1  n=5 rp=1", 5, 0.3, 1.0, 0.0},
        {"cheby1  n=7 rp=0.5", 7, 0.4, 0.5, 0.0},
        {"cheby1  n=9 rp=3", 9, 0.25, 3.0, 0.0},
        {"cheby2  n=5 rs=40", 5, 0.4, 0.0, 40.0},
        {"cheby2  n=7 rs=60", 7, 0.3, 0.0, 60.0},
        {"cheby2  n=9 rs=30", 9, 0.5, 0.0, 30.0},
        {"ellip   n=5 rp=0.5 rs=40", 5, 0.3, 0.5, 40.0},
        {"ellip   n=7 rp=1 rs=60", 7, 0.25, 1.0, 60.0},
        {"ellip   n=9 rp=0.2 rs=40", 9, 0.35, 0.2, 40.0},
    };

    std::cout << "== parallel allpass: magnitude vs analytic reference ==\n";
    for (auto const& c : cases) {
        std::string_view const family = c.name.substr(0, 6);
        if (family == "butter") {
            RunButterworth(c);
        }
        else if (family == "cheby1") {
            RunChebyshev1(c);
        }
        else if (family == "cheby2") {
            RunChebyshev2(c);
        }
        else {
            RunElliptic(c);
        }
    }
    RunFailurePaths();

    if (g_failures == 0) {
        std::cout << "OK\n";
        return 0;
    }
    std::cout << std::format("FAILED: {} check(s)\n", g_failures);
    return 1;
}
