// 椭圆原型「偶数阶修正」的数值测试。
//
// IIRDesign 的阶数恒为偶数(每节一对共轭极点)，而偶数阶椭圆原型的直流电平落在
// 纹波谷(-db_passband dB)、w->inf 的电平落在纹波峰(-db_stopband dB)，于是它的
// 部分分式展开必然带常数项(直接项)。修正把最低反射零点搬到直流、最高传输零点
// 搬到无穷远。这里按**定义式**扫频校验：
//
//   1. 未修正: 直流电平 = -db_passband dB；num_filter 节都有有限零点；极点都在左半平面;
//   2. 修正后: 最后一节零点在无穷远(N-2 个有限零点); 直流增益 1;
//      极点仍在左半平面; 通带 [0,1] 内 |H| 落在 [10^(-rp/20), 1];
//      阻带 (过第一次降到规格的频率之后) |H| <= 10^(-rs/20); w->inf 时 |H| -> 0;
//   3. 修正不搬动通带边沿(w = 1 仍是纹波峰);
//   4. 退化规格(rs <= rp)仍返回 false, 两个实现(Elliptic / EllipticLanden)一致。
//
// 失败返回非零退出码。推导与外部核对见 labs/holters_parker。
#include <qwqdsp/filter/iir_design.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <format>
#include <iostream>
#include <span>
#include <string_view>

namespace {
int g_failures = 0;

void Check(bool ok, std::string_view what, double got, double limit) {
    if (!ok) {
        ++g_failures;
        std::cout << std::format("  FAIL {}: got {:.3e}, limit {:.3e}\n", what, got, limit);
    }
}

using qwqdsp_filter::IIRDesign;
using ZPK = IIRDesign::ZPK;

constexpr size_t kMaxSections = 16;

/// 逐节相乘求 |H(jw)|（每节: k*(s^2+|z|^2)/((s-p)(s-p*))，无零点时分子为 1）
std::complex<double> Evaluate(std::span<ZPK const> sections, double w) {
    std::complex<double> const s{0.0, w};
    std::complex<double> h{1.0, 0.0};
    for (auto const& sec : sections) {
        std::complex<double> numerator{1.0, 0.0};
        if (sec.z) {
            numerator = (s - *sec.z) * (s - std::conj(*sec.z));
        }
        h *= sec.k * numerator / ((s - sec.p) * (s - std::conj(sec.p)));
    }
    return h;
}

struct Scan {
    double passband_min = 1e30;
    double passband_max = 0.0;
    double stopband_max_after_first_crossing = 0.0;
    double at_infinity = 0.0;
};

/// 扫频: 通带取 [1e-4, 1], 阻带取「第一次降到 spec 以下」之后的全部频点
Scan ScanResponse(std::span<ZPK const> sections, double spec_db) {
    Scan out;
    double const spec = std::pow(10.0, -spec_db / 20.0);
    constexpr size_t kNumPoints = 20000;
    bool crossed = false;
    double const lo = 1e-4;
    double const hi = 1e6;
    for (size_t i = 0; i <= kNumPoints; ++i) {
        double const w = lo * std::pow(hi / lo, static_cast<double>(i) / static_cast<double>(kNumPoints));
        double const mag = std::abs(Evaluate(sections, w));
        if (w <= 1.0) {
            out.passband_min = std::min(out.passband_min, mag);
            out.passband_max = std::max(out.passband_max, mag);
        }
        else if (mag < spec) {
            crossed = true;
        }
        if (crossed) {
            out.stopband_max_after_first_crossing = std::max(out.stopband_max_after_first_crossing, mag);
        }
    }
    out.at_infinity = std::abs(Evaluate(sections, 1e9));
    return out;
}

void CheckPolesStable(std::span<ZPK const> sections, std::string_view tag) {
    double worst = -1e30;
    for (auto const& sec : sections) {
        worst = std::max(worst, sec.p.real());
    }
    Check(worst < 0.0, std::format("{}: all poles in LHP (worst Re)", tag), worst, 0.0);
}

void CheckFiniteZeros(std::span<ZPK const> sections, std::string_view tag) {
    // 修正后只有最后一节没有有限零点(零点被搬到无穷远)
    size_t without = 0;
    bool last_without = false;
    for (size_t i = 0; i < sections.size(); ++i) {
        if (!sections[i].z) {
            ++without;
            last_without = (i + 1 == sections.size());
        }
    }
    Check(without == 1 && last_without, std::format("{}: only the last zero goes to infinity", tag),
          static_cast<double>(without), 1.0);
}

void RunCase(size_t num_filter, double rp, double rs) {
    std::string const tag = std::format("n={} rp={} rs={}", 2 * num_filter, rp, rs);

    std::array<ZPK, kMaxSections> plain{};
    std::array<ZPK, kMaxSections> modified{};
    std::span<ZPK> plain_view{plain.data(), num_filter};
    std::span<ZPK> modified_view{modified.data(), num_filter};

    bool const plain_ok = IIRDesign::Elliptic(plain_view, num_filter, rp, rs, false);
    bool const modified_ok = IIRDesign::Elliptic(modified_view, num_filter, rp, rs, true);
    Check(plain_ok && modified_ok, tag + ": design succeeds", plain_ok && modified_ok ? 1.0 : 0.0, 1.0);
    if (!plain_ok || !modified_ok) {
        return;
    }

    // ----- 未修正的既有约定 -----
    CheckPolesStable(plain_view, tag + " plain");
    for (auto const& sec : plain_view) {
        Check(sec.z.has_value(), tag + ": plain keeps every finite zero", 0.0, 0.0);
    }
    double const plain_dc = std::abs(Evaluate(plain_view, 0.0));
    // 直流电平 = -rp dB（旧行为，修正后不再成立）
    Check(std::abs(plain_dc - std::pow(10.0, -rp / 20.0)) < 1e-9, tag + ": plain DC level = -rp dB",
          plain_dc, std::pow(10.0, -rp / 20.0));

    // ----- 修正后的规格 -----
    CheckPolesStable(modified_view, tag + " modified");
    CheckFiniteZeros(modified_view, tag + " modified");
    double const dc = std::abs(Evaluate(modified_view, 0.0));
    Check(std::abs(dc - 1.0) < 1e-9, tag + ": modified DC gain = 1", dc, 1.0);

    auto const scan = ScanResponse(modified_view, rs);
    double const passband_floor = std::pow(10.0, -rp / 20.0);
    Check(scan.passband_min >= passband_floor - 1e-6, tag + ": passband ripple >= -rp dB",
          scan.passband_min, passband_floor);
    Check(scan.passband_max <= 1.0 + 1e-4, tag + ": passband ripple <= 0 dB", scan.passband_max, 1.0);
    double const stopband_ceiling = std::pow(10.0, -rs / 20.0);
    Check(scan.stopband_max_after_first_crossing <= stopband_ceiling * (1.0 + 1e-3),
          tag + ": stopband <= -rs dB after first crossing", scan.stopband_max_after_first_crossing,
          stopband_ceiling);
    Check(scan.at_infinity < 1e-12, tag + ": strictly proper (|H(inf)| = 0)", scan.at_infinity, 1e-12);
    // 修正把直流抬到 0dB, w = 1 处仍是 -rp dB（与 Chebyshev1 的偶数阶修正同一约定）
    double const at_edge = std::abs(Evaluate(modified_view, 1.0));
    Check(std::abs(at_edge - passband_floor) < 1e-9, tag + ": w = 1 still at -rp dB", at_edge, passband_floor);

    // ----- 两个实现（nome 反解 / Landen）在修正上一致 -----
    std::array<ZPK, kMaxSections> landen{};
    std::span<ZPK> landen_view{landen.data(), num_filter};
    if (IIRDesign::EllipticLanden(landen_view, num_filter, rp, rs, true)) {
        double worst = 0.0;
        constexpr size_t kNumPoints = 200;
        for (size_t i = 0; i <= kNumPoints; ++i) {
            double const w = 1e-3 * std::pow(1e3 / 1e-3, static_cast<double>(i) / static_cast<double>(kNumPoints));
            worst = std::max(worst, std::abs(std::abs(Evaluate(modified_view, w)) - std::abs(Evaluate(landen_view, w))));
        }
        Check(worst < 1e-4, tag + ": Elliptic and EllipticLanden agree", worst, 1e-4);
    }
}

void RunFailurePaths() {
    std::array<ZPK, kMaxSections> sections{};
    std::span<ZPK> view{sections.data(), 2};
    // rs <= rp 定义上非法: 两个实现、两种修正都必须拒绝
    for (bool const modify : {false, true}) {
        Check(!IIRDesign::Elliptic(view, 2, 40.0, 20.0, modify), "elliptic rejects rs <= rp", 0.0, 0.0);
        Check(!IIRDesign::EllipticLanden(view, 2, 40.0, 20.0, modify),
              "elliptic landen rejects rs <= rp", 0.0, 0.0);
    }
}
} // namespace

int main() {
    for (auto const& [num_filter, rp, rs] : {std::tuple{1u, 1.0, 40.0}, {2u, 3.0, 30.0}, {3u, 0.1, 80.0},
                                             {4u, 0.5, 60.0}, {8u, 0.01, 100.0}}) {
        RunCase(num_filter, rp, rs);
    }
    RunFailurePaths();

    if (g_failures == 0) {
        std::cout << "iir_design: all checks passed\n";
        return 0;
    }
    std::cout << std::format("iir_design: {} check(s) failed\n", g_failures);
    return 1;
}
