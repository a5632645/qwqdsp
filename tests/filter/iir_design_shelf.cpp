/// @brief 验证 highshelf 原型设计器的三电平契约：直流 0dB、中点 boost/2 dB、幅度极值 boost dB
//
// 每节都是二阶节(一对共轭极点、一对共轭零点)，这里按定义式逐点算 |H(jw)| 后扫频：
//   1. 直流 |H| = 1；
//   2. (1)rad/sec 处 |H| = 10^(boost/40)，即 boost/2 dB；
//   3. 幅度极值(boost > 0 取最大、boost < 0 取最小) = 10^(boost/20)，即 boost dB。
//      未修正的切比雪夫 II / 椭圆高频端本身达不到这个极值，极值出现在棚带的传输零点上，
//      所以先粗扫找局部极值再逐个三分细化，而不是只看 w -> inf；
//   4. 极点全在左半平面(稳定)。
//
// 阈值：巴特沃斯/切比雪夫在 double 下就是精确的(1e-9 dB)；椭圆继承 IIRDesign::Elliptic 形状
// 函数的数值误差(越高阶、阻带越浅越大)，本测试覆盖的规格(<= 16 阶, 阻带 >= 20dB)取 1e-4 dB。
//
// 失败返回非零退出码。
#include <qwqdsp/filter/iir_design_shelf.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <format>
#include <iostream>
#include <span>
#include <string>

namespace {
int g_failures = 0;

using qwqdsp_filter::IIRDesignShelf;
using ZPK = IIRDesignShelf::ZPK;

constexpr size_t kMaxSections = 16;
/// 精确到 double 的组合用的电平偏差上限(dB)
constexpr double kLevelToleranceDb = 1.0e-9;
/// 椭圆组合的电平偏差上限，取决于椭圆原型自身的精度
constexpr double kEllipticToleranceDb = 1.0e-4;
/// 极值只可能落在 (1)rad/sec 以上的棚带里(尺度 lambda <= 1)，这里扫 [1, 1e4]
constexpr double kScanLo = 1.0;
constexpr double kScanHi = 1e4;
constexpr int kCoarsePoints = 1200;
constexpr int kRefineSteps = 60;

/// 逐节相乘求 |H(jw)|^2（无零点节表示零点在无穷远）
double Magnitude2(std::span<ZPK const> sections, double w) {
    std::complex<double> const s{0.0, w};
    std::complex<double> h{1.0, 0.0};
    for (auto const& sec : sections) {
        std::complex<double> numerator{1.0, 0.0};
        if (sec.z) {
            numerator = (s - *sec.z) * (s - std::conj(*sec.z));
        }
        h *= sec.k * numerator / ((s - sec.p) * (s - std::conj(sec.p)));
    }
    return std::norm(h);
}

double LevelDb(std::span<ZPK const> sections, double w) {
    return 10.0 * std::log10(Magnitude2(sections, w));
}

/// 幅度极值(dB)：候选 = 扫频区间两端点 + 所有内部局部极值(三分细化) + w -> inf 的极限
/// @note 单边严格大于(小于)才算局部极值, 免得平坦段里每个点都算一个
double ExtremumDb(std::span<ZPK const> sections, bool want_max) {
    double const log_lo = std::log10(kScanLo);
    double const step = (std::log10(kScanHi) - log_lo) / static_cast<double>(kCoarsePoints);
    std::array<double, kCoarsePoints + 1> sampled{};
    for (int i = 0; i <= kCoarsePoints; ++i) {
        sampled[i] = Magnitude2(sections, std::pow(10.0, log_lo + step * i));
    }
    auto better = [want_max](double a, double b) { return want_max ? (a > b) : (a < b); };
    // 两端点 + w -> inf 的极限(每节都趋于 k, 所以极限是 prod(k)^2)
    double best = better(sampled[0], sampled[kCoarsePoints]) ? sampled[0] : sampled[kCoarsePoints];
    double limit = 1.0;
    for (auto const& sec : sections) {
        limit *= sec.k;
    }
    limit *= limit;
    if (better(limit, best)) {
        best = limit;
    }
    for (int i = 1; i < kCoarsePoints; ++i) {
        bool const extreme = want_max ? ((sampled[i] > sampled[i - 1] && sampled[i] >= sampled[i + 1])
                                         || (sampled[i] >= sampled[i - 1] && sampled[i] > sampled[i + 1]))
                                      : ((sampled[i] < sampled[i - 1] && sampled[i] <= sampled[i + 1])
                                         || (sampled[i] <= sampled[i - 1] && sampled[i] < sampled[i + 1]));
        if (!extreme) {
            continue;
        }
        double a = log_lo + step * (i - 1);
        double b = log_lo + step * (i + 1);
        for (int k = 0; k < kRefineSteps; ++k) {
            double const m1 = a + (b - a) / 3.0;
            double const m2 = b - (b - a) / 3.0;
            double const v1 = Magnitude2(sections, std::pow(10.0, m1));
            double const v2 = Magnitude2(sections, std::pow(10.0, m2));
            if (want_max ? (v1 < v2) : (v1 > v2)) {
                a = m1;
            }
            else {
                b = m2;
            }
        }
        double const v = Magnitude2(sections, std::pow(10.0, 0.5 * (a + b)));
        if (better(v, best)) {
            best = v;
        }
    }
    return 10.0 * std::log10(best);
}

void Check(bool ok, std::string const& what, double got, double limit) {
    if (!ok) {
        ++g_failures;
        std::cout << std::format("  FAIL {}: got {:.6e}, limit {:.6e}\n", what, got, limit);
    }
}

/// 校验一组设计好的节: 直流 0dB、中点 boost/2 dB、极值 boost dB、极点稳定
void CheckLevels(std::string const& tag, std::span<ZPK const> sections, double boost_db, double tolerance_db,
                 bool designed) {
    Check(designed, tag + ": design succeeds", designed ? 1.0 : 0.0, 1.0);
    if (!designed) {
        return;
    }
    double const dc_db = LevelDb(sections, 0.0);
    Check(std::abs(dc_db) < tolerance_db, tag + ": DC = 0dB", dc_db, tolerance_db);
    double const mid_db = LevelDb(sections, 1.0);
    Check(std::abs(mid_db - boost_db / 2.0) < tolerance_db, tag + ": mid = boost/2 dB", mid_db, boost_db / 2.0);
    double const peak_db = ExtremumDb(sections, boost_db > 0.0);
    Check(std::abs(peak_db - boost_db) < tolerance_db, tag + ": extremum = boost dB", peak_db, boost_db);
    for (auto const& sec : sections) {
        Check(sec.p.real() < 0.0, tag + ": pole in left half plane", sec.p.real(), 0.0);
    }
}
} // namespace

int main() {
    std::array<ZPK, kMaxSections> sections{};
    constexpr std::array<size_t, 4> kPairs{1, 2, 4, 8};
    constexpr std::array<double, 6> kBoosts{3.0, 12.0, 24.0, -3.0, -12.0, -24.0};
    constexpr std::array<double, 2> kRipples{0.1, 1.0};
    constexpr std::array<double, 2> kStopbands{20.0, 80.0};

    for (size_t const pairs : kPairs) {
        std::span<ZPK> view{sections.data(), pairs};
        for (double const boost : kBoosts) {
            for (double const ripple : kRipples) {
                std::string const base = std::format("pairs={} boost={} ripple={}", pairs, boost, ripple);
                IIRDesignShelf::Butterworth(view, pairs, boost);
                CheckLevels(base + " butterworth", view, boost, kLevelToleranceDb, true);
                for (bool const modified : {false, true}) {
                    IIRDesignShelf::Chebyshev1(view, pairs, boost, ripple, modified);
                    CheckLevels(base + std::format(" chebyshev1 modify={}", modified), view, boost, kLevelToleranceDb,
                                true);
                    bool const cheb2_ok = IIRDesignShelf::Chebyshev2(view, pairs, boost, ripple, modified);
                    CheckLevels(base + std::format(" chebyshev2 modify={}", modified), view, boost, kLevelToleranceDb,
                                cheb2_ok);
                    for (double const stopband : kStopbands) {
                        bool const ellip_ok = IIRDesignShelf::Elliptic(view, pairs, boost, ripple, stopband, modified);
                        CheckLevels(base + std::format(" elliptic modify={} stop={}", modified, stopband), view, boost,
                                    kEllipticToleranceDb, ellip_ok);
                    }
                }
            }
        }
    }

    // 尺度回退: 当前尺度下三个电平解不出时, 椭圆会把频率尺度推进原型阻带再解
    // 提升方向与衰减方向各覆盖一例
    struct FallbackCase {
        size_t pairs;
        double boost_db;
        double ripple_db;
        double stopband_db;
    };
    for (auto const& c : {
             FallbackCase{2, 24.0,  3.0, 20.0},
             FallbackCase{2, -24.0, 1.0, 10.0},
             FallbackCase{4, 40.0,  1.0, 20.0}
    }) {
        for (bool const modified : {false, true}) {
            std::span<ZPK> view{sections.data(), c.pairs};
            bool const ok = IIRDesignShelf::Elliptic(view, c.pairs, c.boost_db, c.ripple_db, c.stopband_db, modified);
            CheckLevels(std::format("fallback pairs={} boost={} ripple={} stop={} modify={}", c.pairs, c.boost_db,
                                    c.ripple_db, c.stopband_db, modified),
                        view, c.boost_db, kEllipticToleranceDb, ok);
        }
    }

    // 退化输入: boost = 0 输出直通; 阻带 <= 通带纹波 直接失败
    {
        std::span<ZPK> view{sections.data(), 3};
        IIRDesignShelf::Butterworth(view, 3, 0.0);
        Check(std::abs(LevelDb(view, 1.0)) < kLevelToleranceDb, "boost = 0 -> bypass", LevelDb(view, 1.0), 0.0);
        IIRDesignShelf::Chebyshev1(view, 3, 0.0, 1.0, true);
        Check(std::abs(LevelDb(view, 1.0)) < kLevelToleranceDb, "boost = 0 -> bypass (chebyshev1)", LevelDb(view, 1.0),
              0.0);
        bool const bad1 = IIRDesignShelf::Elliptic(view, 3, 12.0, 6.0, 5.0, false);
        bool const bad2 = IIRDesignShelf::Elliptic(view, 3, 12.0, 6.0, 6.0, true);
        Check(!bad1 && !bad2, "elliptic stopband <= ripple -> fail", bad1 || bad2 ? 1.0 : 0.0, 0.0);
    }

    if (g_failures == 0) {
        std::cout << "iir_design_shelf: all checks passed\n";
    }
    return g_failures == 0 ? 0 : 1;
}
