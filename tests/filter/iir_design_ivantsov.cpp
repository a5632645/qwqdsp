// 一般 Ivantsov 变换（IIRDesign::Ivantsov）的数值测试。
//
// 1. sigma = 0 时必须与 Bilinear 逐系数一致（论文 §4 在 sigma = 0 退化）；
// 2. 单个二阶 lowpass / highpass 节必须与 ivantsov.hpp 的 biquad 设计器一致
//    （用同一 sigma，并把模拟节缩放到 biquad 的预畸变 omega）；
// 3. 高阶 LP / HP / BP / BS 与 Chebyshev II shelf：级联稳定、每节有有限零点，
//    且数字响应与模拟原型接近。
//
// 失败返回非零退出码。
#include <qwqdsp/filter/iir_design.hpp>
#include <qwqdsp/filter/iir_design_shelf.hpp>
#include <qwqdsp/filter/ivantsov.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <format>
#include <iostream>
#include <limits>
#include <numbers>
#include <span>
#include <string_view>

namespace {
int g_failures = 0;

using qwqdsp_filter::IIRDesign;
using qwqdsp_filter::IIRDesignShelf;
using BiquadIvantsov = qwqdsp_filter::Ivantsov;
using ZPK = IIRDesign::ZPK;
using BiquadCoeff = qwqdsp_filter::BiquadCoeff;

constexpr double kPi = IIRDesign::pi;
constexpr double kFs = 48000.0;
double const kSigmaOpt = kPi * std::sqrt(2.0 / 3.0);
constexpr size_t kMaxSections = 32;

void Check(bool ok, std::string_view what) {
    if (!ok) {
        ++g_failures;
        std::cout << std::format("  FAIL {}\n", what);
    }
}

void CheckClose(double got, double want, double tolerance, std::string_view what) {
    double const scale = std::max({1.0, std::abs(got), std::abs(want)});
    if (std::abs(got - want) > tolerance * scale) {
        ++g_failures;
        std::cout << std::format("  FAIL {}: got {:.12g}, want {:.12g}, tol {:.1e}\n",
                                 what, got, want, tolerance);
    }
}

double PrewarpedOmega(double wc, double sigma) {
    if (wc > kPi) {
        return 2.0 * kPi / wc;
    }
    double const cot = kPi / std::tan(wc / 2.0);
    return std::sqrt(sigma * sigma + cot * cot);
}

std::complex<double> EvalAnalog(std::span<ZPK const> sections, double w) {
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

/// 数字 ZPK 节在 w(rad/sample) 处的响应; 要求每节都有有限零点
std::complex<double> EvalDigital(std::span<ZPK const> sections, double w) {
    std::complex<double> const zi = std::polar(1.0, -w);
    std::complex<double> h{1.0, 0.0};
    for (auto const& sec : sections) {
        if (!sec.z) {
            return {std::numeric_limits<double>::quiet_NaN(), 0.0};
        }
        std::complex<double> const num = (1.0 - *sec.z * zi) * (1.0 - std::conj(*sec.z) * zi);
        std::complex<double> const den = (1.0 - sec.p * zi) * (1.0 - std::conj(sec.p) * zi);
        h *= sec.k * num / den;
    }
    return h;
}

double MaxDbError(std::span<ZPK const> analog, std::span<ZPK const> digital,
                  double analog_floor_db = -300.0) {
    constexpr size_t kNumPoints = 20000;
    double worst = 0.0;
    double const lo = 20.0;
    double const hi = kFs / 2.0;
    double const floor = std::pow(10.0, analog_floor_db / 20.0);
    for (size_t i = 0; i <= kNumPoints; ++i) {
        double const f = lo * std::pow(hi / lo, static_cast<double>(i) / static_cast<double>(kNumPoints));
        double const w = 2.0 * kPi * f / kFs;
        double const ha = std::abs(EvalAnalog(analog, 2.0 * kPi * f));
        if (ha < floor) {
            continue;
        }
        double const hd = std::abs(EvalDigital(digital, w));
        if (!(ha > 0.0) || !(hd > 0.0)) {
            worst = std::numeric_limits<double>::infinity();
            break;
        }
        worst = std::max(worst, std::abs(20.0 * std::log10(hd) - 20.0 * std::log10(ha)));
    }
    return worst;
}

void CheckStableAndFinite(std::span<ZPK const> sections, std::string_view tag) {
    for (auto const& sec : sections) {
        Check(std::abs(sec.p) < 1.0, std::format("{}: pole inside unit circle", tag));
        Check(sec.z.has_value(), std::format("{}: every section has a finite zero", tag));
    }
}

/// 构造一个 biquad Ivantsov 的原型节, 并缩放到 biquad 的预畸变 omega
ZPK MakeBiquadSection(double wc, double Q, double sigma, bool highpass) {
    double const zeta = 1.0 / (2.0 * Q);
    double const omega = PrewarpedOmega(wc, sigma);
    double const scale = 2.0 * kPi * kFs / omega;
    ZPK sec{};
    sec.p = std::complex<double>{-zeta, std::sqrt(1.0 - zeta * zeta)} * scale;
    if (highpass) {
        sec.z = std::complex<double>{0.0, 0.0};
        sec.k = 1.0;
    }
    else {
        sec.z = std::nullopt;
        sec.k = scale * scale;
    }
    return sec;
}

void CheckDefaultOverload() {
    auto const base = MakeBiquadSection(0.25, 2.0, IIRDesign::kDefaultIvantsovSigma, false);
    std::array<ZPK, 1> implicit{base};
    std::array<ZPK, 1> explicit_{base};
    IIRDesign::Ivantsov(implicit, kFs);
    IIRDesign::Ivantsov(explicit_, kFs, IIRDesign::kDefaultIvantsovSigma);
    CheckClose(std::abs(implicit[0].p - explicit_[0].p), 0.0, 1e-15, "default sigma pole");
    CheckClose(std::abs(*implicit[0].z - *explicit_[0].z), 0.0, 1e-15, "default sigma zero");
    CheckClose(implicit[0].k, explicit_[0].k, 1e-15, "default sigma k");
}

void CheckPrewarp() {
    for (double const freq : {20.0, 100.0, 1000.0, 5000.0, 10000.0, 20000.0}) {
        double const bilinear = IIRDesign::PrewarpBilinear(freq, kFs);
        CheckClose(IIRDesign::PrewarpIvantsov(freq, kFs, 0.0), bilinear, 1e-12,
                   std::format("PrewarpIvantsov sigma=0 @ {} Hz", freq));
        for (double const sigma : {2.0, kSigmaOpt}) {
            double const wa = IIRDesign::PrewarpIvantsov(freq, kFs, sigma);
            double const omega = 2.0 * kPi * kFs / wa;
            double const theta = std::acos((omega * omega - kPi * kPi - sigma * sigma)
                                           / (omega * omega + kPi * kPi - sigma * sigma));
            double const mapped = kFs / (2.0 * kPi) * theta;
            CheckClose(mapped, freq, 1e-9,
                       std::format("PrewarpIvantsov @ {} Hz sigma={}", freq, sigma));
        }
    }
}

void CheckBilinearLimit() {
    std::array<ZPK, kMaxSections> proto{};
    std::array<ZPK, kMaxSections> analytic{};
    std::array<ZPK, kMaxSections> bilinear{};
    std::array<BiquadCoeff, kMaxSections> a_coeff{};
    std::array<BiquadCoeff, kMaxSections> b_coeff{};

    auto run = [&](std::string_view tag, size_t n) {
        std::copy_n(proto.begin(), n, analytic.begin());
        std::copy_n(proto.begin(), n, bilinear.begin());
        IIRDesign::Ivantsov(std::span{analytic}.first(n), kFs, 0.0);
        IIRDesign::Bilinear(std::span{bilinear}.first(n), kFs);
        IIRDesign::TfToBiquad(std::span{analytic}.first(n), std::span{a_coeff}.first(n));
        IIRDesign::TfToBiquad(std::span{bilinear}.first(n), std::span{b_coeff}.first(n));
        for (size_t i = 0; i < n; ++i) {
            auto const& a = a_coeff[i];
            auto const& b = b_coeff[i];
            CheckClose(a.b0, b.b0, 1e-5, std::format("{} sec{} b0 = bilinear", tag, i));
            CheckClose(a.b1, b.b1, 1e-5, std::format("{} sec{} b1 = bilinear", tag, i));
            CheckClose(a.b2, b.b2, 1e-5, std::format("{} sec{} b2 = bilinear", tag, i));
            CheckClose(a.a1, b.a1, 1e-5, std::format("{} sec{} a1 = bilinear", tag, i));
            CheckClose(a.a2, b.a2, 1e-5, std::format("{} sec{} a2 = bilinear", tag, i));
        }
    };

    // 4th-order Butterworth lowpass / highpass
    IIRDesign::Butterworth(std::span{proto}.first(2), 2);
    IIRDesign::ProtyleToLowpass(std::span{proto}.first(2), 2, 2.0 * kPi * 1000.0);
    run("LP", 2);

    std::fill(proto.begin(), proto.end(), ZPK{});
    IIRDesign::Butterworth(std::span{proto}.first(2), 2);
    IIRDesign::ProtyleToHighpass(std::span{proto}.first(2), 2, 2.0 * kPi * 1000.0);
    run("HP", 2);

    // 4th-order bandpass / bandstop: 每节翻倍
    std::fill(proto.begin(), proto.end(), ZPK{});
    IIRDesign::Butterworth(std::span{proto}.first(2), 2);
    IIRDesign::ProtyleToBandpass(std::span{proto}.first(4), 2, 2.0 * kPi * 1000.0, 2.0);
    run("BP", 4);

    std::fill(proto.begin(), proto.end(), ZPK{});
    IIRDesign::Butterworth(std::span{proto}.first(2), 2);
    IIRDesign::ProtyleToBandstop(std::span{proto}.first(4), 2, 2.0 * kPi * 1000.0, 2.0);
    run("BS", 4);
}

void CheckSingleSectionMatchesBiquad() {
    BiquadIvantsov maker;
    double const sigmas[] = {2.0, kSigmaOpt};
    double const wcs[] = {0.1, 0.5, 1.5, 2.5};
    double const qs[] = {0.70710678, 2.0, 10.0};
    for (double sigma : sigmas) {
        for (double wc : wcs) {
            for (double q : qs) {
                for (bool highpass : {false, true}) {
                    std::array<ZPK, 1> analog{MakeBiquadSection(wc, q, sigma, highpass)};
                    IIRDesign::Ivantsov(analog, kFs, sigma);
                    std::array<BiquadCoeff, 1> got{};
                    IIRDesign::TfToBiquad(analog, got);
                    BiquadCoeff const ref = highpass
                        ? maker.Highpass(static_cast<float>(wc), static_cast<float>(q), static_cast<float>(sigma))
                        : maker.Lowpass(static_cast<float>(wc), static_cast<float>(q), static_cast<float>(sigma));
                    std::string const tag = std::format("{} wc={} Q={} sigma={}",
                        highpass ? "HP" : "LP", wc, q, sigma);
                    CheckClose(got[0].b0, ref.b0, 2e-4, tag + " b0");
                    CheckClose(got[0].b1, ref.b1, 2e-4, tag + " b1");
                    CheckClose(got[0].b2, ref.b2, 2e-4, tag + " b2");
                    CheckClose(got[0].a1, ref.a1, 2e-4, tag + " a1");
                    CheckClose(got[0].a2, ref.a2, 2e-4, tag + " a2");
                }
            }
        }
    }
}

void CheckHighOrderTypes() {
    std::array<ZPK, kMaxSections> proto{};
    std::array<ZPK, kMaxSections> analog{};

    auto run = [&](std::string_view tag, size_t n, double sigma, std::string_view what) {
        std::copy_n(proto.begin(), n, analog.begin());
        IIRDesign::Ivantsov(std::span{analog}.first(n), kFs, sigma);
        CheckStableAndFinite(std::span{analog}.first(n), std::format("{} {}", tag, what));
        double const err = MaxDbError(std::span{proto}.first(n), std::span{analog}.first(n), -20.0);
        Check(err < 1.5, std::format("{} {}: passband/transition max dB error {:.3f}", tag, what, err));
    };

    // 4th-order Butterworth lowpass / highpass
    IIRDesign::Butterworth(std::span{proto}.first(2), 2);
    IIRDesign::ProtyleToLowpass(std::span{proto}.first(2), 2, 2.0 * kPi * 5000.0);
    run("LP", 2, kSigmaOpt, "sigma_opt");

    std::fill(proto.begin(), proto.end(), ZPK{});
    IIRDesign::Butterworth(std::span{proto}.first(2), 2);
    IIRDesign::ProtyleToHighpass(std::span{proto}.first(2), 2, 2.0 * kPi * 5000.0);
    run("HP", 2, kSigmaOpt, "sigma_opt");

    // 4th-order bandpass / bandstop
    std::fill(proto.begin(), proto.end(), ZPK{});
    IIRDesign::Butterworth(std::span{proto}.first(2), 2);
    IIRDesign::ProtyleToBandpass(std::span{proto}.first(4), 2, 2.0 * kPi * 5000.0, 2.0);
    run("BP", 4, kSigmaOpt, "sigma_opt");

    std::fill(proto.begin(), proto.end(), ZPK{});
    IIRDesign::Butterworth(std::span{proto}.first(2), 2);
    IIRDesign::ProtyleToBandstop(std::span{proto}.first(4), 2, 2.0 * kPi * 5000.0, 2.0);
    run("BS", 4, kSigmaOpt, "sigma_opt");
}

void CheckChebyshev2Shelf() {
    constexpr size_t kPairs = 4; // order 8
    std::array<ZPK, kMaxSections> proto{};
    std::array<ZPK, kMaxSections> analog{};
    bool const ok = IIRDesignShelf::Chebyshev2(std::span{proto}.first(kPairs), kPairs,
                                               12.0, 1.0, false);
    Check(ok, "Chebyshev2 shelf design succeeds");
    if (!ok) {
        return;
    }
    // shelf 原型在 1 rad/s; 用预畸变尺度把中点搬到 1 kHz
    double const wc = 2.0 * kPi * 1000.0 / kFs;
    double const omega = PrewarpedOmega(wc, kSigmaOpt);
    double const scale = 2.0 * kPi * kFs / omega;
    IIRDesign::ProtyleToLowpass(std::span{proto}.first(kPairs), kPairs, scale);

    std::copy_n(proto.begin(), kPairs, analog.begin());
    IIRDesign::Ivantsov(std::span{analog}.first(kPairs), kFs, kSigmaOpt);
    CheckStableAndFinite(std::span{analog}.first(kPairs), "Chebyshev2 shelf");
    double const err = MaxDbError(std::span{proto}.first(kPairs), std::span{analog}.first(kPairs));
    Check(err < 0.5, std::format("Chebyshev2 shelf: max dB error {:.3f}", err));
    std::cout << std::format("  shelf max dB error = {:.4f}\n", err);
}
} // namespace

int main() {
    CheckDefaultOverload();
    CheckPrewarp();
    CheckBilinearLimit();
    CheckSingleSectionMatchesBiquad();
    CheckHighOrderTypes();
    CheckChebyshev2Shelf();
    if (g_failures == 0) {
        std::cout << "PASS  IIRDesign::Ivantsov\n";
        return 0;
    }
    std::cout << std::format("FAIL  {} check(s)\n", g_failures);
    return 1;
}
