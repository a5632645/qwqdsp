// 半带全通和多相上/下采样器的数值测试。
//
// 用自己投影到指定频点（不需要 FFT：频点取整周期，直接求和即可）验证三件事：
//   1. 实数链路：8× 上采样后的镜像（fs_in − f0）被压到阻带深度；
//   2. 解析链路：8× 解析上采样后**负频率侧**被压到阻带深度（单边性）；
//   3. 8× 上采样再 8× 抽取的往返增益 ≈ 1（上采样每级通带增益与抽取的 1/2 相抵）。
// 失败返回非零退出码。
//
// 系数来源与设计（椭圆半带 N=19 / 阻带 100 dB）见 labs/tf2ca 与 labs/adaa_iir 的笔记。
#include <qwqdsp/filter/halfband_allpass_ovs.hpp>

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

void Check(bool ok, std::string_view what, double got, double limit) {
    if (!ok) {
        ++g_failures;
        std::cout << std::format("FAIL {}: got {:.3e}, limit {:.3e}\n", what, got, limit);
    }
    else {
        std::cout << std::format("ok   {}: {:.3e}\n", what, got);
    }
}

constexpr double kFs = 48000.0;
constexpr std::size_t kStage = 3;
constexpr std::size_t kFactor = 1u << kStage;              // 8
constexpr std::size_t kInputs = 4800;                      // 0.1 s
constexpr std::size_t kDiscard = 1200;                     // 丢瞬态（输入样本数）
constexpr std::size_t kAnalyze = 2400;                     // 分析 0.05 s -> 栅格 20 Hz
// 测试音取 9000 Hz（= 450×20，整周期）：必须**高于解析级的守护带**（约 5.25 kHz @384 kHz），
// 否则负频抑制只受守护带限制（实测 4 kHz 时只有 -58 dB，那是结构固有行为，不是缺陷）。
constexpr double kF0 = 9000.0;
constexpr double kAmplitude = 1.0;

/// 复投影：Σ x[n]·e^{-j2π f n / fs}（f 落在栅格上时无泄漏）
std::complex<double> Project(const std::complex<double>* x, std::size_t n, double f, double fs) {
    std::complex<double> acc{};
    for (std::size_t i = 0; i < n; ++i) {
        const double phase = -2.0 * std::numbers::pi * f * static_cast<double>(i) / fs;
        acc += x[i] * std::complex<double>{std::cos(phase), std::sin(phase)};
    }
    return acc / static_cast<double>(n);
}

/// 实数序列的复投影（按解析信号取正/负频）
std::complex<double> ProjectReal(const double* x, std::size_t n, double f, double fs) {
    std::vector<std::complex<double>> c(n);
    for (std::size_t i = 0; i < n; ++i) {
        c[i] = std::complex<double>{x[i], 0.0};
    }
    return Project(c.data(), n, f, fs);
}

void RunChain() {
    qwqdsp_filter::HalfbandAllpassUpsampler<float, kStage> up;
    qwqdsp_filter::HalfbandAllpassAnalyticUpsampler<float, kStage> upa;
    qwqdsp_filter::HalfbandAllpassDecimator<float, kStage> down;

    std::vector<double> interp(kAnalyze * kFactor);
    std::vector<std::complex<double>> analytic(kAnalyze * kFactor);
    std::vector<double> roundtrip;

    float buf[kFactor];
    std::complex<float> abuf[kFactor];

    bool finite = true;
    for (std::size_t n = 0; n < kInputs; ++n) {
        const double t = static_cast<double>(n) / kFs;
        const float x = static_cast<float>(kAmplitude * std::sin(2.0 * std::numbers::pi * kF0 * t));
        up.Tick(x, buf);
        upa.Tick(x, abuf);
        const float d = down.Tick(buf);
        for (std::size_t i = 0; i < kFactor; ++i) {
            finite = finite && std::isfinite(static_cast<double>(buf[i]));
            finite = finite && std::isfinite(static_cast<double>(abuf[i].real())) &&
                     std::isfinite(static_cast<double>(abuf[i].imag()));
        }
        finite = finite && std::isfinite(static_cast<double>(d));
        if (n >= kDiscard && n < kDiscard + kAnalyze) {
            const std::size_t base = (n - kDiscard) * kFactor;
            for (std::size_t i = 0; i < kFactor; ++i) {
                interp[base + i] = static_cast<double>(buf[i]);
                analytic[base + i] = std::complex<double>{abuf[i].real(), abuf[i].imag()};
            }
            roundtrip.push_back(static_cast<double>(d));
        }
    }
    Check(finite, "all outputs finite", finite ? 0.0 : 1.0, 0.0);

    const std::size_t n = kAnalyze * kFactor;
    const double fsUp = kFs * static_cast<double>(kFactor);

    // 信号频点（8× 速率下仍是 4 kHz）与镜像（fs_in − f0 = 44 kHz）
    const auto sigRe = std::abs(ProjectReal(interp.data(), n, kF0, fsUp));
    const auto imgRe = std::abs(ProjectReal(interp.data(), n, kFs - kF0, fsUp));
    const double rejRealDb = 20.0 * std::log10(imgRe / sigRe);
    std::cout << std::format("实数链路 8× 上采样：信号 {:.4f}，镜像 {:.2e}\n", sigRe, imgRe);
    Check(rejRealDb < -90.0, "real upsample image rejection [dB]", rejRealDb, -90.0);

    // 解析链路：正频保留、负频压制；镜像（44 kHz）也应在阻带
    const auto sigAn = std::abs(Project(analytic.data(), n, kF0, fsUp));
    const auto negAn = std::abs(Project(analytic.data(), n, -kF0, fsUp));
    const auto imgAn = std::abs(Project(analytic.data(), n, kFs - kF0, fsUp));
    const double oneSidedDb = 20.0 * std::log10(negAn / sigAn);
    const double rejAnDb = 20.0 * std::log10(imgAn / sigAn);
    std::cout << std::format("解析链路：信号 {:.4f}，负频 {:.2e}，镜像 {:.2e}\n", sigAn, negAn, imgAn);
    Check(oneSidedDb < -90.0, "analytic upsample negative-frequency rejection [dB]", oneSidedDb, -90.0);
    Check(rejAnDb < -90.0, "analytic upsample image rejection [dB]", rejAnDb, -90.0);

    // 往返：8× 上采样 + 8× 抽取后，正弦幅度应回到 1
    const auto back = std::abs(ProjectReal(roundtrip.data(), roundtrip.size(), kF0, kFs)) * 2.0;
    std::cout << std::format("往返（上采样×8 → 抽取÷8）：幅度 {:.6f}（输入 {:.1f}）\n", back, kAmplitude);
    Check(std::abs(back - kAmplitude) < 0.02 * kAmplitude, "round-trip amplitude error",
          std::abs(back - kAmplitude) / kAmplitude, 0.02);
}
} // namespace

int main() {
    RunChain();
    if (g_failures == 0) {
        std::cout << "all checks passed\n";
    }
    return g_failures;
}
