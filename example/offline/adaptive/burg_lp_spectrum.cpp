// ------------------------------------------------------------
// Burg LPC（标准 + 全通扭曲）的谱估计验证
//
// 对 wormhole.wav 的一段（48 kHz，起点 25013，长度 14217 ≈ 0.30 s）做 Burg LPC，
// 用 matplot++ 把「窗化 DFT 幅度」与几条 **LPC 全极点包络**（标准 Burg + 不同全通
// 扭曲参数的 warped Burg）画在一张图上：
//
//     work_dir/output/burg_lp_spectrum.png
//
// ⚠ matplot++ 靠 gnuplot 落地图片，**运行期需要 PATH 里有 gnuplot 5.2.6+**（见 docs/matplotpp.md）。
//
// 扭曲模型：用一阶全通 A(z) = (a + z^-1)/(1 + a z^-1) 把频率轴扭曲后再做 Burg
// （@ref qwqdsp_adaptive::WarpedBurgLP）。a > 0 压缩低频、a < 0 压缩高频。
// 该模型的频响要把 z^-1 换成 A(z)：H(ω) = 1 / A(A(e^{-jw}))。
//
// 归一化约定：所有曲线都以「DFT 在 20 Hz..20 kHz 内的幅度峰值」为 0 dB。
// 每条 LPC 包络的**绝对增益**在这种比较里不可辨识（激励不是白噪声），所以各自乘一个
// 常数，使其在 DFT 的谐波峰上平均对齐——对齐增益与加权残差都会打印出来，
// 残差才是"形状对不对"的指标。
//
// 读图：看的是**包络形状**——LPC 应贴着强谐波走；在弱谐波（激励造成的凹陷）处
// LPC 偏高是正常的（LPC 拟合共振峰包络，不是单根谐波）。不同 a 的曲线在低频/高频
// 段的疏密不同，就是频率轴被扭曲的结果。
// ------------------------------------------------------------
#include "AudioFile.h"
#include "work_dir.hpp"

#include <qwqdsp/adaptive/burg_lp.hpp>
#include <qwqdsp/adaptive/warped_burg_lp.hpp>
#include <qwqdsp/spectral/real_fft_adv.hpp>
#include <qwqdsp/window/hann.hpp>

#include <matplot/matplot.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <iostream>
#include <numbers>
#include <span>
#include <string>
#include <vector>

namespace {
/// 片段：起点与长度（采样点，48 kHz）
constexpr size_t kStart = 25013;
constexpr size_t kLength = 14217;

/// 画图/统计用的频率范围
constexpr double kFreqMin = 20.0;
constexpr double kFreqMax = 20000.0;
constexpr size_t kNumLogPoints = 2048;

/// 参与对齐拟合的 DFT 峰最低电平（相对带内峰值）——太低会把噪声尖峰也算进来
constexpr double kPeakFloorDb = -30.0;

/// 要叠加的扭曲参数：0 = 标准 Burg（不扭曲），其余为一阶全通系数 a
struct WarpCurve {
    double warp;
    char const* label;
    std::array<float, 3> color;
};
constexpr std::array<WarpCurve, 4> kWarpCurves{{
    {0.0, "LPC (Burg, no warp)", {0.910f, 0.384f, 0.173f}},  // 橙
    {0.3, "LPC (warp a=+0.30)", {0.145f, 0.620f, 0.243f}},   // 绿
    {0.6, "LPC (warp a=+0.60)", {0.760f, 0.243f, 0.713f}},   // 紫
    {-0.3, "LPC (warp a=-0.30)", {0.000f, 0.635f, 0.760f}},  // 青
}};

size_t NextPow2(size_t n) noexcept {
    size_t p = 1;
    while (p < n) {
        p <<= 1;
    }
    return p;
}

double ToDb(double v) noexcept {
    return 20.0 * std::log10(std::max(v, 1e-30));
}
} // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IONBF, 0);

    size_t const order = (argc > 1) ? static_cast<size_t>(std::atoi(argv[1])) : 20;
    if (order == 0 || order >= kLength) {
        std::cout << std::format("order 非法: {}\n", order);
        return 1;
    }

    // ---- 读音频、取片段 ----
    AudioFile<float> audio;
    if (!audio.load(qwqdsp_support::WormholeWav())) {
        std::cout << "FAIL: cannot load wormhole.wav\n";
        return 1;
    }
    float const fs = audio.getSampleRate();
    auto const& samples = audio.samples[0];
    if (kStart + kLength > samples.size()) {
        std::cout << std::format("FAIL: 片段超出音频长度 ({} > {})\n", kStart + kLength, samples.size());
        return 1;
    }
    std::vector<float> x(samples.begin() + kStart, samples.begin() + kStart + kLength);
    double rms = 0.0;
    for (float v : x) {
        rms += static_cast<double>(v) * v;
    }
    rms = std::sqrt(rms / static_cast<double>(x.size()));
    std::cout << std::format("片段: [{}, {}) @ {} Hz, {} 点, RMS {:.4f}, LPC order {}\n", kStart,
                             kStart + kLength, static_cast<int>(fs), x.size(), rms, order);

    // ---- 窗化 DFT ----
    size_t const fft_size = NextPow2(kLength);
    qwqdsp_spectral::RealFftAdv fft;
    fft.Init(fft_size);
    std::vector<float> frame(fft_size, 0.0f);
    std::copy(x.begin(), x.end(), frame.begin());
    qwqdsp_window::Hann::ApplyWindow(std::span{frame}.first(kLength), true);
    std::vector<float> dft_mag(fft.NumBins());
    fft.FFTGainPhase(frame, dft_mag);

    auto bin_freq = [&](size_t b) {
        return static_cast<double>(b) * static_cast<double>(fs) / static_cast<double>(fft_size);
    };
    auto in_band = [](double f) { return f >= kFreqMin && f <= kFreqMax; };

    // ---- 带内 DFT 峰值，作为所有曲线的 0 dB 基准，也是对齐用的样本点 ----
    double dft_peak = 0.0;
    for (size_t b = 1; b + 1 < dft_mag.size(); ++b) {
        if (in_band(bin_freq(b))) {
            dft_peak = std::max(dft_peak, static_cast<double>(dft_mag[b]));
        }
    }
    if (dft_peak <= 0.0) {
        std::cout << "FAIL: 带内 DFT 峰值为 0\n";
        return 1;
    }
    std::vector<std::pair<double, double>> peaks; // {freq, dft_mag}
    double const peak_floor = dft_peak * std::pow(10.0, kPeakFloorDb / 20.0);
    for (size_t b = 1; b + 1 < dft_mag.size(); ++b) {
        if (!in_band(bin_freq(b))) {
            continue;
        }
        if (dft_mag[b] >= dft_mag[b - 1] && dft_mag[b] > dft_mag[b + 1] && dft_mag[b] >= peak_floor) {
            peaks.emplace_back(bin_freq(b), static_cast<double>(dft_mag[b]));
        }
    }

    std::vector<double> freq_grid(kNumLogPoints);
    for (size_t i = 0; i < kNumLogPoints; ++i) {
        double const t = static_cast<double>(i) / static_cast<double>(kNumLogPoints - 1);
        freq_grid[i] = kFreqMin * std::pow(kFreqMax / kFreqMin, t);
    }

    // ---- 对每条曲线：求 k → A(z) → 扭曲轴上的包络 → 对齐增益 ----
    using namespace matplot;
    std::vector<double> dft_x(dft_mag.size());
    std::vector<double> dft_y(dft_mag.size());
    for (size_t b = 0; b < dft_mag.size(); ++b) {
        dft_x[b] = bin_freq(b);
        dft_y[b] = ToDb(static_cast<double>(dft_mag[b]) / dft_peak);
    }
    auto const dft_line = semilogx(dft_x, dft_y, "-");
    dft_line->line_width(0.7f).color({0.231f, 0.490f, 0.847f}); // #3b7dd8
    hold(on);

    std::cout << std::format("\n{:<22s} {:>9s} {:>12s} {:>16s}\n", "curve", "gain(dB)", "resid(dB)", "worst(dB)");
    std::vector<std::string> labels{"DFT (Hann window)"};
    for (auto const& wc : kWarpCurves) {
        // 扭曲 Burg 求反射系数（a=0 时退化为普通 Burg；同一套边界处理便于横向比较）
        qwqdsp_adaptive::WarpedBurgLP burg;
        burg.Init(x.size());
        burg.SetWarp(static_cast<float>(wc.warp));
        std::vector<float> k(order);
        {
            auto x_work = x; // Process 会就地修改 x
            burg.Process(x_work, k);
        }
        std::vector<float> a_coeff(order + 1);
        std::vector<float> a_rev(order + 1);
        qwqdsp_adaptive::BurgLP::Lattice2Tf_KeepK(k, a_coeff, a_rev);

        // 扭曲轴上的包络：H(ω) = 1 / A(A(e^{-jw}))
        auto envelope = [&](double freq) {
            auto const zd = burg.WarpDelay(static_cast<float>(freq), fs); // A(e^{-jw})
            std::complex<double> acc{0.0, 0.0};
            std::complex<double> w{1.0, 0.0};
            for (size_t i = 0; i <= order; ++i) {
                acc += static_cast<double>(a_coeff[i]) * w;
                w *= std::complex<double>{zd.real(), zd.imag()};
            }
            return 1.0 / std::abs(acc);
        };

        // 在 DFT 谐波峰上拟合单一增益常数（按功率加权，强峰主导）
        double weight_sum = 0.0;
        double gain_db = 0.0;
        for (auto const& [f, m] : peaks) {
            double const w = (m / dft_peak) * (m / dft_peak);
            weight_sum += w;
            gain_db += w * ToDb(m / envelope(f));
        }
        gain_db = (weight_sum > 0.0) ? gain_db / weight_sum : 0.0;
        double const gain = std::pow(10.0, gain_db / 20.0);

        double err_mean = 0.0;
        double err_var = 0.0;
        double worst = 0.0;
        for (auto const& [f, m] : peaks) {
            double const w = (m / dft_peak) * (m / dft_peak);
            err_mean += w * ToDb(m / (gain * envelope(f)));
        }
        err_mean = (weight_sum > 0.0) ? err_mean / weight_sum : 0.0;
        for (auto const& [f, m] : peaks) {
            double const w = (m / dft_peak) * (m / dft_peak);
            double const d = ToDb(m / (gain * envelope(f))) - err_mean;
            err_var += w * d * d;
            worst = std::max(worst, std::abs(d + err_mean));
        }
        double const err_std = (weight_sum > 0.0) ? std::sqrt(std::max(0.0, err_var / weight_sum)) : 0.0;
        std::cout << std::format("{:<22s} {:+9.2f} {:7.2f} ± {:<5.2f} {:>15.2f}\n", wc.label, gain_db, err_mean,
                                 err_std, worst);
        labels.emplace_back(wc.label);

        std::vector<double> y(kNumLogPoints);
        for (size_t i = 0; i < kNumLogPoints; ++i) {
            y[i] = ToDb(gain * envelope(freq_grid[i]) / dft_peak);
        }
        auto const line = semilogx(freq_grid, y, "-");
        line->line_width(1.8f).color({wc.color[0], wc.color[1], wc.color[2]});
    }
    hold(off);

    std::cout << std::format("DFT 谐波峰 {} 个（>= {:.0f} dB，按功率加权）\n", peaks.size(), kPeakFloorDb);

    xlim({kFreqMin, kFreqMax});
    ylim({-80.0, 8.0});
    xlabel("Frequency (Hz)");
    ylabel("Magnitude (dB, normalized)");
    title("wormhole.wav [25013, +14217) - DFT vs Burg LPC envelopes (order " + std::to_string(order) + ")");
    // 图例放到绘图区外（gnuplot: set key outside right top），加宽画布给它留位置
    gcf()->size(1000, 450);
    auto const leg = legend(labels);
    leg->inside(false);
    grid(on);

    auto out_dir = qwqdsp_support::GetOutputDir();
    std::filesystem::create_directories(out_dir);
    auto png_path = (out_dir / "burg_lp_spectrum.png").generic_string();
    save(png_path);
    std::cout << std::format("PNG: {}\n", png_path);
    return 0;
}
