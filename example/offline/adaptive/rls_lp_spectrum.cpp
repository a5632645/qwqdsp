// ------------------------------------------------------------
// RLS 线性预测的谱估计验证
//
// 对 wormhole.wav 的一段（48 kHz，起点 25013，长度 14217 ≈ 0.30 s）跑自适应 RLS
// （qwqdsp_adaptive::RLSFIlter），用 matplot++ 把「窗化 DFT 幅度」与 **RLS 估出的
// 全极点包络** 画在一张图上：
//
//     work_dir/output/rls_lp_spectrum.png
//
// ⚠ matplot++ 靠 gnuplot 落地图片，**运行期需要 PATH 里有 gnuplot 5.2.6+**（见 docs/matplotpp.md）。
//
// 喂法：一步预测 `Tick(x[n], x[n+1])`——source/target 取同一段信号、错开一个采样。
// 注意 `RLSFIlter::Tick` 是先把 source 移入延迟线再做预测（回归量含「当前样本」），
// 所以 target 必须给**后一个**样本：若 source 与 target 同值，最小二乘存在零误差解
// w = e0（把当前样本直接抄给输出），RLS 必然收敛到它，系数退化成恒等、不含任何谱信息。
//
// 收敛到的 w 对应预测器 `x̂[m] = Σ_j w[j-1] x[m-j]`，故全极点模型
//
//     H(z) = 1 / (1 - Σ_j w[j-1] z^-j)
//
// 这与 `RLSFIlter::Filter()` 的合成结构一致（该函数用过去的输出做回归量）。
//
// 遗忘因子 λ 决定 RLS 的有效记忆 `1/(1-λ)` 个采样，默认 0.999（≈1000 点 ≈ 21 ms，
// 也是类的默认值）。本段 0.30 s 内是非平稳语音：λ 拉大到 0.9999 / 0.99999（记忆
// ≈1e4 / 1e5 点）会让模型变成整段的折中、包络被抹平，实测一步预测残差从 −39.5 dB
// 退化到 −21.2 / −19.0 dB。λ=1.0（无遗忘=整段 LS）配 P0=1e6 则给出 −48.9 dB，
// 与直接解同一套回归量的法方程（−49.3 dB）一致。
//
// ⚠ 初始逆相关矩阵 P0（`Init(p0)`，类默认 0.01）对本信号太小：P0 相当于"w=0 的先验
// 有多强"，0.01 会让收敛慢到 1.4 万点都到不了 LS 解——实测 order 20、λ=1.0 时尾部
// 残差只有 −18.7 dB（真值 −49.3 dB），画出来的包络是一条平滑低频斜坡、完全不贴共振峰。
// 这里取 1e6（等于几乎没有先验，直接奔 LS 解）。
//
// 归一化约定：所有曲线都以「DFT 在 20 Hz..20 kHz 内的幅度峰值」为 0 dB。包络的
// **绝对增益**在这种比较里不可辨识（激励不是白噪声），所以乘一个常数，使其在 DFT 的
// 谐波峰上按功率加权平均对齐——对齐增益与加权残差都会打印出来，残差才是"形状对不对"
// 的指标。读图：RLS 包络应贴着强谐波走；在弱谐波（激励造成的凹陷）处偏高是正常的。
// ------------------------------------------------------------
#include "AudioFile.h"
#include "matplot_headless.hpp"
#include "work_dir.hpp"

#include <qwqdsp/adaptive/rls_filter.hpp>
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

/// 默认遗忘因子：1/(1-λ) ≈ 1000 点（≈21 ms），类默认值
constexpr double kDefaultForget = 0.999;

/// 初始逆相关矩阵的对角值（`RLSFIlter::Init`）。类默认 0.01 太强，收敛不到 LS 解
constexpr double kRlsP0 = 1e6;

/// 打印残差统计时只看最后这么多点（RLS 已收敛）
constexpr size_t kTailSamples = 2000;

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

/**
 * @brief 跑一遍 RLS 的一步预测，取出收敛后的系数
 * @param x 输入片段，用 `x[n] -> x[n+1]` 的一步预测喂 Tick
 * @param forget 遗忘因子 λ
 * @param w 输出预测系数，`w[j-1] = a_j`，对应 `x̂[m] = Σ_j a_j x[m-j]`
 * @param err_tail_db 输出尾部残差 RMS（相对片段 RMS，dB）
 */
template <int kOrder>
void RunRls(std::vector<float> const& x, double forget, std::array<double, kOrder>& w, double& err_tail_db) {
    double rms = 0.0;
    for (float v : x) {
        rms += static_cast<double>(v) * v;
    }
    rms = std::sqrt(rms / static_cast<double>(x.size()));

    qwqdsp_adaptive::RLSFIlter<kOrder> rls;
    rls.Init(kRlsP0);
    rls.SetForgetParam(static_cast<float>(forget));
    double err2 = 0.0;
    size_t ne = 0;
    for (size_t n = 0; n + 1 < x.size(); ++n) {
        double const pred = rls.Tick(x[n], x[n + 1]);
        if (n + kTailSamples > x.size()) {
            double const e = static_cast<double>(x[n + 1]) - pred;
            err2 += e * e;
            ++ne;
        }
    }
    err_tail_db = ToDb(std::sqrt(err2 / static_cast<double>(ne)) / rms);

    auto const& ww = rls.GetCoeff();
    for (int i = 0; i < kOrder; ++i) {
        w[i] = ww[i];
    }
}
} // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IONBF, 0);

    size_t const order = (argc > 1) ? static_cast<size_t>(std::atoi(argv[1])) : 20;
    double const forget = (argc > 2) ? std::atof(argv[2]) : kDefaultForget;
    if (order == 0 || order > kLength || forget <= 0.0 || forget > 1.0) {
        std::cout << std::format("参数非法: order {} forget {}\n", order, forget);
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
    std::cout << std::format("片段: [{}, {}) @ {} Hz, {} 点, RMS {:.4f}, RLS order {} lambda {}\n", kStart,
                             kStart + kLength, static_cast<int>(fs), x.size(), rms, order, forget);

    // ---- RLS 一步预测（order 是模板参数，用编译期派发） ----
    std::vector<double> w;
    double err_tail_db = 0.0;
    auto run = [&](auto tag) {
        constexpr int kOrder = decltype(tag)::value;
        std::array<double, kOrder> a{};
        RunRls<kOrder>(x, forget, a, err_tail_db);
        w.assign(a.begin(), a.end());
    };
    switch (order) {
    case 8: run(std::integral_constant<int, 8>{}); break;
    case 12: run(std::integral_constant<int, 12>{}); break;
    case 16: run(std::integral_constant<int, 16>{}); break;
    case 20: run(std::integral_constant<int, 20>{}); break;
    case 24: run(std::integral_constant<int, 24>{}); break;
    case 32: run(std::integral_constant<int, 32>{}); break;
    case 48: run(std::integral_constant<int, 48>{}); break;
    case 64: run(std::integral_constant<int, 64>{}); break;
    default:
        std::cout << std::format("order 非法（支持 8/12/16/20/24/32/48/64）: {}\n", order);
        return 1;
    }
    std::cout << std::format("尾部 {} 点一步预测残差: {:.2f} dB（相对片段 RMS），w0 = {:+.4f}\n", kTailSamples,
                             err_tail_db, w[0]);

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

    // ---- 全极点包络 H(ω) = 1 / |1 - Σ_j w[j-1] z^-j| ----
    auto envelope = [&](double freq) {
        double const wr = 2.0 * std::numbers::pi_v<double> * freq / static_cast<double>(fs);
        std::complex<double> const z1 = std::polar(1.0, -wr); // z^-1
        std::complex<double> zp = z1;
        std::complex<double> acc{1.0, 0.0};
        for (size_t j = 1; j <= w.size(); ++j) {
            acc -= w[j - 1] * zp;
            zp *= z1;
        }
        return 1.0 / std::abs(acc);
    };

    // ---- 在 DFT 谐波峰上拟合单一增益常数（按功率加权，强峰主导） ----
    double weight_sum = 0.0;
    double gain_db = 0.0;
    for (auto const& [f, m] : peaks) {
        double const wt = (m / dft_peak) * (m / dft_peak);
        weight_sum += wt;
        gain_db += wt * ToDb(m / envelope(f));
    }
    gain_db = (weight_sum > 0.0) ? gain_db / weight_sum : 0.0;
    double const gain = std::pow(10.0, gain_db / 20.0);

    double err_mean = 0.0;
    double err_var = 0.0;
    double worst = 0.0;
    for (auto const& [f, m] : peaks) {
        double const wt = (m / dft_peak) * (m / dft_peak);
        err_mean += wt * ToDb(m / (gain * envelope(f)));
    }
    err_mean = (weight_sum > 0.0) ? err_mean / weight_sum : 0.0;
    for (auto const& [f, m] : peaks) {
        double const wt = (m / dft_peak) * (m / dft_peak);
        double const d = ToDb(m / (gain * envelope(f))) - err_mean;
        err_var += wt * d * d;
        worst = std::max(worst, std::abs(d + err_mean));
    }
    double const err_std = (weight_sum > 0.0) ? std::sqrt(std::max(0.0, err_var / weight_sum)) : 0.0;

    std::string const label = std::format("RLS all-pole (order {}, lambda {})", order, forget);
    std::cout << std::format("\n{:<22s} {:>9s} {:>12s} {:>16s}\n", "curve", "gain(dB)", "resid(dB)", "worst(dB)");
    std::cout << std::format("{:<22s} {:+9.2f} {:7.2f} ± {:<5.2f} {:>15.2f}\n", label, gain_db, err_mean, err_std,
                             worst);
    std::cout << std::format("DFT 谐波峰 {} 个（>= {:.0f} dB，按功率加权）\n", peaks.size(), kPeakFloorDb);

    // ---- 画图 ----
    using namespace matplot;
    qwqdsp_support::MakeMatplotHeadless(); // 见该函数注释：不切就会每改一次图就重画一遍（闪窗 + 慢）
    std::vector<double> dft_x(dft_mag.size());
    std::vector<double> dft_y(dft_mag.size());
    for (size_t b = 0; b < dft_mag.size(); ++b) {
        dft_x[b] = bin_freq(b);
        dft_y[b] = ToDb(static_cast<double>(dft_mag[b]) / dft_peak);
    }
    auto const dft_line = semilogx(dft_x, dft_y, "-");
    dft_line->line_width(0.7f).color({0.231f, 0.490f, 0.847f}); // #3b7dd8
    hold(on);

    std::vector<double> env_y(kNumLogPoints);
    for (size_t i = 0; i < kNumLogPoints; ++i) {
        env_y[i] = ToDb(gain * envelope(freq_grid[i]) / dft_peak);
    }
    auto const env_line = semilogx(freq_grid, env_y, "-");
    env_line->line_width(1.8f).color({0.910f, 0.384f, 0.173f}); // 橙
    hold(off);

    xlim({kFreqMin, kFreqMax});
    ylim({-80.0, 8.0});
    xlabel("Frequency (Hz)");
    ylabel("Magnitude (dB, normalized)");
    title("wormhole.wav [25013, +14217) - DFT vs RLS all-pole envelope (order " + std::to_string(order) +
          ", lambda " + std::format("{}", forget) + ")");
    gcf()->size(1000, 450);
    legend(std::vector<std::string>{"DFT (Hann window)", label});
    grid(on);

    auto out_dir = qwqdsp_support::GetOutputDir();
    std::filesystem::create_directories(out_dir);
    auto png_path = (out_dir / "rls_lp_spectrum.png").generic_string();
    save(png_path);
    std::cout << std::format("PNG: {}\n", png_path);
    return 0;
}
