#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <vector>
#include <qwqdsp/misc/smoother.hpp>

#include "raylib.h"

/**
 * @brief 论文滑动无窗 DFT 的 NC 方法谱帧
 *
 * 与 nc_reassignment_frame.hpp (线性 FFT bin 对) 不同：本方法以频率网格为基准，
 * y 轴每个像素对应一个 NC bin，bin 中心频率与相应 y 坐标的显示频率相同，
 * 一一映射。每个 bin 采用独立窗长 N_k，并用递归滑动 DFT 逐样本演进，
 * 输出无需窗函数与相位校正。
 *
 * 关键公式（论文 arXiv:2410.07982v3）：
 *   (1)  bin 中心频率:  f(n) = 440·2^((n-69)/12)
 *   (7)  窗长:          N_k = round(round(2·f_c/W_NC)·F_S/(2·f_c))
 *   (5)  左右分量频率:  f_left = f_c − F_S/(2N), f_right = f_c + F_S/(2N)
 *   (8)  NC 幅度:       max(0, −(Re_L·Re_R + Im_L·Im_R))，除以 N 归一化
 *
 * 递归滑动 DFT（相位锚定窗口起点，无相位校正）：
 *   X(n) = W·X(n−1) + x[n] − x[n−N]·W^N,   W = e^{−j·2π·f/F_S}
 *
 * @ref https://arxiv.org/html/2410.07982v3
 */
template <typename Colormap>
struct WindowlessNcFrame {
    /**
     * @brief 每个 NC bin 的硬件描述与滑动 DFT 状态
     */
    struct Bin {
        float f_center{};        // 中心频率(Hz, 仅用于显示/映射)
        double f_left{};         // 左分量参考频率(Hz)
        double f_right{};        // 右分量参考频率(Hz)
        int N{};                 // 窗长(样本数)
        int row{};               // 对应输出行 y(=bin 索引, 用于回填 binDb_)
        // 旋转因子与累加器用 double：float 会随时间累积慢速漂移(见论文 IV-B 节)，
        // double 使滑动 DFT 幅度有界无趋势，消除长期数值误差。
        std::complex<double> Wl{}, WNr_l{};    // 左分量旋转因子 W 与 W^N
        std::complex<double> Wr{}, WNr_r{};    // 右分量旋转因子 W 与 W^N
        std::complex<double> acc_l{}, acc_r{}; // 左右分量滑动 DFT 累加器
        float smoothed_gain{};   // 一阶 IIR 平滑后的线性增益(非 dB)；EMA 在增益域线性作用
        float alpha{};           // 每样本 EMA 系数(抗混叠箱<1, 旁路箱=1)
    };

    /**
     * @brief 初始化 bin 布局与输出行映射
     *
     * @param sampleRate     采样率(Hz)
     * @param fftSize        帧大小(由上层 SpectrogramColumn 提供，本方法只用作数据源)
     * @param hopSize        帧间步进(样本数)；仅第二个及之后的 Process 调用会新增 hopSize 个样本
     * @param zeroPad        零填充因子（与现有 NC 帧接口一致；本方法无零填充，保留占位）
     * @param outputHeight   输出行数（= 画布高度）
     * @param freqMin        最低显示频率(Hz)
     * @param freqMax        最高显示频率(Hz)
     * @param dbFloor        幅度下限(dB)
     * @param bandwidthScale NC bin 带宽缩放系数(默认 1.0)。
     *                       带宽 = 相邻像素频率差 × 此系数：
     *                       > 1 带宽更宽 → 窗长 N 更短 → 频率分辨率更低、时间响应更快；
     *                       < 1 带宽更窄 → 窗长 N 更长 → 频率分辨率更高、时间响应更慢。
     *
     * 说明：不施加统一平滑。滑动 DFT 的窗长 N 本身已是 FIR 时间平滑(分辨率 ≈ N/F_S)，
     * 图像输出采样间隔为 hop/F_S。只有当 N < hop 时箱的时间分辨率才低于图像采样率，
     * 才可能产生时间混叠。因此逐箱自适应：tau_k = max(N_k, hopSize)，使得每个箱的
     * 有效时间分辨率 ≥ 图像采样间隔；N ≥ hop 的箱 tau=N_k，几乎不叠加额外平滑。
     */
    void Init(int sampleRate, int fftSize, int hopSize, int zeroPad, int outputHeight, float freqMin, float freqMax,
              float dbFloor, float bandwidthScale = 1.0f) noexcept {
        sampleRate_ = sampleRate;
        fftSize_ = fftSize;
        hopSize_ = hopSize;
        zeroPad_ = zeroPad;   // 占位，无零填充
        outputHeight_ = outputHeight;
        freqMin_ = freqMin;
        freqMax_ = freqMax;
        dbFloor_ = dbFloor;
        logMin_ = std::log10(freqMin);
        logMax_ = std::log10(freqMax);

        // 幅度归一化：NC 前向公式 gain = sqrt(ncSum)/N 对幅度 A 的纯音，
        // 稳态增益 ≈ A/π（推导：ncSum ≈ (A/2)²·cos(π/N)·N²/2，gain ≈ A/π）。
        // 乘 π 让纯音峰值回到 ≈ 0 dB（与其他算法刻度一致），避免整体偏低 ~10 dB。
        gainNorm_ = static_cast<float>(std::numbers::pi_v<double>);

        // 低音窗长上限(论文 IV-A 节)
        const int maxWindowSamples = std::max(8, static_cast<int>(kMaxWindowS * sampleRate));

        // ── 以频率网格为基准：y 轴每个像素对应一个 NC bin，bin 中心频率 = 该像素频率 ──
        // y=0 在顶部(高频)，y 增大频率递减；bin 数与输出高度一致，一一映射。
        const int n_bins = outputHeight_;
        bins_ema_.clear();
        bins_bypass_.clear();
        int maxN = 8;

        const float logStep = (logMax_ - logMin_) / static_cast<float>(outputHeight_);
        std::vector<float> centers(n_bins);
        for (int y = 0; y < n_bins; ++y)
            centers[y] = std::pow(10.0f, logMax_ - (y + 0.5f) * logStep);

        for (int y = 0; y < n_bins; ++y) {
            Bin b{};
            b.f_center = centers[y];
            b.row = y;

            // 带宽 = 相邻两 bin 中心频率之差 × 缩放系数 (论文: W_NC = f(i+1) - f(i-1))
            // 频率随 y 递减，边界用单侧邻居；取绝对值保证正带宽。
            float f_hi = centers[y];                       // 更高频侧
            float f_lo = centers[y];                       // 更低频侧
            if (y > 0) f_hi = centers[y - 1];              // y-1 在更上方 = 更高频
            if (y + 1 < n_bins) f_lo = centers[y + 1];     // y+1 在更下方 = 更低频
            float w_nc = std::abs(f_hi - f_lo) * bandwidthScale;
            w_nc = std::max(w_nc, 1e-3f);

            // (7) 窗长；clamp 到 [8, maxWindowSamples]
            float q = std::round(2.0f * b.f_center / w_nc);
            float n_float = std::round(q * sampleRate / (2.0f * b.f_center));
            b.N = static_cast<int>(n_float);
            b.N = std::max(8, std::min(b.N, maxWindowSamples));

            // (5) 左右分量频率(用 double 保证旋转因子精度)
            b.f_left = static_cast<double>(b.f_center) - static_cast<double>(sampleRate) / (2.0 * b.N);
            b.f_right = static_cast<double>(b.f_center) + static_cast<double>(sampleRate) / (2.0 * b.N);

            // 旋转因子
            b.Wl = std::polar(1.0, -2.0 * std::numbers::pi_v<double> * b.f_left / sampleRate);
            b.Wr = std::polar(1.0, -2.0 * std::numbers::pi_v<double> * b.f_right / sampleRate);
            b.WNr_l = std::polar(1.0, -2.0 * std::numbers::pi_v<double> * b.f_left * b.N / sampleRate);
            b.WNr_r = std::polar(1.0, -2.0 * std::numbers::pi_v<double> * b.f_right * b.N / sampleRate);

            // 抗混叠 IIR(先带限再采样)：仅对窗长比图像采样间隔短的箱(N < hop)生效。
            // 滤波必须在滑动 DFT 逐样本计算内进行——每 hop 采样一次后再滤波，
            // 采样点已被混叠污染，无法挽回。故 EMA 移到逐样本层：
            //   N >= hop 箱自身窗长已是足够 FIR 带限(分辨率≈N/Fs)，旁路(alpha=1，不做EMA)；
            //   N <  hop 箱时间常数 tau = hop(=图像采样间隔)，每样本系数 = 1-e^(-1/hop)，
            //            在逐样本层低通，保证采样前已带限，消除混叠。
            // 按此分桶存储：bins_ema_(需每样本EMA) 与 bins_bypass_(旁路)，使样本循环
            // 无需 if 分支判断(分离位置在初始化即确定)。
            if (b.N < hopSize_) {
                // b.alpha = 1.0f - std::exp(-1.0f / static_cast<float>(hopSize_));  // 样本级, tau=hop
                // b.alpha = 1.0f - qwqdsp_misc::ExpSmoother::ComputeSmoothFactor2(hopSize_, 4);  // 样本级, tau=hop
                b.alpha = 1.0f;
                b.smoothed_gain = 0.0f;
                bins_ema_.push_back(b);
            } else {
                b.alpha = 1.0f;                                                    // 旁路
                b.smoothed_gain = 0.0f;
                bins_bypass_.push_back(b);
            }

            maxN = std::max(maxN, b.N);
        }

        // ── 共享样本环回缓冲（长度 = 最长窗长）──
        ring_.assign(maxN, 0.0);
        ringSize_ = maxN;
        n_ = 0;

        binDb_.resize(n_bins);
        column_.resize(outputHeight_);
    }

    /**
     * @brief 处理一帧，用递归滑动 DFT 演进所有 bin，输出一列
     *
     * window 与 windowed_frame 传入但被忽略（无窗法）。
     * SpectrogramColumn 提供重叠帧：第一次调用整帧都是新样本；之后每次调用只有
     * 末尾 hopSize_ 个样本是新增的（前面的 overlap 与上一帧重复），故只进给滑动 DFT。
     */
    void Process(std::span<const float> raw_frame, std::span<const float> /*window*/,
                 std::span<const float> /*windowed_frame*/) noexcept {
        // 只有 hopSize_ 个新样本会进入滑动 DFT
        const int new_sample_begin =
            (n_ == 0) ? 0 : std::max(0, static_cast<int>(raw_frame.size()) - hopSize_);
        for (int k = new_sample_begin; k < static_cast<int>(raw_frame.size()); ++k) {
            const double s = static_cast<double>(raw_frame[k]);  // 提升为 double 参与累加

            // ── 抗混叠箱(N < hop)：更新滑动 DFT 累加器 + 每样本 EMA(先带限再采样) ──
            for (Bin& b : bins_ema_) {
                double old = (n_ >= b.N) ? ring_[(n_ - b.N) % ringSize_] : 0.0;
                b.acc_l = b.Wl * b.acc_l + s - old * b.WNr_l;
                b.acc_r = b.Wr * b.acc_r + s - old * b.WNr_r;
                double ncSum = -(b.acc_l.real() * b.acc_r.real() + b.acc_l.imag() * b.acc_r.imag());
                float gain = (ncSum > 0.0) ? static_cast<float>(std::sqrt(ncSum + kEpsSample) / b.N) : 0.0f;
                gain *= gainNorm_;   // 幅度归一化(乘 π)，使纯音峰值回到 ≈0 dB
                b.smoothed_gain += b.alpha * (gain - b.smoothed_gain);   // 每样本 EMA
            }
            // ── 旁路箱(N >= hop)：窗长自身已带限，只更新滑动 DFT 累加器 ──
            for (Bin& b : bins_bypass_) {
                double old = (n_ >= b.N) ? ring_[(n_ - b.N) % ringSize_] : 0.0;
                b.acc_l = b.Wl * b.acc_l + s - old * b.WNr_l;
                b.acc_r = b.Wr * b.acc_r + s - old * b.WNr_r;
            }

            ring_[n_ % ringSize_] = s;
            ++n_;
        }

        // ── 按 row 回填各 bin 的 dB(此时已带限) ──
        constexpr float kEps = 1e-18f;
        // 抗混叠箱: 取每样本 EMA 平滑后的增益
        for (Bin& b : bins_ema_) {
            float gain = b.smoothed_gain;
            float db = (gain > 0.0f) ? 20.0f * std::log10(gain) : dbFloor_;
            binDb_[b.row] = std::clamp(db, dbFloor_, 0.0f);
        }
        // 旁路箱: 输出瞬时增益
        for (Bin& b : bins_bypass_) {
            double ncSum = -(b.acc_l.real() * b.acc_r.real() + b.acc_l.imag() * b.acc_r.imag());
            float gain = (ncSum > 0.0) ? static_cast<float>(std::sqrt(ncSum + kEps) / b.N) : 0.0f;
            gain *= gainNorm_;   // 幅度归一化(乘 π)，使纯音峰值回到 ≈0 dB
            float db = (gain > 0.0f) ? 20.0f * std::log10(gain) : dbFloor_;
            binDb_[b.row] = std::clamp(db, dbFloor_, 0.0f);
        }

        // ── 每行对应一个 bin（一一映射），直接上色 ──
        for (int y = 0; y < outputHeight_; ++y) {
            float db = binDb_[y];
            int idx = static_cast<int>((db - dbFloor_) / (-dbFloor_) * 255.0f);
            idx = std::clamp(idx, 0, 255);
            column_[y] = Colormap::kTable[idx];
        }
    }

    std::span<const Color> GetColumn() const noexcept {
        return {column_.data(), static_cast<size_t>(outputHeight_)};
    }

    int ColumnHeight() const noexcept {
        return outputHeight_;
    }

private:
    static constexpr float kMaxWindowS = 0.075f;   // 低音窗长上限(秒)
    static constexpr double kEpsSample = 1e-18;    // 样本级 NC 增益计算保护(平方根内非负)

    int sampleRate_{}, fftSize_{}, hopSize_{}, zeroPad_{}, outputHeight_{};
    int ringSize_{};
    int n_{};                                          // 已处理的样本计数
    float freqMin_{}, freqMax_{}, logMin_{}, logMax_{}, dbFloor_{};
    float gainNorm_{};                                 // 幅度归一化系数(乘 π，使纯音峰值 ≈0 dB)

    std::vector<Bin> bins_ema_;                       // N<hop 箱：每样本 EMA 抗混叠
    std::vector<Bin> bins_bypass_;                    // N>=hop 箱：旁路(窗长自身带限)
    std::vector<double> ring_;                       // 共享样本环回缓冲(double)
    std::vector<float> binDb_;                         // 每 bin dB
    std::vector<Color> column_;
};
