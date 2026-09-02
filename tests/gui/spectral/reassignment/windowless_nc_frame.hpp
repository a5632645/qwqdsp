#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <vector>

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
        float f_center{};        // 中心频率(Hz)
        float f_left{};          // 左分量参考频率
        float f_right{};         // 右分量参考频率
        int N{};                 // 窗长(样本数)
        std::complex<float> Wl{}, WNr_l{};   // 左分量旋转因子 W 与 W^N
        std::complex<float> Wr{}, WNr_r{};   // 右分量旋转因子 W 与 W^N
        std::complex<float> acc_l{}, acc_r{}; // 左右分量滑动 DFT 累加器
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

        // 低音窗长上限(论文 IV-A 节)
        const int maxWindowSamples = std::max(8, static_cast<int>(kMaxWindowS * sampleRate));

        // ── 以频率网格为基准：y 轴每个像素对应一个 NC bin，bin 中心频率 = 该像素频率 ──
        // y=0 在顶部(高频)，y 增大频率递减；bin 数与输出高度一致，一一映射。
        const int n_bins = outputHeight_;
        bins_.clear();
        bins_.resize(n_bins);
        int maxN = 8;

        const float logStep = (logMax_ - logMin_) / static_cast<float>(outputHeight_);
        std::vector<float> centers(n_bins);
        for (int y = 0; y < n_bins; ++y)
            centers[y] = std::pow(10.0f, logMax_ - (y + 0.5f) * logStep);

        for (int y = 0; y < n_bins; ++y) {
            Bin b{};
            b.f_center = centers[y];

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

            // (5) 左右分量频率
            b.f_left = b.f_center - sampleRate / (2.0f * b.N);
            b.f_right = b.f_center + sampleRate / (2.0f * b.N);

            // 旋转因子
            b.Wl = std::polar(1.0f, -2.0f * std::numbers::pi_v<float> * b.f_left / sampleRate);
            b.Wr = std::polar(1.0f, -2.0f * std::numbers::pi_v<float> * b.f_right / sampleRate);
            b.WNr_l = std::polar(1.0f, -2.0f * std::numbers::pi_v<float> * b.f_left * b.N / sampleRate);
            b.WNr_r = std::polar(1.0f, -2.0f * std::numbers::pi_v<float> * b.f_right * b.N / sampleRate);

            maxN = std::max(maxN, b.N);
            bins_[y] = b;
        }

        // ── 共享样本环回缓冲（长度 = 最长窗长）──
        ring_.assign(maxN, 0.0f);
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
            const float s = raw_frame[k];
            for (Bin& b : bins_) {
                // x[n - N]：从环回缓冲取(样本尚未入窗时视为 0)
                float old = (n_ >= b.N) ? ring_[(n_ - b.N) % ringSize_] : 0.0f;
                b.acc_l = b.Wl * b.acc_l + s - old * b.WNr_l;
                b.acc_r = b.Wr * b.acc_r + s - old * b.WNr_r;
            }
            ring_[n_ % ringSize_] = s;
            ++n_;
        }

        // ── 输出当前时刻各 bin 的 NC dB 值 ──
        constexpr float kEps = 1e-18f;
        const int n_bins = static_cast<int>(bins_.size());
        for (int i = 0; i < n_bins; ++i) {
            const Bin& b = bins_[i];
            float ncSum = -(b.acc_l.real() * b.acc_r.real() + b.acc_l.imag() * b.acc_r.imag());
            if (ncSum < 0.0f) {
                binDb_[i] = dbFloor_;
            } else {
                float nc = std::sqrt(ncSum + kEps) / static_cast<float>(b.N);  // 归一化 ÷N
                binDb_[i] = 20.0f * std::log10(nc);
            }
            binDb_[i] = std::clamp(binDb_[i], dbFloor_, 0.0f);
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
    static constexpr float kMaxWindowS = 0.125f;   // 低音窗长上限(秒)

    int sampleRate_{}, fftSize_{}, hopSize_{}, zeroPad_{}, outputHeight_{};
    int ringSize_{};
    int n_{};                                          // 已处理的样本计数
    float freqMin_{}, freqMax_{}, logMin_{}, logMax_{}, dbFloor_{};

    std::vector<Bin> bins_;
    std::vector<float> ring_;                          // 共享样本环回缓冲
    std::vector<float> binDb_;                         // 每 bin dB
    std::vector<Color> column_;
};
