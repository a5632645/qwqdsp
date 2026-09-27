#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <vector>

#include <qwqdsp/spectral/real_fft_adv.hpp>

#include "raylib.h"

/**
 * @brief 重分配显示的公共累加网格: log 频率行 → 线性子格 + 环形子列
 *
 * 频率轴: 每个 log 行按线性频率细分为 k 个子格
 *   binHz  = sampleRate / fftLen
 *   k[y]   = max(1, lround(行带宽Hz[y] / (binHz · kSubcellScale)))
 * 低频行带宽不足一个子格时 k=1(直接叠加), 高频行 k>1 —— 即"低频在对数网格上
 * 直接叠加, 高频在更细的线性网格上叠加"。bin 落到 (行, 子格) 后:
 *   **子格内求和、行内取 max**。同一主瓣的 bin 重分配后落在同一子格(能量集中),
 *   而宽行内不同频率的 bin 分处不同子格, 不会叠在一起抬高噪声。全程不加权。
 *
 * 时间轴: 群延迟 → 环形子列索引, 环形指针每帧前进一格(免去整块搬移)。
 * subColumns == 1 时退化为纯频率重分配(无时间轴)。
 *
 * 幅度: 标定系数 cal 由 SetWindow 用实际窗函数自动算出 ——
 *   单位幅度、bin 中心的实正弦, 其重分配全谱和 / 谱峰 = 窗相干增益,
 * 使单位纯音在显示上回到 0 dB。
 *
 * @tparam EnableFreqInterp 是否启用频率轴方向的线性分布 (相邻两行加权叠加):
 *                          true  → 纵向平滑、无栅格锯齿, 但落在两行之间的音会衰减
 *                          false → 只落到单一行, 无衰减, 可能有栅格锯齿
 *                          两种情况都只叠加、不做任何归一化
 */
template <bool EnableFreqInterp = true>
struct LogReassignGrid {
    /// @brief 子格宽度相对 FFT bin 宽的比值 (0.25 = bin/4)
    static constexpr float kSubcellScale = 0.25f;

    void Init(int sampleRate, int fftLen, int subColumns, int outputHeight, float freqMin, float freqMax,
              float dbFloor) noexcept {
        sampleRate_ = sampleRate;
        fftLen_ = fftLen;
        subColumns_ = std::max(subColumns, 1);
        outputHeight_ = outputHeight;
        logMin_ = std::log10(freqMin);
        logMax_ = std::log10(freqMax);
        dbFloor_ = dbFloor;
        binSize_ = fftLen / 2 + 1;
        head_ = 0;

        // ── log 频率行 → 频率下界 / 子格宽 / 子格数 ──
        const float bin_hz = static_cast<float>(sampleRate_) / static_cast<float>(fftLen_);
        const float subcell_hz = bin_hz * kSubcellScale;
        const float log_span = logMax_ - logMin_;
        const float rows = static_cast<float>(std::max(outputHeight_ - 1, 1));

        row_f_lo_.resize(outputHeight_);
        row_cell_w_.resize(outputHeight_);
        row_k_.resize(outputHeight_);
        for (int y = 0; y < outputHeight_; ++y) {
            // y_pos = (H-1)·(1-norm), floor(y_pos)=y 对应的 norm 区间
            float norm_lo = static_cast<float>(outputHeight_ - 2 - y) / rows;
            float norm_hi = static_cast<float>(outputHeight_ - 1 - y) / rows;
            float f_lo = std::pow(10.0f, logMin_ + norm_lo * log_span);
            float f_hi = std::pow(10.0f, logMin_ + norm_hi * log_span);
            float row_w = f_hi - f_lo;
            int k = std::max(static_cast<int>(std::lround(row_w / subcell_hz)), 1);
            row_f_lo_[y] = f_lo;
            row_k_[y] = k;
            row_cell_w_[y] = row_w / static_cast<float>(k);
        }
        // 紧凑行偏移 (低频行 k 很小, 避免 H·k_max 的稀疏浪费)
        row_offset_.resize(outputHeight_);
        int total_cells = 0;
        for (int y = 0; y < outputHeight_; ++y) {
            row_offset_[y] = total_cells;
            total_cells += row_k_[y];
        }
        total_cells_ = total_cells;
        col_buf_.assign(static_cast<size_t>(subColumns_) * static_cast<size_t>(total_cells_), 0.0f);

        win_fft_.resize(binSize_);
        cal_in_.resize(fftLen_, 0.0f);
        cal_ = 1.0f;
        cal_ready_ = false;
    }

    /// @brief 把一个 bin 加到网格; inst_freq_hz 需已落在 [freqMin, freqMax]
    /// @param group_delay 群延迟 ∈ [-0.5, 0.5] (无时间重分配时传 0)
    void Add(float inst_freq_hz, float group_delay, float mag) noexcept {
        // ── 频率轴: 按小数位置双线性分配到相邻两行 (EnableFreqInterp), 纵向平滑 ──
        float logF = std::log10(inst_freq_hz);
        float norm = (logF - logMin_) / (logMax_ - logMin_);
        float y_pos = std::clamp(static_cast<float>(outputHeight_ - 1) * (1.0f - norm), 0.0f,
                                 static_cast<float>(outputHeight_ - 1));
        int y0 = static_cast<int>(std::floor(y_pos));
        float y_frac = y_pos - static_cast<float>(y0);
        if (y0 >= outputHeight_ - 1) {
            y0 = outputHeight_ - 1;
            y_frac = 0.0f;
        }
        if constexpr (!EnableFreqInterp) {
            y_frac = 0.0f; // 关闭频率轴线性分布: 只落到 y0
        }
        int const y1 = (y_frac > 0.0f) ? y0 + 1 : y0;

        // ── 时间轴: 群延迟 → 环形子列 + 双线性 ──
        float c_pos = std::clamp((group_delay + 0.5f) * static_cast<float>(subColumns_), 0.0f,
                                 static_cast<float>(subColumns_ - 1));
        int c_idx = static_cast<int>(std::floor(c_pos));
        float c_frac = c_pos - static_cast<float>(c_idx);
        if (c_idx >= subColumns_ - 1) {
            c_idx = subColumns_ - 1;
            c_frac = 0.0f;
        }
        int slot0 = head_ + c_idx;
        if (slot0 >= subColumns_)
            slot0 -= subColumns_;
        int slot1 = slot0 + 1;
        if (slot1 >= subColumns_)
            slot1 -= subColumns_;

        // 行内硬分配子格; 相邻行沿用同一子格序号 (若按"落在该行内"重算, 该频率根本
        // 不在相邻行的频率范围内, 会把所有越界 bin 都夹到边缘子格、叠出假亮边)
        int subcell = static_cast<int>((inst_freq_hz - row_f_lo_[y0]) / row_cell_w_[y0]);
        subcell = std::clamp(subcell, 0, row_k_[y0] - 1);

        auto splat = [&](int row, int cell, float wy) {
            float const v = mag * wy;
            int const base = row_offset_[row] + cell;
            col_buf_[slot0 * total_cells_ + base] += v * (1.0f - c_frac);
            if (c_frac > 0.0f)
                col_buf_[slot1 * total_cells_ + base] += v * c_frac;
        };
        splat(y0, subcell, 1.0f - y_frac);
        if (y1 != y0)
            splat(y1, std::min(subcell, row_k_[y1] - 1), y_frac);
    }

    /// @brief 归约最旧子列(行内取 max) → 标定 → dB → 上色, 然后子列指针前进
    void Emit(std::span<const Color, 256> table, std::span<Color> out) noexcept {
        constexpr float kEps = 1e-12f;
        int const oldest = head_ * total_cells_;
        for (int y = 0; y < outputHeight_; ++y) {
            float v = 0.0f;
            int const base = oldest + row_offset_[y];
            int const nk = row_k_[y];
            for (int j = 0; j < nk; ++j) {
                v = std::max(v, col_buf_[base + j]);
            }
            float dB = 20.0f * std::log10(v / cal_ + kEps);
            dB = std::clamp(dB, dbFloor_, 0.0f);
            int idx = static_cast<int>((dB - dbFloor_) / (-dbFloor_) * 255.0f);
            idx = std::clamp(idx, 0, 255);
            out[y] = table[idx];
        }
        std::fill(col_buf_.begin() + oldest, col_buf_.begin() + oldest + total_cells_, 0.0f);
        head_ = (head_ + 1 == subColumns_) ? 0 : head_ + 1;
    }

    /**
     * @brief 用实际窗函数自动算幅度标定
     *
     * 合成单位幅度、bin 中心的实正弦 (取 fftLen/8, 远离 DC/Nyquist 的负频镜像),
     * 该比值与 bin 位置无关。
     */
    void SetWindow(std::span<const float> window, qwqdsp_spectral::RealFftAdv& fft) noexcept {
        constexpr float two_pi = 2.0f * std::numbers::pi_v<float>;
        int const fft_size = static_cast<int>(window.size());
        int const cal_bin = fftLen_ / 8;
        float const w0 = two_pi * static_cast<float>(cal_bin) / static_cast<float>(fftLen_);

        for (int i = 0; i < fft_size; ++i) {
            cal_in_[i] = window[i] * std::cos(w0 * static_cast<float>(i));
        }
        std::fill(cal_in_.begin() + fft_size, cal_in_.end(), 0.0f);
        fft.FFT(cal_in_, win_fft_);

        float sum = 0.0f;
        float peak = 0.0f;
        for (auto const& c : win_fft_) {
            float m = std::abs(c);
            sum += m;
            peak = std::max(peak, m);
        }
        cal_ = sum / std::max(peak, 1e-20f);
        cal_ready_ = true;
    }

    /// @brief 窗相干增益 (单位纯音的重分配和 / 谱峰), SetWindow 后有效
    float Calibration() const noexcept {
        return cal_;
    }

    bool CalibrationReady() const noexcept {
        return cal_ready_;
    }

private:
    int sampleRate_{}, fftLen_{}, subColumns_{}, outputHeight_{}, binSize_{}, total_cells_{}, head_{};
    float logMin_{}, logMax_{}, dbFloor_{}, cal_{1.0f};
    bool cal_ready_{};

    std::vector<float> row_f_lo_;    // [outputHeight] 该行频率下界 (Hz)
    std::vector<float> row_cell_w_;  // [outputHeight] 该行子格宽 (Hz)
    std::vector<int> row_k_;         // [outputHeight] 该行子格数
    std::vector<int> row_offset_;    // [outputHeight] 该行子格在列块内的起始偏移
    std::vector<float> col_buf_;     // [subColumns · Σk]
    std::vector<float> cal_in_;      // [fftLen] 标定用合成信号
    std::vector<std::complex<float>> win_fft_;
};
