// ------------------------------------------------------------
//  dft_filter_bank.cpp
//  DFT 调制下采样/上采样滤波器组 — 信号分解与还原
//
//  信号链:
//    输入 noise.wav
//    → DftFilterBank::Analysis   (DFT 调制 + 下采样 → M 个子带)
//    → DftFilterBank::Synthesis  (DFT 调制 + 上采样 → 重叠相加还原)
//    → 输出 dft_fb_*.wav
//
//  对比两组参数:
//    · 临界采样 (D = M)      — 子带混叠 (aliasing) 明显 (≈14 dB)
//    · 过采样   (D = M/2)    — 近完美重建 (≈44 dB, Kaiser β=16)
//    · 过采样 + 虚部取反     — 子带取共轭 → 频谱镜像 → 帧内时间反转
//
//  设计说明: 折叠式 WOLA 的重建混叠深度 ≈ 窗函数在滞后 M 处的自相关
//  R(M)/R(0)。L=2M 时 sinc 零点恰好落在窗边界, R(M)/R(0)≈0.04;
//  Kaiser β=16 的低旁瓣进一步压到 ≈0.005, 使过采样重建达到 -44 dB。
// ------------------------------------------------------------

#include <algorithm>
#include <cstddef>
#include <format>
#include <iostream>
#include <numbers>
#include <span>
#include <string>
#include <vector>

#include "AudioFile.h"
#include <qwqdsp/filter/window_fir.hpp>
#include <qwqdsp/spectral/complex_fft.hpp>
#include <qwqdsp/window/kaiser.hpp>
#include <work_dir.hpp>

// ------------------------------------------------------------
//  滤波器组参数
// ------------------------------------------------------------
static constexpr int kNumSubbands = 256;              // M — 子带数 (FFT 长度)
static constexpr int kOverlapFactor = 2;              // K — 重叠因子 (2 使 sinc 零点落在窗边界)
static constexpr int kPrototypeLen = kNumSubbands * kOverlapFactor;  // L = K·M
static constexpr float kKaiserBeta = 16.0f;           // Kaiser 窗 β — 低旁瓣抑制重建混叠

// ════════════════════════════════════════════════════════════
//  DFT 调制滤波器组 (分析 + 合成)
// ════════════════════════════════════════════════════════════

/**
 * @brief DFT 调制滤波器组 — 统一的 WOLA 实现
 *
 * 分析端 (下采样): 加窗 → 时间混叠折叠到 M 点 → M 点 DFT → M 个子带
 * 合成端 (上采样): M 点 IDFT → 展开加窗 → 重叠相加 → COLA 逐样本归一化
 *
 * 抽取/插值因子 R (hop) 由 Init 传入:
 *   · R = M     临界采样 — 子带混叠 (aliasing) 明显
 *   · R = M/2   2x 过采样 — 近完美重建 (NPR)
 *
 * 注意: 分析窗与合成窗相同时, 合成端做逐样本 COLA 归一化
 *       (除以 Σ_m w_a(n-mR)·w_s(n-mR)) 后, 主项 (q=0) 精确还原,
 *       残留误差即子带混叠, 用于对比两种采样率下的重建质量。
 */
struct DftFilterBank {
    /**
     * @brief 初始化滤波器组
     * @param num_subbands 子带数 M (FFT 长度, 需为 2 的幂)
     * @param hop          抽取/插值因子 R
     * @param prototype    原型低通滤波器 (长度 L = K·M, 截止 π/M)
     */
    void Init(int num_subbands, int hop, std::span<const float> prototype) {
        m_ = num_subbands;
        hop_ = hop;
        proto_len_ = static_cast<int>(prototype.size());
        overlap_ = proto_len_ / m_;

        // 分析窗与合成窗 — 演示中使用同一原型窗
        win_a_.assign(prototype.begin(), prototype.end());
        win_s_.assign(prototype.begin(), prototype.end());

        // COLA 归一化系数: norm[r] = Σ_{j≡r (mod R)} w_a(j)·w_s(j)
        norm_.assign(hop_, 0.0f);
        for (int n = 0; n < proto_len_; ++n) {
            norm_[n % hop_] += win_a_[n] * win_s_[n];
        }

        // FFT 工作缓冲
        fft_.Init(m_);
        fold_r_.assign(m_, 0.0f);
        fold_i_.assign(m_, 0.0f);
        ifft_r_.assign(m_, 0.0f);
        ifft_i_.assign(m_, 0.0f);
    }

    /** @brief 帧数 = ceil(N / R) */
    int NumFrames(int num_samples) const noexcept {
        return (num_samples + hop_ - 1) / hop_;
    }

    /**
     * @brief 分析: 信号 → 子带
     * @param x      输入时域信号 (N 样本)
     * @param sub_r  子带实部, 大小 num_frames × M
     * @param sub_i  子带虚部, 大小 num_frames × M
     */
    void Analysis(std::span<const float> x, std::span<float> sub_r, std::span<float> sub_i) noexcept {
        const int num_frames = NumFrames(static_cast<int>(x.size()));
        for (int m = 0; m < num_frames; ++m) {
            // Step 1 — 加窗 + 时间混叠折叠 (fold) 到 M 点
            std::fill(fold_r_.begin(), fold_r_.end(), 0.0f);
            for (int j = 0; j < overlap_; ++j) {
                for (int k = 0; k < m_; ++k) {
                    const int n = m * hop_ + k + j * m_;
                    if (n < static_cast<int>(x.size())) {
                        fold_r_[k] += x[n] * win_a_[k + j * m_];
                    }
                }
            }

            // Step 2 — M 点 DFT → 子带
            float* out_r = sub_r.data() + static_cast<size_t>(m) * m_;
            float* out_i = sub_i.data() + static_cast<size_t>(m) * m_;
            fft_.FFT(fold_r_.data(), fold_i_.data(), out_r, out_i);
        }
    }

    /**
     * @brief 合成: 子带 → 信号 (含 COLA 归一化)
     * @param sub_r 子带实部
     * @param sub_i 子带虚部
     * @param y     输出时域信号, 长度 = 输入长度
     */
    void Synthesis(std::span<const float> sub_r, std::span<const float> sub_i, std::span<float> y) noexcept {
        std::fill(y.begin(), y.end(), 0.0f);
        const int num_frames = NumFrames(static_cast<int>(y.size()));
        const int n_total = static_cast<int>(y.size());

        for (int m = 0; m < num_frames; ++m) {
            // Step 1 — M 点 IDFT
            const float* in_r = sub_r.data() + static_cast<size_t>(m) * m_;
            const float* in_i = sub_i.data() + static_cast<size_t>(m) * m_;
            fft_.IFFT(in_r, in_i, ifft_r_.data(), ifft_i_.data());

            // Step 2 — 展开 + 加窗 + 重叠相加
            for (int j = 0; j < overlap_; ++j) {
                for (int k = 0; k < m_; ++k) {
                    const int n = m * hop_ + k + j * m_;
                    if (n < n_total) {
                        y[n] += ifft_r_[k] * win_s_[k + j * m_];
                    }
                }
            }
        }

        // Step 3 — COLA 逐样本归一化 (主项精确还原, 只残留混叠)
        for (int n = 0; n < n_total; ++n) {
            y[n] /= norm_[n % hop_];
        }
    }

private:
    int m_{};          // 子带数 M
    int hop_{};        // 抽取/插值因子 R
    int proto_len_{};  // 原型长度 L
    int overlap_{};    // 重叠因子 K = L / M

    std::vector<float> win_a_;  // 分析窗 (L)
    std::vector<float> win_s_;  // 合成窗 (L)
    std::vector<float> norm_;   // COLA 归一化系数 (R)

    qwqdsp_spectral::ComplexFFT fft_;
    std::vector<float> fold_r_;  // 折叠缓冲 (M)
    std::vector<float> fold_i_;  // 折叠缓冲虚部 (恒为 0)
    std::vector<float> ifft_r_;  // IDFT 缓冲 (M)
    std::vector<float> ifft_i_;  // IDFT 缓冲虚部
};

// ------------------------------------------------------------
//  工具函数
// ------------------------------------------------------------

/** @brief 设计原型低通: 窗函数法, 截止 π/M, Kaiser 窗 */
static std::vector<float> DesignPrototype() {
    std::vector<float> proto(kPrototypeLen);
    qwqdsp_filter::WindowFIR::Lowpass(proto, std::numbers::pi_v<float> / kNumSubbands);
    qwqdsp_window::Kaiser::ApplyWindow(proto, kKaiserBeta, false);  // 对称窗 (FIR 设计)
    qwqdsp_filter::WindowFIR::Normalize(proto);                     // DC 增益 = 1
    return proto;
}

/** @brief 导出单声道 32-bit float wav */
static void SaveWav(std::span<const float> data, int fs, const std::string& name) {
    AudioFile<float> f;
    f.setNumChannels(1);
    f.setBitDepth(32);
    f.setSampleRate(fs);
    f.setNumSamplesPerChannel(data.size());
    std::copy(data.begin(), data.end(), f.samples[0].begin());
    f.save(qwqdsp_support::OutputFile(name));
}

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------
int main() {
    // 1. 加载输入噪声
    AudioFile<float> file;
    if (!file.load(qwqdsp_support::InputFile("noise.wav"))) {
        std::cout << "failed to load noise.wav\n";
        return 1;
    }
    const int fs = static_cast<int>(file.getSampleRate());
    std::vector<float> x(file.samples[0].begin(), file.samples[0].end());  // 取第一声道
    const int n_samp = static_cast<int>(x.size());
    std::cout << std::format("load noise.wav: fs={} n={}\n", fs, n_samp);

    // 2. 设计原型低通 (M=256, L=2M=512, Kaiser β=16)
    const std::vector<float> proto = DesignPrototype();

    // 3. 临界采样 (R = M) — 子带混叠明显
    {
        DftFilterBank fb;
        fb.Init(kNumSubbands, kNumSubbands, proto);

        std::vector<float> sub_r(static_cast<size_t>(fb.NumFrames(n_samp)) * kNumSubbands);
        std::vector<float> sub_i(sub_r.size());
        fb.Analysis(x, sub_r, sub_i);

        std::vector<float> y(n_samp);
        fb.Synthesis(sub_r, sub_i, y);
        SaveWav(y, fs, "dft_fb_critical_out.wav");
        std::cout << "critical (D=M): saved dft_fb_critical_out.wav\n";
    }

    // 4. 过采样 (R = M/2) — 近完美重建
    {
        DftFilterBank fb;
        fb.Init(kNumSubbands, kNumSubbands / 2, proto);

        std::vector<float> sub_r(static_cast<size_t>(fb.NumFrames(n_samp)) * kNumSubbands);
        std::vector<float> sub_i(sub_r.size());
        fb.Analysis(x, sub_r, sub_i);

        std::vector<float> y(n_samp);
        fb.Synthesis(sub_r, sub_i, y);
        SaveWav(y, fs, "dft_fb_oversampled_out.wav");
        std::cout << "oversampled (D=M/2): saved dft_fb_oversampled_out.wav\n";
    }

    // 5. 过采样 + 虚部取反 (R = M/2) — 子带相位翻转实验
    {
        DftFilterBank fb;
        fb.Init(kNumSubbands, kNumSubbands / 2, proto);

        std::vector<float> sub_r(static_cast<size_t>(fb.NumFrames(n_samp)) * kNumSubbands);
        std::vector<float> sub_i(sub_r.size());
        fb.Analysis(x, sub_r, sub_i);

        // 虚部全部取反 (1i → -1i): 等价于子带取共轭。对实输入, 共轭对称
        // 使得 sub_k → sub_{M-k}, 即频谱镜像 → 合成端帧内时间反转
        for (float& v : sub_i) {
            v = -v;
        }

        std::vector<float> y(n_samp);
        fb.Synthesis(sub_r, sub_i, y);
        SaveWav(y, fs, "dft_fb_oversampled_imneg_out.wav");
        std::cout << "oversampled+imneg (D=M/2): saved dft_fb_oversampled_imneg_out.wav\n";
    }

    std::cout << "done\n";
    return 0;
}
