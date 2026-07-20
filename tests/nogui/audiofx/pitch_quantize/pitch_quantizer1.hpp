#pragma once

#include "scale_helper.hpp"
#include <qwqdsp/window/hann.hpp>
#include <qwqdsp/spectral/real_fft.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <format>
#include <iostream>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

// ------------------------------------------------------------
// PitchQuantizer — 基于朴素 Phase Vocoder 的音高量化器
// ------------------------------------------------------------
/**
 * @brief 基于朴素相位声码器的离线音高量化器。
 *
 * @tparam kClampWhenExceed  超出最大修正范围时的行为
 *   - true （默认）: 修正量 clamp 到 ±fs/(2·hop)
 *   - false         : 不动这个 bin（保持原频率）
 *
 * 原理：
 *   对 STFT 的每个 bin，通过帧间相位差估计瞬时频率，
 *   将瞬时频率拉到最近的调性内 MIDI 音符频率上，
 *   再通过修改综合相位实现频率修正。
 *
 * 约束：
 *   - 朴素 PV 的频率修正上限为 ±fs/(2·hop) Hz
 *   - 分析/综合窗均为 Hann 窗，步长 hop
 */
template <bool kClampWhenExceed = true>
class PitchQuantizer {
public:
    /**
     * @brief 初始化
     * @param sample_rate  采样率 (Hz)
     * @param hop_size     帧移 (样本数)
     * @param fft_size     FFT 大小（须为 2 的幂）
     */
    void Init(float sample_rate, size_t hop_size, size_t fft_size) {
        sample_rate_ = sample_rate;
        hop_size_ = hop_size;
        fft_size_ = fft_size;

        // Hann 窗（分析/综合共用）
        window_.resize(fft_size);
        qwqdsp_window::Hann::Window(window_, true);

        // FFT
        fft_.Init(fft_size);

        // 状态缓冲
        const size_t num_bins = fft_size / 2 + 1;
        prev_phase_.assign(num_bins, 0.0f);
        synth_phase_.assign(num_bins, 0.0f);
        fft_in_.resize(fft_size);
        fft_out_.resize(fft_size + 2);

        first_frame_ = true;

        // 默认调性：C Major
        allowed_mask_ = ScaleHelper::MakeMask(0, ScaleHelper::Type::kMajor);
        root_note_ = 0;
        scale_type_ = ScaleHelper::Type::kMajor;
    }

    /**
     * @brief 以最大修正频率初始化（自动计算满足 Hann COLA 的 hop）
     * @param sample_rate          采样率 (Hz)
     * @param max_correction_freq  期望的最大频率修正范围 (Hz)
     * @param fft_size             FFT 大小（须为 2 的幂）
     *
     * @details
     *   hop 由 max_correction_freq 计算：hop = fs / (2·max_correction_freq)
     *   自动取整到最近的满足 Hann COLA 的 hop（即 fft_size % hop == 0）。
     *   实际 max_correction 可能因取整略小于传入值。
     */
    void Init2(float sample_rate, float max_correction_freq, size_t fft_size) {
        // hop = fs / (2 * max_correction)
        const float hop_desired = sample_rate / (2.0f * max_correction_freq);
        if (hop_desired < 1.0f) {
            std::cout << "  warning: max_correction_freq too large, using hop=1\n" << std::flush;
            Init(sample_rate, 1, fft_size);
            return;
        }

        const size_t hop = FindColaHop(fft_size, static_cast<size_t>(std::round(hop_desired)));
        const float actual_max_correction = sample_rate / (2.0f * static_cast<float>(hop));

        std::cout << std::format("  Init2: desired_hop={:.1f}, cola_hop={}, "
                                 "desired_max_corr={:.1f} Hz, actual_max_corr={:.1f} Hz\n",
                                 hop_desired, hop, max_correction_freq, actual_max_correction)
                  << std::flush;

        Init(sample_rate, hop, fft_size);
    }

    /**
     * @brief 设置调性
     * @param root 根音音级 (0=C, 1=C#, …, 11=B)
     * @param type 音阶类型
     */
    void SetKey(int root, ScaleHelper::Type type) noexcept {
        root_note_ = root % 12;
        scale_type_ = type;
        allowed_mask_ = ScaleHelper::MakeMask(root_note_, type);
    }

    /**
     * @brief 获取当前调性描述
     */
    std::string KeyDescription() const {
        return std::format("{} {}", ScaleHelper::kNoteNames[root_note_],
                           ScaleHelper::TypeName(scale_type_));
    }

    /**
     * @brief 离线处理
     * @param input  输入信号
     * @return 输出信号
     */
    std::vector<float> Process(std::span<const float> input) {
        const size_t N = fft_size_;
        const size_t H = hop_size_;
        const float fs = sample_rate_;
        const size_t num_bins = N / 2 + 1;

        // 帧数
        if (input.size() < N)
            return {};
        const size_t num_frames = (input.size() - N) / H + 1;

        // 输出缓冲（OLA）
        const size_t output_len = H * num_frames + N;
        std::vector<float> output(output_len, 0.0f);

        const float two_pi = 2.0f * std::numbers::pi_v<float>;
        // 最大可修正频率偏移 (Hz)
        const float max_correction = fs / (2.0f * static_cast<float>(H));

        std::cout << std::format("  frames={}, fft={}, hop={}, max_correction={:.1f} Hz\n",
                                 num_frames, N, H, max_correction)
                  << std::flush;

        // 逐帧处理
        for (size_t i = 0; i < num_frames; ++i) {
            // ----------------------------------------------------
            // 分析
            // ----------------------------------------------------
            const float* frame = input.data() + i * H;
            for (size_t j = 0; j < N; ++j) {
                fft_in_[j] = frame[j] * window_[j];
            }
            fft_.FFT(fft_in_.data(), fft_out_.data());

            // 提取幅度 & 相位
            std::vector<float> mag(num_bins);
            std::vector<float> phase(num_bins);
            for (size_t k = 0; k < num_bins; ++k) {
                const float re = fft_out_[2 * k];
                const float im = fft_out_[2 * k + 1];
                mag[k] = std::sqrt(re * re + im * im);
                phase[k] = std::atan2(im, re);
            }

            // ----------------------------------------------------
            // PV 处理：瞬时频率估计 → 音高量化 → 综合相位递推
            // ----------------------------------------------------
            if (!first_frame_) {
                const float inv_N = 1.0f / static_cast<float>(N);
                const float inv_H = 1.0f / static_cast<float>(H);

                for (size_t k = 0; k < num_bins; ++k) {
                    // 跳过 DC / Nyquist（相位信息不可靠）和低幅度 bin
                    if (k == 0 || k == N / 2 || mag[k] < 1e-6f) {
                        synth_phase_[k] += two_pi * static_cast<float>(k)
                                           * static_cast<float>(H) * inv_N;
                        continue;
                    }

                    // ---- 瞬时频率估计 ----
                    float delta = phase[k] - prev_phase_[k];
                    // 解卷绕到 (-π, π]
                    delta -= two_pi * std::round(delta / two_pi);
                    // 减去 bin 中心预期相位差
                    const float expected = two_pi * static_cast<float>(k)
                                           * static_cast<float>(H) * inv_N;
                    const float deviation = delta - expected;

                    // 瞬时频率 (Hz)
                    const float f_inst = static_cast<float>(k) * fs * inv_N
                                         + deviation * fs / (two_pi * static_cast<float>(H));

                    // ---- 音高量化 ----
                    const float f_midi = FindNearestAllowedMidi(f_inst);

                    float df = f_midi - f_inst;
                    if constexpr (kClampWhenExceed) {
                        // 修正量 clamp 到最大范围
                        df = std::clamp(df, -max_correction, max_correction);
                        const float f_target = f_inst + df;
                        const float synth_delta = two_pi * f_target * static_cast<float>(H) / fs;
                        synth_phase_[k] += synth_delta;
                    }
                    else {
                        // 超出最大修正范围时不动这个 bin
                        if (std::abs(df) <= max_correction) {
                            const float f_target = f_inst + df;
                            const float synth_delta = two_pi * f_target * static_cast<float>(H) / fs;
                            synth_phase_[k] += synth_delta;
                        }
                        else {
                            synth_phase_[k] += two_pi * static_cast<float>(k)
                                               * static_cast<float>(H) / static_cast<float>(N);
                        }
                    }
                }
            }
            else {
                // 第一帧：综合相位 = 分析相位
                for (size_t k = 0; k < num_bins; ++k) {
                    synth_phase_[k] = phase[k];
                }
                first_frame_ = false;
            }

            // 保存当前相位
            std::copy(phase.begin(), phase.end(), prev_phase_.begin());

            // ----------------------------------------------------
            // 综合
            // ----------------------------------------------------
            for (size_t k = 0; k < num_bins; ++k) {
                const float c = std::cos(synth_phase_[k]);
                const float s = std::sin(synth_phase_[k]);
                fft_out_[2 * k] = mag[k] * c;
                fft_out_[2 * k + 1] = mag[k] * s;
            }

            fft_.IFFT(fft_out_.data(), fft_in_.data());

            // 加综合窗 + OLA
            for (size_t j = 0; j < N; ++j) {
                output[i * H + j] += fft_in_[j] * window_[j];
            }
        }

        // ---- OLA 增益归一化 ----
        // 对于 Hann 窗 + hop = N/4，OLA 增益理论上均匀
        // 这里用数值方法计算 DC 增益并归一化
        {
            std::vector<float> gain_buf(output_len, 0.0f);
            for (size_t i = 0; i < num_frames; ++i) {
                for (size_t j = 0; j < N; ++j) {
                    gain_buf[i * H + j] += window_[j] * window_[j] / static_cast<float>(N);
                }
            }
            float peak_gain = 0.0f;
            for (float g : gain_buf) {
                peak_gain = std::max(peak_gain, g);
            }
            const float scale = (peak_gain > 1e-10f) ? (1.0f / peak_gain) : 1.0f;
            for (float& x : output) {
                x *= scale;
            }

            std::cout << std::format("  ola_peak_gain={:.4f}, scale={:.4f}\n", peak_gain, scale)
                      << std::flush;
        }

        return output;
    }

    void Reset() noexcept {
        first_frame_ = true;
        std::fill(prev_phase_.begin(), prev_phase_.end(), 0.0f);
        std::fill(synth_phase_.begin(), synth_phase_.end(), 0.0f);
    }

private:
    /**
     * @brief 查找满足 Hann COLA 的最近 hop
     * @details Hann COLA 要求 hop = N/R，R ≥ 4 为整数（即 hop ≤ N/4）。
     *          从 desired 向两侧搜索 N 的约数中满足 hop ≤ N/4 的值。
     */
    static size_t FindColaHop(size_t fft_size, size_t desired) noexcept {
        const size_t kMaxHop = fft_size / 4; // R ≥ 4

        if (desired < 1)
            desired = 1;
        if (desired > kMaxHop)
            desired = kMaxHop;

        // 双向搜索最近的约数（满足 ≤ N/4）
        for (size_t offset = 0; offset <= kMaxHop; ++offset) {
            // 向下搜索
            if (offset < desired) {
                const size_t cand = desired - offset;
                if (cand >= 1 && fft_size % cand == 0)
                    return cand;
            }
            // 向上搜索
            if (desired + offset <= kMaxHop) {
                const size_t cand = desired + offset;
                if (fft_size % cand == 0)
                    return cand;
            }
        }
        // fallback：N/4 本身一定是约数
        return kMaxHop;
    }

    /**
     * @brief 在允许的调性内查找离 freq_hz 最近的 MIDI 音符频率
     */
    float FindNearestAllowedMidi(float freq_hz) const noexcept {
        if (freq_hz <= 0.0f)
            return freq_hz;

        // 粗略最近 MIDI 音符
        const float m_float = 12.0f * std::log2(freq_hz / 440.0f) + 69.0f;
        const int m_round = static_cast<int>(std::round(m_float));

        // 双向扩展搜索最近的允许音级
        for (int offset = 0; offset < 128; ++offset) {
            for (int sign : {1, -1}) {
                const int m = m_round + sign * offset;
                if (m < 0 || m > 127)
                    continue;
                if (ScaleHelper::IsAllowed(m % 12, allowed_mask_)) {
                    return 440.0f * std::pow(2.0f, static_cast<float>(m - 69) / 12.0f);
                }
            }
        }
        return freq_hz; // fallback（不应到达）
    }

    // ---- 参数 ----
    float sample_rate_ = 44100.0f;
    size_t hop_size_ = 512;
    size_t fft_size_ = 2048;
    bool first_frame_ = true;

    // ---- 调性 ----
    uint16_t allowed_mask_ = 0xFFF;
    int root_note_ = 0;
    ScaleHelper::Type scale_type_ = ScaleHelper::Type::kMajor;

    // ---- DSP ----
    qwqdsp_spectral::RealFFT fft_;
    std::vector<float> window_;
    std::vector<float> fft_in_;
    std::vector<float> fft_out_;
    std::vector<float> prev_phase_;
    std::vector<float> synth_phase_;
};
