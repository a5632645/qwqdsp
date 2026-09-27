#pragma once

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

namespace qwqdsp_test {

// ------------------------------------------------------------
// PeakMapPitchShifter — 峰域映射 + Phase Lock 音高移位器
// ------------------------------------------------------------
/**
 * @brief 使用峰域映射和原始分析相位锁定的离线音高移位器。
 *
 * 移植自 pitch_quantize/pitch_quantizer4.hpp 的峰域映射技术，
 * 将"量化到最近允许 MIDI 音符"改为"按 kp 比例映射"：
 *   通过相邻样本 FFT 估计局部谱峰的瞬时频率，
 *   将谱峰按比例 kp 映射到目标位置，
 *   并将峰域内的 bin 整体平移以保留频谱邻域和相位关系。
 *   碰撞时保留幅度最大的 bin，再递推综合相位。
 *
 * 优势：
 *   bin 可自由移动到任意位置，频率修正不再受 ±fs/(2·hop) 约束。
 *
 * @tparam kClampWhenExceed true 时峰可任意移动；false 时频率修正
 *         超过 ±fs/(2·hop) 的峰被跳过（保持原位）。
 */
template <bool kClampWhenExceed = true>
class PeakMapPitchShifter {
public:
    /**
     * @brief 初始化（直接指定 hop）
     */
    void Init(float sample_rate, size_t hop_size, size_t fft_size) {
        sample_rate_ = sample_rate;
        hop_size_ = hop_size;
        fft_size_ = fft_size;
        num_bins_ = fft_size / 2 + 1;

        // Hann 窗（分析/综合共用）
        window_.resize(fft_size);
        qwqdsp_window::Hann::Window(window_, true);

        omega_per_bin_ = 2.0f * std::numbers::pi_v<float> / static_cast<float>(fft_size);
        fs_over_2pi_ = sample_rate / (2.0f * std::numbers::pi_v<float>);

        // FFT
        fft_.Init(fft_size);

        // 缓冲
        // fft_buf_ 以 CCS 格式存 FFT 输入/输出，大小 N+2
        fft_buf_.resize(fft_size + 2);
        fft_buf_delay_.resize(fft_size + 2);
        // 复数格式的分析/综合缓冲
        src_gain_.resize(num_bins_);
        src_omega_.resize(num_bins_);
        src_phase_.resize(num_bins_);
        src_peak_.resize(num_bins_);
        peak_bins_.reserve(num_bins_);
        peak_target_idx_.resize(num_bins_);
        peak_target_omega_.resize(num_bins_);
        synth_gain_.resize(num_bins_);
        synth_omega_.resize(num_bins_);
        synth_phase_.assign(num_bins_, 0.0f);
        synth_source_idx_.resize(num_bins_);
        synth_peak_idx_.resize(num_bins_);
        source_target_idx_.resize(num_bins_);
        first_frame_ = true;
    }

    /**
     * @brief 以最大修正频率初始化（自动计算满足 Hann COLA 的 hop）
     *
     * 即使 bin-mapping 模式无硬性频率修正上限，此接口仍用于
     * 控制 hop 大小以平衡时间/频率分辨率。
     */
    void Init2(float sample_rate, float max_correction_freq, size_t fft_size) {
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
     * @brief 设置音高移位比例（kp>1 升调，kp<1 降调，kp=1 不变）。
     */
    void SetPitchShift(float kp) noexcept {
        kp_ = kp;
    }

    std::vector<float> Process(std::span<const float> input) {
        const size_t N = fft_size_;
        const size_t H = hop_size_;
        const float fs = sample_rate_;
        const size_t nb = num_bins_;

        if (input.size() < N)
            return {};
        const size_t num_frames = (input.size() - N) / H + 1;
        const size_t output_len = H * num_frames + N;
        std::vector<float> output(output_len, 0.0f);

        std::cout << std::format("  frames={}, fft={}, hop={}\n", num_frames, N, H)
                  << std::flush;

        for (size_t i = 0; i < num_frames; ++i) {
            const float* frame = input.data() + i * H;

            // ----------------------------------------------------
            // 分析：相邻样本 FFT 法估计瞬时频率
            // ----------------------------------------------------
            // 标准帧
            for (size_t j = 0; j < N; ++j)
                fft_buf_[j] = frame[j] * window_[j];
            fft_.FFT(fft_buf_.data(), fft_buf_delay_.data());
            // fft_buf_delay_ 现在存 CCS 格式: [re0, im0, re1, im1, ...]

            // 延迟 1 样本帧
            for (size_t j = 1; j < N; ++j)
                fft_buf_[j] = frame[j - 1] * window_[j];
            fft_buf_[0] = 0.0f;
            fft_.FFT(fft_buf_.data(), fft_buf_.data());
            // fft_buf_ 现在存延迟帧的 CCS 格式

            // 提取幅度和瞬时频率（omega = rad/sample）
            for (size_t k = 0; k < nb; ++k) {
                const float re0 = fft_buf_delay_[2 * k];
                const float im0 = fft_buf_delay_[2 * k + 1];
                const float re1 = fft_buf_[2 * k];
                const float im1 = fft_buf_[2 * k + 1];

                const float gain = std::sqrt(re0 * re0 + im0 * im0);
                // X(k) · conj(X_f(k)) — 瞬时频率
                const std::complex<float> c0{re0, im0};
                const std::complex<float> c1{re1, im1};
                const float omega = std::arg(c0 * std::conj(c1));

                src_gain_[k] = gain;
                src_omega_[k] = omega;
                src_phase_[k] = std::atan2(im0, re0);
            }

            // 为每个源 bin 分配最近的局部谱峰，用于 identity phase locking。
            FindSourcePeaks();

            // ----------------------------------------------------
            // 峰域重映射：源峰 → kp 比例目标峰，邻域整体平移
            // ----------------------------------------------------
            std::fill(synth_gain_.begin(), synth_gain_.end(), 0.0f);
            std::fill(synth_omega_.begin(), synth_omega_.end(), 0.0f);
            std::fill(synth_source_idx_.begin(), synth_source_idx_.end(), kNoBin);
            std::fill(synth_peak_idx_.begin(), synth_peak_idx_.end(), kNoBin);
            std::fill(source_target_idx_.begin(), source_target_idx_.end(), kNoBin);
            std::fill(peak_target_idx_.begin(), peak_target_idx_.end(), kNoBin);

            const float inv_omega_bin = 1.0f / omega_per_bin_;

            // 仅映射局部谱峰；峰域内的其他 bin 跟随峰整体移动。
            for (const size_t peak : peak_bins_) {
                if (src_gain_[peak] < 1e-6f)
                    continue;

                // 瞬时频率 (Hz)
                const float f_inst = src_omega_[peak] * fs_over_2pi_;
                // 目标频率 = 源频率 × kp
                const float f_target = f_inst * kp_;
                // 目标 omega (rad/sample)
                const float omega_target = f_target / fs_over_2pi_;

                if constexpr (!kClampWhenExceed) {
                    // 检查是否需要跳过
                    const float df = f_target - f_inst;
                    const float max_corr = fs / (2.0f * static_cast<float>(H));
                    if (std::abs(df) > max_corr)
                        continue;
                }

                // 目标 bin 索引
                size_t target_idx = static_cast<size_t>(omega_target * inv_omega_bin + 0.5f);
                if (target_idx == 0 || target_idx >= nb - 1)
                    continue;

                peak_target_idx_[peak] = target_idx;
                peak_target_omega_[peak] = omega_target;
            }

            for (size_t k = 1; k + 1 < nb; ++k) {
                if (src_gain_[k] < 1e-6f)
                    continue;

                const size_t peak = src_peak_[k];
                const size_t peak_target = peak_target_idx_[peak];
                if (peak_target == kNoBin)
                    continue;

                const std::ptrdiff_t offset = static_cast<std::ptrdiff_t>(k)
                                              - static_cast<std::ptrdiff_t>(peak);
                const std::ptrdiff_t target_signed = static_cast<std::ptrdiff_t>(peak_target) + offset;
                if (target_signed <= 0 || target_signed >= static_cast<std::ptrdiff_t>(nb - 1))
                    continue;

                const size_t target_idx = static_cast<size_t>(target_signed);
                source_target_idx_[k] = target_idx;

                // 碰撞处理：保留幅度最大的源
                if (src_gain_[k] > synth_gain_[target_idx]) {
                    synth_gain_[target_idx] = src_gain_[k];
                    synth_omega_[target_idx] = peak_target_omega_[peak]
                                               + omega_per_bin_ * static_cast<float>(offset);
                    synth_source_idx_[target_idx] = k;
                }
            }

            // 只有仍保留在输出频谱中的峰才能作为相位锁定锚点。
            for (size_t k = 1; k + 1 < nb; ++k) {
                if (synth_source_idx_[k] == kNoBin)
                    continue;

                const size_t peak_source = src_peak_[synth_source_idx_[k]];
                const size_t peak_target = source_target_idx_[peak_source];
                if (peak_target < nb && synth_source_idx_[peak_target] == peak_source)
                    synth_peak_idx_[k] = peak_target;
                else
                    synth_peak_idx_[k] = k;
            }

            // ----------------------------------------------------
            // 综合相位递推
            // ----------------------------------------------------
            if (first_frame_) {
                // 第一帧：用分析相位初始化综合相位
                for (size_t k = 0; k < nb; ++k) {
                    if (synth_source_idx_[k] != kNoBin)
                        synth_phase_[k] = src_phase_[synth_source_idx_[k]];
                }
                first_frame_ = false;
            }
            else {
                for (size_t k = 0; k < nb; ++k) {
                    if (synth_gain_[k] > 0.0f) {
                        synth_phase_[k] += synth_omega_[k] * static_cast<float>(H);
                        synth_phase_[k] = WrapPi(synth_phase_[k]);
                    }
                    else {
                        // 无内容的 bin：按 bin 中心频率递推
                        synth_phase_[k] += omega_per_bin_ * static_cast<float>(k)
                                           * static_cast<float>(H);
                        synth_phase_[k] = WrapPi(synth_phase_[k]);
                    }
                }
            }

            // Identity phase locking：保持每个谱峰邻域内的分析相对相位。
            for (size_t k = 1; k + 1 < nb; ++k) {
                const size_t peak_target = synth_peak_idx_[k];
                if (peak_target == kNoBin || peak_target == k)
                    continue;

                const size_t source = synth_source_idx_[k];
                const size_t peak_source = synth_source_idx_[peak_target];
                const float relative_phase = WrapPi(src_phase_[source] - src_phase_[peak_source]);
                synth_phase_[k] = WrapPi(synth_phase_[peak_target] + relative_phase);
            }

            // ----------------------------------------------------
            // 综合：重构频谱 → IFFT → 加窗 → OLA
            // ----------------------------------------------------
            // 用 synth_gain_ / synth_phase_ 填充 CCS 格式
            for (size_t k = 0; k < nb; ++k) {
                const float g = synth_gain_[k];
                const float p = synth_phase_[k];
                fft_buf_[2 * k] = g * std::cos(p);
                fft_buf_[2 * k + 1] = g * std::sin(p);
            }

            fft_.IFFT(fft_buf_.data(), fft_buf_delay_.data());

            for (size_t j = 0; j < N; ++j) {
                output[i * H + j] += fft_buf_delay_[j] * window_[j];
            }

        }

        // ---- OLA 增益归一化 ----
        {
            std::vector<float> gain_buf(output_len, 0.0f);
            for (size_t i = 0; i < num_frames; ++i) {
                for (size_t j = 0; j < N; ++j) {
                    gain_buf[i * H + j] += window_[j] * window_[j] / static_cast<float>(N);
                }
            }
            float peak_gain = 0.0f;
            for (float g : gain_buf)
                peak_gain = std::max(peak_gain, g);
            const float scale = (peak_gain > 1e-10f) ? (1.0f / peak_gain) : 1.0f;
            for (float& x : output)
                x *= scale;

            std::cout << std::format("  ola_peak_gain={:.4f}, scale={:.4f}\n", peak_gain, scale)
                      << std::flush;
        }

        return output;
    }

    void Reset() noexcept {
        first_frame_ = true;
        std::fill(synth_phase_.begin(), synth_phase_.end(), 0.0f);
    }

private:
    static constexpr size_t kNoBin = std::numeric_limits<size_t>::max();

    /**
     * @brief 查找局部谱峰，并为每个源 bin 分配最近的峰。
     *
     * 仅使用内部 bin 作为锚点，直流和 Nyquist bin 维持其自身索引。
     */
    void FindSourcePeaks() {
        peak_bins_.clear();
        for (size_t k = 1; k + 1 < num_bins_; ++k) {
            const float gain = src_gain_[k];
            if (gain >= src_gain_[k - 1] && gain >= src_gain_[k + 1]
                && (gain > src_gain_[k - 1] || gain > src_gain_[k + 1])) {
                peak_bins_.push_back(k);
            }
        }

        if (peak_bins_.empty()) {
            for (size_t k = 0; k < num_bins_; ++k)
                src_peak_[k] = k;
            return;
        }

        size_t peak_index = 0;
        for (size_t k = 0; k < num_bins_; ++k) {
            while (peak_index + 1 < peak_bins_.size()
                   && k >= peak_bins_[peak_index]
                               + (peak_bins_[peak_index + 1] - peak_bins_[peak_index] + 1) / 2) {
                ++peak_index;
            }
            src_peak_[k] = peak_bins_[peak_index];
        }
    }

    static float WrapPi(float x) noexcept {
        constexpr float two_pi = 2.0f * std::numbers::pi_v<float>;
        while (x > std::numbers::pi_v<float>)
            x -= two_pi;
        while (x < -std::numbers::pi_v<float>)
            x += two_pi;
        return x;
    }

    static size_t FindColaHop(size_t fft_size, size_t desired) noexcept {
        const size_t kMaxHop = fft_size / 4;
        if (desired < 1)
            desired = 1;
        if (desired > kMaxHop)
            desired = kMaxHop;

        for (size_t offset = 0; offset <= kMaxHop; ++offset) {
            if (offset < desired) {
                const size_t cand = desired - offset;
                if (cand >= 1 && fft_size % cand == 0)
                    return cand;
            }
            if (desired + offset <= kMaxHop) {
                const size_t cand = desired + offset;
                if (fft_size % cand == 0)
                    return cand;
            }
        }
        return kMaxHop;
    }

    // ---- 参数 ----
    float sample_rate_ = 44100.0f;
    size_t hop_size_ = 256;
    size_t fft_size_ = 2048;
    size_t num_bins_ = 1025;
    float kp_ = 1.0f;
    bool first_frame_ = true;
    float omega_per_bin_ = 0.0f;  // 2π / N
    float fs_over_2pi_ = 0.0f;    // fs / (2π)

    // ---- DSP ----
    qwqdsp_spectral::RealFFT fft_;
    std::vector<float> window_;

    // FFT 缓冲 (CCS 格式，大小 N+2)
    std::vector<float> fft_buf_;
    std::vector<float> fft_buf_delay_;

    // 分析结果
    std::vector<float> src_gain_;
    std::vector<float> src_omega_;
    std::vector<float> src_phase_;
    std::vector<size_t> src_peak_;
    std::vector<size_t> peak_bins_;
    std::vector<size_t> peak_target_idx_;
    std::vector<float> peak_target_omega_;

    // 综合结果
    std::vector<float> synth_gain_;
    std::vector<float> synth_omega_;
    std::vector<float> synth_phase_;
    std::vector<size_t> synth_source_idx_;
    std::vector<size_t> synth_peak_idx_;
    std::vector<size_t> source_target_idx_;
};

} // namespace qwqdsp_test
