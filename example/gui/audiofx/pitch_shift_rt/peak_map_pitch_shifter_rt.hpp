#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>

#include <qwqdsp/window/hann.hpp>
#include <qwqdsp/spectral/real_fft.hpp>

namespace pitch_shift_rt {

// ------------------------------------------------------------
// RealtimePeakMapPitchShifter — 实时峰域映射音高移位器
// ------------------------------------------------------------
/**
 * @brief 使用峰域映射和原始分析相位锁定的实时音高移位器。
 *
 * 离线版 peak_map_pitch_shifter.hpp 的流式化：
 *   通过相邻样本 FFT 估计局部谱峰的瞬时频率，
 *   将谱峰按比例 kp 映射到目标位置，峰域内的 bin 整体平移以保留频谱邻域，
 *   碰撞时保留幅度最大的 bin，再递推综合相位并做 identity phase locking。
 *
 * 实时化改造（与 RealtimePitchShifter 的输入环形缓冲模式同构）：
 *   - 输入逐采样写入 2048 点环形缓冲，每收满 512 个新样本触发一帧处理
 *   - 合成帧 OLA 到固定容量合成环形缓冲，输出每样本读出
 *   - Hann 四倍重叠下 Σw² 为常数 1.5，OLA 归一化因子固定 N/1.5
 *   - 全部固定容量 std::array，音频路径无动态内存分配
 *
 * @tparam kClampWhenExceed true 时峰可任意移动；false 时频率修正
 *         超过 ±fs/(2·hop) 的峰被跳过（保持原位）。
 */
template <bool kClampWhenExceed = true>
class RealtimePeakMapPitchShifter {
public:
    static constexpr std::size_t kFrameSize = 2048;
    static constexpr std::size_t kHopSize = 512;
    static constexpr std::size_t kNumBins = kFrameSize / 2 + 1;  // 1025
    // 合成环形缓冲：帧长 2048 / hop 512，16 hop 容量，够覆盖帧跨度
    static constexpr std::size_t kSynthesisRingSize = kFrameSize * 4;

    /**
     * @brief 初始化 FFT、窗函数和处理状态。
     * @param[in] sample_rate 采样率，单位 Hz（用于瞬时频率换算）
     */
    void init(float sample_rate) {
        sample_rate_ = sample_rate;
        omega_per_bin_ = 2.0f * std::numbers::pi_v<float> / static_cast<float>(kFrameSize);
        fs_over_2pi_ = sample_rate / (2.0f * std::numbers::pi_v<float>);

        fft_.Init(kFrameSize);
        qwqdsp_window::Hann::Window(window_, true);
        reset();
    }

    /**
     * @brief 设置音高移动量（半音语义，正升负降）。
     * @param[in] semitones 半音数，内部限制在 ±12
     */
    void setPitchShift(float semitones) {
        kp_ = std::exp2(std::clamp(semitones, -12.0f, 12.0f) / 12.0f);
    }

    /**
     * @brief 清空流式缓冲和全部算法状态。
     */
    void reset() noexcept {
        input_ring_.fill(0.0f);
        synthesis_ring_.fill(0.0f);
        input_write_ = 0;
        collected_samples_ = 0;
        samples_since_frame_ = 0;
        analysis_started_ = false;
        synthesis_frame_start_ = 0;
        synthesis_available_ = 0;
        synthesis_read_ = 0;
        first_frame_ = true;
        synth_phase_.fill(0.0f);
    }

    /**
     * @brief 处理一段单声道音频（输出与输入等长，稳态延迟 kFrameSize）。
     * @param[in] input 单声道浮点输入缓冲
     * @param[out] output 单声道浮点输出缓冲
     * @param[in] frame_count 缓冲中的采样数
     */
    void process(const float* input, float* output, std::size_t frame_count) {
        if (output == nullptr) {
            return;
        }
        if (input == nullptr) {
            std::fill_n(output, frame_count, 0.0f);
            return;
        }

        for (std::size_t i = 0; i < frame_count; ++i) {
            input_ring_[input_write_] = input[i];
            input_write_ = (input_write_ + 1) % kFrameSize;
            collected_samples_ = std::min(collected_samples_ + 1, kFrameSize);

            if (!analysis_started_ && collected_samples_ == kFrameSize) {
                processFrame();
                analysis_started_ = true;
                samples_since_frame_ = 0;
            }
            else if (analysis_started_) {
                ++samples_since_frame_;
                if (samples_since_frame_ >= kHopSize) {
                    processFrame();
                    samples_since_frame_ = 0;
                }
            }

            output[i] = renderOutputSample();
        }
    }

    /**
     * @brief 返回稳态处理延迟。
     * @return 延迟采样数
     */
    static constexpr std::size_t latencySamples() noexcept {
        return kFrameSize;
    }

private:
    static constexpr std::size_t kNoBin = std::numeric_limits<std::size_t>::max();
    // 分析/综合共用 periodic Hann，2048/512 四倍重叠下任意点的 Σw ≡ 2.0，
    // OLA 归一化因子取 1/Σw，保证输出幅度与输入一致
    // （离线版沿用 quantizer4 的 N/Σw² 缩放，会放大千倍，仅靠 Normalize 兜底，
    //   实时输出直通声卡，必须用正确的 1/Σw）。
    static constexpr float kOlaScale = 0.5f;

    /**
     * @brief 处理一帧：分析瞬时频率 → 峰域映射 → 综合相位 → IFFT → OLA。
     */
    void processFrame() noexcept {
        // ----------------------------------------------------
        // 分析：相邻样本 FFT 法估计瞬时频率
        // ----------------------------------------------------
        // 标准帧（input_write_ 指向环形缓冲中最旧样本）
        for (std::size_t j = 0; j < kFrameSize; ++j)
            fft_buf_[j] = input_ring_[(input_write_ + j) % kFrameSize] * window_[j];
        fft_.FFT(fft_buf_.data(), fft_buf_delay_.data());
        // fft_buf_delay_ 现在存 CCS 格式: [re0, im0, re1, im1, ...]

        // 延迟 1 样本帧（窗首样本为 0，j=0 处自然为零）
        for (std::size_t j = 0; j < kFrameSize; ++j)
            fft_buf_[j] = input_ring_[(input_write_ + j - 1 + kFrameSize) % kFrameSize] * window_[j];
        fft_.FFT(fft_buf_.data(), fft_buf_.data());
        // fft_buf_ 现在存延迟帧的 CCS 格式

        // 提取幅度和瞬时频率（omega = rad/sample）
        for (std::size_t k = 0; k < kNumBins; ++k) {
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
        findSourcePeaks();

        // ----------------------------------------------------
        // 峰域重映射：源峰 → kp 比例目标峰，邻域整体平移
        // ----------------------------------------------------
        synth_gain_.fill(0.0f);
        synth_omega_.fill(0.0f);
        synth_source_idx_.fill(kNoBin);
        synth_peak_idx_.fill(kNoBin);
        source_target_idx_.fill(kNoBin);
        peak_target_idx_.fill(kNoBin);

        const float inv_omega_bin = 1.0f / omega_per_bin_;

        // 仅映射局部谱峰；峰域内的其他 bin 跟随峰整体移动。
        for (std::size_t peak_index = 0; peak_index < peak_count_; ++peak_index) {
            const std::size_t peak = peak_bins_[peak_index];
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
                const float max_corr = sample_rate_ / (2.0f * static_cast<float>(kHopSize));
                if (std::fabs(df) > max_corr)
                    continue;
            }

            // 目标 bin 索引
            std::size_t target_idx = static_cast<std::size_t>(omega_target * inv_omega_bin + 0.5f);
            if (target_idx == 0 || target_idx >= kNumBins - 1)
                continue;

            peak_target_idx_[peak] = target_idx;
            peak_target_omega_[peak] = omega_target;
        }

        for (std::size_t k = 1; k + 1 < kNumBins; ++k) {
            if (src_gain_[k] < 1e-6f)
                continue;

            const std::size_t peak = src_peak_[k];
            const std::size_t peak_target = peak_target_idx_[peak];
            if (peak_target == kNoBin)
                continue;

            const std::ptrdiff_t offset = static_cast<std::ptrdiff_t>(k)
                                          - static_cast<std::ptrdiff_t>(peak);
            const std::ptrdiff_t target_signed = static_cast<std::ptrdiff_t>(peak_target) + offset;
            if (target_signed <= 0 || target_signed >= static_cast<std::ptrdiff_t>(kNumBins - 1))
                continue;

            const std::size_t target_idx = static_cast<std::size_t>(target_signed);
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
        for (std::size_t k = 1; k + 1 < kNumBins; ++k) {
            if (synth_source_idx_[k] == kNoBin)
                continue;

            const std::size_t peak_source = src_peak_[synth_source_idx_[k]];
            const std::size_t peak_target = source_target_idx_[peak_source];
            if (peak_target < kNumBins && synth_source_idx_[peak_target] == peak_source)
                synth_peak_idx_[k] = peak_target;
            else
                synth_peak_idx_[k] = k;
        }

        // ----------------------------------------------------
        // 综合相位递推
        // ----------------------------------------------------
        if (first_frame_) {
            // 第一帧：用分析相位初始化综合相位
            for (std::size_t k = 0; k < kNumBins; ++k) {
                if (synth_source_idx_[k] != kNoBin)
                    synth_phase_[k] = src_phase_[synth_source_idx_[k]];
            }
            first_frame_ = false;
        }
        else {
            for (std::size_t k = 0; k < kNumBins; ++k) {
                if (synth_gain_[k] > 0.0f) {
                    synth_phase_[k] += synth_omega_[k] * static_cast<float>(kHopSize);
                    synth_phase_[k] = WrapPi(synth_phase_[k]);
                }
                else {
                    // 无内容的 bin：按 bin 中心频率递推
                    synth_phase_[k] += omega_per_bin_ * static_cast<float>(k)
                                       * static_cast<float>(kHopSize);
                    synth_phase_[k] = WrapPi(synth_phase_[k]);
                }
            }
        }

        // Identity phase locking：保持每个谱峰邻域内的分析相对相位。
        for (std::size_t k = 1; k + 1 < kNumBins; ++k) {
            const std::size_t peak_target = synth_peak_idx_[k];
            if (peak_target == kNoBin || peak_target == k)
                continue;

            const std::size_t source = synth_source_idx_[k];
            const std::size_t peak_source = synth_source_idx_[peak_target];
            const float relative_phase = WrapPi(src_phase_[source] - src_phase_[peak_source]);
            synth_phase_[k] = WrapPi(synth_phase_[peak_target] + relative_phase);
        }

        // ----------------------------------------------------
        // 综合：重构频谱 → IFFT → 加窗 → OLA
        // ----------------------------------------------------
        for (std::size_t k = 0; k < kNumBins; ++k) {
            const float g = synth_gain_[k];
            const float p = synth_phase_[k];
            fft_buf_[2 * k] = g * std::cos(p);
            fft_buf_[2 * k + 1] = g * std::sin(p);
        }

        fft_.IFFT(fft_buf_.data(), fft_buf_delay_.data());

        for (std::size_t j = 0; j < kFrameSize; ++j) {
            const std::size_t index = static_cast<std::size_t>(
                (synthesis_frame_start_ + j) % kSynthesisRingSize);
            synthesis_ring_[index] += fft_buf_delay_[j] * window_[j];
        }
        // 新合成帧在 hop 边界处的窗值为零，因此该边界样本已完整可读。
        synthesis_available_ = synthesis_frame_start_ + kHopSize + 1;
        synthesis_frame_start_ += kHopSize;
    }

    /**
     * @brief 查找局部谱峰，并为每个源 bin 分配最近的峰。
     *
     * 仅使用内部 bin 作为锚点，直流和 Nyquist bin 维持其自身索引。
     */
    void findSourcePeaks() noexcept {
        peak_count_ = 0;
        for (std::size_t k = 1; k + 1 < kNumBins; ++k) {
            const float gain = src_gain_[k];
            if (gain >= src_gain_[k - 1] && gain >= src_gain_[k + 1]
                && (gain > src_gain_[k - 1] || gain > src_gain_[k + 1])) {
                peak_bins_[peak_count_++] = k;
            }
        }

        if (peak_count_ == 0) {
            for (std::size_t k = 0; k < kNumBins; ++k)
                src_peak_[k] = k;
            return;
        }

        std::size_t peak_index = 0;
        for (std::size_t k = 0; k < kNumBins; ++k) {
            while (peak_index + 1 < peak_count_
                   && k >= peak_bins_[peak_index]
                               + (peak_bins_[peak_index + 1] - peak_bins_[peak_index] + 1) / 2) {
                ++peak_index;
            }
            src_peak_[k] = peak_bins_[peak_index];
        }
    }

    static float WrapPi(float x) noexcept {
        constexpr float kTwoPi = 2.0f * std::numbers::pi_v<float>;
        while (x > std::numbers::pi_v<float>)
            x -= kTwoPi;
        while (x < -std::numbers::pi_v<float>)
            x += kTwoPi;
        return x;
    }

    /**
     * @brief 渲染一个输出样本：从合成环形缓冲读出并做 OLA 归一化。
     */
    float renderOutputSample() noexcept {
        if (synthesis_read_ + 1 > synthesis_available_)
            return 0.0f;
        const std::size_t index = static_cast<std::size_t>(synthesis_read_ % kSynthesisRingSize);
        const float x = synthesis_ring_[index];
        synthesis_ring_[index] = 0.0f;
        ++synthesis_read_;
        return x * kOlaScale;
    }

    // ---- 参数 ----
    float sample_rate_ = 48000.0f;
    float kp_ = 1.0f;
    bool first_frame_ = true;
    float omega_per_bin_ = 0.0f;  // 2π / N
    float fs_over_2pi_ = 0.0f;    // fs / (2π)

    // ---- DSP ----
    qwqdsp_spectral::RealFFT fft_;
    std::array<float, kFrameSize> window_{};

    // 输入环形缓冲（分析窗的滑动窗口）
    std::array<float, kFrameSize> input_ring_{};
    // 合成环形缓冲（OLA 累计）
    std::array<float, kSynthesisRingSize> synthesis_ring_{};

    // FFT 缓冲 (CCS 格式，大小 N+2)
    std::array<float, kFrameSize + 2> fft_buf_{};
    std::array<float, kFrameSize + 2> fft_buf_delay_{};

    // 分析结果
    std::array<float, kNumBins> src_gain_{};
    std::array<float, kNumBins> src_omega_{};
    std::array<float, kNumBins> src_phase_{};
    std::array<std::size_t, kNumBins> src_peak_{};
    std::array<std::size_t, kNumBins> peak_bins_{};
    std::size_t peak_count_ = 0;
    std::array<std::size_t, kNumBins> peak_target_idx_{};
    std::array<float, kNumBins> peak_target_omega_{};

    // 综合结果
    std::array<float, kNumBins> synth_gain_{};
    std::array<float, kNumBins> synth_omega_{};
    std::array<float, kNumBins> synth_phase_{};
    std::array<std::size_t, kNumBins> synth_source_idx_{};
    std::array<std::size_t, kNumBins> synth_peak_idx_{};
    std::array<std::size_t, kNumBins> source_target_idx_{};

    // 流式状态
    std::size_t input_write_ = 0;
    std::size_t collected_samples_ = 0;
    std::size_t samples_since_frame_ = 0;
    bool analysis_started_ = false;
    std::uint64_t synthesis_frame_start_ = 0;
    std::uint64_t synthesis_available_ = 0;
    std::uint64_t synthesis_read_ = 0;
};

} // namespace pitch_shift_rt
