#pragma once

#include "../../../nogui/audiofx/pitch_quantize/scale_helper.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <qwqdsp/filter/int_delay.hpp>
#include <qwqdsp/segement/analyze_synthsis_online.hpp>
#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/window/hann.hpp>
#include <vector>

// ------------------------------------------------------------
// PitchQuantizerRTV2 — 实时流式音高量化器 (v2)
// ------------------------------------------------------------
/**
 * @brief 基于 PitchQuantizerRT v1 的改进版本，color 参数扩展为 [0, 200%]。
 *
 * color 行为：
 *   0%   — 频率不改变
 *   100% — 将 bin 频率按 src_midi → target_midi 的比值缩放
 *   200% — 与 v1 100% 一致（完全量化到最近允许 MIDI 频率）
 *
 * 原理同 v1：相邻样本 FFT 估计局部谱峰的瞬时频率，
 *       将谱峰量化到目标 MIDI 音符对应的 bin 位置，
 *       峰域整体平移 + identity phase locking。
 */
class PitchQuantizerRTV2 {
public:
    void Init(float sample_rate, size_t hop_size, size_t fft_size) {
        sr_ = sample_rate;
        hop_ = hop_size;
        fft_ = fft_size;
        nb_ = fft_size / 2 + 1;

        win_.resize(fft_size);
        qwqdsp_window::Hann::Window(win_, true);

        omega_per_bin_ = 2.0f * std::numbers::pi_v<float> / static_cast<float>(fft_size);
        fs_over_2pi_ = sample_rate / (2.0f * std::numbers::pi_v<float>);

        fft_obj_.Init(fft_size);

        fft_buf_.resize(fft_size + 2);
        fft_buf_delay_.resize(fft_size + 2);

        src_gain_.resize(nb_);
        src_omega_.resize(nb_);
        src_phase_.resize(nb_);
        src_peak_.resize(nb_);
        peak_bins_.reserve(nb_);
        peak_target_idx_.resize(nb_);
        peak_target_omega_.resize(nb_);
        synth_gain_.resize(nb_);
        synth_omega_.resize(nb_);
        synth_phase_.assign(nb_, 0.0f);
        synth_source_idx_.resize(nb_);
        synth_peak_idx_.resize(nb_);
        source_target_idx_.resize(nb_);

        prev_last_sample_ = 0.0f;
        first_frame_ = true;

        aso_.SetSize(fft_size);
        aso_.SetHop(hop_size);
        dry_delay_.Init(fft_size + 1);

        SetKey(0, ScaleHelper::Type::kMajor);
    }

    void SetKey(int root, ScaleHelper::Type type) noexcept {
        root_note_ = root % 12;
        scale_type_ = type;
        allowed_mask_ = ScaleHelper::MakeMask(root_note_, type);
    }

    /**
     * @brief 处理一整个 block 的音频
     * @param block 输入/输出缓冲（in-place 处理）
     */
    void Process(std::span<float> block) noexcept {
        // 1. 保存干信号副本
        if (mix_ < 1.0f) {
            dry_buf_.assign(block.begin(), block.end());
        }

        // 2. 处理 block（in-place，湿信号覆写）
        aso_.Process(block, [this](std::span<const float> in, std::span<float> out) { ProcessFrame(in, out); });

        // 3. 干湿比混合：干信号延迟 fft_size 对齐湿信号的分析延迟
        if (mix_ < 1.0f && !dry_buf_.empty()) {
            for (size_t i = 0; i < block.size(); ++i) {
                dry_delay_.Push(dry_buf_[i]);
                float dry = dry_delay_.GetAfterPush(fft_);
                block[i] = mix_ * block[i] + (1.0f - mix_) * dry;
            }
        }
    }

    void Reset() noexcept {
        first_frame_ = true;
        prev_last_sample_ = 0.0f;
        std::fill(synth_phase_.begin(), synth_phase_.end(), 0.0f);
        aso_.Reset();
        dry_delay_.Reset();
    }

    float mix_ = 1.0f;
    float color_ = 1.0f; // [0, 2], 0=off, 1=ratio, 2=full quantize
    float low_cut_hz_ = 0.0f;
    float high_cut_hz_ = 20000.0f;

    void SetColor(float color) noexcept {
        color_ = std::clamp(color, 0.0f, 2.0f);
    }
    void SetLowCut(float hz) noexcept {
        low_cut_hz_ = hz;
    }
    void SetHighCut(float hz) noexcept {
        high_cut_hz_ = hz;
    }
private:
    static constexpr size_t kNoBin = SIZE_MAX;

    void ProcessFrame(std::span<const float> in, std::span<float> out) noexcept {
        const size_t N = fft_;
        const size_t H = hop_;
        const size_t nb = nb_;

        // ---- 标准帧（加窗 + FFT）----
        for (size_t j = 0; j < N; ++j)
            fft_buf_[j] = in[j] * win_[j];
        fft_obj_.FFT(fft_buf_.data(), fft_buf_delay_.data());

        // ---- 延迟 1 样本帧（加窗 + FFT）----
        fft_buf_[0] = prev_last_sample_ * win_[0];
        for (size_t j = 1; j < N; ++j)
            fft_buf_[j] = in[j - 1] * win_[j];
        prev_last_sample_ = in[N - 1];
        fft_obj_.FFT(fft_buf_.data(), fft_buf_.data());

        // ---- 提取幅度和瞬时频率 ----
        for (size_t k = 0; k < nb; ++k) {
            const float re0 = fft_buf_delay_[2 * k];
            const float im0 = fft_buf_delay_[2 * k + 1];
            const float re1 = fft_buf_[2 * k];
            const float im1 = fft_buf_[2 * k + 1];
            const float gain = std::sqrt(re0 * re0 + im0 * im0);
            const std::complex<float> c0{re0, im0};
            const std::complex<float> c1{re1, im1};
            const float omega = std::arg(c0 * std::conj(c1));
            src_gain_[k] = gain;
            src_omega_[k] = omega;
            src_phase_[k] = std::atan2(im0, re0);
        }

        FindSourcePeaks();

        // ---- 峰域重映射 ----
        std::fill(synth_gain_.begin(), synth_gain_.end(), 0.0f);
        std::fill(synth_omega_.begin(), synth_omega_.end(), 0.0f);
        std::fill(synth_source_idx_.begin(), synth_source_idx_.end(), kNoBin);
        std::fill(synth_peak_idx_.begin(), synth_peak_idx_.end(), kNoBin);
        std::fill(source_target_idx_.begin(), source_target_idx_.end(), kNoBin);
        std::fill(peak_target_idx_.begin(), peak_target_idx_.end(), kNoBin);

        const float inv_omega_bin = 1.0f / omega_per_bin_;

        for (const size_t peak : peak_bins_) {
            if (src_gain_[peak] < 1e-6f)
                continue;
            const float f_inst = src_omega_[peak] * fs_over_2pi_;

            // 带通滤波
            const float eff_low = (low_cut_hz_ <= 21.0f) ? 0.0f : low_cut_hz_;
            const float eff_high = (high_cut_hz_ >= 19999.0f) ? 1e9f : high_cut_hz_;
            float f_blend;
            if (f_inst >= eff_low && f_inst <= eff_high) {
                const float f_midi = FindNearestAllowedMidi(f_inst);

                // v2 color 逻辑：两段插值
                //   [0%, 100%]: f_inst → f_inst * ratio
                //   [100%, 200%]: f_inst * ratio → f_midi
                const float midi_float = FreqToMidi(f_inst);
                const float midi_src = std::round(midi_float);
                const float f_midi_src = MidiToFreq(midi_src);
                const float ratio = f_midi / f_midi_src;

                if (color_ <= 1.0f) {
                    // [0%, 100%]: lerp between f_inst and f_inst * ratio
                    f_blend = f_inst * (1.0f + (ratio - 1.0f) * color_);
                }
                else {
                    // [100%, 200%]: lerp between f_inst * ratio and f_midi
                    const float c = color_ - 1.0f; // [0, 1]
                    f_blend = (f_inst * ratio) * (1.0f - c) + f_midi * c;
                }
            }
            else {
                f_blend = f_inst; // 直通，不量化
            }

            const float omega_target = f_blend / fs_over_2pi_;
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
            const std::ptrdiff_t offset = static_cast<std::ptrdiff_t>(k) - static_cast<std::ptrdiff_t>(peak);
            const std::ptrdiff_t target_signed = static_cast<std::ptrdiff_t>(peak_target) + offset;
            if (target_signed <= 0 || target_signed >= static_cast<std::ptrdiff_t>(nb - 1))
                continue;
            const size_t target_idx = static_cast<size_t>(target_signed);
            source_target_idx_[k] = target_idx;
            if (src_gain_[k] > synth_gain_[target_idx]) {
                synth_gain_[target_idx] = src_gain_[k];
                synth_omega_[target_idx] = peak_target_omega_[peak] + omega_per_bin_ * static_cast<float>(offset);
                synth_source_idx_[target_idx] = k;
            }
        }

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

        // ---- 综合相位递推 ----
        if (first_frame_) {
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
                    synth_phase_[k] += omega_per_bin_ * static_cast<float>(k) * static_cast<float>(H);
                    synth_phase_[k] = WrapPi(synth_phase_[k]);
                }
            }
        }

        for (size_t k = 1; k + 1 < nb; ++k) {
            const size_t peak_target = synth_peak_idx_[k];
            if (peak_target == kNoBin || peak_target == k)
                continue;
            const size_t source = synth_source_idx_[k];
            const size_t peak_source = synth_source_idx_[peak_target];
            const float relative_phase = WrapPi(src_phase_[source] - src_phase_[peak_source]);
            synth_phase_[k] = WrapPi(synth_phase_[peak_target] + relative_phase);
        }

        // ---- 合成频谱 → IFFT → 加窗输出 ----
        for (size_t k = 0; k < nb; ++k) {
            const float g = synth_gain_[k];
            const float p = synth_phase_[k];
            fft_buf_[2 * k] = g * std::cos(p);
            fft_buf_[2 * k + 1] = g * std::sin(p);
        }
        fft_obj_.IFFT(fft_buf_.data(), fft_buf_delay_.data());

        constexpr float kWolaScale = 8.0f / 3.0f;
        for (size_t j = 0; j < N; ++j) {
            out[j] = fft_buf_delay_[j] * win_[j] * kWolaScale;
        }
    }

    // ---- 内部工具 ----

    void FindSourcePeaks() noexcept {
        peak_bins_.clear();
        for (size_t k = 1; k + 1 < nb_; ++k) {
            const float gain = src_gain_[k];
            if (gain >= src_gain_[k - 1] && gain >= src_gain_[k + 1]
                && (gain > src_gain_[k - 1] || gain > src_gain_[k + 1])) {
                peak_bins_.push_back(k);
            }
        }
        if (peak_bins_.empty()) {
            for (size_t k = 0; k < nb_; ++k)
                src_peak_[k] = k;
            return;
        }
        size_t peak_index = 0;
        for (size_t k = 0; k < nb_; ++k) {
            while (peak_index + 1 < peak_bins_.size()
                   && k >= peak_bins_[peak_index] + (peak_bins_[peak_index + 1] - peak_bins_[peak_index] + 1) / 2) {
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

    /** @brief 频率 → MIDI 序号（浮点，不取整） */
    static float FreqToMidi(float freq_hz) noexcept {
        return 12.0f * std::log2(freq_hz / 440.0f) + 69.0f;
    }

    /** @brief MIDI 序号 → 频率 */
    static float MidiToFreq(float midi) noexcept {
        return 440.0f * std::pow(2.0f, (midi - 69.0f) / 12.0f);
    }

    float FindNearestAllowedMidi(float freq_hz) const noexcept {
        if (freq_hz <= 0.0f)
            return freq_hz;
        const float m_float = 12.0f * std::log2(freq_hz / 440.0f) + 69.0f;
        const int m_round = static_cast<int>(std::round(m_float));
        for (int offset = 0; offset < 128; ++offset) {
            for (int sign : {1, -1}) {
                const int m = m_round + sign * offset;
                if (m < 0 || m > 127)
                    continue;
                if (ScaleHelper::IsAllowed(m % 12, allowed_mask_))
                    return 440.0f * std::pow(2.0f, static_cast<float>(m - 69) / 12.0f);
            }
        }
        return freq_hz;
    }

    float sr_ = 48000.0f;
    size_t hop_ = 512;
    size_t fft_ = 2048;
    size_t nb_ = 1025;
    bool first_frame_ = true;
    float omega_per_bin_ = 0.0f;
    float fs_over_2pi_ = 0.0f;

    uint16_t allowed_mask_ = 0xFFF;
    int root_note_ = 0;
    ScaleHelper::Type scale_type_ = ScaleHelper::Type::kMajor;

    qwqdsp_spectral::RealFFT fft_obj_;
    std::vector<float> win_;
    qwqdsp_segement::AnalyzeSynthsisOnline aso_;
    qwqdsp_filter::IntDelay dry_delay_;

    std::vector<float> fft_buf_;
    std::vector<float> fft_buf_delay_;
    float prev_last_sample_ = 0.0f;
    std::vector<float> dry_buf_;

    std::vector<float> src_gain_;
    std::vector<float> src_omega_;
    std::vector<float> src_phase_;
    std::vector<size_t> src_peak_;
    std::vector<size_t> peak_bins_;
    std::vector<size_t> peak_target_idx_;
    std::vector<float> peak_target_omega_;

    std::vector<float> synth_gain_;
    std::vector<float> synth_omega_;
    std::vector<float> synth_phase_;
    std::vector<size_t> synth_source_idx_;
    std::vector<size_t> synth_peak_idx_;
    std::vector<size_t> source_target_idx_;
};
