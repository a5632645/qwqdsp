#include <AudioFile.h>
#include <audio_ops.hpp>
#include <qwqdsp/window/hann.hpp>
#include <qwqdsp/spectral/real_fft.hpp>
#include <work_dir.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <format>
#include <iostream>
#include <numbers>
#include <span>
#include <vector>

// ------------------------------------------------------------
// ScaleHelper — 调性辅助工具
// ------------------------------------------------------------
/**
 * @brief 提供音阶级别上的调性判断与 MIDI 音符查找。
 *
 * 音级 (note class) 映射: 0=C, 1=C#, 2=D, ..., 11=B。
 * 通过位掩码 (uint16_t, 低 12 位) 表示哪些音级被允许。
 */
struct ScaleHelper {
    enum class Type : uint8_t {
        kMajor,         ///< 大调音阶
        kMinorNatural,  ///< 自然小调
        kMinorHarmonic, ///< 和声小调
        kMinorMelodic,  ///< 旋律小调（上行）
        kChromatic,     ///< 半音阶（全部 12 个音）
    };

    static constexpr std::array<int, 7> kMajorIntervals{0, 2, 4, 5, 7, 9, 11};
    static constexpr std::array<int, 7> kMinorNaturalIntervals{0, 2, 3, 5, 7, 8, 10};
    static constexpr std::array<int, 7> kMinorHarmonicIntervals{0, 2, 3, 5, 7, 8, 11};
    static constexpr std::array<int, 7> kMinorMelodicIntervals{0, 2, 3, 5, 7, 9, 11};
    static constexpr std::array<int, 12> kChromaticIntervals{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    static uint16_t MakeMask(int root, Type type) noexcept {
        uint16_t mask = 0;
        std::span<const int> intervals;
        switch (type) {
        case Type::kMajor:
            intervals = kMajorIntervals;
            break;
        case Type::kMinorNatural:
            intervals = kMinorNaturalIntervals;
            break;
        case Type::kMinorHarmonic:
            intervals = kMinorHarmonicIntervals;
            break;
        case Type::kMinorMelodic:
            intervals = kMinorMelodicIntervals;
            break;
        case Type::kChromatic:
            intervals = kChromaticIntervals;
            break;
        }
        for (int iv : intervals) {
            mask |= static_cast<uint16_t>(1u << ((root + iv) % 12));
        }
        return mask;
    }

    static bool IsAllowed(int note_class, uint16_t mask) noexcept {
        return (mask >> static_cast<uint16_t>(note_class % 12)) & 1u;
    }

    static const char* TypeName(Type type) noexcept {
        switch (type) {
        case Type::kMajor:
            return "Major";
        case Type::kMinorNatural:
            return "Minor(Natural)";
        case Type::kMinorHarmonic:
            return "Minor(Harmonic)";
        case Type::kMinorMelodic:
            return "Minor(Melodic)";
        case Type::kChromatic:
            return "Chromatic";
        }
        return "Unknown";
    }

    static constexpr const char* kNoteNames[12] = {"C", "C#", "D", "D#", "E", "F",
                                                    "F#", "G", "G#", "A", "A#", "B"};
};

// ------------------------------------------------------------
// PitchQuantizer2 — 基于 Bin-mapping Phase Vocoder 的音高量化器
// ------------------------------------------------------------
/**
 * @brief 基于 bin-mapping 相位声码器的离线音高量化器。
 *
 * 原理（与 pitch_shifter2.hpp 的 PhaseVocoder 相同）：
 *   通过相邻样本 FFT 估计每个 bin 的瞬时频率，
 *   将幅度和相位重新映射到目标 MIDI 音符对应的 bin 位置，
 *   碰撞时保留幅度最大的 bin，再递推综合相位。
 *
 * 优势：
 *   bin 可自由移动到任意位置，频率修正不再受 ±fs/(2·hop) 约束。
 */
template <bool kClampWhenExceed = true>
class PitchQuantizer2 {
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
        max_correction_ = 0.0f; // bin-mapping 模式无硬上限

        // FFT
        fft_.Init(fft_size);

        // 缓冲
        // fft_buf_ 以 CCS 格式存 FFT 输入/输出，大小 N+2
        fft_buf_.resize(fft_size + 2);
        fft_buf_delay_.resize(fft_size + 2);
        // 复数格式的分析/综合缓冲
        src_gain_.resize(num_bins_);
        src_omega_.resize(num_bins_);
        synth_gain_.resize(num_bins_);
        synth_omega_.resize(num_bins_);
        synth_phase_.assign(num_bins_, 0.0f);
        first_frame_ = true;

        // 默认调性：C Major
        allowed_mask_ = ScaleHelper::MakeMask(0, ScaleHelper::Type::kMajor);
        root_note_ = 0;
        scale_type_ = ScaleHelper::Type::kMajor;
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

    void SetKey(int root, ScaleHelper::Type type) noexcept {
        root_note_ = root % 12;
        scale_type_ = type;
        allowed_mask_ = ScaleHelper::MakeMask(root_note_, type);
    }

    std::string KeyDescription() const {
        return std::format("{} {}", ScaleHelper::kNoteNames[root_note_],
                           ScaleHelper::TypeName(scale_type_));
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
            }

            // ----------------------------------------------------
            // Bin 重映射：源 bin → 目标 MIDI 音符位置
            // ----------------------------------------------------
            std::fill(synth_gain_.begin(), synth_gain_.end(), 0.0f);
            std::fill(synth_omega_.begin(), synth_omega_.end(), 0.0f);

            const float inv_omega_bin = 1.0f / omega_per_bin_;

            for (size_t k = 0; k < nb; ++k) {
                if (k == 0 || k == N / 2 || src_gain_[k] < 1e-6f)
                    continue;

                // 瞬时频率 (Hz)
                const float f_inst = src_omega_[k] * fs_over_2pi_;
                // 最近允许 MIDI 频率
                const float f_midi = FindNearestAllowedMidi(f_inst);
                // 目标 omega (rad/sample)
                float omega_target = f_midi / fs_over_2pi_;

                if constexpr (!kClampWhenExceed) {
                    // 检查是否需要跳过
                    const float df = f_midi - f_inst;
                    const float max_corr = fs / (2.0f * static_cast<float>(H));
                    if (std::abs(df) > max_corr)
                        continue;
                }

                // 目标 bin 索引
                size_t target_idx = static_cast<size_t>(omega_target * inv_omega_bin + 0.5f);
                if (target_idx >= nb)
                    target_idx = nb - 1;
                if (target_idx == 0)
                    target_idx = 1;

                // 碰撞处理：保留幅度最大的源
                if (src_gain_[k] > synth_gain_[target_idx]) {
                    synth_gain_[target_idx] = src_gain_[k];
                    synth_omega_[target_idx] = omega_target;
                }
            }

            // ----------------------------------------------------
            // 综合相位递推
            // ----------------------------------------------------
            if (first_frame_) {
                // 第一帧：用分析相位初始化综合相位
                for (size_t k = 0; k < nb; ++k) {
                    const float re = fft_buf_delay_[2 * k];
                    const float im = fft_buf_delay_[2 * k + 1];
                    synth_phase_[k] = std::atan2(im, re);
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
                if (ScaleHelper::IsAllowed(m % 12, allowed_mask_)) {
                    return 440.0f * std::pow(2.0f, static_cast<float>(m - 69) / 12.0f);
                }
            }
        }
        return freq_hz;
    }

    // ---- 参数 ----
    float sample_rate_ = 44100.0f;
    size_t hop_size_ = 256;
    size_t fft_size_ = 2048;
    size_t num_bins_ = 1025;
    bool first_frame_ = true;
    float omega_per_bin_ = 0.0f;     // 2π / N
    float fs_over_2pi_ = 0.0f;       // fs / (2π)
    float max_correction_ = 0.0f;

    // ---- 调性 ----
    uint16_t allowed_mask_ = 0xFFF;
    int root_note_ = 0;
    ScaleHelper::Type scale_type_ = ScaleHelper::Type::kMajor;

    // ---- DSP ----
    qwqdsp_spectral::RealFFT fft_;
    std::vector<float> window_;

    // FFT 缓冲 (CCS 格式，大小 N+2)
    std::vector<float> fft_buf_;
    std::vector<float> fft_buf_delay_;

    // 分析结果
    std::vector<float> src_gain_;
    std::vector<float> src_omega_;

    // 综合结果
    std::vector<float> synth_gain_;
    std::vector<float> synth_omega_;
    std::vector<float> synth_phase_;
};

// ------------------------------------------------------------
// main
// ------------------------------------------------------------
int main() {
    const auto wav_path = qwqdsp_support::WormholeWav();
    std::cout << std::format("loading {}\n", wav_path) << std::flush;

    AudioFile<float> file{wav_path};
    auto& x_vec = file.samples.front();
    const float fs = file.getSampleRate();
    std::cout << std::format("sample_rate={}, len={}\n", fs, x_vec.size()) << std::flush;

    constexpr size_t kFftSize = 2048;

    // 使用 Init2：由最大修正频率自动计算 hop
    // bin-mapping 模式无硬性上限，这里设一个合理的值控制 hop
    PitchQuantizer2<true> pq;
    pq.Init(fs, kFftSize / 4, kFftSize);

    // 调性: C Major
    pq.SetKey(0, ScaleHelper::Type::kMajor);
    std::cout << std::format("key: {}\n", pq.KeyDescription()) << std::flush;

    auto out = pq.Process(x_vec);
    std::cout << std::format("output len={}\n", out.size()) << std::flush;

    qwqdsp_support::AudioOps::Normalize(out);

    file.setNumSamplesPerChannel(out.size());
    file.samples[0] = out;
    file.setNumChannels(1);
    file.save(qwqdsp_support::OutputFile("pitch_quantizer2.wav"));

    std::cout << "saved pitch_quantizer2.wav\n" << std::flush;
    return 0;
}
