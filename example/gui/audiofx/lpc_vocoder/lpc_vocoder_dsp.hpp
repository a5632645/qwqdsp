#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>

#include "qwqdsp/adaptive/burg_lp.hpp"
#include "qwqdsp/oscillator/noise.hpp"
#include "qwqdsp/pitch/yin.hpp"
#include "qwqdsp/window/hann.hpp"

namespace lpc_vocoder::detail {

/**
 * @brief 单帧分析结果
 */
struct VocoderStats {
    float f0_hz{};          ///< 有声帧的基频，清音帧为零
    float voicing{};        ///< 脉冲激励权重，范围为零到一
    float input_peak{};     ///< 当前缓冲的输入峰值
    float output_peak{};    ///< 当前缓冲的输出峰值
};

/**
 * @brief 用 Burg LPC 分析与基频同步脉冲激励合成语音
 *
 * 每 256 个采样更新一次模型：Burg 求反射系数，YIN 求基频与清浊音；
 * 浊音用 Rosenberg 声门流量导数激励，清音用白噪声，两者按等功率混合。
 * 激励幅度由输入帧 RMS 与完整合成链的稳态增益换算；合成端用反射系数
 * 格型滤波器实现全极点模型，并在一个 hop 内线性插值反射系数、基频与增益。
 */
class LpcVocoder {
public:
    static constexpr std::size_t kPitchFrame = 2048;
    static constexpr std::size_t kLpcFrame = 1024;
    static constexpr std::size_t kHop = 256;
    static constexpr std::size_t kMaxOrder = 80;
    static constexpr std::size_t kImpulseLength = 1024;
    static constexpr float kMinPitchHz = 70.0f;
    static constexpr float kMaxPitchHz = 500.0f;

    /// @brief 初始化分析器、合成器和窗函数
    /// @param[in] sample_rate 采样率，单位为 Hz
    void Init(float sample_rate) {
        fs_ = sample_rate;
        yin_.Init(fs_, static_cast<int>(kPitchFrame));
        yin_.SetMinPitch(kMinPitchHz);
        yin_.SetMaxPitch(kMaxPitchHz);
        yin_.SetThreshold(0.2f);
        burg_.Init(kLpcFrame);
        qwqdsp_window::Hann::Window(window_, true);

        ring_.fill(0.0f);
        lpc_block_.fill(0.0f);
        pitch_block_.fill(0.0f);
        k_target_.fill(0.0f);
        k_current_.fill(0.0f);
        k_step_.fill(0.0f);
        lattice_.fill(0.0f);

        ring_pos_ = 0;
        samples_seen_ = 0;
        hop_counter_ = 0;
        phase_ = 0.0f;
        deemphasis_ = 0.0f;
        dc_input_ = 0.0f;
        dc_output_ = 0.0f;
        pulse_previous_ = 0.0f;
        f0_current_ = 200.0f;
        f0_target_ = 200.0f;
        f0_step_ = 0.0f;
        pulse_amp_current_ = 0.0f;
        pulse_amp_target_ = 0.0f;
        pulse_amp_step_ = 0.0f;
        noise_amp_current_ = 0.0f;
        noise_amp_target_ = 0.0f;
        noise_amp_step_ = 0.0f;
        has_model_ = false;
        stats_ = {};
        noise_.SetSeed(0x9E3779B9u);
    }

    /// @brief 处理一段单声道音频
    /// @param[in] input 输入采样缓冲
    /// @param[out] output 输出采样缓冲
    /// @param[in] count 缓冲中的采样数
    /// @param[in] order LPC 阶数，内部会夹到有效范围
    /// @param[in] mix 干湿混合，一为纯合成
    /// @param[in] noise_mix 噪声激励比例，零为纯脉冲、一为纯噪声
    /// @param[in] voicing_threshold 清浊音判定阈值，YIN 非周期比低于它才判为浊音
    /// @return 本段音频的基频、清浊音与峰值状态
    VocoderStats Process(const float* input, float* output, std::size_t count, int order, float mix, float noise_mix,
                         float voicing_threshold) noexcept {
        order = std::clamp(order, 4, static_cast<int>(kMaxOrder));
        mix = std::clamp(mix, 0.0f, 1.0f);
        noise_mix_ = std::clamp(noise_mix, 0.0f, 1.0f);
        voicing_threshold_ = std::clamp(voicing_threshold, 0.05f, 0.95f);

        float input_peak = 0.0f;
        float output_peak = 0.0f;
        const bool ready = samples_seen_ >= kPitchFrame;
        for (std::size_t i = 0; i < count; ++i) {
            const float x = input[i];
            input_peak = std::max(input_peak, std::abs(x));

            ring_[ring_pos_] = x;
            ring_pos_ = (ring_pos_ + 1) % kPitchFrame;
            ++samples_seen_;

            if (samples_seen_ >= kPitchFrame) {
                ++hop_counter_;
                if (hop_counter_ >= kHop) {
                    hop_counter_ = 0;
                    Analyze(order);
                }
            }

            const float wet = ready ? Synthesize() : x;
            const float y = std::clamp(mix * wet + (1.0f - mix) * x, -1.0f, 1.0f);
            output[i] = y;
            output_peak = std::max(output_peak, std::abs(y));
        }

        stats_.input_peak = input_peak;
        stats_.output_peak = output_peak;
        return stats_;
    }

private:
    /// @brief 对最近一帧做基频和 LPC 分析并设置下一 hop 的目标参数
    void Analyze(int order) noexcept {
        for (std::size_t i = 0; i < kPitchFrame; ++i) {
            pitch_block_[i] = ring_[(ring_pos_ + i) % kPitchFrame];
        }
        yin_.Process(pitch_block_);
        const auto pitch = yin_.GetPitch();

        float voicing = 0.0f;
        if (pitch.pitch_hz >= kMinPitchHz && pitch.pitch_hz <= kMaxPitchHz) {
            voicing = std::clamp((voicing_threshold_ - pitch.non_period_ratio) / kVoiceRampRange, 0.0f, 1.0f);
        }
        // 等功率交叉淡化：旋钮为零出纯脉冲，为一出纯噪声
        const float pulse_share = voicing * std::cos(noise_mix_ * std::numbers::pi_v<float> * 0.5f);

        // 预加重只用于分析；合成端由去加重恢复低频
        const std::size_t lpc_start = kPitchFrame - kLpcFrame;
        float previous = ring_[(ring_pos_ + lpc_start - 1) % kPitchFrame];
        float frame_energy = 0.0f;
        for (std::size_t i = 0; i < kLpcFrame; ++i) {
            const float x = ring_[(ring_pos_ + lpc_start + i) % kPitchFrame];
            const float pre = x - kPreEmphasis * previous;
            previous = x;
            frame_energy += x * x;
            lpc_block_[i] = pre * window_[i];
        }
        const float frame_rms = std::sqrt(frame_energy / static_cast<float>(kLpcFrame));

        std::array<float, kMaxOrder> reflection{};
        burg_.Process(std::span<float>{lpc_block_}, std::span<float>{reflection}.first(order));
        // 静音帧会让 Burg 的归一化除零，得到 NaN；非有限值退化为透明级
        for (std::size_t i = 0; i < static_cast<std::size_t>(order); ++i) {
            reflection[i] = std::isfinite(reflection[i]) ? std::clamp(reflection[i], -0.999f, 0.999f) : 0.0f;
        }

        float f0 = f0_target_;
        if (voicing > 0.0f) {
            f0 = std::clamp(pitch.pitch_hz, kMinPitchHz, kMaxPitchHz);
        }

        for (std::size_t i = 0; i < kMaxOrder; ++i) {
            k_target_[i] = i < static_cast<std::size_t>(order) ? reflection[i] : 0.0f;
        }

        // 仿真直接使用与合成相同的格型结构，避免高阶直接型多项式的数值病态
        auto lattice_tick = [order, &reflection](float excitation,
                                                std::array<float, kMaxOrder + 1>& state) noexcept {
            float forward = excitation;
            for (int m = static_cast<int>(kMaxOrder); m >= 1; --m) {
                const std::size_t index = static_cast<std::size_t>(m - 1);
                const float k = m <= order ? reflection[index] : 0.0f;
                const float backward_previous = state[index];
                const float forward_new = forward - k * backward_previous;
                state[static_cast<std::size_t>(m)] = backward_previous + k * forward_new;
                forward = forward_new;
            }
            state[0] = forward;
            return forward;
        };

        // 噪声链路的脉冲响应能量（含去加重与隔直），用于噪声激励定标
        std::array<float, kMaxOrder + 1> impulse_state{};
        float impulse_deemphasis = 0.0f;
        float impulse_dc_input = 0.0f;
        float impulse_dc_output = 0.0f;
        float noise_filter_energy = 0.0f;
        for (std::size_t n = 0; n < kImpulseLength; ++n) {
            const float y = lattice_tick(n == 0 ? 1.0f : 0.0f, impulse_state);
            const float deemphasized = y + kPreEmphasis * impulse_deemphasis;
            impulse_deemphasis = deemphasized;
            const float dc_blocked = deemphasized - impulse_dc_input + kDcBlockPole * impulse_dc_output;
            impulse_dc_input = deemphasized;
            impulse_dc_output = dc_blocked;
            noise_filter_energy += dc_blocked * dc_blocked;
        }
        noise_filter_energy = std::max(noise_filter_energy, 1e-6f);

        // 用四个周期的稳态响应估计单位幅度声门脉冲串（含倾斜补偿）的最终输出 RMS
        const std::size_t period = std::max<std::size_t>(1, static_cast<std::size_t>(std::lround(fs_ / f0)));
        std::array<float, kMaxOrder + 1> pulse_state{};
        float pulse_previous = 0.0f;
        float pulse_deemphasis = 0.0f;
        float pulse_dc_input = 0.0f;
        float pulse_dc_output = 0.0f;
        float pulse_energy = 0.0f;
        for (int cycle = 0; cycle < 4; ++cycle) {
            for (std::size_t n = 0; n < period; ++n) {
                const float phase = std::fmod(static_cast<float>(cycle * period + n) * f0 / fs_, 1.0f);
                const float raw_pulse = RosenbergDerivative(phase);
                const float shaped_pulse = raw_pulse - kPulseTilt * pulse_previous;
                pulse_previous = raw_pulse;
                const float y = lattice_tick(shaped_pulse, pulse_state);
                const float deemphasized = y + kPreEmphasis * pulse_deemphasis;
                pulse_deemphasis = deemphasized;
                const float dc_blocked = deemphasized - pulse_dc_input + kDcBlockPole * pulse_dc_output;
                pulse_dc_input = deemphasized;
                pulse_dc_output = dc_blocked;
                if (cycle == 3) {
                    pulse_energy += dc_blocked * dc_blocked;
                }
            }
        }
        const float pulse_rms = std::sqrt(pulse_energy / static_cast<float>(period));

        pulse_amp_target_ = frame_rms * pulse_share / std::max(pulse_rms, 1e-6f);
        noise_amp_target_ = frame_rms * std::sqrt(std::max(1.0f - pulse_share * pulse_share, 0.0f))
                            * std::sqrt(3.0f / noise_filter_energy);
        f0_target_ = f0;

        if (!has_model_) {
            k_current_ = k_target_;
            k_step_.fill(0.0f);
            f0_current_ = f0_target_;
            pulse_amp_current_ = pulse_amp_target_;
            noise_amp_current_ = noise_amp_target_;
            f0_step_ = 0.0f;
            pulse_amp_step_ = 0.0f;
            noise_amp_step_ = 0.0f;
            has_model_ = true;
        }
        else {
            constexpr float kInvHop = 1.0f / static_cast<float>(kHop);
            for (std::size_t i = 0; i < kMaxOrder; ++i) {
                k_step_[i] = (k_target_[i] - k_current_[i]) * kInvHop;
            }
            f0_step_ = (f0_target_ - f0_current_) * kInvHop;
            pulse_amp_step_ = (pulse_amp_target_ - pulse_amp_current_) * kInvHop;
            noise_amp_step_ = (noise_amp_target_ - noise_amp_current_) * kInvHop;
        }

        stats_.f0_hz = voicing > 0.0f ? f0 : 0.0f;
        stats_.voicing = voicing;
    }

    /// @brief 合成一个采样并推进参数插值
    /// @return 去加重后的合成采样
    float Synthesize() noexcept {
        phase_ += f0_current_ / fs_;
        if (phase_ >= 1.0f) {
            phase_ -= std::floor(phase_);
        }
        const float raw_pulse = RosenbergDerivative(phase_);
        const float pulse = raw_pulse - kPulseTilt * pulse_previous_;
        pulse_previous_ = raw_pulse;

        const float excitation = pulse_amp_current_ * pulse + noise_amp_current_ * noise_.Next();

        // 全极点格型合成：注入第 p 级前向误差，反向递推到第零级
        float forward = excitation;
        for (int m = static_cast<int>(kMaxOrder); m >= 1; --m) {
            const std::size_t index = static_cast<std::size_t>(m - 1);
            const float backward_previous = lattice_[index];
            const float forward_new = forward - k_current_[index] * backward_previous;
            lattice_[static_cast<std::size_t>(m)] = backward_previous + k_current_[index] * forward_new;
            forward = forward_new;
        }
        lattice_[0] = forward;

        const float y = forward + kPreEmphasis * deemphasis_;
        deemphasis_ = y;

        // 去加重对直流有二十倍增益，单极性脉冲串会留下明显偏置
        const float dc_blocked = y - dc_input_ + kDcBlockPole * dc_output_;
        dc_input_ = y;
        dc_output_ = dc_blocked;

        for (std::size_t i = 0; i < kMaxOrder; ++i) {
            k_current_[i] += k_step_[i];
        }
        f0_current_ += f0_step_;
        pulse_amp_current_ += pulse_amp_step_;
        noise_amp_current_ += noise_amp_step_;
        return dc_blocked;
    }

    /// @brief Rosenberg 声门流量导数，相位范围为零到一
    /// @note 开启段占 40% 周期，关闭段占 16%，其余为零；调用方再做一阶微分补偿声门源频谱倾斜
    static float RosenbergDerivative(float phase) noexcept {
        if (phase < kOpenPhase) {
            const float u = phase / kOpenPhase;
            return (6.0f * u - 6.0f * u * u) / kOpenPhase;
        }
        if (phase < kOpenPhase + kClosePhase) {
            const float u = (phase - kOpenPhase) / kClosePhase;
            return -2.0f * u / kClosePhase;
        }
        return 0.0f;
    }

    static constexpr float kPreEmphasis = 0.95f;
    static constexpr float kDcBlockPole = 0.995f;
    /// @brief 一阶微分，把声门源的频谱倾斜抵消到接近平坦
    static constexpr float kPulseTilt = 1.0f;
    static constexpr float kOpenPhase = 0.4f;
    static constexpr float kClosePhase = 0.16f;
    static constexpr float kVoiceRampRange = 0.25f;

    qwqdsp_pitch::Yin yin_;
    qwqdsp_adaptive::BurgLP burg_;
    qwqdsp_oscillator::WhiteNoise noise_;

    std::array<float, kPitchFrame> ring_{};
    std::array<float, kLpcFrame> lpc_block_{};
    std::array<float, kPitchFrame> pitch_block_{};
    std::array<float, kLpcFrame> window_{};
    std::array<float, kMaxOrder> k_target_{};
    std::array<float, kMaxOrder> k_current_{};
    std::array<float, kMaxOrder> k_step_{};
    std::array<float, kMaxOrder + 1> lattice_{};

    std::size_t ring_pos_{};
    std::size_t samples_seen_{};
    std::size_t hop_counter_{};
    float fs_{48000.0f};
    float phase_{};
    float deemphasis_{};
    float dc_input_{};
    float dc_output_{};
    float pulse_previous_{};
    float f0_current_{200.0f};
    float f0_target_{200.0f};
    float f0_step_{};
    float pulse_amp_current_{};
    float pulse_amp_target_{};
    float pulse_amp_step_{};
    float noise_amp_current_{};
    float noise_amp_target_{};
    float noise_amp_step_{};
    float noise_mix_{0.0f};
    float voicing_threshold_{0.45f};
    bool has_model_{};
    VocoderStats stats_{};
};

} // namespace lpc_vocoder::detail
