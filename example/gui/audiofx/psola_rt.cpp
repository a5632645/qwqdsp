#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <iostream>
#include <numbers>
#include <span>
#include <string>

#include "miniaudio.h"
#include "raylib.h"
#include "slider.hpp"

#include "qwqdsp/pitch/yin.hpp"

namespace psola_rt2 {

static constexpr float kSampleRate = 48000.0f;
static constexpr int kWindowWidth = 560;
static constexpr int kWindowHeight = 300;
static constexpr const char* kWindowTitle = "PSOLA RT2";
static constexpr float kMinPitchSemitones = -12.0f;
static constexpr float kMaxPitchSemitones = 12.0f;
static constexpr float kMinFormantSemitones = -8.0f;
static constexpr float kMaxFormantSemitones = 8.0f;

namespace detail {

// ------------------------------------------------------------
// 公共数值工具
// ------------------------------------------------------------

struct PsolaMath {
    /**
     * @brief 将半音偏移量换算为频率倍率
     * @param[in] semitones 半音偏移量
     * @return 对应的频率倍率
     */
    static float semitonesToRatio(float semitones) {
        return std::exp2(semitones / 12.0f);
    }

    /**
     * @brief 计算以脉冲标记为中心的 Hann 窗权重
     * @param[in] offset 相对脉冲标记的采样偏移
     * @param[in] half_width 窗函数的半宽
     * @return 范围为零到一的窗权重
     */
    static float hannWindow(int offset, int half_width) {
        if (half_width <= 0 || std::abs(offset) >= half_width) {
            return 0.0f;
        }

        const float phase = static_cast<float>(offset) / static_cast<float>(half_width);
        return 0.5f + 0.5f * std::cos(std::numbers::pi_v<float> * phase);
    }
};

// ------------------------------------------------------------
// 基频检测
// ------------------------------------------------------------

struct PitchEstimate {
    float period_samples = kSampleRate / 200.0f;
    float pitch_hz = 0.0f;
    float confidence = 0.0f;
    bool voiced = false;
};

struct PitchDetector {
    static constexpr std::size_t kFrameSize = 2048;
    static constexpr std::size_t kHopSize = 256;
    static constexpr float kMinPitchHz = 70.0f;
    static constexpr float kMaxPitchHz = 500.0f;
    static constexpr float kYinThreshold = 0.25f;
    static constexpr float kRmsThreshold = 0.0025f;
    // static constexpr float kConfidenceThreshold = 0.62f;
    static constexpr float kConfidenceThreshold = 1 - kYinThreshold;

    /**
     * @brief 初始化 YIN 检测器及分析缓冲
     * @param[in] sample_rate 音频采样率，单位为 Hz
     */
    void init(float sample_rate) {
        sample_rate_ = sample_rate;
        yin_.Init(sample_rate, static_cast<int>(kFrameSize));
        yin_.SetMinPitch(kMinPitchHz);
        yin_.SetMaxPitch(kMaxPitchHz);
        yin_.SetThreshold(kYinThreshold);

        input_ring_.fill(0.0f);
        analysis_frame_.fill(0.0f);
        write_position_ = 0;
        samples_received_ = 0;
        samples_until_analysis_ = kHopSize;
        estimate_ = PitchEstimate{};
    }

    /**
     * @brief 接收一个输入采样并按固定步长更新基频估计
     * @param[in] input_sample 当前单声道输入采样
     * @return 本次调用完成新一帧分析时返回 true
     */
    bool processSample(float input_sample) {
        input_ring_[write_position_] = input_sample;
        write_position_ = (write_position_ + 1) % kFrameSize;
        samples_received_ = std::min(samples_received_ + 1, kFrameSize);

        if (samples_received_ < kFrameSize) {
            return false;
        }

        if (--samples_until_analysis_ > 0) {
            return false;
        }
        samples_until_analysis_ = kHopSize;

        float squared_sum = 0.0f;
        for (std::size_t i = 0; i < kFrameSize; ++i) {
            const float sample = input_ring_[(write_position_ + i) % kFrameSize];
            analysis_frame_[i] = sample;
            squared_sum += sample * sample;
        }

        yin_.Process(std::span<const float>{analysis_frame_});
        const qwqdsp_pitch::Pitch pitch = yin_.GetPitch();
        const float rms = std::sqrt(squared_sum / static_cast<float>(kFrameSize));
        const float confidence = std::clamp(1.0f - pitch.non_period_ratio, 0.0f, 1.0f);
        const bool pitch_in_range = pitch.pitch_hz >= kMinPitchHz && pitch.pitch_hz <= kMaxPitchHz;

        estimate_.pitch_hz = pitch_in_range ? pitch.pitch_hz : 0.0f;
        estimate_.confidence = confidence;
        estimate_.voiced = pitch_in_range && rms >= kRmsThreshold && confidence >= kConfidenceThreshold;
        if (estimate_.voiced) {
            estimate_.period_samples = sample_rate_ / pitch.pitch_hz;
        }
        return true;
    }

    /**
     * @brief 查询检测器是否已接收完整分析帧
     * @return 可以提供稳定分析结果时返回 true
     */
    bool ready() const noexcept {
        return samples_received_ >= kFrameSize;
    }

    /**
     * @brief 获取最近一次基频估计
     * @return 当前基频周期、频率、置信度和清浊音状态
     */
    PitchEstimate estimate() const noexcept {
        return estimate_;
    }
private:
    qwqdsp_pitch::Yin yin_;
    std::array<float, kFrameSize> input_ring_{};
    std::array<float, kFrameSize> analysis_frame_{};
    std::size_t write_position_ = 0;
    std::size_t samples_received_ = 0;
    std::size_t samples_until_analysis_ = kHopSize;
    float sample_rate_ = kSampleRate;
    PitchEstimate estimate_{};
};

// ------------------------------------------------------------
// 实时 TD-PSOLA
// ------------------------------------------------------------

struct RealtimePsola {
    static constexpr std::size_t kRingSize = 32768;
    static constexpr std::size_t kRingMask = kRingSize - 1;
    static constexpr std::size_t kVoicedHangoverUpdates = 3;
    static constexpr float kVoicedReleaseConfidence = 0.45f;
    static constexpr float kPeriodSmoothing = 0.18f;
    static constexpr float kUnvoicedMaxJitterRatio = 0.30f;
    static constexpr float kPitchMarkSearchHalfPeriods = 0.55f;
    static constexpr float kSourceGrainHalfPeriods = 1.15f;
    static constexpr float kTargetOverlapHalfPeriods = 0.55f;
    static constexpr float kMinimumOlaNormalization = 1.0f;
    static constexpr std::uint64_t kLatency = 6144;
    static constexpr double kScheduleLookahead = 3200.0;
    static constexpr float kMinPeriod = kSampleRate / PitchDetector::kMaxPitchHz;
    static constexpr float kMaxPeriod = kSampleRate / PitchDetector::kMinPitchHz;
    static_assert((kRingSize & (kRingSize - 1)) == 0, "环形缓冲长度必须为二的幂");

    /**
     * @brief 初始化 PSOLA 状态和固定容量环形缓冲
     * @param[in] sample_rate 音频采样率，单位为 Hz
     */
    void init(float sample_rate) {
        pitch_detector_.init(sample_rate);
        input_ring_.fill(0.0f);
        overlap_ring_.fill(0.0f);
        weight_ring_.fill(0.0f);
        voiced_ring_.fill(0.0f);

        input_position_ = 0;
        next_target_mark_ = 0.0;
        smoothed_period_ = sample_rate / 200.0f;
        voiced_target_ = 0.0f;
        voiced_hangover_updates_ = 0;
        mark_polarity_ = 1.0f;
        noise_state_ = 0x9E3779B9U;
        scheduler_started_ = false;
        polarity_initialized_ = false;
    }

    /**
     * @brief 处理一段单声道实时音频
     * @param[in] input 输入采样缓冲
     * @param[out] output 输出采样缓冲
     * @param[in] frame_count 缓冲中的采样数
     * @param[in] pitch_ratio 目标音高倍率
     * @param[in] formant_ratio 目标共振峰倍率
     */
    void process(const float* input, float* output, std::size_t frame_count, float pitch_ratio, float formant_ratio) {
        const float safe_pitch_ratio = std::clamp(pitch_ratio, 0.5f, 2.0f);
        const float safe_formant_ratio = std::clamp(formant_ratio, 0.5f, 2.0f);

        for (std::size_t i = 0; i < frame_count; ++i) {
            output[i] = processSample(input[i], safe_pitch_ratio, safe_formant_ratio);
        }
    }

    /**
     * @brief 获取最近检测到的基频
     * @return 有声输入的基频，无法可靠检测时返回零
     */
    float detectedPitchHz() const noexcept {
        return pitch_detector_.estimate().voiced ? pitch_detector_.estimate().pitch_hz : 0.0f;
    }

    /**
     * @brief 获取最近基频检测的置信度
     * @return 范围为零到一的置信度
     */
    float pitchConfidence() const noexcept {
        return pitch_detector_.estimate().confidence;
    }
private:
    /**
     * @brief 将绝对采样位置映射到环形缓冲索引
     * @param[in] position 绝对采样位置
     * @return 环形缓冲中的数组索引
     */
    static std::size_t ringIndex(std::uint64_t position) noexcept {
        return static_cast<std::size_t>(position) & kRingMask;
    }

    /**
     * @brief 按清浊音类型处理一个采样并输出延迟结果
     * @param[in] input_sample 当前输入采样
     * @param[in] pitch_ratio 目标音高倍率
     * @param[in] formant_ratio 目标共振峰倍率
     * @return 当前输出采样
     * @note 低权重 OLA 区域不做满幅归一化，防止降调时恢复被跳过的源脉冲
     */
    float processSample(float input_sample, float pitch_ratio, float formant_ratio) {
        const std::size_t current_index = ringIndex(input_position_);
        input_ring_[current_index] = input_sample;

        if (pitch_detector_.processSample(input_sample)) {
            updatePitchEstimate(pitch_detector_.estimate());
        }
        voiced_ring_[current_index] = voiced_target_;

        if (!scheduler_started_ && pitch_detector_.ready()
            && static_cast<double>(input_position_) > kScheduleLookahead) {
            next_target_mark_ = static_cast<double>(input_position_) - kScheduleLookahead;
            scheduler_started_ = true;
        }
        scheduleAvailableGrains(pitch_ratio, formant_ratio);

        float dry_sample = 0.0f;
        if (input_position_ >= kLatency) {
            const std::size_t delayed_index = ringIndex(input_position_ - kLatency);
            dry_sample = input_ring_[delayed_index];
        }

        float wet_sample = dry_sample;
        if (weight_ring_[current_index] > 0.02f) {
            const float normalization = std::max(weight_ring_[current_index], kMinimumOlaNormalization);
            wet_sample = overlap_ring_[current_index] / normalization;
        }

        overlap_ring_[current_index] = 0.0f;
        weight_ring_[current_index] = 0.0f;
        ++input_position_;
        return std::clamp(wet_sample, -1.0f, 1.0f);
    }

    /**
     * @brief 平滑基频周期并以迟滞方式更新二值浊音门控
     * @param[in] estimate 最新 YIN 基频估计
     * @note 置信度只参与清浊音状态判定，不作为干湿混合比例
     */
    void updatePitchEstimate(const PitchEstimate& estimate) {
        if (!estimate.voiced) {
            const bool weakly_voiced =
                voiced_target_ > 0.0f && estimate.pitch_hz > 0.0f && estimate.confidence >= kVoicedReleaseConfidence;
            if (weakly_voiced) {
                voiced_hangover_updates_ = kVoicedHangoverUpdates;
                return;
            }

            if (voiced_hangover_updates_ > 0) {
                --voiced_hangover_updates_;
                voiced_target_ = 1.0f;
            }
            else {
                voiced_target_ = 0.0f;
            }
            return;
        }

        voiced_target_ = 1.0f;
        voiced_hangover_updates_ = kVoicedHangoverUpdates;

        float candidate_period = std::clamp(estimate.period_samples, kMinPeriod, kMaxPeriod);
        const float period_ratio = candidate_period / smoothed_period_;
        if (period_ratio > 1.80f && period_ratio < 2.20f) {
            candidate_period *= 0.5f;
        }
        else if (period_ratio > 0.45f && period_ratio < 0.56f) {
            candidate_period *= 2.0f;
        }

        smoothed_period_ += kPeriodSmoothing * (candidate_period - smoothed_period_);
        smoothed_period_ = std::clamp(smoothed_period_, kMinPeriod, kMaxPeriod);
    }

    /**
     * @brief 按清浊音类型安排具备输入前瞻的 PSOLA 颗粒
     * @param[in] pitch_ratio 目标音高倍率
     * @param[in] formant_ratio 目标共振峰倍率
     */
    void scheduleAvailableGrains(float pitch_ratio, float formant_ratio) {
        if (!scheduler_started_) {
            return;
        }

        const double latest_schedulable = static_cast<double>(input_position_) - kScheduleLookahead;
        while (next_target_mark_ <= latest_schedulable) {
            if (readAlignedVoicing(next_target_mark_) >= 0.5f) {
                const float target_period = smoothed_period_ / pitch_ratio;
                addVoicedGrain(next_target_mark_, smoothed_period_, target_period, formant_ratio);
                next_target_mark_ += static_cast<double>(target_period);
            }
            else {
                addUnvoicedGrain(next_target_mark_, smoothed_period_, formant_ratio);
                next_target_mark_ += static_cast<double>(smoothed_period_);
            }
        }
    }

    /**
     * @brief 读取与输入帧中心对齐的浊音状态
     * @param[in] target_position 输入时间轴上的目标位置
     * @return 浊音返回一，清音或越界位置返回零
     */
    float readAlignedVoicing(double target_position) const noexcept {
        const std::int64_t target = static_cast<std::int64_t>(target_position);
        const std::int64_t aligned_position = target + static_cast<std::int64_t>(PitchDetector::kFrameSize / 2);
        if (aligned_position < 0 || static_cast<std::uint64_t>(aligned_position) > input_position_) {
            return 0.0f;
        }
        return voiced_ring_[ringIndex(static_cast<std::uint64_t>(aligned_position))];
    }

    /**
     * @brief 定位目标时刻附近同极性的输入脉冲
     * @param[in] target_position 目标时刻对应的输入位置
     * @param[in] source_period 当前输入基音周期
     * @return 对齐后的整数输入脉冲位置
     */
    std::int64_t findSourceMark(double target_position, float source_period) {
        const std::int64_t target = static_cast<std::int64_t>(std::llround(target_position));
        const int search_half = std::max(1, static_cast<int>(std::lround(source_period * kPitchMarkSearchHalfPeriods)));
        std::int64_t best_position = target;

        if (!polarity_initialized_) {
            float largest_magnitude = 0.0f;
            for (int offset = -search_half; offset <= search_half; ++offset) {
                const float sample = readInput(target + offset);
                if (std::abs(sample) > largest_magnitude) {
                    largest_magnitude = std::abs(sample);
                    best_position = target + offset;
                    mark_polarity_ = sample < 0.0f ? -1.0f : 1.0f;
                }
            }
            polarity_initialized_ = true;
            return best_position;
        }

        float best_score = -1.0f;
        for (int offset = -search_half; offset <= search_half; ++offset) {
            const float score = mark_polarity_ * readInput(target + offset);
            if (score > best_score) {
                best_score = score;
                best_position = target + offset;
            }
        }
        return best_position;
    }

    /**
     * @brief 将一个经共振峰重采样的脉冲片段叠加到输出环形缓冲
     * @param[in] target_position 合成标记在输入时间轴上的位置
     * @param[in] source_period 输入基音周期
     * @param[in] target_period 输出基音周期
     * @param[in] formant_ratio 共振峰频率倍率
     * @note 窗宽同时限制源域脉冲数量并保证目标标记之间具有稳定重叠
     */
    void addVoicedGrain(double target_position, float source_period, float target_period, float formant_ratio) {
        const std::int64_t source_mark = findSourceMark(target_position, source_period);
        const std::uint64_t output_mark =
            static_cast<std::uint64_t>(std::llround(target_position + static_cast<double>(kLatency)));
        const float source_limited_half_width = kSourceGrainHalfPeriods * source_period / formant_ratio;
        const float overlap_limited_half_width = kTargetOverlapHalfPeriods * target_period;
        const int half_width =
            std::max(2, static_cast<int>(std::ceil(std::max(source_limited_half_width, overlap_limited_half_width))));

        overlapGrain(source_mark, output_mark, half_width, formant_ratio);
    }

    /**
     * @brief 使用最近浊音周期合成非周期清音颗粒
     * @param[in] target_position 合成标记在输入时间轴上的位置
     * @param[in] held_period 最近可靠浊音的基音周期
     * @param[in] formant_ratio 共振峰频率倍率
     * @note Formant 非零时随机扰动分析标记，避免稳定的周期性相消
     */
    void addUnvoicedGrain(double target_position, float held_period, float formant_ratio) {
        const float formant_octaves = std::abs(std::log2(formant_ratio));
        const float jitter_ratio = std::min(kUnvoicedMaxJitterRatio, formant_octaves * 0.45f);
        const double source_position =
            target_position + static_cast<double>(nextNoiseSample() * jitter_ratio * held_period);
        const std::int64_t source_mark = static_cast<std::int64_t>(std::llround(source_position));
        const std::uint64_t output_mark =
            static_cast<std::uint64_t>(std::llround(target_position + static_cast<double>(kLatency)));
        const int half_width = std::max(2, static_cast<int>(std::ceil(held_period)));

        overlapGrain(source_mark, output_mark, half_width, formant_ratio);
    }

    /**
     * @brief 生成用于清音分析标记扰动的双极性伪随机数
     * @return 范围约为负一到一的伪随机值
     */
    float nextNoiseSample() noexcept {
        noise_state_ ^= noise_state_ << 13U;
        noise_state_ ^= noise_state_ >> 17U;
        noise_state_ ^= noise_state_ << 5U;
        constexpr float kScale = 1.0f / 16777215.0f;
        const float normalized = static_cast<float>(noise_state_ >> 8U) * kScale;
        return 2.0f * normalized - 1.0f;
    }

    /**
     * @brief 对颗粒内部重采样并执行带权重归一化的重叠相加
     * @param[in] source_mark 输入颗粒中心的绝对采样位置
     * @param[in] output_mark 输出颗粒中心的绝对采样位置
     * @param[in] half_width Hann 窗半宽
     * @param[in] formant_ratio 共振峰频率倍率
     */
    void overlapGrain(std::int64_t source_mark, std::uint64_t output_mark, int half_width, float formant_ratio) {
        for (int offset = -half_width; offset <= half_width; ++offset) {
            const float window = PsolaMath::hannWindow(offset, half_width);
            if (window <= 0.0f) {
                continue;
            }

            const double source_position =
                static_cast<double>(source_mark) + static_cast<double>(offset) * static_cast<double>(formant_ratio);
            const float sample = readInputLinear(source_position);
            const std::uint64_t output_position =
                static_cast<std::uint64_t>(static_cast<std::int64_t>(output_mark) + offset);
            const std::size_t output_index = ringIndex(output_position);
            overlap_ring_[output_index] += sample * window;
            weight_ring_[output_index] += window;
        }
    }

    /**
     * @brief 读取指定绝对位置的输入采样
     * @param[in] position 有符号绝对采样位置
     * @return 有效历史范围内的采样，越界时返回零
     */
    float readInput(std::int64_t position) const noexcept {
        if (position < 0 || static_cast<std::uint64_t>(position) > input_position_) {
            return 0.0f;
        }
        return input_ring_[ringIndex(static_cast<std::uint64_t>(position))];
    }

    /**
     * @brief 以线性插值读取非整数位置的输入采样
     * @param[in] position 浮点绝对采样位置
     * @return 插值后的输入采样
     */
    float readInputLinear(double position) const {
        const std::int64_t first_position = static_cast<std::int64_t>(std::floor(position));
        const float fraction = static_cast<float>(position - static_cast<double>(first_position));
        return std::lerp(readInput(first_position), readInput(first_position + 1), fraction);
    }

    PitchDetector pitch_detector_;
    std::array<float, kRingSize> input_ring_{};
    std::array<float, kRingSize> overlap_ring_{};
    std::array<float, kRingSize> weight_ring_{};
    std::array<float, kRingSize> voiced_ring_{};

    std::uint64_t input_position_ = 0;
    double next_target_mark_ = 0.0;
    float smoothed_period_ = kSampleRate / 200.0f;
    float voiced_target_ = 0.0f;
    std::size_t voiced_hangover_updates_ = 0;
    std::uint32_t noise_state_ = 0x9E3779B9U;
    float mark_polarity_ = 1.0f;
    bool scheduler_started_ = false;
    bool polarity_initialized_ = false;
};

} // namespace detail

// ------------------------------------------------------------
// 音频与界面
// ------------------------------------------------------------

struct PsolaDemo {
    /**
     * @brief 初始化 DSP 与界面旋钮
     * @param[in] sample_rate 音频设备采样率，单位为 Hz
     */
    void init(float sample_rate) {
        processor_.init(sample_rate);
        createControls();
    }

    /**
     * @brief 创建相互独立的音高和共振峰控制旋钮
     */
    void createControls() {
        pitch_knob_.set_bound(118, 104, 130, 132)
            .set_title("Pitch")
            .set_range(kMinPitchSemitones, kMaxPitchSemitones, 0.1f, 0.0f)
            .set_value(0.0f)
            .set_sensitivity(2)
            .set_name_font_size(17)
            .set_number_font_size(13)
            .set_fore_color(Color{235, 235, 238, 255})
            .set_bg_color(Color{24, 26, 29, 255});
        pitch_knob_.value_to_text_function = [](float value) { return std::format("{:+.1f} st", value); };
        pitch_knob_.on_value_change = [this](float value) { pitch_semitones_.store(value, std::memory_order_relaxed); };

        formant_knob_.set_bound(312, 104, 130, 132)
            .set_title("Formant")
            .set_range(kMinFormantSemitones, kMaxFormantSemitones, 0.1f, 0.0f)
            .set_value(0.0f)
            .set_sensitivity(2)
            .set_name_font_size(17)
            .set_number_font_size(13)
            .set_fore_color(Color{235, 235, 238, 255})
            .set_bg_color(Color{24, 26, 29, 255});
        formant_knob_.value_to_text_function = [](float value) { return std::format("{:+.1f} st", value); };
        formant_knob_.on_value_change = [this](float value) {
            formant_semitones_.store(value, std::memory_order_relaxed);
        };
    }

    /**
     * @brief 处理音频设备提供的一段单声道采样
     * @param[in] input 单声道输入缓冲
     * @param[out] output 单声道输出缓冲
     * @param[in] frame_count 缓冲中的采样数
     */
    void processAudio(const float* input, float* output, std::size_t frame_count) {
        const float pitch_ratio = detail::PsolaMath::semitonesToRatio(pitch_semitones_.load(std::memory_order_relaxed));
        const float formant_ratio =
            detail::PsolaMath::semitonesToRatio(formant_semitones_.load(std::memory_order_relaxed));

        processor_.process(input, output, frame_count, pitch_ratio, formant_ratio);

        float input_peak = 0.0f;
        float output_peak = 0.0f;
        for (std::size_t i = 0; i < frame_count; ++i) {
            input_peak = std::max(input_peak, std::abs(input[i]));
            output_peak = std::max(output_peak, std::abs(output[i]));
        }

        input_peak_.store(std::max(input_peak, input_peak_.load(std::memory_order_relaxed) * 0.82f),
                          std::memory_order_relaxed);
        output_peak_.store(std::max(output_peak, output_peak_.load(std::memory_order_relaxed) * 0.82f),
                           std::memory_order_relaxed);
        detected_pitch_hz_.store(processor_.detectedPitchHz(), std::memory_order_relaxed);
        pitch_confidence_.store(processor_.pitchConfidence(), std::memory_order_relaxed);
    }

    /**
     * @brief 绘制参数控制、检测状态与输入输出电平
     * @param[in] audio_running 音频设备是否成功启动
     */
    void draw(bool audio_running) {
        static constexpr Color kBackground{18, 20, 23, 255};
        static constexpr Color kPanel{24, 26, 29, 255};
        static constexpr Color kText{232, 233, 235, 255};
        static constexpr Color kMuted{126, 132, 139, 255};
        static constexpr Color kAccent{76, 201, 151, 255};
        static constexpr Color kWarning{225, 172, 74, 255};

        ClearBackground(kBackground);
        DrawRectangle(0, 0, kWindowWidth, 72, kPanel);
        DrawRectangle(0, 71, kWindowWidth, 1, Color{48, 52, 57, 255});
        DrawText("TD-PSOLA", 24, 18, 25, kText);
        DrawText("REAL-TIME VOICE PROCESSOR", 25, 48, 10, kMuted);

        const Color status_color = audio_running ? kAccent : kWarning;
        DrawCircle(445, 30, 5.0f, status_color);
        DrawText(audio_running ? "AUDIO ONLINE" : "AUDIO OFFLINE", 458, 24, 12, status_color);

        DrawRectangle(279, 105, 1, 128, Color{48, 52, 57, 255});
        pitch_knob_.display();
        formant_knob_.display();

        const float detected_pitch = detected_pitch_hz_.load(std::memory_order_relaxed);
        const float confidence = pitch_confidence_.load(std::memory_order_relaxed);
        const std::string pitch_text =
            detected_pitch > 0.0f ? std::format("F0  {:6.1f} Hz", detected_pitch) : std::string{"F0       -- Hz"};
        DrawText(pitch_text.c_str(), 24, 264, 13, detected_pitch > 0.0f ? kText : kMuted);

        drawMeter(184, 266, 116, input_peak_.load(std::memory_order_relaxed), "IN", kAccent);
        drawMeter(326, 266, 116, output_peak_.load(std::memory_order_relaxed), "OUT", Color{91, 160, 225, 255});

        const int confidence_width = static_cast<int>(72.0f * std::clamp(confidence, 0.0f, 1.0f));
        DrawText("TRACK", 462, 264, 10, kMuted);
        DrawRectangle(462, 280, 72, 3, Color{45, 48, 52, 255});
        DrawRectangle(462, 280, confidence_width, 3, status_color);

        input_peak_.store(input_peak_.load(std::memory_order_relaxed) * 0.96f, std::memory_order_relaxed);
        output_peak_.store(output_peak_.load(std::memory_order_relaxed) * 0.96f, std::memory_order_relaxed);
    }

    /**
     * @brief 响应 miniaudio 的全双工单声道数据请求
     * @param[in] device 发起回调的音频设备
     * @param[out] output_buffer 单声道浮点输出缓冲
     * @param[in] input_buffer 单声道浮点输入缓冲
     * @param[in] frame_count 本次请求的音频帧数
     */
    static void audioCallback(ma_device* device, void* output_buffer, const void* input_buffer, ma_uint32 frame_count) {
        auto* output = static_cast<float*>(output_buffer);
        const auto* input = static_cast<const float*>(input_buffer);
        auto* demo = static_cast<PsolaDemo*>(device->pUserData);

        if (output == nullptr) {
            return;
        }
        if (input == nullptr || demo == nullptr) {
            for (ma_uint32 i = 0; i < frame_count; ++i) {
                output[i] = 0.0f;
            }
            return;
        }

        demo->processAudio(input, output, static_cast<std::size_t>(frame_count));
    }
private:
    /**
     * @brief 绘制一个紧凑的水平峰值电平表
     * @param[in] x 电平表左侧坐标
     * @param[in] y 电平表顶部坐标
     * @param[in] width 电平表宽度
     * @param[in] value 线性峰值幅度
     * @param[in] label 电平表标签
     * @param[in] color 填充颜色
     */
    static void drawMeter(int x, int y, int width, float value, const char* label, Color color) {
        const float normalized = std::clamp(value, 0.0f, 1.0f);
        const int bar_width = width - 24;
        const int fill_width = static_cast<int>(static_cast<float>(bar_width) * normalized);
        DrawText(label, x, y - 2, 10, Color{126, 132, 139, 255});
        DrawRectangle(x + 24, y, bar_width, 7, Color{45, 48, 52, 255});
        DrawRectangle(x + 24, y, fill_width, 7, color);
    }

    detail::RealtimePsola processor_;
    Knob pitch_knob_;
    Knob formant_knob_;
    std::atomic<float> pitch_semitones_{0.0f};
    std::atomic<float> formant_semitones_{0.0f};
    std::atomic<float> detected_pitch_hz_{0.0f};
    std::atomic<float> pitch_confidence_{0.0f};
    std::atomic<float> input_peak_{0.0f};
    std::atomic<float> output_peak_{0.0f};
};

} // namespace psola_rt2

/**
 * @brief 启动实时单声道 PSOLA 演示程序
 * @return 正常退出返回零，音频设备初始化失败时仍保留界面并返回零
 */
int main() {
    using namespace psola_rt2;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, kWindowTitle);
    SetTargetFPS(60);

    static PsolaDemo demo;
    demo.init(kSampleRate);

    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.playback.format = ma_format_f32;
    config.playback.channels = 1;
    config.sampleRate = static_cast<ma_uint32>(kSampleRate);
    config.dataCallback = PsolaDemo::audioCallback;
    config.pUserData = &demo;
    config.periodSizeInMilliseconds = 5;

    ma_device device{};
    const ma_result init_result = ma_device_init(nullptr, &config, &device);
    ma_result audio_result = init_result;
    bool audio_running = false;
    if (init_result == MA_SUCCESS) {
        audio_result = ma_device_start(&device);
        audio_running = audio_result == MA_SUCCESS;
    }

    if (!audio_running) {
        std::cout << std::format("miniaudio: device start failed ({})\n", static_cast<int>(audio_result));
    }

    while (!WindowShouldClose()) {
        BeginDrawing();
        demo.draw(audio_running);
        EndDrawing();
    }

    if (init_result == MA_SUCCESS) {
        if (audio_running) {
            ma_device_stop(&device);
        }
        ma_device_uninit(&device);
    }
    CloseWindow();
    return 0;
}
