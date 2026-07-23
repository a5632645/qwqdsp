#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>

#include <qwqdsp/pitch/yin.hpp>

#include "miniaudio.h"
#include "raylib.h"
#include "slider.hpp"

// ------------------------------------------------------------
// 常量
// ------------------------------------------------------------
static constexpr int kWindowWidth = 480;
static constexpr int kWindowHeight = 240;
static constexpr float kSampleRate = 48000.0f;

static constexpr float kMinF0 = 75.0f;
static constexpr float kMaxF0 = 500.0f;

static constexpr int kInputBufSize = 8192;     // 须为 2 的幂
static constexpr int kOutputBufSize = 32768;   // 须为 2 的幂
static constexpr int kAlgorithmLatency = 2048; // 合成输出和清音直通的公共延迟
static constexpr int kGrainMaxLen = 4096;
static constexpr int kPendingGrainCount = 8;

// YIN 音高检测参数
static constexpr int kYinFrameSize = 2048;
static constexpr int kYinInterval = static_cast<int>(0.01f * kSampleRate);
static constexpr float kYinThreshold = 0.18f;

// 窗函数 LUT
static constexpr int kWinTableSize = 2048;

// ------------------------------------------------------------
// HannWindow
// ------------------------------------------------------------
struct HannWindow {
    /**
     * @brief 初始化 Hann 窗查找表。
     */
    static void init() noexcept {
        for (int i = 0; i < kWinTableSize; ++i) {
            table_[i] = 0.5f
                      * (1.0f
                         - std::cos(2.0f * std::numbers::pi_v<float>
                                    * static_cast<float>(i) / static_cast<float>(kWinTableSize - 1)));
        }
    }

    /**
     * @brief 返回指定长度 Hann 窗中的一个采样值。
     * @param[in] index 窗内采样索引。
     * @param[in] length 窗的总采样数。
     * @return 范围为 [0, 1] 的窗权重。
     */
    static float value(int index, int length) noexcept {
        if (length <= 1) {
            return 1.0f;
        }

        const float pos =
            static_cast<float>(index) / static_cast<float>(length - 1) * static_cast<float>(kWinTableSize - 1);
        const int index0 = static_cast<int>(pos);
        const int index1 = std::min(index0 + 1, kWinTableSize - 1);
        const float fraction = pos - static_cast<float>(index0);
        return table_[index0] + fraction * (table_[index1] - table_[index0]);
    }
private:
    inline static std::array<float, kWinTableSize> table_{};
};

// ------------------------------------------------------------
// FormantGrain
// ------------------------------------------------------------
struct FormantGrain {
    std::array<float, kGrainMaxLen> samples{};
    int length = 0;
    int t0 = 0;
};

struct PendingGrain {
    std::int64_t epoch_time = 0;
    int t0 = 0;
    int radius = 0;
    bool active = false;
};

// ------------------------------------------------------------
// 实时 PSOLA
// ------------------------------------------------------------
/**
 * @brief 支持清浊音切换、完整 Grain 提取及归一化 OLA 的实时 TD-PSOLA。
 *
 * 浊音由 YIN 周期估计和残差信号中的周期峰值共同驱动。清音不参与
 * PSOLA 合成，而是通过同一输出时间轴上的延迟直通路径输出，并在
 * 清浊音边界进行平滑交叉淡化。
 */
class RealtimePsola {
public:
    /**
     * @brief 重置处理器并设置采样率。
     * @param[in] sample_rate 采样率，单位为 Hz。
     */
    void init(float sample_rate) {
        sample_rate_ = sample_rate;
        min_period_ = static_cast<int>(sample_rate_ / kMaxF0);
        max_period_ = static_cast<int>(sample_rate_ / kMinF0);
        default_period_ = (min_period_ + max_period_) / 2;

        std::fill(input_buf_.begin(), input_buf_.end(), 0.0f);
        std::fill(output_buf_.begin(), output_buf_.end(), 0.0f);
        std::fill(output_weight_buf_.begin(), output_weight_buf_.end(), 0.0f);
        std::fill(dry_buf_.begin(), dry_buf_.end(), 0.0f);
        std::fill(voicing_buf_.begin(), voicing_buf_.end(), 0.0f);
        pending_grains_.fill(PendingGrain{});

        latest_grain_ = FormantGrain{};
        synth_grain_ = FormantGrain{};
        input_wpos_ = 0;
        output_rpos_ = 0;
        input_time_ = 0;
        phasor_phase_ = 0.0f;
        synthesis_active_ = false;
        pitch_scale_smoothed_ = 1.0f;
        voiced_mix_ = 0.0f;

        yin_counter_ = 0;
        yin_period_ = 0.0f;
        yin_previous_period_ = 0.0f;
        yin_clarity_ = 0.0f;
        yin_voiced_ = false;
        voiced_candidate_count_ = 0;
        unvoiced_candidate_count_ = 0;
        yin_.Init(sample_rate_, kYinFrameSize);
        yin_.SetMinPitch(kMinF0);
        yin_.SetMaxPitch(kMaxF0);
        yin_.SetThreshold(kYinThreshold);

        previous_epoch_time_ = -1;
        previous_period_ = static_cast<float>(default_period_);
        gci_since_last_ = 0;
        gci_searching_ = false;
        gci_peak_ = 0.0f;
        gci_peak_time_ = 0;
        gci_peak_age_ = 0;
        residual_previous_input_ = 0.0f;
    }

    std::atomic<float> pitch_shift_semitones_{0.0f};
    std::atomic<float> formant_shift_{1.0f};
    std::atomic<float> pitch_scale_target_{1.0f};

    /**
     * @brief 处理一段单声道实时音频。
     * @param[in] input 输入采样缓冲区。
     * @param[out] output 输出采样缓冲区。
     * @param[in] num_samples 缓冲区中的采样数。
     */
    void process(const float* input, float* output, int num_samples) noexcept {
        const float pitch_target = pitch_scale_target_.load(std::memory_order_relaxed);
        const float formant_scale = std::clamp(formant_shift_.load(std::memory_order_relaxed), 0.05f, 4.0f);
        constexpr float kParameterSmooth = 0.001f;
        constexpr float kVoicingFade = 0.004f;

        for (int i = 0; i < num_samples; ++i) {
            const float x = input[i];
            pitch_scale_smoothed_ += kParameterSmooth * (pitch_target - pitch_scale_smoothed_);
            const float pitch_scale = std::max(pitch_scale_smoothed_, 0.05f);

            input_buf_[input_wpos_] = x;
            runAnalyzer(input_time_, x);

            const int synthesis_center =
                wrapIndex(static_cast<std::int64_t>(output_rpos_) + kAlgorithmLatency, kOutputBufSize);
            dry_buf_[synthesis_center] = x;
            voicing_buf_[synthesis_center] = yin_voiced_ ? 1.0f : 0.0f;

            const bool can_synthesize = yin_voiced_ && latest_grain_.length > 0;
            if (can_synthesize && !synthesis_active_) {
                phasor_phase_ = 1.0f;
            }

            if (can_synthesize) {
                const float t0 = static_cast<float>(latest_grain_.t0);
                phasor_phase_ += pitch_scale / std::max(t0, 1.0f);
                if (phasor_phase_ >= 1.0f) {
                    phasor_phase_ -= std::floor(phasor_phase_);
                    synth_grain_ = latest_grain_;
                    if (std::abs(formant_scale - 1.0f) > 1e-4f) {
                        synth_grain_ = resampleGrain(synth_grain_, formant_scale);
                    }
                    renderGrain(synth_grain_, synthesis_center);
                }
            }
            else {
                phasor_phase_ = 0.0f;
            }
            synthesis_active_ = can_synthesize;

            const float weight = output_weight_buf_[output_rpos_];
            const bool has_psola = weight > 0.02f;
            const float mix_target = has_psola ? voicing_buf_[output_rpos_] : 0.0f;
            voiced_mix_ += kVoicingFade * (mix_target - voiced_mix_);

            const float dry = dry_buf_[output_rpos_];
            const float wet = has_psola ? output_buf_[output_rpos_] / weight : dry;
            output[i] = dry + voiced_mix_ * (wet - dry);

            output_buf_[output_rpos_] = 0.0f;
            output_weight_buf_[output_rpos_] = 0.0f;
            dry_buf_[output_rpos_] = 0.0f;
            voicing_buf_[output_rpos_] = 0.0f;

            input_wpos_ = wrapIndex(static_cast<std::int64_t>(input_wpos_) + 1, kInputBufSize);
            output_rpos_ = wrapIndex(static_cast<std::int64_t>(output_rpos_) + 1, kOutputBufSize);
            ++input_time_;
        }
    }
private:
    /**
     * @brief 将有符号位置映射到 2 的幂大小的环形缓冲区。
     * @param[in] index 单调采样位置或相对位置。
     * @param[in] size 环形缓冲区大小。
     * @return 缓冲区内的有效索引。
     */
    static int wrapIndex(std::int64_t index, int size) noexcept {
        return static_cast<int>(index & static_cast<std::int64_t>(size - 1));
    }

    /**
     * @brief 运行周期估计、GCI 跟踪和延迟 Grain 发布。
     * @param[in] time 当前输入的单调采样位置。
     * @param[in] x 当前输入采样。
     */
    void runAnalyzer(std::int64_t time, float x) noexcept {
        if (++yin_counter_ >= kYinInterval) {
            yin_counter_ = 0;
            runYin(time);
        }

        runGci(time, x);
        publishReadyGrains(time);
    }

    /**
     * @brief 使用 FFT 加速 YIN 估计基音周期和周期清晰度。
     * @param[in] time 当前输入的单调采样位置。
     */
    void runYin(std::int64_t time) noexcept {
        float raw_energy = 0.0f;
        const std::int64_t frame_start = time - kYinFrameSize + 1;
        for (int i = 0; i < kYinFrameSize; ++i) {
            const float sample = input_buf_[wrapIndex(frame_start + i, kInputBufSize)];
            yin_frame_[i] = sample;
            raw_energy += sample * sample;
        }
        const float rms = std::sqrt(raw_energy / static_cast<float>(kYinFrameSize));

        yin_.Process(yin_frame_);
        const qwqdsp_pitch::Pitch result = yin_.GetPitch();
        const float clarity = std::clamp(1.0f - result.non_period_ratio, 0.0f, 1.0f);

        float detected_period = 0.0f;
        if (result.pitch_hz >= kMinF0 && result.pitch_hz <= kMaxF0 && result.non_period_ratio <= 0.45f) {
            detected_period = sample_rate_ / result.pitch_hz;
            if (yin_previous_period_ > 0.0f) {
                const float ratio = detected_period / yin_previous_period_;
                if (ratio > 1.55f && detected_period * 0.5f >= static_cast<float>(min_period_)) {
                    detected_period *= 0.5f;
                }
                else if (ratio < 0.65f && detected_period * 2.0f <= static_cast<float>(max_period_)) {
                    detected_period *= 2.0f;
                }
            }
        }

        if (detected_period > 0.0f && clarity >= 0.35f) {
            if (yin_previous_period_ > 0.0f) {
                detected_period = 0.55f * yin_previous_period_ + 0.45f * detected_period;
            }
            yin_period_ = detected_period;
            yin_previous_period_ = detected_period;
        }
        yin_clarity_ = clarity;
        updateVoicing(detected_period > 0.0f, clarity, rms);
    }

    /**
     * @brief 通过双阈值和连续帧确认更新清浊音状态。
     * @param[in] has_period 当前帧是否检测到有效周期。
     * @param[in] clarity 当前帧的 YIN 周期清晰度。
     * @param[in] rms 当前分析帧的线性 RMS。
     */
    void updateVoicing(bool has_period, float clarity, float rms) noexcept {
        constexpr float kSilenceRms = 0.001f;
        constexpr float kVoicedEnterClarity = 0.62f;
        constexpr float kVoicedExitClarity = 0.45f;
        const bool previous_state = yin_voiced_;

        if (yin_voiced_) {
            const bool lost_voice = rms < kSilenceRms || !has_period || clarity < kVoicedExitClarity;
            unvoiced_candidate_count_ = lost_voice ? unvoiced_candidate_count_ + 1 : 0;
            if (unvoiced_candidate_count_ >= 2) {
                yin_voiced_ = false;
            }
        }
        else {
            const bool found_voice = rms >= kSilenceRms && has_period && clarity >= kVoicedEnterClarity;
            voiced_candidate_count_ = found_voice ? voiced_candidate_count_ + 1 : 0;
            if (voiced_candidate_count_ >= 2) {
                yin_voiced_ = true;
            }
        }

        if (yin_voiced_ != previous_state) {
            voiced_candidate_count_ = 0;
            unvoiced_candidate_count_ = 0;
            resetPitchMarks();
            if (!yin_voiced_) {
                yin_period_ = 0.0f;
                yin_previous_period_ = 0.0f;
            }
        }
    }

    /**
     * @brief 清除 GCI、待发布 Grain 和当前合成 Grain 状态。
     */
    void resetPitchMarks() noexcept {
        previous_epoch_time_ = -1;
        previous_period_ = yin_period_ > 0.0f ? yin_period_ : static_cast<float>(default_period_);
        gci_since_last_ = 0;
        gci_searching_ = false;
        gci_peak_ = 0.0f;
        gci_peak_age_ = 0;
        latest_grain_.length = 0;
        pending_grains_.fill(PendingGrain{});
        synthesis_active_ = false;
    }

    /**
     * @brief 在预测周期附近搜索残差峰值并生成 GCI 标记。
     * @param[in] time 当前输入的单调采样位置。
     * @param[in] x 当前输入采样。
     */
    void runGci(std::int64_t time, float x) noexcept {
        const float residual = 0.5f * (x - residual_previous_input_);
        residual_previous_input_ = x;

        if (!yin_voiced_ || yin_period_ <= 0.0f) {
            return;
        }

        ++gci_since_last_;
        const float expected_period = previous_epoch_time_ >= 0 ? previous_period_ : yin_period_;
        const int search_start = std::max(1, static_cast<int>(0.55f * expected_period));
        const int earliest_finish = std::max(search_start + 1, static_cast<int>(0.78f * expected_period));
        const int latest_finish = std::max(earliest_finish + 1, static_cast<int>(1.35f * expected_period));

        if (!gci_searching_ && gci_since_last_ >= search_start) {
            gci_searching_ = true;
            gci_peak_ = std::abs(residual);
            gci_peak_time_ = time;
            gci_peak_age_ = 0;
        }
        else if (gci_searching_) {
            ++gci_peak_age_;
            const float magnitude = std::abs(residual);
            if (magnitude > gci_peak_) {
                gci_peak_ = magnitude;
                gci_peak_time_ = time;
                gci_peak_age_ = 0;
            }
        }

        const int settling_samples = std::max(4, static_cast<int>(0.04f * expected_period));
        const bool peak_settled =
            gci_searching_ && gci_since_last_ >= earliest_finish && gci_peak_age_ >= settling_samples;
        const bool search_expired = gci_searching_ && gci_since_last_ >= latest_finish;
        if (!peak_settled && !search_expired) {
            return;
        }

        int period = static_cast<int>(std::round(yin_period_));
        if (previous_epoch_time_ >= 0) {
            const std::int64_t measured = gci_peak_time_ - previous_epoch_time_;
            const float ratio = static_cast<float>(measured) / std::max(yin_period_, 1.0f);
            if (ratio >= 0.70f && ratio <= 1.40f) {
                period = static_cast<int>(measured);
            }
        }
        period = std::clamp(period, min_period_, max_period_);

        queueGrain(gci_peak_time_, period);
        previous_epoch_time_ = gci_peak_time_;
        previous_period_ = 0.7f * previous_period_ + 0.3f * static_cast<float>(period);
        gci_since_last_ = static_cast<int>(time - gci_peak_time_);
        gci_searching_ = false;
        gci_peak_ = 0.0f;
        gci_peak_age_ = 0;
    }

    /**
     * @brief 将 GCI 加入延迟提取队列，等待完整右半窗到达。
     * @param[in] epoch_time GCI 的单调采样位置。
     * @param[in] t0 当前基音周期，单位为采样数。
     */
    void queueGrain(std::int64_t epoch_time, int t0) noexcept {
        const float pitch_scale = std::max(pitch_scale_smoothed_, 0.05f);
        const float formant_scale = std::clamp(formant_shift_.load(std::memory_order_relaxed), 0.05f, 4.0f);
        const float support_ratio = std::max(1.0f, formant_scale / pitch_scale);
        const int radius =
            std::clamp(static_cast<int>(std::ceil(static_cast<float>(t0) * support_ratio)), t0, kGrainMaxLen / 2 - 1);

        PendingGrain* slot = nullptr;
        for (PendingGrain& pending : pending_grains_) {
            if (!pending.active) {
                slot = &pending;
                break;
            }
            if (slot == nullptr || pending.epoch_time < slot->epoch_time) {
                slot = &pending;
            }
        }

        *slot = PendingGrain{epoch_time, t0, radius, true};
    }

    /**
     * @brief 发布所有右半窗已经完整写入输入缓冲区的 Grain。
     * @param[in] time 当前输入的单调采样位置。
     */
    void publishReadyGrains(std::int64_t time) noexcept {
        while (true) {
            PendingGrain* ready = nullptr;
            for (PendingGrain& pending : pending_grains_) {
                if (!pending.active || time - pending.epoch_time < pending.radius) {
                    continue;
                }
                if (ready == nullptr || pending.epoch_time < ready->epoch_time) {
                    ready = &pending;
                }
            }
            if (ready == nullptr) {
                break;
            }

            publishGrain(*ready);
            ready->active = false;
        }
    }

    /**
     * @brief 从输入环形缓冲区提取一个未加窗的完整 Grain。
     * @param[in] pending 待发布 Grain 的时间和尺寸信息。
     */
    void publishGrain(const PendingGrain& pending) noexcept {
        const int length = std::min(2 * pending.radius, kGrainMaxLen);
        for (int i = 0; i < length; ++i) {
            const std::int64_t source_time = pending.epoch_time + i - pending.radius;
            latest_grain_.samples[i] = input_buf_[wrapIndex(source_time, kInputBufSize)];
        }
        latest_grain_.length = length;
        latest_grain_.t0 = pending.t0;
    }

    /**
     * @brief 围绕 Grain 中心进行线性插值重采样。
     * @param[in] source 原始未加窗 Grain。
     * @param[in] ratio 重采样倍率，大于 1 时缩短 Grain。
     * @return 重采样后的 Grain；尺寸越界时返回原 Grain。
     */
    static FormantGrain resampleGrain(const FormantGrain& source, float ratio) noexcept {
        FormantGrain destination;
        const int output_length = static_cast<int>(std::round(static_cast<float>(source.length) / ratio));
        if (output_length < 2 || output_length > kGrainMaxLen) {
            return source;
        }

        destination.length = output_length;
        destination.t0 = source.t0;
        const float source_center = 0.5f * static_cast<float>(source.length - 1);
        const float output_center = 0.5f * static_cast<float>(output_length - 1);
        for (int i = 0; i < output_length; ++i) {
            const float source_position = source_center + (static_cast<float>(i) - output_center) * ratio;
            const int index0 = std::clamp(static_cast<int>(std::floor(source_position)), 0, source.length - 1);
            const int index1 = std::min(index0 + 1, source.length - 1);
            const float fraction = std::clamp(source_position - static_cast<float>(index0), 0.0f, 1.0f);
            destination.samples[i] =
                source.samples[index0] + fraction * (source.samples[index1] - source.samples[index0]);
        }
        return destination;
    }

    /**
     * @brief 将一个 Grain 加窗并叠加到输出及归一化权重缓冲区。
     * @param[in] grain 待渲染的未加窗 Grain。
     * @param[in] center Grain 在输出环形缓冲区中的中心位置。
     */
    void renderGrain(const FormantGrain& grain, int center) noexcept {
        const int half = grain.length / 2;
        for (int i = 0; i < grain.length; ++i) {
            const float window = HannWindow::value(i, grain.length);
            const int output_index = wrapIndex(static_cast<std::int64_t>(center) + i - half, kOutputBufSize);
            output_buf_[output_index] += grain.samples[i] * window;
            output_weight_buf_[output_index] += window;
        }
    }

    float sample_rate_ = kSampleRate;
    int min_period_ = 0;
    int max_period_ = 0;
    int default_period_ = 0;

    std::array<float, kInputBufSize> input_buf_{};
    int input_wpos_ = 0;
    std::int64_t input_time_ = 0;

    FormantGrain latest_grain_{};
    FormantGrain synth_grain_{};
    std::array<PendingGrain, kPendingGrainCount> pending_grains_{};
    float phasor_phase_ = 0.0f;
    bool synthesis_active_ = false;
    float pitch_scale_smoothed_ = 1.0f;

    std::array<float, kOutputBufSize> output_buf_{};
    std::array<float, kOutputBufSize> output_weight_buf_{};
    std::array<float, kOutputBufSize> dry_buf_{};
    std::array<float, kOutputBufSize> voicing_buf_{};
    int output_rpos_ = 0;
    float voiced_mix_ = 0.0f;

    qwqdsp_pitch::Yin yin_{};
    std::array<float, kYinFrameSize> yin_frame_{};
    int yin_counter_ = 0;
    float yin_period_ = 0.0f;
    float yin_previous_period_ = 0.0f;
    float yin_clarity_ = 0.0f;
    bool yin_voiced_ = false;
    int voiced_candidate_count_ = 0;
    int unvoiced_candidate_count_ = 0;

    std::int64_t previous_epoch_time_ = -1;
    float previous_period_ = 0.0f;
    int gci_since_last_ = 0;
    bool gci_searching_ = false;
    float gci_peak_ = 0.0f;
    std::int64_t gci_peak_time_ = 0;
    int gci_peak_age_ = 0;
    float residual_previous_input_ = 0.0f;
};

static RealtimePsola s_psola;

// ------------------------------------------------------------
// miniaudio 回调
// ------------------------------------------------------------
/**
 * @brief 将 miniaudio 双工输入交给实时 PSOLA 处理器。
 * @param[in] device miniaudio 设备指针。
 * @param[out] output_buffer 输出缓冲区。
 * @param[in] input_buffer 输入缓冲区。
 * @param[in] frame_count 单声道帧数。
 */
extern "C" void maCallback(ma_device* device, void* output_buffer, const void* input_buffer, ma_uint32 frame_count) {
    (void)device;
    const auto* input = static_cast<const float*>(input_buffer);
    auto* output = static_cast<float*>(output_buffer);
    if (input == nullptr || output == nullptr) {
        return;
    }
    s_psola.process(input, output, static_cast<int>(frame_count));
}

// ------------------------------------------------------------
// 旋钮工具
// ------------------------------------------------------------
/**
 * @brief 创建并排列一个参数旋钮。
 * @param[in,out] x 旋钮横坐标，返回时移动到下一个位置。
 * @param[in] y 旋钮纵坐标。
 * @param[in] width 旋钮宽度。
 * @param[in] height 旋钮高度。
 * @param[in] title 显示标题。
 * @param[in] minimum 参数最小值。
 * @param[in] maximum 参数最大值。
 * @param[in] step 参数步进值。
 * @param[in] default_value 参数默认值。
 * @return 配置完成的旋钮。
 */
static Knob makeKnob(int& x, int y, int width, int height, const char* title, float minimum, float maximum, float step,
                     float default_value) {
    Knob knob;
    knob.set_bound(x, y, width, height)
        .set_title(title)
        .set_range(minimum, maximum, step, default_value)
        .set_value(default_value)
        .set_name_font_size(9)
        .set_number_font_size(9)
        .set_fore_color(Color{200, 200, 210, 255})
        .set_bg_color(Color{30, 30, 35, 255});
    x += width + 14;
    return knob;
}

// ------------------------------------------------------------
// main
// ------------------------------------------------------------
/**
 * @brief 启动实时 PSOLA 演示程序。
 * @return 正常退出时返回 0。
 */
int main() {
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, "Realtime PSOLA - Voiced/Unvoiced");
    SetTargetFPS(60);

    HannWindow::init();
    s_psola.init(kSampleRate);
    s_psola.pitch_scale_target_.store(1.0f, std::memory_order_relaxed);
    s_psola.pitch_shift_semitones_.store(0.0f, std::memory_order_relaxed);
    s_psola.formant_shift_.store(1.0f, std::memory_order_relaxed);

    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.playback.format = ma_format_f32;
    config.playback.channels = 1;
    config.sampleRate = static_cast<ma_uint32>(kSampleRate);
    config.dataCallback = maCallback;
    config.pUserData = nullptr;
    config.periodSizeInMilliseconds = 10;

    ma_device device;
    const bool audio_ok = ma_device_init(nullptr, &config, &device) == MA_SUCCESS;
    if (audio_ok) {
        ma_device_start(&device);
    }

    constexpr int kKnobWidth = 80;
    constexpr int kKnobHeight = 82;
    constexpr int kGap = 30;
    int knob_x = (kWindowWidth - (2 * kKnobWidth + kGap)) / 2;
    const int knob_y = (kWindowHeight - kKnobHeight) / 2 - 10;

    Knob pitch_knob = makeKnob(knob_x, knob_y, kKnobWidth, kKnobHeight, "Pitch", -12.0f, 12.0f, 0.5f, 0.0f);
    Knob formant_knob = makeKnob(knob_x, knob_y, kKnobWidth, kKnobHeight, "Formant", 0.5f, 2.0f, 0.05f, 1.0f);

    pitch_knob.on_value_change = [](float value) {
        s_psola.pitch_scale_target_.store(std::pow(2.0f, value / 12.0f), std::memory_order_relaxed);
        s_psola.pitch_shift_semitones_.store(value, std::memory_order_relaxed);
    };
    pitch_knob.value_to_text_function = [](float value) -> std::string {
        return value >= 0.0f ? TextFormat("+%.1f st", value) : TextFormat("%.1f st", value);
    };
    formant_knob.on_value_change = [](float value) { s_psola.formant_shift_.store(value, std::memory_order_relaxed); };
    formant_knob.value_to_text_function = [](float value) -> std::string { return TextFormat("%.2fx", value); };

    const float pitch_value = pitch_knob.get_value();
    s_psola.pitch_scale_target_.store(std::pow(2.0f, pitch_value / 12.0f), std::memory_order_relaxed);
    s_psola.pitch_shift_semitones_.store(pitch_value, std::memory_order_relaxed);
    s_psola.formant_shift_.store(formant_knob.get_value(), std::memory_order_relaxed);

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(Color{10, 10, 12, 255});
        pitch_knob.display();
        formant_knob.display();
        DrawText("Realtime PSOLA - normalized V/UV synthesis", 10, kWindowHeight - 24, 10, Color{80, 80, 90, 255});
        EndDrawing();
    }

    if (audio_ok) {
        ma_device_stop(&device);
        ma_device_uninit(&device);
    }
    CloseWindow();
    return 0;
}
