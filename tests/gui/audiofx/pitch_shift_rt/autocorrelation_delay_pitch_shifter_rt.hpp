#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>

namespace pitch_shift_rt {

/**
 * @brief 基于自相关周期估计和可变延迟线的实时移调器。
 *
 * 输入逐采样写入环形延迟线，读头以移调比例倍速前进产生音高移动。
 * 处理器周期性搜索归一化自相关并估计基频周期；当读头与写头的距离越出
 * 上下限时，用双读头交叉淡化切换到自相关对齐的最佳位置，避免跳变咔哒声。
 * 所有缓冲均为固定容量，适合音频回调线程。
 */
class AutocorrelationDelayPitchShifter {
public:
    static constexpr std::size_t kBufferSize = 16384;      ///< 延迟线容量（采样）
    static constexpr std::size_t kAnalysisSize = 1024;     ///< 周期分析窗口
    static constexpr std::size_t kAnalysisHop = 128;       ///< 周期分析间隔
    static constexpr std::size_t kMinLag = 20;             ///< 自相关最小滞后
    static constexpr std::size_t kMaxLag = 1000;           ///< 自相关最大滞后
    static constexpr std::size_t kJumpCorrSize = 128;      ///< 跳转相关窗上限
    static constexpr std::size_t kFadeMax = 1024;          ///< 交叉淡化最大长度
    static constexpr std::size_t kFadeSafety = 16;         ///< 淡化期间读头安全余量
    static constexpr std::size_t kMaxJumpDistance = 4096;  ///< 跳转距离上限
    static constexpr std::size_t kInitialDelay = 1800;     ///< 升调启动延迟（首个跳转须在分析就绪后）
    static constexpr float kPeriodThreshold = 0.5f;        ///< 周期峰显著阈值

    /** @brief 初始化采样率并清空处理状态。 */
    void init(float sample_rate) noexcept {
        sample_rate_ = sample_rate;
        reset();
    }

    /** @brief 设置移调半音数，范围限制为 +/-12。 */
    void setPitchShift(float semitones) noexcept {
        if (!std::isfinite(semitones)) return;
        target_ratio_ = std::exp2(std::clamp(semitones, -12.0f, 12.0f) / 12.0f);
    }

    /** @brief 清空延迟线、自相关窗口和交叉淡化状态（保留移调设置）。 */
    void reset() noexcept {
        buffer_.fill(0.0f);
        write_pos_ = 0;
        written_ = 0;
        sample_count_ = 0;
        analysis_counter_ = 0;
        period_ = std::clamp(sample_rate_ / 300.0f,
                             static_cast<float>(kMinLag), static_cast<float>(kMaxLag));
        current_ratio_ = target_ratio_;
        read_b_active_ = false;
        fade_progress_ = 1.0f;
        fade_len_ = 1.0f;
        read_b_ = 0.0f;
        // 初始延迟按移调方向选择：升调预留延迟收缩空间，保证首个跳转发生在
        // 周期分析就绪之后；降调从最小延迟开始增长。
        const float initial_age = (target_ratio_ >= 1.0f)
            ? static_cast<float>(kInitialDelay)
            : minimumDelay();
        read_a_ = std::fmod(static_cast<float>(kBufferSize) - initial_age,
                            static_cast<float>(kBufferSize));
    }

    /**
     * @brief 处理一段单声道音频。
     * @param[in] input 输入采样缓冲
     * @param[out] output 输出采样缓冲
     * @param[in] frame_count 采样数
     */
    void process(const float* input, float* output, std::size_t frame_count) noexcept {
        if (output == nullptr) return;
        if (input == nullptr) {
            std::fill_n(output, frame_count, 0.0f);
            return;
        }
        for (std::size_t i = 0; i < frame_count; ++i) {
            buffer_[write_pos_] = input[i];
            write_pos_ = (write_pos_ + 1) % kBufferSize;
            written_ = std::min(written_ + 1, kBufferSize);
            sample_count_ = std::min(sample_count_ + 1, kAnalysisSize);
            if (++analysis_counter_ >= kAnalysisHop && sample_count_ == kAnalysisSize) {
                analysis_counter_ = 0;
                estimatePeriod();
            }
            current_ratio_ += 0.02f * (target_ratio_ - current_ratio_);
            const float ratio = current_ratio_;

            read_a_ = std::fmod(read_a_ + ratio, static_cast<float>(kBufferSize));
            if (read_b_active_) {
                read_b_ = std::fmod(read_b_ + ratio, static_cast<float>(kBufferSize));
                fade_progress_ = std::min(1.0f, fade_progress_ + 1.0f / fade_len_);
            }

            // 读头到写头（下一写槽）的距离，即当前延迟
            const float age = std::fmod(
                static_cast<float>(write_pos_) + static_cast<float>(kBufferSize) - read_a_,
                static_cast<float>(kBufferSize));

            if (!read_b_active_) {
                if (ratio > 1.0f && age < minimumDelay()) {
                    if (canStartJump(age, true)) startJump(true, age);
                } else if (ratio < 1.0f && age > maximumDelay()) {
                    if (canStartJump(age, false)) startJump(false, age);
                }
            }

            const float a = interpolate(read_a_);
            if (!read_b_active_) {
                output[i] = a;
                continue;
            }
            const float b = interpolate(read_b_);
            // 等功率（余弦）交叉淡化：A/B 相位不完全一致时线性淡化会产生
            // 幅度凹陷与边带；余弦曲线保持总能量恒定，拼接伪影更小。
            const float t = std::min(1.0f, fade_progress_);
            const float w = 0.5f - 0.5f * std::cos(t * static_cast<float>(std::numbers::pi));
            output[i] = a * (1.0f - w) + b * w;
            if (fade_progress_ >= 1.0f) {
                read_a_ = read_b_;
                read_b_active_ = false;
            }
        }
    }

    /** @brief 返回近似稳态延迟（采样）。 */
    static constexpr std::size_t latencySamples() noexcept { return kInitialDelay; }

private:
    // ------------------------------------------------------------
    // 延迟上下限
    // ------------------------------------------------------------

    /** @brief 相关窗长度：高音区随周期缩小，保证窗口不越过写头。 */
    float corrWindowSize() const noexcept {
        return std::clamp(period_ * 1.5f, 64.0f, static_cast<float>(kJumpCorrSize));
    }

    /**
     * @brief 触发升调跳转的最小延迟。
     *
     * 必须为交叉淡化预留空间：淡化期间读头延迟以 (ratio-1) 速率收缩，
     * 结束时的延迟不能低于相关窗下限，否则下一次跳转无法启动。因此
     * 触发延迟 = 窗下限 + 约半个拼接间隔（淡化时长）。
     */
    float minimumDelay() const noexcept {
        const float win_floor = corrWindowSize() + static_cast<float>(kFadeSafety) + 1.0f;
        const float fade_room = std::min(period_ * 0.5f, 96.0f);
        return std::max(std::clamp(period_ * 1.25f, 48.0f, 600.0f),
                        win_floor + fade_room);
    }

    /** @brief 触发降调跳转的最大延迟。 */
    float maximumDelay() const noexcept {
        return std::clamp(period_ * 3.0f, 160.0f, 3000.0f);
    }

    /**
     * @brief 判断当前延迟是否满足跳转条件。
     *
     * 要求周期分析已就绪、源窗全部位于已写入历史内。
     * 注意触发条件 age < minimumDelay() 与这里的下限不能互相矛盾：
     * minimumDelay() 恒大于此处下限，保证触发区间非空。
     */
    bool canStartJump(float age, bool pitch_up) const noexcept {
        if (sample_count_ != kAnalysisSize) return false;
        // 源相关窗 [current, current+win) 必须全部为已写入数据
        const float win_floor = corrWindowSize() + 1.0f;
        if (age < win_floor) return false;
        if (static_cast<float>(written_) < age + 1.0f) return false;
        const float min_delay = minimumDelay();
        const float max_delay = maximumDelay();
        // 落点延迟上限（升调取带顶，降调取最大延迟）必须在已写历史内
        const float max_candidate_age = pitch_up ? min_delay + period_ : max_delay;
        return static_cast<float>(written_) >=
               max_candidate_age + static_cast<float>(kFadeSafety);
    }

    /** @brief 开始一次双读头交叉淡化跳转。 */
    void startJump(bool pitch_up, float age) noexcept {
        const std::size_t current = static_cast<std::size_t>(read_a_);
        // 跳转距离 = 平滑周期估计（固定重叠长度自相关 + 正确符号的抛物线细化，
        // 峰位置无偏差，收敛到真实输入周期）。逐拼接距离偏差会在输出相位上
        // 累积，表现为频率偏移与低频 FM 边带，因此直接用周期估计而非相关细化
        // 的抖动值。
        const float jump_d = std::clamp(period_, static_cast<float>(kMinLag),
                                        static_cast<float>(kMaxJumpDistance));
        // 升调向旧数据跳（增大延迟），降调向写头方向跳（减小延迟）。
        // 基准必须用带小数部分的 read_a_：相关搜索在整数网格上进行，若以
        // 截断后的整数位置为基准，落点与读头实际相位会差一个分数采样，
        // 每次拼接都会引入系统性相位偏差（频率漂移 + 调制伪影）。
        const float offset = pitch_up ? -jump_d : jump_d;
        read_b_ = std::fmod(read_a_ + offset + static_cast<float>(kBufferSize),
                            static_cast<float>(kBufferSize));
        read_b_active_ = true;
        fade_progress_ = 0.0f;
        fade_len_ = computeFadeLength(pitch_up, age, jump_d);
    }

    /**
     * @brief 计算交叉淡化长度。
     *
     * 淡化时长覆盖整个拼接间隔（period/(ratio-1)），充分平滑跳转的残余
     * 相位步进；升调时还必须保证淡出读头在淡化结束前不低于相关窗下限
     * （由触发延迟预留保证）。
     */
    float computeFadeLength(bool pitch_up, float age, float jump_d) const noexcept {
        const float ratio = current_ratio_;
        const float speed = std::max(0.05f, std::fabs(ratio - 1.0f));
        float len = period_ / speed;
        len = std::clamp(len, 16.0f, static_cast<float>(kFadeMax));
        if (pitch_up && ratio > 1.0f) {
            const float win_floor = corrWindowSize() + static_cast<float>(kFadeSafety) + 1.0f;
            const float max_a = (age - win_floor) / speed;
            const float max_b = (age + jump_d - static_cast<float>(kFadeSafety)) / speed;
            len = std::min(len, std::min(max_a, max_b));
        }
        return std::max(1.0f, len);
    }

    // ------------------------------------------------------------
    // 周期估计
    // ------------------------------------------------------------

    /**
     * @brief 计算归一化自相关（窗口为最近 kAnalysisSize 个采样）。
     *
     * 所有滞后共用相同的重叠长度（最近 512 个采样），避免 (N-lag) 因子
     * 使抛物线峰位置偏向低滞后——否则周期估计会有约 -0.3 采样偏差，导致
     * 跳转距离偏离真实周期，产生频率漂移与 FM 边带。
     */
    float autocorrScore(std::size_t lag) const noexcept {
        constexpr std::size_t kOverlap = 512;
        float corr = 0.0f, e0 = 0.0f, e1 = 0.0f;
        const std::size_t end = kAnalysisSize;
        for (std::size_t i = end - kOverlap; i < end; ++i) {
            const float a = buffer_[(write_pos_ + kBufferSize - kAnalysisSize + i) % kBufferSize];
            const float b = buffer_[(write_pos_ + kBufferSize - kAnalysisSize + i - lag) % kBufferSize];
            corr += a * b;
            e0 += a * a;
            e1 += b * b;
        }
        return corr / std::sqrt(e0 * e1 + 1.0e-12f);
    }

    /**
     * @brief 升序扫描自相关，取首个超过阈值且为局部极大的滞后作为基频周期。
     *
     * 优先选择第一个显著峰，避免将正弦波锁到二次谐波；无显著峰（静音/噪声）
     * 时保持原估计。峰位置用抛物线拟合细化到亚采样精度。
     */
    void estimatePeriod() noexcept {
        float prev = autocorrScore(kMinLag - 1);
        for (std::size_t lag = kMinLag; lag <= kMaxLag; ++lag) {
            const float cur = autocorrScore(lag);
            const float next = autocorrScore(lag + 1);
            if (cur > kPeriodThreshold && cur >= prev && cur >= next) {
                const float denom = prev - 2.0f * cur + next;
                float frac = 0.0f;
                if (std::fabs(denom) > 1.0e-9f) {
                    frac = std::clamp(0.5f * (prev - next) / denom, -0.5f, 0.5f);
                }
                const float refined = static_cast<float>(lag) + frac;
                period_ = std::clamp(0.5f * period_ + 0.5f * refined,
                                     static_cast<float>(kMinLag), static_cast<float>(kMaxLag));
                return;
            }
            prev = cur;
        }
    }

    // ------------------------------------------------------------
    // 插值与状态
    // ------------------------------------------------------------

    /** @brief 4 点 Catmull-Rom 插值读取延迟线。 */
    float interpolate(float position) const noexcept {
        const auto i = static_cast<std::size_t>(position);
        const std::size_t i0 = (i + kBufferSize - 1) % kBufferSize;
        const std::size_t i2 = (i + 1) % kBufferSize;
        const std::size_t i3 = (i + 2) % kBufferSize;
        const float t = position - std::floor(position);
        const float t2 = t * t;
        const float t3 = t2 * t;
        const float a0 = -0.5f * t3 + t2 - 0.5f * t;
        const float a1 = 1.5f * t3 - 2.5f * t2 + 1.0f;
        const float a2 = -1.5f * t3 + 2.0f * t2 + 0.5f * t;
        const float a3 = 0.5f * t3 - 0.5f * t2;
        return a0 * buffer_[i0] + a1 * buffer_[i] + a2 * buffer_[i2] + a3 * buffer_[i3];
    }

    float sample_rate_ = 48000.0f;
    std::array<float, kBufferSize> buffer_{};
    std::size_t write_pos_ = 0;
    std::size_t written_ = 0;
    std::size_t sample_count_ = 0;
    std::size_t analysis_counter_ = 0;
    float period_ = 160.0f;
    float target_ratio_ = 1.0f;
    float current_ratio_ = 1.0f;
    float read_a_ = 0.0f;
    float read_b_ = 0.0f;
    bool read_b_active_ = false;
    float fade_progress_ = 1.0f;
    float fade_len_ = 1.0f;
};

} // namespace pitch_shift_rt
