#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>

namespace pitch_shift_rt {

/**
 * @brief 基于互相关对齐的延迟线 WSOLA 实时移调器。
 *
 * 机制与 PitcherWsola 相同：
 *   输入逐采样写入环形延迟线；单个读头以 (1+speed) 倍速前进产生音高移动
 *   （read_head1_ 即当前延迟，每样本 -= speed）。
 *   升调时延迟收缩，降至下界（交叉淡化长度）触发跳转，读头跳到延迟线上界
 *   （固定目标位置）；降调时延迟增长，超过上界触发跳转，读头跳回目标位置。
 *   跳转前将源段（src）与目标段（target）各 512 采样做互相关，找到两段波形
 *   的最佳对齐偏移 cmax，新读头置于 target - cmax，再以 512 采样交叉淡化
 *   从旧读头平滑过渡到新读头，避免跳变咔哒与内容重播。
 *   固定容量环形缓冲，适合音频回调线程。延迟线容量决定可处理的最低频率
 *   （约采样率/延迟线大小）。
 */
class AutocorrelationDelayPitchShifter {
public:
    static constexpr std::size_t kBufferSize = 4096; ///< 环形延迟线容量（采样）
    static constexpr std::size_t kCrossFade = 512;   ///< 交叉淡化长度（也即互相关块大小）
    static constexpr std::size_t kTargetDelay = 2048; ///< 跳转目标延迟（缓冲中心，淡化两侧空间充足）
    static constexpr std::size_t kMinDelay = kCrossFade;      ///< 升调触发下界
    static constexpr std::size_t kMaxDelay = kBufferSize - kCrossFade; ///< 降调触发上界

    /** @brief 初始化采样率并清空处理状态。 */
    void init(float sample_rate) noexcept {
        sample_rate_ = sample_rate;
        reset();
    }

    /** @brief 设置移调半音数，范围限制为 +/-12。 */
    void setPitchShift(float semitones) noexcept {
        if (!std::isfinite(semitones)) return;
        target_speed_ = std::exp2(std::clamp(semitones, -12.0f, 12.0f) / 12.0f) - 1.0f;
    }

    /** @brief 清空延迟线与交叉淡化状态（保留移调设置）。 */
    void reset() noexcept {
        buffer_.fill(0.0f);
        write_pos_ = 0;
        read_head1_ = static_cast<float>(kTargetDelay);
        read_head2_ = -1.0f;
        speed_ = target_speed_;
        fade_ = false;
        fade_progress_ = 1.0f;
        fade_len_ = static_cast<float>(kCrossFade);
        out_sample_ = 0.0f;
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
            update(input[i]);
            output[i] = out_sample_;
        }
    }

    /** @brief 返回近似稳态延迟（采样）。 */
    static constexpr std::size_t latencySamples() noexcept { return kBufferSize / 2; }

private:
    // ------------------------------------------------------------
    // 读头与跳转
    // ------------------------------------------------------------

    void update(float sample) noexcept {
        // 速度平滑（慢 EMA，与 PitcherWsola 一致，避免速度抖动影响拼接）
        speed_ += 0.001f * (target_speed_ - speed_);

        // 读头前进（延迟变化）
        read_head1_ -= speed_;

        // 写入新样本
        buffer_[write_pos_] = sample;
        write_pos_ = (write_pos_ + 1) % kBufferSize;

        const float l1 = readInterp(read_head1_);

        if (fade_) {
            // 交叉淡化中：两路读头都前进（对照 PitcherWsola），
            // fade 结束读头1 接管读头2 位置。
            const float l2 = readInterp(read_head2_);
            read_head2_ -= speed_;
            fade_progress_ += 1.0f / fade_len_;
            const float t = std::min(1.0f, fade_progress_);
            // 余弦权重：w 从 0→1，从 l2（新位置）淡到 l1（旧位置），
            // fade 结束读头1 接管读头2（与 PitcherWsola 一致）。
            const float w = 0.5f - 0.5f * std::cos(t * static_cast<float>(std::numbers::pi));
            out_sample_ = l1 * w + l2 * (1.0f - w);
            if (fade_progress_ >= 1.0f) {
                read_head1_ = read_head2_;
                fade_ = false;
            }
        }
        else {
            bool crit;
            if (speed_ <= 0.0f) {
                // 降调：延迟增长，越过上界时跳回中心
                crit = read_head1_ >= static_cast<float>(kMaxDelay);
            }
            else {
                // 升调：延迟收缩，降至下界时跳回中心。触发点取 2×kMinDelay，
                // 保证 fade 期间读头1 从触发点继续收缩 kCrossFade 采样后仍
                // 高于 0（不读到未来数据）。
                crit = read_head1_ <= static_cast<float>(kMinDelay) * 2.0f;
            }
            if (crit) {
                // 互相关找目标段相对源段的最佳对齐偏移（高精度）
                copySegment(mem1_.data(), static_cast<int>(read_head1_));
                copySegment(mem2_.data(), static_cast<int>(kTargetDelay));
                const float cmax = computeMaxAcfPosition();
                // 新读头置于目标位置减去对齐偏移，使拼接处波形连续
                read_head2_ = static_cast<float>(kTargetDelay) - cmax;
                // 淡化时长 = kCrossFade/|speed|（与 PitcherWsola 一致）：
                // fade 期间读头移动恰好 kCrossFade 采样，之后有稳态段。
                // 升调 fade 长度同时受读头1 剩余空间限制（不越 0）。
                const float speed_abs = std::max(0.05f, std::fabs(speed_));
                fade_len_ = std::floor(static_cast<float>(kCrossFade) / speed_abs);
                if (speed_ > 0.0f) {
                    // fade 期间读头1 从触发点(1024)收缩，不能越过 0
                    fade_len_ = std::min(fade_len_,
                                         std::floor(static_cast<float>(kMinDelay) * 2.0f / speed_abs));
                }
                fade_len_ = std::clamp(fade_len_, 16.0f, static_cast<float>(kBufferSize) / 2.0f);
                fade_progress_ = 0.0f;
                fade_ = true;
            }
            out_sample_ = l1;
        }
    }

    /** @brief 将延迟线上延迟为 delay 的 kCrossFade 个采样复制到 dst。 */
    void copySegment(float* dst, int delay) const noexcept {
        int pos = write_pos_ - delay;
        while (pos < 0) pos += static_cast<int>(kBufferSize);
        pos %= static_cast<int>(kBufferSize);
        for (std::size_t i = 0; i < kCrossFade; ++i) {
            dst[i] = buffer_[(pos + static_cast<int>(i)) % kBufferSize];
        }
    }

    /**
     * @brief 高精度互相关：计算跳转对齐偏移 cmax。
     *
     * 对偏移 o ∈ [0, kCrossFade)，计算源段 mem1_ 与目标段 mem2_ 平移 o 的
     * 归一化互相关。取"第一个超过阈值 0.9 的显著局部峰"（避免谐波），
     * 抛物线细化到亚采样（限制 ±0.5 采样）。
     */
    float computeMaxAcfPosition() const noexcept {
        constexpr std::size_t kSearch = kCrossFade;
        constexpr float kPeakThresh = 0.9f;
        std::array<float, kCrossFade> scores{};
        for (std::size_t o = 0; o < kSearch; ++o) {
            float acc = 0.0f, e1 = 1.0e-9f, e2 = 1.0e-9f;
            for (std::size_t i = 0; i + o < kCrossFade; ++i) {
                const float a = mem1_[i];
                const float b = mem2_[i + o];
                acc += a * b;
                e1 += a * a;
                e2 += b * b;
            }
            scores[o] = acc / std::sqrt(e1 * e2);
        }
        // 从 0 向上扫，找第一个显著局部极大峰（互相关峰 o）
        std::size_t best_o = 0;
        for (std::size_t o = 1; o + 1 < kSearch; ++o) {
            const float s = scores[o];
            if (s >= kPeakThresh && s >= scores[o - 1] && s >= scores[o + 1]) {
                best_o = o;
                break;
            }
        }
        if (best_o == 0) {
            // 无显著峰（非周期段）：用全局最高分
            for (std::size_t o = 1; o < kSearch; ++o) {
                if (scores[o] > scores[best_o]) best_o = o;
            }
        }
        // 抛物线亚采样细化（限制 ±0.5 采样）
        float best_f = static_cast<float>(best_o);
        if (best_o > 0 && best_o + 1 < kSearch) {
            const float s_lo = scores[best_o - 1];
            const float s_hi = scores[best_o + 1];
            const float s_mid = scores[best_o];
            const float denom = s_lo - 2.0f * s_mid + s_hi;
            if (std::fabs(denom) > 1.0e-9f) {
                const float frac = 0.5f * (s_lo - s_hi) / denom;
                if (std::fabs(frac) < 0.5f) {
                    best_f = static_cast<float>(best_o) + frac;
                }
            }
        }
        return best_f;
    }

    // ------------------------------------------------------------
    // 插值
    // ------------------------------------------------------------

    /** @brief 4 点 Catmull-Rom 插值读取延迟线（dtime 为延迟）。 */
    float readInterp(float dtime) const noexcept {
        float read_pos = static_cast<float>(write_pos_) - dtime;
        while (read_pos < 0.0f) read_pos += static_cast<float>(kBufferSize);
        while (read_pos >= static_cast<float>(kBufferSize)) read_pos -= static_cast<float>(kBufferSize);
        const auto i1 = static_cast<std::size_t>(std::floor(read_pos));
        const float t = read_pos - std::floor(read_pos);
        const auto i0 = (i1 + kBufferSize - 1) % kBufferSize;
        const auto i2 = (i1 + 1) % kBufferSize;
        const auto i3 = (i1 + 2) % kBufferSize;
        const float t2 = t * t;
        const float t3 = t2 * t;
        const float a0 = -0.5f * t3 + t2 - 0.5f * t;
        const float a1 = 1.5f * t3 - 2.5f * t2 + 1.0f;
        const float a2 = -1.5f * t3 + 2.0f * t2 + 0.5f * t;
        const float a3 = 0.5f * t3 - 0.5f * t2;
        return a0 * buffer_[i0] + a1 * buffer_[i1] + a2 * buffer_[i2] + a3 * buffer_[i3];
    }

    float sample_rate_ = 48000.0f;
    std::array<float, kBufferSize> buffer_{};
    std::array<float, kCrossFade> mem1_{};
    std::array<float, kCrossFade> mem2_{};
    std::size_t write_pos_ = 0;
    float read_head1_ = 0.0f;  ///< 主读头（当前延迟，采样）
    float read_head2_ = -1.0f; ///< 交叉淡化读头（新位置，采样）
    float speed_ = 0.0f;       ///< 读头速度偏移（ratio - 1）
    float target_speed_ = 0.0f;///< 目标速度偏移
    bool fade_ = false;
    float fade_progress_ = 1.0f;
    float fade_len_ = static_cast<float>(kCrossFade);
    float out_sample_ = 0.0f;
};

} // namespace pitch_shift_rt
