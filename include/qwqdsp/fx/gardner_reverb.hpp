#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <memory>
#include <numbers>

//------------------------------------------------------------------------------
// Gardner 扩散混响器（嵌套全通网络）。结构与参数取自
//
// Gardner, W. G. 1992. "The Virtual Acoustic Room." M.S. thesis, MIT Media Lab.
// §4.7 嵌套全通定义（式 4.9）、§4.9 通用全通混响器（图 4.10）、
// §4.10 三种房间尺寸的扩散混响器（图 4.11，取其中小房间与大厅两组）。
//
// 与论文的差异：
//   - 论文的晚期扩散混响器是单声道，此处左右声道各跑一份独立实例；论文的声道去相关
//     来自各声道不同的早期反射（§4.11），此处不含早期反射。
//   - 论文用实测插值把混响时间映射到环路增益，此处按环路总延迟解析换算（见 SetDecay）。
//------------------------------------------------------------------------------

namespace qwqdsp_fx {

namespace detail::gardner {

// ------------------------------------------------------------
// Delay — 整数抽头循环延迟线
// ------------------------------------------------------------
/// @brief 整数抽头循环延迟线，缓冲区取 2 的幂以便用掩码取模
class Delay {
public:
    /// @brief 分配缓冲区
    /// @param num_samples 延迟长度（样本），不足 1 按 1 处理
    void Init(size_t num_samples) noexcept {
        num_samples_ = std::max<size_t>(num_samples, 1);
        size_t buffer_size = 1;
        while (buffer_size < num_samples_ + 1) {
            buffer_size *= 2;
        }
        buffer_.reset(new float[buffer_size]);
        std::memset(buffer_.get(), 0, buffer_size * sizeof(float));
        mask_ = buffer_size - 1;
        write_idx_ = 0;
    }

    /// @brief 返回 num_samples 个样本前写入的值，并写入 x
    float TapAndPush(float x) noexcept {
        float const out = buffer_[(write_idx_ - num_samples_) & mask_];
        buffer_[write_idx_] = x;
        write_idx_ = (write_idx_ + 1) & mask_;
        return out;
    }

    void Reset() noexcept {
        std::memset(buffer_.get(), 0, (mask_ + 1) * sizeof(float));
        write_idx_ = 0;
    }
private:
    std::unique_ptr<float[]> buffer_;
    size_t mask_ = 0;
    size_t num_samples_ = 1;
    size_t write_idx_ = 0;
};

// ------------------------------------------------------------
// OnePole — 环路空气吸收低通
// ------------------------------------------------------------
/// @brief 直流增益为 1 的一阶低通
class OnePole {
public:
    void Init(float sample_rate) noexcept {
        sample_rate_ = sample_rate;
        Recalc();
    }

    void SetCutoff(float cutoff /* Hz */) noexcept {
        cutoff_ = cutoff;
        Recalc();
    }

    float Process(float x) noexcept {
        z_ = x * a_ + z_ * b_;
        return z_;
    }

    void Reset() noexcept {
        z_ = 0.0f;
    }
private:
    void Recalc() noexcept {
        b_ = std::exp(-2 * std::numbers::pi_v<float> * cutoff_ / sample_rate_);
        a_ = 1 - b_;
    }

    float sample_rate_ = 1.0f;
    float cutoff_ = 0.0f;
    float a_ = 0.0f;
    float b_ = 0.0f;
    float z_ = 0.0f;
};

// ------------------------------------------------------------
// Allpass — 式 4.9 中 G 为纯延迟的情况
// ------------------------------------------------------------
/// @brief Schroeder 全通，传递函数 (z^-D - g) / (1 - g z^-D)
class Allpass {
public:
    void Init(size_t num_samples, float gain) noexcept {
        delay_.Init(num_samples);
        gain_ = gain;
        pending_ = 0.0f;
    }

    float Tick(float x) noexcept {
        float const d = delay_.TapAndPush(pending_);
        float const w = x + gain_ * d;
        pending_ = w;
        return d - gain_ * w;
    }

    void Reset() noexcept {
        delay_.Reset();
        pending_ = 0.0f;
    }
private:
    Delay delay_;
    float gain_ = 0.0f;
    // 延迟线写入比抽头晚一拍，等效多一拍延迟，对混响无影响
    float pending_ = 0.0f;
};

// ------------------------------------------------------------
// NestedAllpass — 式 4.9，延迟元件换成「纯延迟 + 内层全通级联」
// ------------------------------------------------------------
/// @brief 嵌套全通，延迟元件 = 纯延迟(p) 接 kNumInner 级内层全通（p = 外圈延迟 − 内层延迟之和）
template <size_t kNumInner>
class NestedAllpass {
public:
    /// @brief 初始化，延迟长度按 ms 换算为样本
    /// @param sample_rate 采样率（Hz）
    /// @param outer_ms 外圈延迟长度（ms），须不小于 inner_ms 之和
    /// @param outer_gain 外圈增益
    /// @param inner_ms 各级内层全通的延迟长度（ms），对应图 4.11 括号内自左至右的顺序
    /// @param inner_gain 各级内层全通的增益
    void Init(float sample_rate, float outer_ms, float outer_gain, std::array<float, kNumInner> inner_ms,
              std::array<float, kNumInner> inner_gain) noexcept {
        float inner_total_ms = 0.0f;
        for (size_t i = 0; i < kNumInner; ++i) {
            inner_total_ms += inner_ms[i];
            inner_[i].Init(MsToSamples(inner_ms[i], sample_rate), inner_gain[i]);
        }
        plain_ms_ = outer_ms - inner_total_ms;
        outer_gain_ = outer_gain;
        plain_.Init(MsToSamples(plain_ms_, sample_rate));
        pending_ = 0.0f;
    }

    float Tick(float x) noexcept {
        float d = plain_.TapAndPush(pending_);
        for (size_t i = 0; i < kNumInner; ++i) {
            d = inner_[i].Tick(d);
        }
        float const w = x + outer_gain_ * d;
        pending_ = w;
        return d - outer_gain_ * w;
    }

    void Reset() noexcept {
        plain_.Reset();
        for (auto& ap : inner_) {
            ap.Reset();
        }
        pending_ = 0.0f;
    }
private:
    static size_t MsToSamples(float ms, float sample_rate) noexcept {
        return static_cast<size_t>(std::lround(ms * sample_rate * 0.001f));
    }

    Delay plain_;
    std::array<Allpass, kNumInner> inner_{};
    float plain_ms_ = 0.0f;
    float outer_gain_ = 0.0f;
    float pending_ = 0.0f;
};

} // namespace detail::gardner

//------------------------------------------------------------------------------
// 环路增益与混响时间的关系：每绕环路一圈幅度乘以环路增益，一圈耗时 kLoopDelayMs
//------------------------------------------------------------------------------
/// @brief 由目标 RT60 求环路增益
/// @param rt60_ms RT60 时间（ms）
/// @param loop_delay_ms 环路总延迟（ms）
inline float GardnerLoopGain(float rt60_ms, float loop_delay_ms) noexcept {
    if (rt60_ms <= 0.0f) {
        return 0.0f;
    }
    return std::exp(loop_delay_ms / rt60_ms * std::log(0.001f));
}

//------------------------------------------------------------------------------
// GardnerRoomReverb（小房间，图 4.11 上）
//
// 输入 → 24ms 延迟 → 双嵌套全通 35(0.3)[内层 22(0.4) + 8.3(0.6)]
//      → 嵌套全通 66(0.1)[内层 30(0.4)] → 环路低通 4.2kHz → 环路增益 → 输入
// 输出 = 0.5 ×（双嵌套全通输出）+ 0.5 ×（环路末端）
//
// 括号内为全通增益，数字为延迟长度（ms）。
//------------------------------------------------------------------------------

//==============================================================================
/** Gardner 小房间扩散混响器。
 */
class GardnerRoomReverb {
public:
    static constexpr float kLoopDelayMs = 24.0f + 35.0f + 66.0f;

    // 实测标定的有效环路延迟：嵌套全通的内层谐振使包络衰减慢于环路总延迟的解析预期（见 SetDecay）
    static constexpr float kEffectiveLoopDelayMs = 170.0f;

    // 论文给出的可用区间为 380~570ms；下限即无环路反馈时全通链自身的衰减（实测 T30 375ms）
    static constexpr float kMinDecayMs = 380.0f;
    static constexpr float kMaxDecayMs = 1200.0f;
    static constexpr float kDefaultDecayMs = 500.0f;
    static constexpr float kDefaultDamping = 4200.0f;
    static constexpr float kDefaultMix = 0.35f;

    /// @brief 初始化采样率，预分配所有延迟线
    /// @param sample_rate 采样率（Hz）
    void Init(float sample_rate) noexcept {
        for (auto& ch : channels_) {
            ch.in_delay_.Init(SamplesOf(24.0f, sample_rate));
            ch.ap1_.Init(sample_rate, 35.0f, 0.3f, {22.0f, 8.3f}, {0.4f, 0.6f});
            ch.ap2_.Init(sample_rate, 66.0f, 0.1f, {30.0f}, {0.4f});
            ch.loop_lpf_.Init(sample_rate);
            ch.Reset();
        }
        SetDamping(kDefaultDamping);
        SetDecay(kDefaultDecayMs);
    }

    /// @brief 设置干湿比
    /// @param m 干湿比 [0, 1]
    void SetMix(float m) noexcept {
        mix_ = std::clamp(m, 0.0f, 1.0f);
    }

    /// @brief 设置环路低通截止频率，控制尾部高频衰减
    /// @param cutoff 截止频率（Hz），范围 [16, 20000]
    void SetDamping(float cutoff) noexcept {
        cutoff = std::clamp(cutoff, 16.0f, 20000.0f);
        for (auto& ch : channels_) {
            ch.loop_lpf_.SetCutoff(cutoff);
        }
    }

    /// @brief 设置 RT60 混响衰减时间，按有效环路延迟换算环路增益
    /// @param rt60_ms RT60 时间（ms），范围 [0, 20000]
    /// @note 有效环路延迟为实测标定值（kEffectiveLoopDelayMs），需要精确控制衰减时用 SetLoopGain。
    void SetDecay(float rt60_ms) noexcept {
        SetLoopGain(GardnerLoopGain(std::clamp(rt60_ms, 0.0f, 20000.0f), kEffectiveLoopDelayMs));
    }

    /// @brief 设置环路增益，即论文图中环路上的 gain
    /// @param gain 环路增益，范围 [0, 1)
    void SetLoopGain(float gain) noexcept {
        loop_gain_ = std::clamp(gain, 0.0f, 0.999f);
    }

    /// @brief 处理一对立体声样本
    /// @param dry_left 干信号左声道
    /// @param dry_right 干信号右声道
    /// @param[out] left_out 处理后左声道
    /// @param[out] right_out 处理后右声道
    void Process(float dry_left, float dry_right, float* left_out, float* right_out) noexcept {
        float const wet_left = channels_[0].Tick(dry_left, loop_gain_);
        float const wet_right = channels_[1].Tick(dry_right, loop_gain_);
        *left_out = dry_left * (1 - mix_) + wet_left * mix_;
        *right_out = dry_right * (1 - mix_) + wet_right * mix_;
    }

    /// @brief 批量处理立体声样本（in-place）
    /// @param l 左声道输入/输出缓冲
    /// @param r 右声道输入/输出缓冲
    /// @param num 样本帧数
    void Process(float* l, float* r, size_t num) noexcept {
        for (size_t i = 0; i < num; ++i) {
            Process(l[i], r[i], &l[i], &r[i]);
        }
    }

    /// @brief 复位所有延迟线与滤波器状态
    void Reset() noexcept {
        for (auto& ch : channels_) {
            ch.Reset();
        }
    }
private:
    struct Channel {
        float Tick(float x, float loop_gain) noexcept {
            float y = in_delay_.TapAndPush(x + feedback_);
            y = ap1_.Tick(y);
            float const tap1 = y;
            y = ap2_.Tick(y);
            feedback_ = loop_lpf_.Process(y) * loop_gain;
            return 0.5f * tap1 + 0.5f * y;
        }

        void Reset() noexcept {
            in_delay_.Reset();
            ap1_.Reset();
            ap2_.Reset();
            loop_lpf_.Reset();
            feedback_ = 0.0f;
        }

        detail::gardner::Delay in_delay_;
        detail::gardner::NestedAllpass<2> ap1_;
        detail::gardner::NestedAllpass<1> ap2_;
        detail::gardner::OnePole loop_lpf_;
        float feedback_ = 0.0f;
    };

    static size_t SamplesOf(float ms, float sample_rate) noexcept {
        return static_cast<size_t>(std::lround(ms * sample_rate * 0.001f));
    }

    float mix_ = 0.0f;
    float loop_gain_ = 0.0f;
    std::array<Channel, 2> channels_;
};

//------------------------------------------------------------------------------
// GardnerHallReverb（大厅，图 4.11 下）
//
// 输入 → 全通 8(0.3) → 全通 12(0.3) → 4ms → 17ms → 嵌套全通 87(0.5)[内层 62(0.25)]
//      → 31ms → 3ms → 双嵌套全通 120(0.5)[内层 76(0.25) + 30(0.25)]
//      → 环路低通 2.6kHz → 环路增益 → 输入
// 输出 = 0.34 ×（17ms 后）+ 0.14 ×（31ms 后）+ 0.14 ×（环路末端）
//
// 括号内为全通增益，数字为延迟长度（ms）。
//------------------------------------------------------------------------------

//==============================================================================
/** Gardner 大厅扩散混响器。
 */
class GardnerHallReverb {
public:
    static constexpr float kLoopDelayMs = 8.0f + 12.0f + 4.0f + 17.0f + 87.0f + 31.0f + 3.0f + 120.0f;

    // 实测标定的有效环路延迟：嵌套全通的内层谐振使包络衰减慢于环路总延迟的解析预期（见 SetDecay）
    static constexpr float kEffectiveLoopDelayMs = 450.0f;

    // 论文给出的可用区间为 1.3s 以上；下限附近全通链自身的衰减（实测 T30 1607ms）会主导
    static constexpr float kMinDecayMs = 1500.0f;
    static constexpr float kMaxDecayMs = 8000.0f;
    static constexpr float kDefaultDecayMs = 2500.0f;
    static constexpr float kDefaultDamping = 2600.0f;
    static constexpr float kDefaultMix = 0.35f;

    /// @brief 初始化采样率，预分配所有延迟线
    /// @param sample_rate 采样率（Hz）
    void Init(float sample_rate) noexcept {
        for (auto& ch : channels_) {
            ch.ap1_.Init(sample_rate, 8.0f, 0.3f, {}, {});
            ch.ap2_.Init(sample_rate, 12.0f, 0.3f, {}, {});
            ch.d1_.Init(SamplesOf(4.0f, sample_rate));
            ch.d2_.Init(SamplesOf(17.0f, sample_rate));
            ch.ap3_.Init(sample_rate, 87.0f, 0.5f, {62.0f}, {0.25f});
            ch.d3_.Init(SamplesOf(31.0f, sample_rate));
            ch.d4_.Init(SamplesOf(3.0f, sample_rate));
            ch.ap4_.Init(sample_rate, 120.0f, 0.5f, {76.0f, 30.0f}, {0.25f, 0.25f});
            ch.loop_lpf_.Init(sample_rate);
            ch.Reset();
        }
        SetDamping(kDefaultDamping);
        SetDecay(kDefaultDecayMs);
    }

    /// @brief 设置干湿比
    /// @param m 干湿比 [0, 1]
    void SetMix(float m) noexcept {
        mix_ = std::clamp(m, 0.0f, 1.0f);
    }

    /// @brief 设置环路低通截止频率，控制尾部高频衰减
    /// @param cutoff 截止频率（Hz），范围 [16, 20000]
    void SetDamping(float cutoff) noexcept {
        cutoff = std::clamp(cutoff, 16.0f, 20000.0f);
        for (auto& ch : channels_) {
            ch.loop_lpf_.SetCutoff(cutoff);
        }
    }

    /// @brief 设置 RT60 混响衰减时间，按有效环路延迟换算环路增益
    /// @param rt60_ms RT60 时间（ms），范围 [0, 20000]
    /// @note 有效环路延迟为实测标定值（kEffectiveLoopDelayMs），需要精确控制衰减时用 SetLoopGain。
    void SetDecay(float rt60_ms) noexcept {
        SetLoopGain(GardnerLoopGain(std::clamp(rt60_ms, 0.0f, 20000.0f), kEffectiveLoopDelayMs));
    }

    /// @brief 设置环路增益，即论文图中环路上的 gain
    /// @param gain 环路增益，范围 [0, 1)
    void SetLoopGain(float gain) noexcept {
        loop_gain_ = std::clamp(gain, 0.0f, 0.999f);
    }

    /// @brief 处理一对立体声样本
    /// @param dry_left 干信号左声道
    /// @param dry_right 干信号右声道
    /// @param[out] left_out 处理后左声道
    /// @param[out] right_out 处理后右声道
    void Process(float dry_left, float dry_right, float* left_out, float* right_out) noexcept {
        float const wet_left = channels_[0].Tick(dry_left, loop_gain_);
        float const wet_right = channels_[1].Tick(dry_right, loop_gain_);
        *left_out = dry_left * (1 - mix_) + wet_left * mix_;
        *right_out = dry_right * (1 - mix_) + wet_right * mix_;
    }

    /// @brief 批量处理立体声样本（in-place）
    /// @param l 左声道输入/输出缓冲
    /// @param r 右声道输入/输出缓冲
    /// @param num 样本帧数
    void Process(float* l, float* r, size_t num) noexcept {
        for (size_t i = 0; i < num; ++i) {
            Process(l[i], r[i], &l[i], &r[i]);
        }
    }

    /// @brief 复位所有延迟线与滤波器状态
    void Reset() noexcept {
        for (auto& ch : channels_) {
            ch.Reset();
        }
    }
private:
    struct Channel {
        float Tick(float x, float loop_gain) noexcept {
            float y = ap1_.Tick(x + feedback_);
            y = ap2_.Tick(y);
            y = d1_.TapAndPush(y);
            y = d2_.TapAndPush(y);
            float const tap1 = y;
            y = ap3_.Tick(y);
            y = d3_.TapAndPush(y);
            float const tap2 = y;
            y = d4_.TapAndPush(y);
            y = ap4_.Tick(y);
            feedback_ = loop_lpf_.Process(y) * loop_gain;
            return 0.34f * tap1 + 0.14f * tap2 + 0.14f * y;
        }

        void Reset() noexcept {
            ap1_.Reset();
            ap2_.Reset();
            ap3_.Reset();
            ap4_.Reset();
            d1_.Reset();
            d2_.Reset();
            d3_.Reset();
            d4_.Reset();
            loop_lpf_.Reset();
            feedback_ = 0.0f;
        }

        detail::gardner::NestedAllpass<0> ap1_;
        detail::gardner::NestedAllpass<0> ap2_;
        detail::gardner::NestedAllpass<1> ap3_;
        detail::gardner::NestedAllpass<2> ap4_;
        detail::gardner::Delay d1_;
        detail::gardner::Delay d2_;
        detail::gardner::Delay d3_;
        detail::gardner::Delay d4_;
        detail::gardner::OnePole loop_lpf_;
        float feedback_ = 0.0f;
    };

    static size_t SamplesOf(float ms, float sample_rate) noexcept {
        return static_cast<size_t>(std::lround(ms * sample_rate * 0.001f));
    }

    float mix_ = 0.0f;
    float loop_gain_ = 0.0f;
    std::array<Channel, 2> channels_;
};

} // namespace qwqdsp_fx
