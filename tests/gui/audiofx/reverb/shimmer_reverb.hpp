#pragma once
#include "qwqdsp/fx/pitch_shifter.hpp"
#include "qwqdsp/polymath.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <memory>
#include <numbers>

namespace detail {

// ------------------------------------------------------------
//  OnePoleFilter
// ------------------------------------------------------------
class OnePoleFilter {
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
        a_ = 0;
        b_ = 0;
        z_ = 0;
    }
private:
    float sample_rate_ = 1;
    float cutoff_ = 0;

    float a_ = 0;
    float b_ = 0;
    float z_ = 0;

    void Recalc() noexcept {
        b_ = std::exp(-2 * std::numbers::pi_v<float> * cutoff_ / sample_rate_);
        a_ = 1 - b_;
    }
};

// ------------------------------------------------------------
//  DelayLine
// ------------------------------------------------------------
class DelayLine {
public:
    void Init(size_t size) noexcept {
        size_ = size;
        // For speed, create a bigger buffer than we really need.
        size_t buffer_size = CeilPowerOfTwo(size_);
        buffer_.reset(new float[size_t(buffer_size)]);
        std::memset(buffer_.get(), 0, size_t(buffer_size) * sizeof(float));

        mask_ = buffer_size - 1;
        write_idx_ = 0;
    }

    void Push(float val) noexcept {
        buffer_[size_t(write_idx_++)] = val;
        write_idx_ &= mask_;
    }

    float Tap(float delay /* samples */) const noexcept {
        assert(delay <= static_cast<float>(size_));

        size_t d = (size_t)delay;
        float frac = 1 - (delay - static_cast<float>(d));

        size_t read_idx = (write_idx_ - 1) - d;
        float a = buffer_[size_t((read_idx - 1) & mask_)];
        float b = buffer_[size_t(read_idx & mask_)];

        return a + (b - a) * frac;
    }

    // This does read-before-write.
    float TapAndPush(float delay, float val) noexcept {
        float out = Tap(delay);
        Push(val);
        return out;
    }

    void Reset() noexcept {
        std::memset(buffer_.get(), 0, size_t(CeilPowerOfTwo(size_)) * sizeof(float));
        write_idx_ = 0;
    }

    size_t GetSize() const noexcept {
        return size_;
    }
private:
    size_t size_ = 0;

    std::unique_ptr<float[]> buffer_;
    size_t mask_ = 0;

    size_t write_idx_ = 0;

    static size_t CeilPowerOfTwo(size_t n) noexcept {
        size_t r = 1;
        while (r < n) {
            r *= 2;
        }
        return r;
    }
};

//------------------------------------------
// DelayAllpass
//------------------------------------------

class DelayAllpass {
public:
    void Init(size_t size, float gain) noexcept {
        delay_line_.Init(size);
        gain_ = gain;
    }

    float Process(float x, float delay) noexcept {
        float wd = delay_line_.Tap(delay);
        float w = x + gain_ * wd;
        float y = -gain_ * w + wd;
        delay_line_.Push(w);
        return y;
    }

    void SetGain(float gain) noexcept {
        gain_ = gain;
    }

    float Tap(float delay) const noexcept {
        return delay_line_.Tap(delay);
    }

    size_t GetSize() const noexcept {
        return delay_line_.GetSize();
    }

    void Reset() noexcept {
        delay_line_.Reset();
    }
private:
    DelayLine delay_line_;
    float gain_ = 0.0f;
};

//--------------------------------------------------------------
// Lfo
//--------------------------------------------------------------

class Lfo {
public:
    void Init(float sample_rate) noexcept {
        sample_rate_ = sample_rate;
        Recalc();
    }

    void SetFrequency(float freq) noexcept {
        freq_ = freq;
        Recalc();
    }

    float Process() noexcept {
        float out = -qwqdsp::polymath::SinParabola(phase_);

        phase_ += phase_inc_;
        if (phase_ > std::numbers::pi_v<float>)
            phase_ = -std::numbers::pi_v<float>;

        return out;
    }

    void Reset() noexcept {
        phase_ = 0;
    }
private:
    float sample_rate_ = 1;
    float freq_ = 0;

    float phase_inc_ = 0;
    float phase_ = -std::numbers::pi_v<float>;

    void Recalc() noexcept {
        phase_inc_ = freq_ / sample_rate_;
        phase_inc_ *= 2 * std::numbers::pi_v<float>;
    }
};

//------------------------------------------
// Tank
//------------------------------------------

class Tank {
public:
    void Init(size_t apf1_size, float apf1_gain, // First APF
              float max_mod_depth,
              size_t delay1_size,                // First delay
              size_t apf2_size, float apf2_gain, // Second APF
              size_t delay2_size                 // Second delay
              ) noexcept {
        apf1_size_ = apf1_size;
        max_mod_depth_ = max_mod_depth;
        float max_apf1_size = static_cast<float>(apf1_size_) + max_mod_depth_ + 1;
        apf1_.Init(size_t(max_apf1_size), apf1_gain);

        del1_.Init(delay1_size);
        apf2_.Init(apf2_size, apf2_gain);
        del2_.Init(delay2_size);

        RecalcSizeRatio();
    }

    void SetDecay(float rt60_ms, float sample_rate) noexcept {
        rt60_ms_ = rt60_ms;
        sample_rate_ = sample_rate;
        RecalcDecay();
    }

    float GetDecayGain() const noexcept {
        return decay_gain_;
    }

    void SetSizeRatio(float size_ratio) noexcept {
        size_ratio_ = size_ratio;
        RecalcSizeRatio();
        RecalcDecay();
    }

    void Process(float val) noexcept {
        // APF1: "Controls density of tail."
        val = apf1_.Process(val, apf1_delay_ + lfo_.Process() * mod_depth_);
        val = del1_.TapAndPush(del1_delay_, val);

        val = damping_.Process(val);
        val *= decay_gain_;

        // APF2: "Decorrelates tank signals."
        val = apf2_.Process(val, apf2_delay_);
        val = del2_.TapAndPush(del2_delay_, val);

        out_ = val;

        // Shimmer: 移调 + 软限幅，供下一帧交叉反馈使用
        shifted_out_ = SoftClip(shifter_.Tick(out_));
    }

    void SetPitchShift(float semitones) noexcept {
        shifter_.SetPitchShift(semitones);
    }

    void Reset() noexcept {
        apf1_.Reset();
        apf2_.Reset();
        del1_.Reset();
        del2_.Reset();
        damping_.Reset();
        lfo_.Reset();

        shifter_ = {};
        shifted_out_ = 0.0f;
    }

    float GetTankLength() const noexcept {
        return apf1_delay_ + apf2_delay_ + del1_delay_ + del2_delay_
               + static_cast<float>(qwqdsp_fx::PitchShifter::kLatency);
    }

    // --- public members (accessed directly by ShimmerReverb) ---
    float out_ = 0.0f;
    float shifted_out_ = 0.0f;  // 移调 + 限幅后的输出，供下一帧反馈

    qwqdsp_fx::PitchShifter shifter_;
    DelayAllpass apf1_;
    DelayAllpass apf2_;
    DelayLine del1_;
    DelayLine del2_;
    detail::OnePoleFilter damping_;
    detail::Lfo lfo_;
private:
    void RecalcDecay() noexcept {
        float rt60_samples = rt60_ms_ / 1000.0f * sample_rate_;
        float length = GetTankLength();
        if (length > 0.0f && rt60_samples > 0.0f) [[likely]] {
            decay_gain_ = std::exp(length * std::log(0.001f) / rt60_samples);
        }
        else {
            decay_gain_ = 0.0f;
        }
        apf2_.SetGain(std::clamp(decay_gain_ * 1.3f, 0.25f, 0.5f));
    }

    size_t apf1_size_ = 0;
    float max_mod_depth_ = 0;
    float mod_depth_ = 0;

    float apf1_delay_ = 0;
    float apf2_delay_ = 0;
    float del1_delay_ = 0;
    float del2_delay_ = 0;

    float decay_gain_ = 0.0f;
    float rt60_ms_ = 2000.0f;
    float sample_rate_ = 48000.0f;
    float size_ratio_ = 0;

    void RecalcSizeRatio() noexcept {
        apf1_delay_ = static_cast<float>(apf1_size_) * size_ratio_;
        mod_depth_ = max_mod_depth_ * size_ratio_;

        apf2_delay_ = static_cast<float>(apf2_.GetSize()) * size_ratio_;
        del1_delay_ = static_cast<float>(del1_.GetSize()) * size_ratio_;
        del2_delay_ = static_cast<float>(del2_.GetSize()) * size_ratio_;
    }

    /** 三次软限幅：线性通过小信号，±1.5 处饱和到 ±1.0 */
    static float SoftClip(float x) noexcept {
        constexpr float kThresh = 1.5f;
        float ax = std::abs(x);
        if (ax <= kThresh)
            return x - (4.0f / 27.0f) * x * x * x;
        return (x > 0 ? 1.0f : -1.0f);
    }
};

} // namespace detail

//------------------------------------------------------------------------------
// PlateReverb is an implementation of the classic plate reverb algorithm
// described by Jon Dattorro.
//
// Dattorro, J. 1997. "Effect Design Part 1: Reverberators and Other Filters."
// Journal of the Audio Engineering Society, Vol. 45, No. 9
//
// https://ccrma.stanford.edu/~dattorro/EffectDesignPart1.pdf
//
// Parameters:
//
//    mix:        Dry/wet mix.
//    predelay:   Delay before reverb.
//    lowpass:    Apply a lowpass filter before reverb.
//    decay:      RT60 混响衰减时间（ms）。每个 Tank 根据自身环路长度独立计算每采样点衰减增益。
//    size:       The size of our imaginary plate.
//    damping:    How much high frequencies are filtered during reverb.
//
//------------------------------------------------------------------------------

//==============================================================================
/** Shimmer reverb — Plate reverb with pitch shifters in the tank feedback path.
 */
class ShimmerReverb {
public:
    static constexpr float kMaxPredelay = 0.1f; // s
    static constexpr float kMaxSize = 2.0f;

    /**
     * @brief 初始化采样率，预分配所有延迟线。
     * @param sample_rate 采样率（Hz）
     */
    void Init(float sample_rate) noexcept {
        sample_rate_ = sample_rate;

        // Ratio of our sample rate to the sample rate that is used in
        // Dattorro's paper.
        float r = float(sample_rate_ / 29761.0f);

        // Predelay
        predelay_line_.Init((size_t)std::ceil(sample_rate_ * kMaxPredelay));

        // Lowpass filters
        lowpass_.Init(sample_rate_);
        left_tank_.damping_.Init(sample_rate_);
        right_tank_.damping_.Init(sample_rate_);

        // Diffusers
        diffusers_[0].Init((size_t)std::ceil(142 * r), 0.75f);
        diffusers_[1].Init((size_t)std::ceil(107 * r), 0.75f);
        diffusers_[2].Init((size_t)std::ceil(379 * r), 0.625f);
        diffusers_[3].Init((size_t)std::ceil(277 * r), 0.625f);

        // Tanks
        float max_mod_depth = float(8.0f * kMaxSize * r);
        left_tank_.Init((size_t)std::ceil(kMaxSize * 672 * r), -0.7f, // apf1
                        max_mod_depth,
                        (size_t)std::ceil(kMaxSize * 4453 * r),       // del1
                        (size_t)std::ceil(kMaxSize * 1800 * r), 0.5f, // apf2
                        (size_t)std::ceil(kMaxSize * 3720 * r)        // del2
        );
        right_tank_.Init((size_t)std::ceil(kMaxSize * 908 * r), -0.7f, // apf1
                         max_mod_depth,
                         (size_t)std::ceil(kMaxSize * 4217 * r),       // del1
                         (size_t)std::ceil(kMaxSize * 2656 * r), 0.5f, // apf2
                         (size_t)std::ceil(kMaxSize * 3163 * r)        // del2
        );

        left_tank_.lfo_.Init(sample_rate_);
        right_tank_.lfo_.Init(sample_rate_);
        left_tank_.lfo_.SetFrequency(1.0f);
        right_tank_.lfo_.SetFrequency(0.95f);

        // Tap points
        base_left_taps_ = {
            266.0f * r,  // right_tank_.del1_
            2974.0f * r, // right_tank_.del1_
            1913.0f * r, // right_tank_.apf2_
            1996.0f * r, // right_tank_.del2_
            1990.0f * r, // left_tank_.del1_
            187.0f * r,  // left_tank_.apf2_
            1066.0f * r, // left_tank_.del2_
        };
        base_right_taps_ = {
            353.0f * r,  // left_tank_.del1_
            3627.0f * r, // left_tank_.del1_
            1228.0f * r, // left_tank_.apf2_
            2673.0f * r, // left_tank_.del2_
            2111.0f * r, // right_tank_.del1_
            335.0f * r,  // right_tank_.apf2_
            121.0f * r,  // right_tank_.del2_
        };

        // 默认 RT60 2000ms
        SetDecay(2000.0f);

        // Shimmer: 初始化 Tank 内移调器
        left_tank_.SetPitchShift(pitch_shift_);
        right_tank_.SetPitchShift(pitch_shift_);
    }

    /**
     * @brief 设置干湿比。
     * @param m 干湿比 [0, 1]
     */
    void SetMix(float m) noexcept {
        mix_ = std::clamp(m, 0.0f, 1.0f);
    }

    /**
     * @brief 设置预延迟时间。
     * @param pd 预延迟（s），范围 [0, 0.1]
     */
    void SetPredelay(float pd) noexcept {
        predelay_ = std::clamp(pd, 0.0f, kMaxPredelay) * sample_rate_;
    }

    /**
     * @brief 设置输入低通滤波器截止频率。
     * @param cutoff 截止频率（Hz），范围 [16, 20000]
     */
    void SetLowpass(float cutoff) noexcept {
        cutoff = std::clamp(cutoff, 16.0f, 20000.0f);
        lowpass_.SetCutoff(cutoff);
    }

    /**
     * @brief 设置 RT60 混响衰减时间。每个 Tank 根据自身环路长度独立计算每采样点衰减增益。
     * @param rt60_ms RT60 时间（ms），范围 [0, 10000]
     */
    void SetDecay(float rt60_ms) noexcept {
        rt60_ms = std::clamp(rt60_ms, 0.0f, 10000.0f);
        left_tank_.SetDecay(rt60_ms, sample_rate_);
        right_tank_.SetDecay(rt60_ms, sample_rate_);
    }

    /**
     * @brief 设置板尺寸。缩放 Tank 内所有延迟线、APF 以及 Tap 点的延迟时间。
     * @note Dattorro 原论文中无此参数，此为扩展。
     * @param sz 板尺寸 [0, 2]
     */
    void SetSize(float sz) noexcept {
        float size_ratio = std::clamp(sz, 0.0f, kMaxSize) / kMaxSize;

        // Scale the tank delays and APFs in each tank
        left_tank_.SetSizeRatio(size_ratio);
        right_tank_.SetSizeRatio(size_ratio);

        // Scale the taps
        for (size_t i = 0; i < kNumTaps; i++) {
            left_taps_[size_t(i)] = base_left_taps_[size_t(i)] * size_ratio;
            right_taps_[size_t(i)] = base_right_taps_[size_t(i)] * size_ratio;
        }
    }

    /**
     * @brief 设置移调量。
     * @param semitones 半音数（如 +12 = 一个八度）
     */
    void SetPitchShift(float semitones) noexcept {
        pitch_shift_ = semitones;
        left_tank_.SetPitchShift(semitones);
        right_tank_.SetPitchShift(semitones);
    }

    /**
     * @brief 设置混响阻尼低通截止频率。控制混响尾部高频衰减量。
     * @param cutoff 阻尼截止频率（Hz），范围 [16, 20000]
     */
    void SetDamping(float cutoff) noexcept {
        cutoff = std::clamp(cutoff, 16.0f, 20000.0f);

        left_tank_.damping_.SetCutoff(cutoff);
        right_tank_.damping_.SetCutoff(cutoff);
    }

    /**
     * @brief 处理一对立体声样本。
     * @param dry_left  干信号左声道
     * @param dry_right 干信号右声道
     * @param[out] left_out  处理后左声道
     * @param[out] right_out 处理后右声道
     * @note 使用合成立体声（synthetic stereo）：左右声道求和后送入混响，
     *       再通过交叉 Tap 生成立体声输出。
     */
    void Process(float dry_left, float dry_right, float* left_out, float* right_out) noexcept {
        // Note that this is "synthetic stereo".  We produce a stereo pair
        // of output samples based on the summed input.
        float sum = dry_left + dry_right;

        // Predelay
        sum = predelay_line_.TapAndPush(predelay_, sum);

        // Input lowpass
        sum = lowpass_.Process(sum);

        // Diffusers
        sum = diffusers_[0].Process(sum, (float)diffusers_[0].GetSize());
        sum = diffusers_[1].Process(sum, (float)diffusers_[1].GetSize());
        sum = diffusers_[2].Process(sum, (float)diffusers_[2].GetSize());
        sum = diffusers_[3].Process(sum, (float)diffusers_[3].GetSize());

        // Tanks (cross-feedback with pitch-shifted signals)
        // 每个 Tank 内部已做移调 + 软限幅，结果存在 shifted_out_
        float left_in = sum + right_tank_.shifted_out_ * right_tank_.GetDecayGain();
        float right_in = sum + left_tank_.shifted_out_ * left_tank_.GetDecayGain();
        left_tank_.Process(left_in);
        right_tank_.Process(right_in);

        // Tap for output
        float wet_left = right_tank_.del1_.Tap(left_taps_[0]) //  266
                       + right_tank_.del1_.Tap(left_taps_[1]) // 2974
                       - right_tank_.apf2_.Tap(left_taps_[2]) // 1913
                       + right_tank_.del2_.Tap(left_taps_[3]) // 1996
                       - left_tank_.del1_.Tap(left_taps_[4])  // 1990
                       - left_tank_.apf2_.Tap(left_taps_[5])  //  187
                       - left_tank_.del2_.Tap(left_taps_[6]); // 1066

        float wet_right = left_tank_.del1_.Tap(right_taps_[0])   //  353
                        + left_tank_.del1_.Tap(right_taps_[1])   // 3627
                        - left_tank_.apf2_.Tap(right_taps_[2])   // 1228
                        + left_tank_.del2_.Tap(right_taps_[3])   // 2673
                        - right_tank_.del1_.Tap(right_taps_[4])  // 2111
                        - right_tank_.apf2_.Tap(right_taps_[5])  //  335
                        - right_tank_.del2_.Tap(right_taps_[6]); //  121

        // Mix
        *left_out = dry_left * (1 - mix_) + wet_left * mix_;
        *right_out = dry_right * (1 - mix_) + wet_right * mix_;
    }

    /**
     * @brief 批量处理立体声样本（in-place）。
     * @param l   左声道输入/输出缓冲
     * @param r   右声道输入/输出缓冲
     * @param num 样本帧数
     */
    void Process(float* l, float* r, size_t num) noexcept {
        for (size_t i = 0; i < num; i++)
            Process(l[i], r[i], &l[i], &r[i]);
    }

    /**
     * @brief 复位所有延迟线、滤波器及 PitchShifter 状态。
     */
    void Reset() noexcept {
        predelay_line_.Reset();
        lowpass_.Reset();

        for (auto& d : diffusers_)
            d.Reset();

        left_tank_.Reset();
        right_tank_.Reset();
    }
private:
    float sample_rate_ = 1.0f;

    float mix_ = 0.0f;
    float predelay_ = 0.0f;
    detail::DelayLine predelay_line_;
    detail::OnePoleFilter lowpass_;
    std::array<detail::DelayAllpass, 4> diffusers_;

    detail::Tank left_tank_;
    detail::Tank right_tank_;

    float pitch_shift_ = 12.0f;  // 默认 +1 八度

    static constexpr size_t kNumTaps = 7;
    std::array<float, kNumTaps> base_left_taps_ = {};
    std::array<float, kNumTaps> base_right_taps_ = {};
    std::array<float, kNumTaps> left_taps_ = {};
    std::array<float, kNumTaps> right_taps_ = {};
};
