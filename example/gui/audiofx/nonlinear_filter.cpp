#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <vector>

#include "raylib.h"
#include "slider.hpp"

#include "qwqdsp/convert.hpp"
#include "qwqdsp/fx/oversample.hpp"
#include "qwqdsp/oscillator/polyblep.hpp"

// ------------------------------------------------------------
// copy filter implement (used by nonlinear filter)
// ------------------------------------------------------------

static float Nonlinear(float x) noexcept;

namespace qwqdsp_filter {
class SvfTPT {
public:
    void Reset() noexcept {
        s1_ = 0;
        s2_ = 0;
    }

    /**
     * @param w [0, pi]
     * @param r2 <0:不稳定 =0:无阻尼 0~1:复共轭极点 >1:分裂成两个单极点
     */
    void SetCoeffSVF(float w, float r2) noexcept {
        g_ = std::tan(w / 2);
        R2_ = r2;
        g1_ = r2 + g_;
        d_ = 1 / (1 + r2 * g_ + g_ * g_);
    }

    void SetCoeffQ(float w, float Q) noexcept {
        SetCoeffSVF(w, 1 / Q);
    }

    /**
     * @return [hp, bp, lp]
     */
    std::array<float, 3> TickMultiMode(float x) noexcept {
        float hp = (x - g1_ * s1_ - s2_) * d_;
        float v1 = g_ * hp;
        float bp = v1 + s1_;
        float v2 = g_ * bp;
        float lp = v2 + s2_;
        s1_ = Nonlinear(bp + v1);
        s2_ = lp + v2;
        return {hp, bp, lp};
    }
private:
    float R2_{};
    float g_{};
    float g1_{};
    float d_{};
    float s1_{};
    float s2_{};
};

class Ladder {
public:
    void Reset() noexcept {
        s1_ = 0;
        s2_ = 0;
        s3_ = 0;
        s4_ = 0;
    }

    auto TickMultiMode(float x) noexcept {
        struct Output {
            float hp;
            float lp1;
            float lp2;
            float lp3;
            float lp4;
        };
        Output r;

        float S = g2_ * s2_ + glp_ * (s3_ + s1_ * g2_) + s4_;
        S /= (1 + g_);
        r.hp = (x - k_ * S) / (1 + k_ * g4_);
        r.lp1 = Nonlinear(TickLpTPT(r.hp, s1_, glp_));
        r.lp2 = Nonlinear(TickLpTPT(r.lp1, s2_, glp_));
        r.lp3 = Nonlinear(TickLpTPT(r.lp2, s3_, glp_));
        r.lp4 = Nonlinear(TickLpTPT(r.lp3, s4_, glp_));
        return r;
    }

    /**
     * @param w [0, pi-0.1]
     * @param k [0, 3.99] => [实数极点, 自震荡]
     */
    void Set(float w, float k, float Q) noexcept {
        R_ = 0.5f / Q;
        g_ = std::tan(w / 2);
        glp_ = g_ / (1 + g_);
        g2_ = glp_ * glp_;
        g4_ = g2_ * g2_;
        k_ = k;
    }
private:
    static float TickLpTPT(float x, float& s, float geps) noexcept {
        float const delta = geps * (x - s);
        s += delta;
        float const y = s;
        s += delta;
        return y;
    }

    float R_{};
    float k_{};
    float g_{};
    float g2_{};
    float g4_{};
    float glp_{};
    float s1_{};
    float s2_{};
    float s3_{};
    float s4_{};
};

class Ladder8Pole {
public:
    void Reset() noexcept {
        std::fill_n(s_, 8, 0.0f);
    }

    /**
     * @param k [0,1.884]
     */
    void Set(float w, float k) noexcept {
        g_ = std::tan(w / 2) * (1 + std::numbers::sqrt2_v<float>);
        k_ = k;
    }

    auto Tick(float x) noexcept {
        struct Output {
            float lp1;
            float lp2;
            float lp3;
            float lp4;
            float lp5;
            float lp6;
            float lp7;
            float lp8;
        } output;

        float B = 1 / (1 + g_);
        float A = g_ * B;
        float S = s_[7] + A * (s_[6] + A * (s_[5] + A * (s_[4] + A * (s_[3] + A * (s_[2] + A * (s_[1] + A * s_[0]))))));
        float u = x - B * k_ * S;
        u /= (1 + k_ * A * A * A * A * A * A * A * A);

        output.lp1 = Nonlinear(TickLpTPT(u, s_[0], A));
        output.lp2 = Nonlinear(TickLpTPT(output.lp1, s_[1], A));
        output.lp3 = Nonlinear(TickLpTPT(output.lp2, s_[2], A));
        output.lp4 = Nonlinear(TickLpTPT(output.lp3, s_[3], A));
        output.lp5 = Nonlinear(TickLpTPT(output.lp4, s_[4], A));
        output.lp6 = Nonlinear(TickLpTPT(output.lp5, s_[5], A));
        output.lp7 = Nonlinear(TickLpTPT(output.lp6, s_[6], A));
        output.lp8 = Nonlinear(TickLpTPT(output.lp7, s_[7], A));

        return output;
    }
private:
    static float TickLpTPT(float x, float& s, float geps) noexcept {
        float const delta = geps * (x - s);
        s += delta;
        float const y = s;
        s += delta;
        return y;
    }

    float s_[8]{};
    float g_{};
    float k_{};
};

class TransposeSallenKey {
public:
    void Reset() noexcept {
        s1_ = 0;
        s2_ = 0;
    }

    auto Tick(float x) noexcept {
        struct Output {
            float lp2;
            float bp2;
            float hp2;
        } output;

        float G = (1 + g_);
        float u = x * G * G - G * k_ * s2_ + k_ * s1_;
        u /= (G * G - g_ * k_);
        u = Nonlinear(u);
        float lp1 = TickLpTPT(u, s1_, glp_);
        float hp1 = u - lp1;
        output.lp2 = TickLpTPT(lp1, s2_, glp_);
        output.bp2 = lp1 - output.lp2;
        output.hp2 = hp1 - output.bp2;
        return output;
    }

    /**
     * @param w [0,pi]
     * @param k [0,2] => [实数极点,自震荡]
     */
    void Set(float w, float k) noexcept {
        g_ = std::tan(w / 2);
        glp_ = g_ / (1 + g_);
        k_ = k;
    }
private:
    static float TickLpTPT(float x, float& s, float geps) noexcept {
        float const delta = geps * (x - s);
        s += delta;
        float const y = s;
        s += delta;
        return y;
    }

    float g_{};
    float glp_{};
    float k_{};
    float s1_{};
    float s2_{};
};

class TransposeSallenKeyTrackingHp {
public:
    void Set(float lp_w, float hp_w, float k) noexcept {
        float const lp_g = std::tan(lp_w * 0.5f);
        float const hp_g = std::tan(hp_w * 0.5f);
        lp_coeff_ = lp_g / (1.0f + lp_g);
        hp_lp_coeff_ = hp_g / (1.0f + hp_g);
        float const gain = k / (1.0f + hp_g);
        // 保留代数环分母的余量，避免高 LP 截止频率和低 HP 截止频率下出现奇点。
        feedback_gain_ = gain / std::max(1.0f, gain * lp_coeff_ / 0.9f);
        inverse_ = 1.0f / (1.0f - feedback_gain_ * lp_coeff_);
    }

    auto Tick(float x) noexcept {
        struct Output {
            float lp2;
            float bp2;
            float hp2;
        };

        // 反馈取自第一级 LP，经独立的一阶 HP 返回输入。
        float const lp1_offset = (1.0f - lp_coeff_) * lp1_state_;
        float const u = Nonlinear((x + feedback_gain_ * (lp1_offset - hp_state_)) * inverse_);
        float const lp1 = TickLp(u, lp1_state_, lp_coeff_);
        TickLp(lp1, hp_state_, hp_lp_coeff_);
        float const lp2 = TickLp(lp1, lp2_state_, lp_coeff_);
        float const bp2 = lp1 - lp2;
        float const hp2 = (u - lp1) - bp2;
        return Output{lp2, bp2, hp2};
    }

private:
    static float TickLp(float x, float& state, float coeff) noexcept {
        float const delta = coeff * (x - state);
        float const y = state + delta;
        state = y + delta;
        return y;
    }

    float lp_coeff_{};
    float hp_lp_coeff_{};
    float feedback_gain_{};
    float inverse_{};
    float lp1_state_{};
    float lp2_state_{};
    float hp_state_{};
};

namespace detail {
struct SallenKey4PoleOutput {
    float lp4;
    float bp4;
    float hp4;
};

static inline double TptLp(double x, double state, double coeff) noexcept {
    return state + coeff * (x - state);
}

static inline double TickTptLp(double x, double& state, double coeff) noexcept {
    double const y = TptLp(x, state, coeff);
    state = 2.0 * y - state;
    return y;
}
} // namespace detail

// TODO: looks like that feedback network 1-bp2 and 2-hp2 are not hp_w? 
class TransposeSallenKey4Pole {
public:
    void Set(float w, float zeta) noexcept {
        double const g = std::tan(static_cast<double>(w) * 0.5);
        double const t = std::sqrt(1.0 + zeta);
        bp_gain_ = 4.0 - 2.0 * t;
        hp_gain_ = -2.0 * (t - 1.0) * (t - 1.0);
        coeff_ = g / (1.0 + g);

        double const direct = coeff_ * (1.0 - coeff_);
        inverse_ = 1.0 / (1.0 - bp_gain_ * direct - hp_gain_ * direct * direct);
    }

    detail::SallenKey4PoleOutput Tick(float x) noexcept {
        double const lp1_offset = detail::TptLp(0.0, state_[0], coeff_);
        double const lp2_offset = detail::TptLp(lp1_offset, state_[1], coeff_);
        double const lp3_offset = detail::TptLp(lp2_offset, state_[2], coeff_);
        double const lp4_offset = detail::TptLp(lp3_offset, state_[3], coeff_);
        double const bp_offset = lp1_offset - lp2_offset;
        double const hp_offset = (lp2_offset - lp3_offset) - (lp3_offset - lp4_offset);
        double const feedback = bp_gain_ * bp_offset + hp_gain_ * hp_offset;
        double const u = Nonlinear(static_cast<float>((x + feedback) * inverse_));

        double const lp1 = detail::TickTptLp(u, state_[0], coeff_);
        double const lp2 = detail::TickTptLp(lp1, state_[1], coeff_);
        double const lp3 = detail::TickTptLp(lp2, state_[2], coeff_);
        double const lp4 = detail::TickTptLp(lp3, state_[3], coeff_);
        double const bp4 = (lp2 - lp3) - (lp3 - lp4);
        double const hp4 = (u - lp1) - 3.0 * (lp1 - lp2) + (3.0 * (lp2 - lp3) - (lp3 - lp4));
        return {static_cast<float>(lp4), static_cast<float>(bp4), static_cast<float>(hp4)};
    }

private:
    std::array<double, 4> state_{};
    double bp_gain_{};
    double hp_gain_{};
    double coeff_{};
    double inverse_{};
};

class TransposeSallenKey4PoleTrackingHp {
public:
    void Set(float lp_w, float hp_w, float zeta) noexcept {
        double const lp_g = std::tan(static_cast<double>(lp_w) * 0.5);
        double const hp_g = std::tan(static_cast<double>(hp_w) * 0.5);
        double const hp_ratio = hp_g / lp_g;
        double const t = std::sqrt(1.0 + zeta);
        bp_gain_ = 4.0 - 2.0 * t;
        hp_gain_ = -2.0 * (t - 1.0) * (t - 1.0);
        lp_coeff_ = lp_g / (1.0 + lp_g);
        hp_coeff_ = hp_g / (1.0 + hp_g);

        // u = x + k1·BP2(u) + k2·HP2(LP2(u))，用两条支路的直通量解 ZDF。
        double const bp_direct = hp_coeff_ * (1.0 - hp_coeff_);
        double const hp_direct = lp_coeff_ * (1.0 - hp_coeff_);
        inverse_ = 1.0 / (1.0 - bp_gain_ * bp_direct - hp_gain_ * hp_direct * hp_direct);

        hp_mix_ = 1.0 + 2.0 * hp_ratio;
        hp_ratio_squared_ = hp_ratio * hp_ratio;
    }

    detail::SallenKey4PoleOutput Tick(float x) noexcept {
        double const lp1_offset = detail::TptLp(0.0, lp_state_[0], lp_coeff_);
        double const lp2_offset = detail::TptLp(lp1_offset, lp_state_[1], lp_coeff_);
        double const hp_lp1_offset = detail::TptLp(lp2_offset, hp_state_[0], hp_coeff_);
        double const hp_lp2_offset = detail::TptLp(hp_lp1_offset, hp_state_[1], hp_coeff_);
        double const bp_lp1_offset = detail::TptLp(0.0, bp_state_[0], hp_coeff_);
        double const bp_lp2_offset = detail::TptLp(bp_lp1_offset, bp_state_[1], hp_coeff_);
        double const bp_offset = bp_lp1_offset - bp_lp2_offset;
        double const hp_offset = (lp2_offset - hp_lp1_offset) - (hp_lp1_offset - hp_lp2_offset);
        double const feedback = bp_gain_ * bp_offset + hp_gain_ * hp_offset;
        double const u = Nonlinear(static_cast<float>((x + feedback) * inverse_));
        double const bp_lp1 = detail::TickTptLp(u, bp_state_[0], hp_coeff_);
        detail::TickTptLp(bp_lp1, bp_state_[1], hp_coeff_);
        double const lp1 = detail::TickTptLp(u, lp_state_[0], lp_coeff_);
        double const lp2 = detail::TickTptLp(lp1, lp_state_[1], lp_coeff_);
        double const hp_branch_lp1 = detail::TickTptLp(lp2, hp_state_[0], hp_coeff_);
        double const hp_branch_lp2 = detail::TickTptLp(hp_branch_lp1, hp_state_[1], hp_coeff_);
        double const lp3 = detail::TickTptLp(lp2, lp_state_[2], lp_coeff_);
        double const lp4 = detail::TickTptLp(lp3, lp_state_[3], lp_coeff_);

        double const bp4 = (lp2 - lp3) - (lp3 - lp4);
        // 使用反馈支路的状态合成 HP4，对应模拟分子 p⁴。
        double const hp4 = (u - lp1) - hp_mix_ * (lp1 - lp2)
                         + hp_ratio_squared_ * (3.0 * (lp2 - hp_branch_lp1)
                                              - (hp_branch_lp1 - hp_branch_lp2));
        return {static_cast<float>(lp4), static_cast<float>(bp4), static_cast<float>(hp4)};
    }

private:
    std::array<double, 4> lp_state_{};
    std::array<double, 2> bp_state_{};
    std::array<double, 2> hp_state_{};
    double bp_gain_{};
    double hp_gain_{};
    double hp_mix_{};
    double hp_ratio_squared_{};
    double lp_coeff_{};
    double hp_coeff_{};
    double inverse_{};
};

/**
 * NI Massive "Scream": 低通滤波器 + 内部饱和反馈.
 * 反馈路径把低通输出放大后送进一个高通, 饱和后再与输入相加, 形成强烈的谐波和共振峰.
 * 参考开源复刻 https://github.com/Cure-Audio/Scream
 *  - 低通固定 Q = 1/sqrt(2)
 *  - 高通 Q 随 resonance 从 1/sqrt(2) 变到 sqrt(2)
 *  - 反馈增益随 resonance 从 -12dB 变到 +12dB
 */
class Scream {
public:
    void Reset() noexcept {
        lp_ic1_ = 0;
        lp_ic2_ = 0;
        hp_ic1_ = 0;
        hp_ic2_ = 0;
        feedback_ = 0;
    }

    /**
     * @param lp_w 低通截止频率, [0, pi]
     * @param hp_w 反馈路径高通截止频率, [0, pi]
     * @param res  [0,1]
     */
    void Set(float lp_w, float hp_w, float res) noexcept {
        constexpr float kSqrt2 = std::numbers::sqrt2_v<float>;
        float const lp_q = kSqrt2 * 0.5f;
        float const hp_q = lp_q + res * (kSqrt2 - lp_q);
        float const lp_k = 1.0f / lp_q;
        float const hp_k = 1.0f / hp_q;

        float const lp_g = std::tan(lp_w * 0.5f);
        float const hp_g = std::tan(hp_w * 0.5f);

        lp_a1_ = 1.0f / (1.0f + lp_g * (lp_g + lp_k));
        lp_a2_ = lp_g * lp_a1_;
        lp_a3_ = lp_g * lp_a2_;

        hp_a1_ = 1.0f / (1.0f + hp_g * (hp_g + hp_k));
        hp_a2_ = hp_g * hp_a1_;
        hp_a3_ = hp_g * hp_a2_;
        hp_k_ = hp_k;

        feedback_gain_ = std::pow(10.0f, (res * 24.0f - 12.0f) / 20.0f);
    }

    float Tick(float x) noexcept {
        // 前馈: 输入 + 反馈, 饱和后过低通
        float const u = std::tanh(x + feedback_);
        float const y = ProcessLp(u);

        // 反馈: 低通输出放大 -> 高通 -> 饱和
        feedback_ = std::tanh(ProcessHp(y * feedback_gain_));
        return y;
    }
private:
    float ProcessLp(float x) noexcept {
        float const v3 = x - lp_ic2_;
        float const v1 = lp_a1_ * lp_ic1_ + lp_a2_ * v3;
        float const v2 = lp_ic2_ + lp_a2_ * lp_ic1_ + lp_a3_ * v3;
        lp_ic1_ = 2.0f * v1 - lp_ic1_;
        lp_ic2_ = 2.0f * v2 - lp_ic2_;
        return v2;
    }

    float ProcessHp(float x) noexcept {
        float const v3 = x - hp_ic2_;
        float const v1 = hp_a1_ * hp_ic1_ + hp_a2_ * v3;
        float const v2 = hp_ic2_ + hp_a2_ * hp_ic1_ + hp_a3_ * v3;
        hp_ic1_ = 2.0f * v1 - hp_ic1_;
        hp_ic2_ = 2.0f * v2 - hp_ic2_;
        return x - hp_k_ * v1 - v2;
    }

    float lp_a1_{};
    float lp_a2_{};
    float lp_a3_{};
    float hp_a1_{};
    float hp_a2_{};
    float hp_a3_{};
    float hp_k_{};
    float lp_ic1_{};
    float lp_ic2_{};
    float hp_ic1_{};
    float hp_ic2_{};
    float feedback_{};
    float feedback_gain_{};
};
} // namespace qwqdsp_filter

// ------------------------------------------------------------
// nonlinear
// ------------------------------------------------------------

enum class NonlinearType {
    Tanh,
    HardClip,
    SoftClip,
    ArcTan,
    Cubic,
    NumTypes
};

static constexpr const char* kNonlinearTypeNames[]{"tanh", "hard clip", "soft clip", "atan", "cubic"};

static std::atomic<NonlinearType> g_nonlinear_type{NonlinearType::Tanh};
static std::atomic<float> g_nonlinear_drive{0.0f};
static std::atomic<float> g_nonlinear_bias{0.0f};

static float Nonlinear(float x) noexcept {
    float const drive = g_nonlinear_drive.load(std::memory_order_relaxed);
    float const bias = g_nonlinear_bias.load(std::memory_order_relaxed);
    auto const type = g_nonlinear_type.load(std::memory_order_relaxed);
    float const linear_drive = std::pow(10.0f, drive / 20.0f);
    float y = linear_drive * x + bias;

    switch (type) {
        case NonlinearType::Tanh:
            return std::tanh(y);
        case NonlinearType::HardClip:
            return std::clamp(y, -1.0f, 1.0f);
        case NonlinearType::SoftClip:
            return y / (1.0f + std::abs(y));
        case NonlinearType::ArcTan:
            return std::atan(y);
        case NonlinearType::Cubic:
            return std::clamp(y - y * y * y / 3.0f, -1.0f, 1.0f);
        default:
            return y;
    }
}

// ------------------------------------------------------------
// app
// ------------------------------------------------------------

static constexpr int kWidth = 760;
static constexpr int kHeight = 384;
static constexpr Color kBackground{18, 22, 28, 255};
static constexpr Color kPanel{28, 34, 42, 255};
static constexpr Color kText{231, 237, 245, 255};
static constexpr Color kMutedText{153, 167, 184, 255};
static constexpr float kFs = 48000.0f;
static constexpr float kPi = std::numbers::pi_v<float>;

enum FilterType {
    SvfTPT = 0,
    Ladder,
    Ladder8Pole,
    Scream,
    SallenKey,
    SallenKeyTrackingHp,
    SallenKey4Pole,
    SallenKey4PoleTrackingHp,
    NumFilterTypes
};
static constexpr const char* kFilterNames[]{
    "SVF TPT", "Ladder", "Ladder 8P", "Scream",
    "Sallen-Key", "Sallen-Key HP", "Sallen-Key 4P", "Sallen-Key 4P HP"
};

enum OutputMode {
    Lowpass = 0,
    Bandpass,
    Highpass,
    NumOutputModes
};
static constexpr const char* kOutputModeNames[]{"lp", "bp", "hp"};

using SawOsc = qwqdsp_oscillator::PolyBlep<qwqdsp_oscillator::blep_coeff::BlackmanNutallApprox>;

static std::atomic<FilterType> filter_type{FilterType::SvfTPT};
static std::atomic<OutputMode> output_mode{OutputMode::Lowpass};
static std::atomic<float> pitch{48.0f};
static std::atomic<float> cutoff_pitch{72.0f};
static std::atomic<float> resonance{0.2f};
static std::atomic<float> output_gain{0.35f};
static std::atomic<float> osc_level{0.5f};
static std::atomic<float> scream_amount{0.5f};

static SawOsc saw;
static qwqdsp_filter::SvfTPT svf_tpt;
static qwqdsp_filter::Ladder ladder;
static qwqdsp_filter::Ladder8Pole ladder_8pole;
static qwqdsp_filter::TransposeSallenKey sallen_key;
static qwqdsp_filter::Scream scream;
static qwqdsp_filter::TransposeSallenKeyTrackingHp sallen_key_tracking_hp;
static qwqdsp_filter::TransposeSallenKey4Pole sallen_key_4pole;
static qwqdsp_filter::TransposeSallenKey4PoleTrackingHp sallen_key_4pole_tracking_hp;

static qwqdsp_fx::Oversample g_oversample;

static float scale(float value, float min, float max) noexcept {
    return min + value * (max - min);
}

static bool isOutputModeAvailable(FilterType type, OutputMode mode) noexcept {
    if (mode == OutputMode::Lowpass) {
        return true;
    }

    return type == FilterType::SvfTPT || type == FilterType::Ladder ||
           type == FilterType::SallenKey || type == FilterType::SallenKeyTrackingHp ||
           type == FilterType::SallenKey4Pole || type == FilterType::SallenKey4PoleTrackingHp;
}

static OutputMode effectiveOutputMode(FilterType type, OutputMode mode) noexcept {
    return isOutputModeAvailable(type, mode) ? mode : OutputMode::Lowpass;
}

static float TrackingHpRatio(float amount) noexcept {
    constexpr float kScreamOctaveSpan = 131.5892f / 12.0f;
    return std::exp2((amount - 1.0f) * kScreamOctaveSpan);
}

static void updateFilterCoeff(FilterType type, float w, float res, float scream_amt) noexcept {
    float const zeta = scale(res, 1.0f, -0.5f);
    switch (type) {
        case FilterType::SvfTPT:
            svf_tpt.SetCoeffSVF(w, scale(res, 1.0f, -0.5f));
            break;
        case FilterType::Ladder:
            ladder.Set(w, scale(res, 0.0f, 4.5f), 1.0f);
            break;
        case FilterType::Ladder8Pole:
            ladder_8pole.Set(w, scale(res, 0.0f, 2.5f));
            break;
        case FilterType::SallenKey:
            sallen_key.Set(w, scale(res, 0.0f, 2.6f));
            break;
        case FilterType::Scream: {
            // scream=0 时高通在低通下方约 11 个八度, scream=1 时两者重合
            float const lp_freq = w * kFs / kPi;
            float const ratio = TrackingHpRatio(scream_amt);
            float const hp_freq = std::clamp(lp_freq * ratio, 5.0f, kFs * 0.9f);
            scream.Set(w, qwqdsp::convert::Freq2W(hp_freq, kFs * 2), res);
            break;
        }
        case FilterType::SallenKeyTrackingHp:
            sallen_key_tracking_hp.Set(w, w * TrackingHpRatio(scream_amt), scale(res, 0.0f, 2.6f));
            break;
        case FilterType::SallenKey4Pole:
            sallen_key_4pole.Set(w, zeta);
            break;
        case FilterType::SallenKey4PoleTrackingHp:
            sallen_key_4pole_tracking_hp.Set(w, w * TrackingHpRatio(scream_amt), zeta);
            break;
        default:
            break;
    }
}

static float SelectSallenKeyOutput(const auto& output, OutputMode mode) noexcept {
    if (mode == OutputMode::Bandpass) {
        return output.bp2;
    }
    if (mode == OutputMode::Highpass) {
        return output.hp2;
    }
    return output.lp2;
}

static float SelectSallenKey4PoleOutput(const qwqdsp_filter::detail::SallenKey4PoleOutput& output,
                                       OutputMode mode) noexcept {
    if (mode == OutputMode::Bandpass) {
        return output.bp4;
    }
    if (mode == OutputMode::Highpass) {
        return output.hp4;
    }
    return output.lp4;
}

static float tickFilter(FilterType type, OutputMode mode, float x) noexcept {
    mode = effectiveOutputMode(type, mode);

    switch (type) {
        case FilterType::SvfTPT: {
            auto const r = svf_tpt.TickMultiMode(x);
            if (mode == OutputMode::Bandpass)
                return r[1];
            if (mode == OutputMode::Highpass)
                return r[0];
            return r[2];
        }
        case FilterType::Ladder: {
            auto const r = ladder.TickMultiMode(x);
            if (mode == OutputMode::Bandpass)
                return r.lp2 + r.lp4 - 2.0f * r.lp3;
            if (mode == OutputMode::Highpass)
                return r.hp - 4.0f * (r.lp1 + r.lp3) + r.lp4 + 6.0f * r.lp2;
            return r.lp4;
        }
        case FilterType::Ladder8Pole:
            return ladder_8pole.Tick(x).lp8;
        case FilterType::SallenKey:
            return SelectSallenKeyOutput(sallen_key.Tick(x), mode);
        case FilterType::Scream:
            return scream.Tick(x);
        case FilterType::SallenKeyTrackingHp:
            return SelectSallenKeyOutput(sallen_key_tracking_hp.Tick(x), mode);
        case FilterType::SallenKey4Pole:
            return SelectSallenKey4PoleOutput(sallen_key_4pole.Tick(x), mode);
        case FilterType::SallenKey4PoleTrackingHp:
            return SelectSallenKey4PoleOutput(sallen_key_4pole_tracking_hp.Tick(x), mode);
        default:
            return 0.0f;
    }
}

static void audioInputCallback(void* _buffer, unsigned int frames) {
    struct T {
        float l;
        float r;
    };
    std::span buffer{reinterpret_cast<T*>(_buffer), frames};

    auto const type = filter_type.load(std::memory_order_relaxed);
    auto const mode = output_mode.load(std::memory_order_relaxed);
    auto const gain = output_gain.load(std::memory_order_relaxed);
    auto const osc_lvl = osc_level.load(std::memory_order_relaxed);
    auto const scream_lvl = scream_amount.load(std::memory_order_relaxed);

    // 2x 超采样: 振荡器与滤波器以 96kHz 运行
    saw.SetFreq(qwqdsp::convert::Pitch2Freq(pitch.load(std::memory_order_relaxed)), kFs * 2);
    {
        float const freq = qwqdsp::convert::Pitch2Freq(cutoff_pitch.load(std::memory_order_relaxed));
        float const safe_freq = std::clamp(freq, 10.0f, kFs * 2 * 0.45f);
        float const w = qwqdsp::convert::Freq2W(safe_freq, kFs * 2);
        updateFilterCoeff(type, w, resonance.load(std::memory_order_relaxed), scream_lvl);
    }

    static std::vector<float> up;
    static std::vector<float> down;
    up.resize(frames * 2);
    down.resize(frames);

    for (size_t i = 0; i < up.size(); ++i) {
        float const x = saw.Sawtooth() * osc_lvl;
        up[i] = tickFilter(type, mode, x);
        // up[i] = x;
    }

    // 多相下采样 2→1
    g_oversample.Downsample(up, down);

    for (size_t i = 0; i < down.size(); ++i) {
        buffer[i].l = down[i] * gain;
        buffer[i].r = down[i] * gain;
    }
}

static void setupKnob(Knob& knob, Rectangle bound, const char* title) {
    knob.set_bound(bound);
    knob.set_bg_color(kPanel);
    knob.set_fore_color(kText);
    knob.set_title(title);
}

template <size_t kCount>
static size_t DrawSelector(Rectangle first_row, size_t columns, const char* const (&names)[kCount],
                           size_t selected, auto is_enabled, auto on_select) {
    constexpr float kGap = 8.0f;
    constexpr int kFontSize = 16;
    constexpr Color kSelected{154, 221, 205, 255};
    constexpr Color kHovered{48, 61, 72, 255};
    auto const mouse_pos = GetMousePosition();
    float const width = (first_row.width + kGap) / static_cast<float>(columns) - kGap;
    auto item_bounds = [&](size_t i) {
        return Rectangle{
            first_row.x + (width + kGap) * static_cast<float>(i % columns),
            first_row.y + (first_row.height + kGap) * static_cast<float>(i / columns),
            width,
            first_row.height
        };
    };

    if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
        for (size_t i = 0; i < kCount; ++i) {
            if (is_enabled(i) && CheckCollisionPointRec(mouse_pos, item_bounds(i))) {
                selected = i;
                on_select(i);
                break;
            }
        }
    }
    for (size_t i = 0; i < kCount; ++i) {
        Rectangle const item = item_bounds(i);
        bool const enabled = is_enabled(i);
        bool const hovered = enabled && CheckCollisionPointRec(mouse_pos, item);
        Color const background = i == selected ? kSelected : (hovered ? kHovered : kPanel);
        Color const foreground = i == selected ? kBackground : (enabled ? kText : GRAY);
        DrawRectangleRounded(item, 0.24f, 6, background);
        int const text_width = MeasureText(names[i], kFontSize);
        DrawText(names[i], static_cast<int>(item.x + (item.width - text_width) * 0.5f),
                 static_cast<int>(item.y + (item.height - kFontSize) * 0.5f), kFontSize, foreground);
    }
    return selected;
}

static constexpr std::array<float, 65> kHalfbandCoeffs = {
    -0.0000003856f, -0.0003044707f, 0.0000010085f,  0.0005326581f,  -0.0000017348f, -0.0009739876f, 0.0000027020f,
    0.0016310629f,  -0.0000039780f, -0.0025679886f, 0.0000054613f,  0.0038614894f,  -0.0000072091f, -0.0056054843f,
    0.0000090709f,  0.0079210177f,  -0.0000111919f, -0.0109748173f, 0.0000132931f,  0.0150179768f,  -0.0000152761f,
    -0.0204685818f, 0.0000171944f,  0.0281074181f,  -0.0000190393f, -0.0396091882f, 0.0000203860f,  0.0593530469f,
    -0.0000213413f, -0.1034683564f, 0.0000220592f,  0.3174232112f,  0.4999779612f,  0.3174232112f,  0.0000220592f,
    -0.1034683564f, -0.0000213413f, 0.0593530469f,  0.0000203860f,  -0.0396091882f, -0.0000190393f, 0.0281074181f,
    0.0000171944f,  -0.0204685818f, -0.0000152761f, 0.0150179768f,  0.0000132931f,  -0.0109748173f, -0.0000111919f,
    0.0079210177f,  0.0000090709f,  -0.0056054843f, -0.0000072091f, 0.0038614894f,  0.0000054613f,  -0.0025679886f,
    -0.0000039780f, 0.0016310629f,  0.0000027020f,  -0.0009739876f, -0.0000017348f, 0.0005326581f,  0.0000010085f,
    -0.0003044707f, -0.0000003856f,
};

int main(void) {
    InitWindow(kWidth, kHeight, "nonlinear filter");

    InitAudioDevice();
    SetAudioStreamBufferSizeDefault(512);
    AudioStream stream = LoadAudioStream(48000, 32, 2);
    SetAudioStreamCallback(stream, audioInputCallback);
    PlayAudioStream(stream);

    // 初始化 2x 超采样器（半带 FIR 系数）
    g_oversample.Init(kHalfbandCoeffs, 1);

    constexpr float kMargin = 24.0f;
    constexpr float kContentWidth = kWidth - 2.0f * kMargin;
    constexpr Rectangle kKnobPanel{kMargin, 244.0f, kContentWidth, 104.0f};
    Rectangle bound{kMargin, 258.0f, kContentWidth / 8.0f, 76.0f};

    Knob pitch_knob;
    pitch_knob.on_value_change = [](float value) { pitch.store(value, std::memory_order_relaxed); };
    setupKnob(pitch_knob, bound, "pitch");
    pitch_knob.set_range(0.0f, 127.0f, 0.1f, pitch.load(std::memory_order_relaxed));

    bound.x += bound.width;
    Knob cutoff_knob;
    cutoff_knob.on_value_change = [](float value) { cutoff_pitch.store(value, std::memory_order_relaxed); };
    setupKnob(cutoff_knob, bound, "cutoff");
    cutoff_knob.set_range(0.0f, 127.0f, 0.1f, cutoff_pitch.load(std::memory_order_relaxed));

    bound.x += bound.width;
    Knob resonance_knob;
    resonance_knob.on_value_change = [](float value) { resonance.store(value, std::memory_order_relaxed); };
    setupKnob(resonance_knob, bound, "res");
    resonance_knob.set_range(0.0f, 1.0f, 0.001f, resonance.load(std::memory_order_relaxed));

    bound.x += bound.width;
    Knob scream_knob;
    scream_knob.on_value_change = [](float value) { scream_amount.store(value, std::memory_order_relaxed); };
    setupKnob(scream_knob, bound, "scream");
    scream_knob.set_range(0.0f, 1.0f, 0.001f, scream_amount.load(std::memory_order_relaxed));

    bound.x += bound.width;
    Knob output_knob;
    output_knob.on_value_change = [](float value) { output_gain.store(value, std::memory_order_relaxed); };
    setupKnob(output_knob, bound, "output");
    output_knob.set_range(0.0f, 1.0f, 0.001f, output_gain.load(std::memory_order_relaxed));

    bound.x += bound.width;
    Knob osc_level_knob;
    osc_level_knob.on_value_change = [](float value) { osc_level.store(value, std::memory_order_relaxed); };
    setupKnob(osc_level_knob, bound, "osc lvl");
    osc_level_knob.set_range(0.0f, 4.0f, 0.01f, osc_level.load(std::memory_order_relaxed));

    bound.x += bound.width;
    Knob drive_knob;
    drive_knob.on_value_change = [](float value) { g_nonlinear_drive.store(value, std::memory_order_relaxed); };
    setupKnob(drive_knob, bound, "drive");
    drive_knob.set_range(0.0f, 1.0f, 0.01f, g_nonlinear_drive.load(std::memory_order_relaxed));

    bound.x += bound.width;
    Knob bias_knob;
    bias_knob.on_value_change = [](float value) { g_nonlinear_bias.store(value, std::memory_order_relaxed); };
    setupKnob(bias_knob, bound, "bias");
    bias_knob.set_range(-0.3f, 0.3f, 0.001f, g_nonlinear_bias.load(std::memory_order_relaxed));

    constexpr Rectangle kFilterRow{kMargin, 52.0f, kContentWidth, 42.0f};
    constexpr Rectangle kModeRow{kMargin, 178.0f, 220.0f, 38.0f};
    constexpr Rectangle kNonlinearRow{kMargin + 244.0f, 178.0f, kContentWidth - 244.0f, 38.0f};

    SetTargetFPS(30);
    while (!WindowShouldClose()) {
        BeginDrawing();
        {
            ClearBackground(kBackground);

            DrawText("Filter topology", kMargin, 22, 20, kText);
            DrawText("Output mode", kModeRow.x, 154, 16, kMutedText);
            DrawText("Nonlinearity", kNonlinearRow.x, 154, 16, kMutedText);
            DrawRectangleRounded(kKnobPanel, 0.16f, 6, kPanel);

            pitch_knob.display();
            cutoff_knob.display();
            resonance_knob.display();
            scream_knob.display();
            output_knob.display();
            osc_level_knob.display();
            drive_knob.display();
            bias_knob.display();

            auto const always_enabled = [](size_t) { return true; };
            auto const current_filter = static_cast<FilterType>(DrawSelector(
                kFilterRow, 4, kFilterNames, filter_type.load(std::memory_order_relaxed), always_enabled, [](size_t i) {
                    filter_type.store(static_cast<FilterType>(i), std::memory_order_relaxed);
                }));

            auto const current_mode = effectiveOutputMode(current_filter, output_mode.load(std::memory_order_relaxed));
            DrawSelector(kModeRow, NumOutputModes, kOutputModeNames, current_mode, [current_filter](size_t i) {
                    return isOutputModeAvailable(current_filter, static_cast<OutputMode>(i));
                }, [](size_t i) {
                    output_mode.store(static_cast<OutputMode>(i), std::memory_order_relaxed);
                });

            DrawSelector(kNonlinearRow, static_cast<size_t>(NonlinearType::NumTypes), kNonlinearTypeNames,
                         static_cast<size_t>(g_nonlinear_type.load(std::memory_order_relaxed)),
                         always_enabled, [](size_t i) {
                             g_nonlinear_type.store(static_cast<NonlinearType>(i), std::memory_order_relaxed);
                         });

            DrawText("Right-click a knob to reset", kMargin, 360, 13, kMutedText);
        }
        EndDrawing();
    }

    UnloadAudioStream(stream);
    CloseAudioDevice();
    CloseWindow();
}
