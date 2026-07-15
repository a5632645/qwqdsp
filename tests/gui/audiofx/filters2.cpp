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
// copy filter implement
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

        float S = g2_ * s2_ + glp_ * (s3_ + s1_ * glp_) + s4_;
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
            return y - y * y * y / 3.0f;
        default:
            return y;
    }
}

// ------------------------------------------------------------
// app
// ------------------------------------------------------------

static constexpr int kWidth = 620;
static constexpr int kHeight = 400;
static constexpr float kFs = 48000.0f;
static constexpr float kPi = std::numbers::pi_v<float>;

enum FilterType {
    SvfTPT = 0,
    Ladder,
    Ladder8Pole,
    SallenKey,
    NumFilterTypes
};
static constexpr const char* kFilterNames[]{"svf tpt", "ladder", "ladder 8p", "sallen key"};

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

static SawOsc saw;
static qwqdsp_filter::SvfTPT svf_tpt;
static qwqdsp_filter::Ladder ladder;
static qwqdsp_filter::Ladder8Pole ladder_8pole;
static qwqdsp_filter::TransposeSallenKey sallen_key;

static qwqdsp_fx::Oversample g_oversample;

static float scale(float value, float min, float max) noexcept {
    return min + value * (max - min);
}

static bool isOutputModeAvailable(FilterType type, OutputMode mode) noexcept {
    if (mode == OutputMode::Lowpass) {
        return true;
    }

    return type == FilterType::SvfTPT || type == FilterType::Ladder || type == FilterType::SallenKey;
}

static OutputMode effectiveOutputMode(FilterType type, OutputMode mode) noexcept {
    return isOutputModeAvailable(type, mode) ? mode : OutputMode::Lowpass;
}

static void updateFilterCoeff(FilterType type, float w, float res) noexcept {
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
        default:
            break;
    }
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
        case FilterType::SallenKey: {
            auto const output = sallen_key.Tick(x);
            if (mode == OutputMode::Bandpass) {
                return output.bp2;
            }
            if (mode == OutputMode::Highpass) {
                return output.hp2;
            }
            return output.lp2;
        }
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

    // 2x 超采样: 振荡器与滤波器以 96kHz 运行
    saw.SetFreq(qwqdsp::convert::Pitch2Freq(pitch.load(std::memory_order_relaxed)), kFs * 2);
    {
        float const freq = qwqdsp::convert::Pitch2Freq(cutoff_pitch.load(std::memory_order_relaxed));
        float const safe_freq = std::clamp(freq, 10.0f, kFs * 2 * 0.45f);
        float const w = qwqdsp::convert::Freq2W(safe_freq, kFs * 2);
        updateFilterCoeff(type, w, resonance.load(std::memory_order_relaxed));
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
    knob.set_bg_color(BLACK);
    knob.set_fore_color(RAYWHITE);
    knob.set_title(title);
}

static void drawSelector(Rectangle bound, size_t num_items, const char* const* names, size_t selected,
                         bool (*is_enabled)(size_t)) {
    auto const mouse_pos = GetMousePosition();
    float const each_width = bound.width / static_cast<float>(num_items);

    for (size_t i = 0; i < num_items; ++i) {
        float const x = bound.x + static_cast<float>(i) * each_width;
        float const y = bound.y;
        float const w = each_width - 2.0f;
        float const h = bound.height;
        bool const enabled = is_enabled == nullptr || is_enabled(i);
        Color const fore = enabled ? WHITE : DARKGRAY;

        if (i == selected) {
            DrawRectangle(x, y, w, h, WHITE);
            DrawText(names[i], x + 4.0f, y + 8.0f, 12, BLACK);
        }
        else {
            DrawRectangleLines(x, y, w, h, fore);
            DrawText(names[i], x + 4.0f, y + 8.0f, 12, fore);
        }
    }
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
    InitWindow(kWidth, kHeight, "saw filter");

    InitAudioDevice();
    SetAudioStreamBufferSizeDefault(512);
    AudioStream stream = LoadAudioStream(48000, 32, 2);
    SetAudioStreamCallback(stream, audioInputCallback);
    PlayAudioStream(stream);

    // 初始化 2x 超采样器（半带 FIR 系数）
    g_oversample.Init(kHalfbandCoeffs, 1);

    Rectangle bound;
    bound.x = 0;
    bound.y = 0;
    bound.width = 70;
    bound.height = 60;

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
    bias_knob.set_range(-0.1f, 0.1f, 0.001f, g_nonlinear_bias.load(std::memory_order_relaxed));

    Rectangle filter_bound;
    filter_bound.x = 0;
    filter_bound.y = 90;
    filter_bound.width = static_cast<float>(GetScreenWidth());
    filter_bound.height = 30;

    Rectangle mode_bound = filter_bound;
    mode_bound.y += mode_bound.height + 8.0f;

    Rectangle nonlinear_bound;
    nonlinear_bound.x = 0;
    nonlinear_bound.y = mode_bound.y + mode_bound.height + 8.0f;
    nonlinear_bound.width = static_cast<float>(GetScreenWidth());
    nonlinear_bound.height = 30;

    SetTargetFPS(30);
    while (!WindowShouldClose()) {
        BeginDrawing();
        {
            ClearBackground(BLACK);

            pitch_knob.display();
            cutoff_knob.display();
            resonance_knob.display();
            output_knob.display();
            osc_level_knob.display();
            drive_knob.display();
            bias_knob.display();

            auto const mouse_pos = GetMousePosition();

            float const filter_width = filter_bound.width / static_cast<float>(NumFilterTypes);
            for (size_t i = 0; i < static_cast<size_t>(NumFilterTypes); ++i) {
                Rectangle const item{filter_bound.x + static_cast<float>(i) * filter_width, filter_bound.y,
                                     filter_width - 2.0f, filter_bound.height};
                if (CheckCollisionPointRec(mouse_pos, item) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
                    filter_type.store(static_cast<FilterType>(i), std::memory_order_relaxed);
                }
            }

            auto const current_filter = filter_type.load(std::memory_order_relaxed);
            auto const mode_width = mode_bound.width / static_cast<float>(NumOutputModes);
            for (size_t i = 0; i < static_cast<size_t>(NumOutputModes); ++i) {
                auto const mode = static_cast<OutputMode>(i);
                Rectangle const item{mode_bound.x + static_cast<float>(i) * mode_width, mode_bound.y, mode_width - 2.0f,
                                     mode_bound.height};
                if (CheckCollisionPointRec(mouse_pos, item) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)
                    && isOutputModeAvailable(current_filter, mode)) {
                    output_mode.store(mode, std::memory_order_relaxed);
                }
            }

            auto const current_mode = effectiveOutputMode(current_filter, output_mode.load(std::memory_order_relaxed));
            drawSelector(filter_bound, static_cast<size_t>(NumFilterTypes), kFilterNames,
                         static_cast<size_t>(current_filter), nullptr);
            drawSelector(mode_bound, static_cast<size_t>(NumOutputModes), kOutputModeNames,
                         static_cast<size_t>(current_mode), [](size_t i) {
                             return isOutputModeAvailable(filter_type.load(std::memory_order_relaxed),
                                                          static_cast<OutputMode>(i));
                         });

            auto const current_nl_type = static_cast<size_t>(g_nonlinear_type.load(std::memory_order_relaxed));
            float const nl_width = nonlinear_bound.width / static_cast<float>(NonlinearType::NumTypes);
            for (size_t i = 0; i < static_cast<size_t>(NonlinearType::NumTypes); ++i) {
                Rectangle const item{nonlinear_bound.x + static_cast<float>(i) * nl_width, nonlinear_bound.y,
                                     nl_width - 2.0f, nonlinear_bound.height};
                if (CheckCollisionPointRec(mouse_pos, item) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
                    g_nonlinear_type.store(static_cast<NonlinearType>(i), std::memory_order_relaxed);
                }
            }
            drawSelector(nonlinear_bound, static_cast<size_t>(NonlinearType::NumTypes), kNonlinearTypeNames,
                         current_nl_type, nullptr);

            DrawText("right click knob: reset", 8, kHeight - 22, 12, DARKGRAY);
        }
        EndDrawing();
    }

    UnloadAudioStream(stream);
    CloseAudioDevice();
    CloseWindow();
}
