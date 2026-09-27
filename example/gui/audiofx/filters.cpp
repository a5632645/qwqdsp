#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>

#include "raylib.h"
#include "slider.hpp"

#include "qwqdsp/convert.hpp"
#include "qwqdsp/filter/ladder.hpp"
#include "qwqdsp/filter/ladder_8pole.hpp"
#include "qwqdsp/filter/ota_one_pole.hpp"
#include "qwqdsp/filter/svf_tpt.hpp"
#include "qwqdsp/filter/transpose_sallen_key.hpp"
#include "qwqdsp/oscillator/polyblep.hpp"

static constexpr int kWidth = 620;
static constexpr int kHeight = 400;
static constexpr float kFs = 48000.0f;
static constexpr float kPi = std::numbers::pi_v<float>;

enum FilterType {
    SvfTPT = 0,
    Ladder,
    Ladder8Pole,
    SallenKey,
    OTAOnePole,
    NumFilterTypes
};
static constexpr const char* kFilterNames[]{
    "svf tpt",
    "ladder",
    "ladder 8p",
    "sallen key",
    "ota 1p"
};

enum OutputMode {
    Lowpass = 0,
    Bandpass,
    Highpass,
    NumOutputModes
};
static constexpr const char* kOutputModeNames[]{
    "lp",
    "bp",
    "hp"
};

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
static qwqdsp_filter::OTAOnePole ota_one_pole;

static float scale(float value, float min, float max) noexcept {
    return min + value * (max - min);
}

static float cutoffPitchToW(float value) noexcept {
    float const freq = qwqdsp::convert::Pitch2Freq(value);
    float const safe_freq = std::clamp(freq, 10.0f, kFs * 0.45f);
    return std::clamp(2.0f * kPi * safe_freq / kFs, 0.0001f, kPi - 0.1f);
}

static bool isOutputModeAvailable(FilterType type, OutputMode mode) noexcept {
    if (mode == OutputMode::Lowpass) {
        return true;
    }

    return type == FilterType::SvfTPT
        || type == FilterType::Ladder
        || type == FilterType::SallenKey;
}

static OutputMode effectiveOutputMode(FilterType type, OutputMode mode) noexcept {
    return isOutputModeAvailable(type, mode) ? mode : OutputMode::Lowpass;
}

static void updateFilterCoeff(FilterType type, float w, float res) noexcept {
    switch (type) {
    case FilterType::SvfTPT:
        svf_tpt.SetCoeffQ(w, scale(res, 0.5f, 20.0f));
        break;
    case FilterType::Ladder:
        ladder.Set(w, scale(res, 0.0f, 3.9f), 1.0f);
        break;
    case FilterType::Ladder8Pole:
        ladder_8pole.Set(w, scale(res, 0.0f, 1.8f));
        break;
    case FilterType::SallenKey:
        sallen_key.Set(w, scale(res, 0.0f, 1.95f));
        break;
    case FilterType::OTAOnePole:
        ota_one_pole.Set(w);
        break;
    default:
        break;
    }
}

static float tickFilter(FilterType type, OutputMode mode, float x) noexcept {
    mode = effectiveOutputMode(type, mode);

    switch (type) {
    case FilterType::SvfTPT:
        if (mode == OutputMode::Bandpass) {
            return svf_tpt.TickBandpass(x);
        }
        if (mode == OutputMode::Highpass) {
            return svf_tpt.TickHighpass(x);
        }
        return svf_tpt.TickLowpass(x);
    case FilterType::Ladder:
        if (mode == OutputMode::Bandpass) {
            return ladder.TickBandpass(x);
        }
        if (mode == OutputMode::Highpass) {
            return ladder.TickHighpass(x);
        }
        return ladder.Tick(x);
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
    case FilterType::OTAOnePole:
        return ota_one_pole.Tick(x);
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

    saw.SetFreq(qwqdsp::convert::Pitch2Freq(pitch.load(std::memory_order_relaxed)), kFs);
    updateFilterCoeff(
        type,
        cutoffPitchToW(cutoff_pitch.load(std::memory_order_relaxed)),
        resonance.load(std::memory_order_relaxed));

    for (auto& s : buffer) {
        float const x = saw.Sawtooth() * osc_level.load(std::memory_order_relaxed);
        s.l = tickFilter(type, mode, x) * gain;
        s.r = s.l;
    }
}

static void setupKnob(Knob& knob, Rectangle bound, const char* title) {
    knob.set_bound(bound);
    knob.set_bg_color(BLACK);
    knob.set_fore_color(RAYWHITE);
    knob.set_title(title);
}

static void drawSelector(
    Rectangle bound,
    size_t num_items,
    const char* const* names,
    size_t selected,
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

int main(void) {
    InitWindow(kWidth, kHeight, "saw filter");

    InitAudioDevice();
    SetAudioStreamBufferSizeDefault(512);
    AudioStream stream = LoadAudioStream(48000, 32, 2);
    SetAudioStreamCallback(stream, audioInputCallback);
    PlayAudioStream(stream);

    Rectangle bound;
    bound.x = 0;
    bound.y = 0;
    bound.width = 70;
    bound.height = 60;

    Knob pitch_knob;
    pitch_knob.on_value_change = [](float value) {
        pitch.store(value, std::memory_order_relaxed);
    };
    setupKnob(pitch_knob, bound, "pitch");
    pitch_knob.set_range(0.0f, 127.0f, 0.1f, pitch.load(std::memory_order_relaxed));

    bound.x += bound.width;
    Knob cutoff_knob;
    cutoff_knob.on_value_change = [](float value) {
        cutoff_pitch.store(value, std::memory_order_relaxed);
    };
    setupKnob(cutoff_knob, bound, "cutoff");
    cutoff_knob.set_range(0.0f, 127.0f, 0.1f, cutoff_pitch.load(std::memory_order_relaxed));

    bound.x += bound.width;
    Knob resonance_knob;
    resonance_knob.on_value_change = [](float value) {
        resonance.store(value, std::memory_order_relaxed);
    };
    setupKnob(resonance_knob, bound, "res");
    resonance_knob.set_range(0.0f, 1.0f, 0.001f, resonance.load(std::memory_order_relaxed));

    bound.x += bound.width;
    Knob output_knob;
    output_knob.on_value_change = [](float value) {
        output_gain.store(value, std::memory_order_relaxed);
    };
    setupKnob(output_knob, bound, "output");
    output_knob.set_range(0.0f, 1.0f, 0.001f, output_gain.load(std::memory_order_relaxed));

    bound.x += bound.width;
    Knob osc_level_knob;
    osc_level_knob.on_value_change = [](float value) {
        osc_level.store(value, std::memory_order_relaxed);
    };
    setupKnob(osc_level_knob, bound, "osc lvl");
    osc_level_knob.set_range(0.0f, 4.0f, 0.01f, osc_level.load(std::memory_order_relaxed));

    Rectangle filter_bound;
    filter_bound.x = 0;
    filter_bound.y = 90;
    filter_bound.width = static_cast<float>(GetScreenWidth());
    filter_bound.height = 30;

    Rectangle mode_bound = filter_bound;
    mode_bound.y += mode_bound.height + 8.0f;

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

            auto const mouse_pos = GetMousePosition();

            float const filter_width = filter_bound.width / static_cast<float>(NumFilterTypes);
            for (size_t i = 0; i < static_cast<size_t>(NumFilterTypes); ++i) {
                Rectangle const item{
                    filter_bound.x + static_cast<float>(i) * filter_width,
                    filter_bound.y,
                    filter_width - 2.0f,
                    filter_bound.height
                };
                if (CheckCollisionPointRec(mouse_pos, item) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
                    filter_type.store(static_cast<FilterType>(i), std::memory_order_relaxed);
                }
            }

            auto const current_filter = filter_type.load(std::memory_order_relaxed);
            auto const mode_width = mode_bound.width / static_cast<float>(NumOutputModes);
            for (size_t i = 0; i < static_cast<size_t>(NumOutputModes); ++i) {
                auto const mode = static_cast<OutputMode>(i);
                Rectangle const item{
                    mode_bound.x + static_cast<float>(i) * mode_width,
                    mode_bound.y,
                    mode_width - 2.0f,
                    mode_bound.height
                };
                if (CheckCollisionPointRec(mouse_pos, item)
                    && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)
                    && isOutputModeAvailable(current_filter, mode)) {
                    output_mode.store(mode, std::memory_order_relaxed);
                }
            }

            auto const current_mode = effectiveOutputMode(
                current_filter,
                output_mode.load(std::memory_order_relaxed));
            drawSelector(
                filter_bound,
                static_cast<size_t>(NumFilterTypes),
                kFilterNames,
                static_cast<size_t>(current_filter),
                nullptr);
            drawSelector(
                mode_bound,
                static_cast<size_t>(NumOutputModes),
                kOutputModeNames,
                static_cast<size_t>(current_mode),
                [](size_t i) {
                    return isOutputModeAvailable(
                        filter_type.load(std::memory_order_relaxed),
                        static_cast<OutputMode>(i));
                });

            DrawText("right click knob: reset", 8, kHeight - 22, 12, DARKGRAY);
        }
        EndDrawing();
    }

    UnloadAudioStream(stream);
    CloseAudioDevice();
    CloseWindow();
}
