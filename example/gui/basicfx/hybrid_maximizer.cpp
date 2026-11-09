// ------------------------------------------------------------
//  Hybrid Maximizer — feedback limiter + feedforward brickwall
//
//  feedback 级没有 look-ahead，会 overshoot；
//  后面串一个 feedforward look-ahead limiter 保证最终 ceiling。
// ------------------------------------------------------------
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <format>
#include <functional>
#include <iostream>
#include <string>
#include <utility>

#include "miniaudio.h"
#include "raylib.h"
#include "slider.hpp"

#include "hybrid_maximizer_dsp.hpp"

namespace hybrid_maximizer_demo {

constexpr int kWindowWidth = 780;
constexpr int kWindowHeight = 460;
constexpr float kSampleRate = 48000.0f;
constexpr float kMeterFloorDb = -60.0f;
constexpr float kMeterCeilDb = 0.0f;

static std::atomic<float> s_input_gain_db{6.0f};
static std::atomic<float> s_threshold_db{-6.0f};
static std::atomic<float> s_ratio{20.0f};
static std::atomic<float> s_attack_ms{1.0f};
static std::atomic<float> s_release_ms{100.0f};
static std::atomic<float> s_ceiling_db{-0.1f};
static std::atomic<float> s_lookahead_ms{5.0f};

static std::atomic<float> s_input_peak{0.0f};
static std::atomic<float> s_output_peak{0.0f};
static std::atomic<float> s_feedback_gain_reduction_db{0.0f};
static std::atomic<float> s_final_gain_reduction_db{0.0f};

static HybridMaximizer s_processor;

// ------------------------------------------------------------
//  增益历史曲线（UI 线程按帧采样）
// ------------------------------------------------------------
static constexpr size_t kHistorySize = 240;

struct HistoryPlot {
    std::array<float, kHistorySize> feedback{};
    std::array<float, kHistorySize> final{};
    size_t pos{};

    void Push(float feedback_db, float final_db) noexcept {
        feedback[pos] = feedback_db;
        final[pos] = final_db;
        pos = (pos + 1) % kHistorySize;
    }

    float LastFeedback() const noexcept {
        return feedback[(pos + kHistorySize - 1) % kHistorySize];
    }

    float LastFinal() const noexcept {
        return final[(pos + kHistorySize - 1) % kHistorySize];
    }

    void Draw(float x, float y, float width, float height) const {
        constexpr float kRangeDb = 24.0f;

        DrawRectangle(static_cast<int>(x), static_cast<int>(y),
                      static_cast<int>(width), static_cast<int>(height),
                      Color{20, 22, 25, 255});
        DrawRectangleLinesEx(Rectangle{x, y, width, height}, 1.0f, Color{60, 60, 70, 255});

        for (int db = 0; db >= -24; db -= 6) {
            float const gy = y + (-static_cast<float>(db) / kRangeDb) * height;
            DrawLineEx(Vector2{x, gy}, Vector2{x + width, gy}, 1.0f, Color{40, 43, 48, 255});
            std::string const label = std::format("{}", db);
            DrawText(label.c_str(), static_cast<int>(x) + 4, static_cast<int>(gy) + 2, 9,
                     Color{110, 115, 122, 255});
        }

        auto draw_curve = [&](std::array<float, kHistorySize> const& data, Color color) {
            for (size_t i = 1; i < kHistorySize; ++i) {
                float const v0 = data[(pos + i - 1) % kHistorySize];
                float const v1 = data[(pos + i) % kHistorySize];
                float const x0 = x + static_cast<float>(i - 1) * width / static_cast<float>(kHistorySize - 1);
                float const x1 = x + static_cast<float>(i) * width / static_cast<float>(kHistorySize - 1);
                float const y0 = y + std::clamp(-v0 / kRangeDb, 0.0f, 1.0f) * height;
                float const y1 = y + std::clamp(-v1 / kRangeDb, 0.0f, 1.0f) * height;
                DrawLineEx(Vector2{x0, y0}, Vector2{x1, y1}, 1.5f, color);
            }
        };

        draw_curve(feedback, Color{225, 172, 74, 255});
        draw_curve(final, Color{220, 226, 235, 255});

        DrawLineEx(Vector2{x + 8, y + 12}, Vector2{x + 24, y + 12}, 2.0f, Color{225, 172, 74, 255});
        DrawText("FEEDBACK", static_cast<int>(x) + 30, static_cast<int>(y) + 6, 9,
                 Color{180, 185, 192, 255});
        DrawLineEx(Vector2{x + 110, y + 12}, Vector2{x + 126, y + 12}, 2.0f, Color{220, 226, 235, 255});
        DrawText("FINAL", static_cast<int>(x) + 132, static_cast<int>(y) + 6, 9,
                 Color{180, 185, 192, 255});

        std::string const current = std::format("FB {:.1f} dB   FINAL {:.1f} dB",
                                                LastFeedback(), LastFinal());
        int const text_width = MeasureText(current.c_str(), 9);
        DrawText(current.c_str(), static_cast<int>(x + width - static_cast<float>(text_width) - 8.0f),
                 static_cast<int>(y) + 6, 9, Color{160, 165, 172, 255});
    }
};

static HistoryPlot s_history;

static std::array<Knob, 7> s_knobs;

constexpr float kKnobX = 20.0f;
constexpr float kKnobY = 76.0f;
constexpr float kKnobWidth = 96.0f;
constexpr float kKnobHeight = 130.0f;
constexpr float kKnobGap = 6.0f;

static void SetupKnob(Knob& knob, size_t index, const char* title, float min, float max, float step,
                      float value, int sensitivity, std::function<std::string(float)> format,
                      std::function<void(float)> on_change) {
    knob.set_bound(
                static_cast<int>(kKnobX + static_cast<float>(index) * (kKnobWidth + kKnobGap)),
                static_cast<int>(kKnobY),
                static_cast<int>(kKnobWidth),
                static_cast<int>(kKnobHeight))
        .set_title(title)
        .set_range(min, max, step, value)
        .set_value(value)
        .set_sensitivity(sensitivity)
        .set_name_font_size(12)
        .set_number_font_size(11)
        .set_fore_color(Color{235, 235, 238, 255})
        .set_bg_color(Color{24, 26, 29, 255});
    knob.value_to_text_function = std::move(format);
    knob.on_value_change = std::move(on_change);
}

static void CreateKnobs() {
    SetupKnob(s_knobs[0], 0, "Input", 0.0f, 24.0f, 0.5f, s_input_gain_db.load(),
              2, [](float v) { return std::format("{:+.1f} dB", v); },
              [](float v) { s_input_gain_db.store(v, std::memory_order_relaxed); });

    SetupKnob(s_knobs[1], 1, "Threshold", -48.0f, 0.0f, 0.5f, s_threshold_db.load(),
              2, [](float v) { return std::format("{:+.1f} dB", v); },
              [](float v) { s_threshold_db.store(v, std::memory_order_relaxed); });

    SetupKnob(s_knobs[2], 2, "Ratio", 1.0f, 50.0f, 1.0f, s_ratio.load(),
              2, [](float v) { return std::format("{:.0f}:1", v); },
              [](float v) { s_ratio.store(v, std::memory_order_relaxed); });

    SetupKnob(s_knobs[3], 3, "Attack", 0.1f, 20.0f, 0.1f, s_attack_ms.load(),
              2, [](float v) { return std::format("{:.1f} ms", v); },
              [](float v) { s_attack_ms.store(v, std::memory_order_relaxed); });

    SetupKnob(s_knobs[4], 4, "Release", 10.0f, 1000.0f, 5.0f, s_release_ms.load(),
              2, [](float v) { return std::format("{:.0f} ms", v); },
              [](float v) { s_release_ms.store(v, std::memory_order_relaxed); });

    SetupKnob(s_knobs[5], 5, "Ceiling", -12.0f, 0.0f, 0.1f, s_ceiling_db.load(),
              1, [](float v) { return std::format("{:+.1f} dB", v); },
              [](float v) { s_ceiling_db.store(v, std::memory_order_relaxed); });

    SetupKnob(s_knobs[6], 6, "Look-ahead", 0.5f, 10.0f, 0.1f, s_lookahead_ms.load(),
              2, [](float v) { return std::format("{:.1f} ms", v); },
              [](float v) { s_lookahead_ms.store(v, std::memory_order_relaxed); });

    s_input_gain_db.store(s_knobs[0].get_value(), std::memory_order_relaxed);
    s_threshold_db.store(s_knobs[1].get_value(), std::memory_order_relaxed);
    s_ratio.store(s_knobs[2].get_value(), std::memory_order_relaxed);
    s_attack_ms.store(s_knobs[3].get_value(), std::memory_order_relaxed);
    s_release_ms.store(s_knobs[4].get_value(), std::memory_order_relaxed);
    s_ceiling_db.store(s_knobs[5].get_value(), std::memory_order_relaxed);
    s_lookahead_ms.store(s_knobs[6].get_value(), std::memory_order_relaxed);
}

static float AmplitudeToNorm(float amplitude) noexcept {
    float const db = 20.0f * std::log10(std::max(amplitude, 1.0e-6f));
    return std::clamp((db - kMeterFloorDb) / (kMeterCeilDb - kMeterFloorDb), 0.0f, 1.0f);
}

static void DrawVerticalMeter(float x, float y, float width, float height, float value,
                              const char* label, Color fill_color) {
    DrawRectangle(static_cast<int>(x), static_cast<int>(y), static_cast<int>(width), static_cast<int>(height),
                  Color{25, 25, 30, 255});
    DrawRectangleLines(static_cast<int>(x), static_cast<int>(y), static_cast<int>(width), static_cast<int>(height),
                       Color{60, 60, 70, 255});

    float const norm = AmplitudeToNorm(value);
    if (norm > 0.0f) {
        int const fill_height = static_cast<int>(norm * (height - 2.0f));
        DrawRectangle(static_cast<int>(x) + 1, static_cast<int>(y + height - 1.0f) - fill_height,
                      static_cast<int>(width) - 2, fill_height, fill_color);
    }

    int const label_width = MeasureText(label, 10);
    DrawText(label, static_cast<int>(x + (width - static_cast<float>(label_width)) * 0.5f),
             static_cast<int>(y + height + 6.0f), 10, Color{140, 140, 150, 255});
}

} // namespace hybrid_maximizer_demo

// ------------------------------------------------------------
//  miniaudio 回调（双声道入 → 双声道出）
// ------------------------------------------------------------
extern "C" void MaCallback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;

    auto* input = static_cast<const float*>(pInput);
    auto* output = static_cast<float*>(pOutput);

    if (input == nullptr || output == nullptr) {
        return;
    }

    using namespace hybrid_maximizer_demo;

    s_processor.SetParameters(
        s_input_gain_db.load(std::memory_order_relaxed),
        s_threshold_db.load(std::memory_order_relaxed),
        s_ratio.load(std::memory_order_relaxed),
        s_attack_ms.load(std::memory_order_relaxed),
        s_release_ms.load(std::memory_order_relaxed),
        s_ceiling_db.load(std::memory_order_relaxed),
        s_lookahead_ms.load(std::memory_order_relaxed));

    std::copy_n(input, static_cast<size_t>(frameCount) * 2, output);
    s_processor.Process(output, static_cast<size_t>(frameCount));

    float input_peak = 0.0f;
    float output_peak = 0.0f;
    for (size_t i = 0; i < static_cast<size_t>(frameCount) * 2; ++i) {
        input_peak = std::max(input_peak, std::abs(input[i]));
        output_peak = std::max(output_peak, std::abs(output[i]));
    }

    s_input_peak.store(input_peak, std::memory_order_relaxed);
    s_output_peak.store(output_peak, std::memory_order_relaxed);
    s_feedback_gain_reduction_db.store(s_processor.GetFeedbackGainReductionDb(),
                                       std::memory_order_relaxed);
    s_final_gain_reduction_db.store(s_processor.GetFinalGainReductionDb(),
                                    std::memory_order_relaxed);
}

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------
int main(void) {
    using namespace hybrid_maximizer_demo;

    s_processor.Init(kSampleRate);
    CreateKnobs();

    // ── miniaudio full-duplex（双声道入 → 双声道出）──
    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 2;
    config.playback.format = ma_format_f32;
    config.playback.channels = 2;
    config.sampleRate = static_cast<ma_uint32>(kSampleRate);
    config.dataCallback = MaCallback;
    config.pUserData = nullptr;
    config.periodSizeInMilliseconds = 10;

    ma_device device;
    ma_result result = ma_device_init(nullptr, &config, &device);
    if (result != MA_SUCCESS) {
        std::cerr << "failed to initialize audio device\n";
        return 1;
    }
    result = ma_device_start(&device);
    if (result != MA_SUCCESS) {
        std::cerr << "failed to start audio device\n";
        ma_device_uninit(&device);
        return 1;
    }

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, "Hybrid Maximizer");
    SetTargetFPS(60);

    float display_input = 0.0f;
    float display_output = 0.0f;

    while (!WindowShouldClose()) {
        constexpr float kDecay = 0.86f;
        display_input = std::max(s_input_peak.load(std::memory_order_relaxed), display_input * kDecay);
        display_output = std::max(s_output_peak.load(std::memory_order_relaxed), display_output * kDecay);
        s_history.Push(s_feedback_gain_reduction_db.load(std::memory_order_relaxed),
                       s_final_gain_reduction_db.load(std::memory_order_relaxed));

        BeginDrawing();
        ClearBackground(Color{18, 20, 23, 255});

        // ── 顶栏 ──
        DrawRectangle(0, 0, kWindowWidth, 64, Color{24, 26, 29, 255});
        DrawRectangle(0, 63, kWindowWidth, 1, Color{48, 52, 57, 255});
        DrawText("HYBRID MAXIMIZER", 20, 14, 24, Color{232, 233, 235, 255});
        DrawText("FEEDBACK LIMITER -> FEEDFORWARD BRICKWALL LIMITER", 21, 43, 10,
                 Color{126, 132, 139, 255});

        // ── 旋钮 ──
        for (auto& knob : s_knobs) {
            knob.display();
        }

        // ── 电平表 ──
        constexpr float kMeterY = 250.0f;
        constexpr float kMeterWidth = 28.0f;
        constexpr float kMeterHeight = 130.0f;
        DrawVerticalMeter(24.0f, kMeterY, kMeterWidth, kMeterHeight, display_input, "IN",
                          Color{76, 201, 151, 255});
        DrawVerticalMeter(72.0f, kMeterY, kMeterWidth, kMeterHeight, display_output, "OUT",
                          Color{76, 201, 151, 255});

        DrawText("GAIN HISTORY", 150, 232, 10, Color{140, 140, 150, 255});
        s_history.Draw(150.0f, 250.0f, 600.0f, 140.0f);

        // ── 底部提示 ──
        DrawText("drag knob: adjust   shift/alt: fine   right click: reset", 20,
                 kWindowHeight - 20, 10, Color{80, 80, 90, 255});

        EndDrawing();
    }

    ma_device_stop(&device);
    ma_device_uninit(&device);
    CloseWindow();
    return 0;
}
