// ------------------------------------------------------------
//  Multiband Maximizer — 经典三段多段最大化器
//
//  120 Hz / 2.5 kHz LR4 分频
//  每段 Drive -> look-ahead limiter -> 段 Ceiling
//  三段相加后再过一个最终 look-ahead limiter，保证总输出 ceiling。
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

#include "multiband_maximizer_dsp.hpp"

namespace multiband_maximizer_demo {

constexpr int kWindowWidth = 780;
constexpr int kWindowHeight = 640;
constexpr float kSampleRate = 48000.0f;
constexpr float kMeterFloorDb = -60.0f;
constexpr float kMeterCeilDb = 0.0f;

static constexpr std::array<const char*, kNumBands> kBandNames{
    "LOW",
    "MID",
    "HIGH",
};

static constexpr std::array<const char*, kNumBands> kBandRanges{
    "< 120 Hz",
    "120 Hz - 2.5 kHz",
    "> 2.5 kHz",
};

static constexpr std::array<Color, kNumBands> kBandColors{{
    {76, 201, 151, 255},
    {64, 170, 220, 255},
    {225, 172, 74, 255},
}};

static constexpr std::array<float, kNumBands> kDefaultDriveDb{6.0f, 6.0f, 6.0f};
static constexpr std::array<float, kNumBands> kDefaultBandCeilingDb{-0.1f, -0.1f, -0.1f};

// ------------------------------------------------------------
//  参数（UI 线程写入，音频线程读取）
// ------------------------------------------------------------
static std::array<std::atomic<float>, kNumBands> s_drive_db;
static std::array<std::atomic<float>, kNumBands> s_band_ceiling_db;
static std::atomic<float> s_release_ms{100.0f};
static std::atomic<float> s_lookahead_ms{5.0f};
static std::atomic<float> s_output_ceiling_db{-0.1f};

// ------------------------------------------------------------
//  电平表（音频线程写入，UI 线程读取）
// ------------------------------------------------------------
static std::atomic<float> s_input_peak{0.0f};
static std::atomic<float> s_output_peak{0.0f};
static std::array<std::atomic<float>, kNumBands> s_band_gain_reduction_db;
static std::atomic<float> s_final_gain_reduction_db{0.0f};

static MultibandMaximizer s_processor;

// ------------------------------------------------------------
//  控件布局
// ------------------------------------------------------------
static std::array<std::array<Knob, 2>, kNumBands> s_band_knobs;
static Knob s_release_knob;
static Knob s_lookahead_knob;
static Knob s_output_knob;

constexpr float kGlobalKnobX = 20.0f;
constexpr float kGlobalKnobY = 76.0f;
constexpr float kGlobalKnobWidth = 96.0f;
constexpr float kGlobalKnobHeight = 112.0f;
constexpr float kGlobalKnobGap = 6.0f;

constexpr float kBandKnobX = 100.0f;
constexpr float kBandKnobY = 210.0f;
constexpr float kBandKnobWidth = 96.0f;
constexpr float kBandKnobHeight = 112.0f;
constexpr float kBandKnobGap = 6.0f;
constexpr float kBandRowStep = 122.0f;

static void SetupKnob(Knob& knob, float x, float y, float width, float height, const char* title,
                      float min, float max, float step, float value, int sensitivity,
                      std::function<std::string(float)> format, std::function<void(float)> on_change) {
    knob.set_bound(static_cast<int>(x), static_cast<int>(y), static_cast<int>(width),
                   static_cast<int>(height))
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
    for (size_t band = 0; band < kNumBands; ++band) {
        float const y = kBandKnobY + static_cast<float>(band) * kBandRowStep;

        SetupKnob(s_band_knobs[band][0], kBandKnobX, y, kBandKnobWidth, kBandKnobHeight,
                  "Drive", 0.0f, 24.0f, 0.5f, kDefaultDriveDb[band], 2,
                  [](float v) { return std::format("{:+.1f} dB", v); },
                  [band](float v) { s_drive_db[band].store(v, std::memory_order_relaxed); });

        SetupKnob(s_band_knobs[band][1], kBandKnobX + kBandKnobWidth + kBandKnobGap, y,
                  kBandKnobWidth, kBandKnobHeight,
                  "Ceiling", -12.0f, 0.0f, 0.1f, kDefaultBandCeilingDb[band], 1,
                  [](float v) { return std::format("{:+.1f} dB", v); },
                  [band](float v) {
                      s_band_ceiling_db[band].store(v, std::memory_order_relaxed);
                  });

        s_drive_db[band].store(kDefaultDriveDb[band], std::memory_order_relaxed);
        s_band_ceiling_db[band].store(kDefaultBandCeilingDb[band], std::memory_order_relaxed);
    }

    SetupKnob(s_release_knob, kGlobalKnobX, kGlobalKnobY, kGlobalKnobWidth, kGlobalKnobHeight,
              "Release", 10.0f, 1000.0f, 5.0f, 100.0f, 2,
              [](float v) { return std::format("{:.0f} ms", v); },
              [](float v) { s_release_ms.store(v, std::memory_order_relaxed); });

    SetupKnob(s_lookahead_knob, kGlobalKnobX + kGlobalKnobWidth + kGlobalKnobGap, kGlobalKnobY,
              kGlobalKnobWidth, kGlobalKnobHeight,
              "Look-ahead", 0.5f, 10.0f, 0.1f, 5.0f, 2,
              [](float v) { return std::format("{:.1f} ms", v); },
              [](float v) { s_lookahead_ms.store(v, std::memory_order_relaxed); });

    SetupKnob(s_output_knob, kGlobalKnobX + 2.0f * (kGlobalKnobWidth + kGlobalKnobGap), kGlobalKnobY,
              kGlobalKnobWidth, kGlobalKnobHeight,
              "Output", -12.0f, 0.0f, 0.1f, -0.1f, 1,
              [](float v) { return std::format("{:+.1f} dB", v); },
              [](float v) { s_output_ceiling_db.store(v, std::memory_order_relaxed); });

    s_release_ms.store(s_release_knob.get_value(), std::memory_order_relaxed);
    s_lookahead_ms.store(s_lookahead_knob.get_value(), std::memory_order_relaxed);
    s_output_ceiling_db.store(s_output_knob.get_value(), std::memory_order_relaxed);
}

// ------------------------------------------------------------
//  电平表与 band 绘制
// ------------------------------------------------------------
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

static void DrawGainReductionMeter(float x, float y, float width, float height,
                                   float gain_reduction_db, const char* label, Color color) {
    constexpr float kMaxReductionDb = 24.0f;

    DrawText(label, static_cast<int>(x), static_cast<int>(y - 18.0f), 10,
             Color{140, 140, 150, 255});

    DrawRectangle(static_cast<int>(x), static_cast<int>(y), static_cast<int>(width), static_cast<int>(height),
                  Color{25, 25, 30, 255});
    DrawRectangleLines(static_cast<int>(x), static_cast<int>(y), static_cast<int>(width), static_cast<int>(height),
                       Color{60, 60, 70, 255});

    float const norm = std::clamp(-gain_reduction_db / kMaxReductionDb, 0.0f, 1.0f);
    if (norm > 0.0f) {
        int const fill_width = static_cast<int>(norm * (width - 2.0f));
        DrawRectangle(static_cast<int>(x) + 1, static_cast<int>(y) + 1, fill_width,
                      static_cast<int>(height) - 2, color);
    }

    std::string const text = std::format("{:.1f} dB", gain_reduction_db);
    int const text_width = MeasureText(text.c_str(), 11);
    DrawText(text.c_str(), static_cast<int>(x + width - static_cast<float>(text_width) - 8.0f),
             static_cast<int>(y + (height - 11.0f) * 0.5f), 11, Color{235, 235, 238, 255});
}

static void DrawBandRow(size_t band, float display_gain_reduction_db) {
    float const y = kBandKnobY + static_cast<float>(band) * kBandRowStep;
    float const panel_y = y - 8.0f;
    constexpr float kPanelHeight = 120.0f;

    DrawRectangle(10, static_cast<int>(panel_y), kWindowWidth - 20, static_cast<int>(kPanelHeight),
                  Color{22, 24, 27, 255});
    DrawRectangleLinesEx(Rectangle{10.0f, panel_y, static_cast<float>(kWindowWidth - 20), kPanelHeight},
                         1.0f, Color{38, 41, 46, 255});

    DrawText(kBandNames[band], 20, static_cast<int>(y + 32.0f), 15, kBandColors[band]);
    DrawText(kBandRanges[band], 20, static_cast<int>(y + 56.0f), 9, Color{130, 135, 142, 255});

    for (auto& knob : s_band_knobs[band]) {
        knob.display();
    }

    DrawGainReductionMeter(330.0f, y + 46.0f, 430.0f, 18.0f, display_gain_reduction_db,
                           "BAND GR", kBandColors[band]);
}

} // namespace multiband_maximizer_demo

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

    using namespace multiband_maximizer_demo;

    std::array<float, kNumBands> drives;
    std::array<float, kNumBands> ceilings;
    for (size_t band = 0; band < kNumBands; ++band) {
        drives[band] = s_drive_db[band].load(std::memory_order_relaxed);
        ceilings[band] = s_band_ceiling_db[band].load(std::memory_order_relaxed);
    }

    s_processor.SetParameters(
        drives, ceilings,
        s_release_ms.load(std::memory_order_relaxed),
        s_lookahead_ms.load(std::memory_order_relaxed),
        s_output_ceiling_db.load(std::memory_order_relaxed));

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
    for (size_t band = 0; band < kNumBands; ++band) {
        s_band_gain_reduction_db[band].store(
            s_processor.GetBandGainReductionDb(band), std::memory_order_relaxed);
    }
    s_final_gain_reduction_db.store(s_processor.GetFinalGainReductionDb(),
                                    std::memory_order_relaxed);
}

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------
int main(void) {
    using namespace multiband_maximizer_demo;

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
    InitWindow(kWindowWidth, kWindowHeight, "Multiband Maximizer");
    SetTargetFPS(60);

    float display_input = 0.0f;
    float display_output = 0.0f;
    std::array<float, kNumBands> display_band_gain_reduction{};
    float display_final_gain_reduction = 0.0f;

    while (!WindowShouldClose()) {
        constexpr float kDecay = 0.86f;
        display_input = std::max(s_input_peak.load(std::memory_order_relaxed), display_input * kDecay);
        display_output = std::max(s_output_peak.load(std::memory_order_relaxed), display_output * kDecay);
        for (size_t band = 0; band < kNumBands; ++band) {
            display_band_gain_reduction[band] = std::min(
                s_band_gain_reduction_db[band].load(std::memory_order_relaxed),
                display_band_gain_reduction[band] * kDecay);
        }
        display_final_gain_reduction = std::min(
            s_final_gain_reduction_db.load(std::memory_order_relaxed),
            display_final_gain_reduction * kDecay);

        BeginDrawing();
        ClearBackground(Color{18, 20, 23, 255});

        // ── 顶栏 ──
        DrawRectangle(0, 0, kWindowWidth, 64, Color{24, 26, 29, 255});
        DrawRectangle(0, 63, kWindowWidth, 1, Color{48, 52, 57, 255});
        DrawText("MULTIBAND MAXIMIZER", 20, 14, 24, Color{232, 233, 235, 255});
        DrawText("3-BAND LR4 LIMITERS -> SUM -> FINAL LIMITER", 21, 43, 10,
                 Color{126, 132, 139, 255});

        // ── 全局参数与主电平表 ──
        s_release_knob.display();
        s_lookahead_knob.display();
        s_output_knob.display();

        DrawText("MASTER", 350, 76, 9, Color{140, 140, 150, 255});
        DrawVerticalMeter(350.0f, 90.0f, 24.0f, 88.0f, display_input, "IN", Color{76, 201, 151, 255});
        DrawVerticalMeter(386.0f, 90.0f, 24.0f, 88.0f, display_output, "OUT", Color{76, 201, 151, 255});

        DrawText("CROSSOVERS", 450, 88, 9, Color{140, 140, 150, 255});
        DrawText("120 Hz  /  2.5 kHz", 450, 104, 12, Color{220, 223, 228, 255});
        DrawText("3 BAND LIMITERS + FINAL LIMITER", 450, 124, 9, Color{110, 115, 122, 255});
        DrawGainReductionMeter(450.0f, 162.0f, 300.0f, 18.0f, display_final_gain_reduction,
                               "FINAL GR", Color{235, 235, 238, 255});

        // ── 三个 band ──
        for (size_t band = 0; band < kNumBands; ++band) {
            DrawBandRow(band, display_band_gain_reduction[band]);
        }

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
