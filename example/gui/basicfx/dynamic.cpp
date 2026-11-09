// ------------------------------------------------------------
//  Dynamics — 双声道 real-time gate / expander / compressor
//
//  feed-forward、双声道联动 block 峰值检测：
//    block 内先取峰值，再按当前模式计算目标增益（hard knee），
//    最后用 attack/release 平滑逐样本施加增益。
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

#include "dynamic_dsp.hpp"

namespace dynamic_demo {

constexpr int kWindowWidth = 740;
constexpr int kWindowHeight = 460;
constexpr float kSampleRate = 48000.0f;
constexpr float kMeterFloorDb = -60.0f;
constexpr float kMeterCeilDb = 0.0f;

static constexpr const char* kModeNames[] = {
    "COMPRESSOR",
    "EXPANDER",
    "GATE",
};

static constexpr const char* kKnobTitles[] = {
    "Threshold",
    "Ratio",
    "Range",
    "Attack",
    "Release",
    "Makeup",
};

// ------------------------------------------------------------
//  参数（UI 线程写入，音频线程读取）
// ------------------------------------------------------------
static std::atomic<int> s_mode{static_cast<int>(DynamicsMode::Compressor)};
static std::atomic<float> s_threshold_db{-24.0f};
static std::atomic<float> s_ratio{4.0f};
static std::atomic<float> s_range_db{60.0f};
static std::atomic<float> s_attack_ms{10.0f};
static std::atomic<float> s_release_ms{150.0f};
static std::atomic<float> s_makeup_db{0.0f};

// ------------------------------------------------------------
//  电平表（音频线程写入，UI 线程读取）
// ------------------------------------------------------------
static std::atomic<float> s_input_peak{0.0f};
static std::atomic<float> s_output_peak{0.0f};
static std::atomic<float> s_gain_reduction_db{0.0f};

static DynamicsProcessor s_processor;

// ------------------------------------------------------------
//  旋钮
// ------------------------------------------------------------
static std::array<Knob, 6> s_knobs;

static constexpr float kKnobX = 20.0f;
static constexpr float kKnobY = 118.0f;
static constexpr float kKnobWidth = 108.0f;
static constexpr float kKnobHeight = 130.0f;
static constexpr float kKnobGap = 8.0f;

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
        .set_name_font_size(13)
        .set_number_font_size(12)
        .set_fore_color(Color{235, 235, 238, 255})
        .set_bg_color(Color{24, 26, 29, 255});
    knob.value_to_text_function = std::move(format);
    knob.on_value_change = std::move(on_change);
}

static void CreateKnobs() {
    SetupKnob(s_knobs[0], 0, "Threshold", -60.0f, 0.0f, 0.5f, s_threshold_db.load(),
              2, [](float v) { return std::format("{:+.1f} dB", v); },
              [](float v) { s_threshold_db.store(v, std::memory_order_relaxed); });
    SetupKnob(s_knobs[1], 1, "Ratio", 1.0f, 20.0f, 0.1f, s_ratio.load(),
              1, [](float v) { return std::format("{:.1f}:1", v); },
              [](float v) { s_ratio.store(v, std::memory_order_relaxed); });
    SetupKnob(s_knobs[2], 2, "Range", 0.0f, 100.0f, 1.0f, s_range_db.load(),
              2, [](float v) { return std::format("{:.0f} dB", v); },
              [](float v) { s_range_db.store(v, std::memory_order_relaxed); });
    SetupKnob(s_knobs[3], 3, "Attack", 0.1f, 200.0f, 0.1f, s_attack_ms.load(),
              10, [](float v) { return std::format("{:.1f} ms", v); },
              [](float v) { s_attack_ms.store(v, std::memory_order_relaxed); });
    SetupKnob(s_knobs[4], 4, "Release", 5.0f, 1000.0f, 5.0f, s_release_ms.load(),
              2, [](float v) { return std::format("{:.0f} ms", v); },
              [](float v) { s_release_ms.store(v, std::memory_order_relaxed); });
    SetupKnob(s_knobs[5], 5, "Makeup", 0.0f, 24.0f, 0.5f, s_makeup_db.load(),
              2, [](float v) { return std::format("{:+.1f} dB", v); },
              [](float v) { s_makeup_db.store(v, std::memory_order_relaxed); });

    // 同步一次，防止旋钮初始值和音频线程参数不一致
    s_threshold_db.store(s_knobs[0].get_value(), std::memory_order_relaxed);
    s_ratio.store(s_knobs[1].get_value(), std::memory_order_relaxed);
    s_range_db.store(s_knobs[2].get_value(), std::memory_order_relaxed);
    s_attack_ms.store(s_knobs[3].get_value(), std::memory_order_relaxed);
    s_release_ms.store(s_knobs[4].get_value(), std::memory_order_relaxed);
    s_makeup_db.store(s_knobs[5].get_value(), std::memory_order_relaxed);
}

static bool IsKnobEnabled(size_t index, DynamicsMode mode) noexcept {
    switch (index) {
    case 1: // Ratio
        return mode != DynamicsMode::Gate;
    case 2: // Range
        return mode != DynamicsMode::Compressor;
    case 5: // Makeup
        return mode != DynamicsMode::Gate;
    default:
        return true;
    }
}

static void DrawDisabledKnob(size_t index) {
    constexpr int kNameFont = 13;
    constexpr int kNumberFont = 12;

    float const x = kKnobX + static_cast<float>(index) * (kKnobWidth + kKnobGap);
    char const* title = kKnobTitles[index];

    int const title_width = MeasureText(title, kNameFont);
    DrawText(title, static_cast<int>(x + (kKnobWidth - static_cast<float>(title_width)) * 0.5f),
             static_cast<int>(kKnobY), kNameFont, Color{90, 94, 100, 255});

    float const slider_y = kKnobY + static_cast<float>(kNameFont);
    float const slider_height = kKnobHeight - static_cast<float>(kNameFont + kNumberFont);
    float const center_x = x + kKnobWidth * 0.5f;
    float const center_y = slider_y + slider_height * 0.5f;
    float const radius = std::min(kKnobWidth, slider_height) * 0.5f;

    DrawCircle(static_cast<int>(center_x), static_cast<int>(center_y), radius, Color{20, 22, 25, 255});
    DrawCircleLines(static_cast<int>(center_x), static_cast<int>(center_y), radius, Color{52, 55, 60, 255});

    char const* dash = "--";
    int const dash_width = MeasureText(dash, kNumberFont);
    DrawText(dash, static_cast<int>(x + (kKnobWidth - static_cast<float>(dash_width)) * 0.5f),
             static_cast<int>(kKnobY + kKnobHeight - static_cast<float>(kNumberFont)), kNumberFont,
             Color{90, 94, 100, 255});
}

static void DisplayKnobs(DynamicsMode mode) {
    for (size_t i = 0; i < s_knobs.size(); ++i) {
        bool const enabled = IsKnobEnabled(i, mode);
        s_knobs[i].SetEnable(enabled);
        if (enabled) {
            s_knobs[i].display();
        }
        else {
            DrawDisabledKnob(i);
        }
    }
}

// ------------------------------------------------------------
//  模式选择
// ------------------------------------------------------------
namespace mode_bar {

constexpr float kX = 20.0f;
constexpr float kY = 72.0f;
constexpr float kWidth = 120.0f;
constexpr float kHeight = 28.0f;
constexpr float kGap = 8.0f;
constexpr int kFontSize = 12;

static Rectangle Bounds(int index) {
    return Rectangle{
        kX + static_cast<float>(index) * (kWidth + kGap),
        kY,
        kWidth,
        kHeight,
    };
}

static int HitTest(Vector2 mouse) {
    for (int i = 0; i < static_cast<int>(DynamicsMode::kNumModes); ++i) {
        if (CheckCollisionPointRec(mouse, Bounds(i))) {
            return i;
        }
    }
    return -1;
}

static void Draw(int active, Vector2 mouse) {
    for (int i = 0; i < static_cast<int>(DynamicsMode::kNumModes); ++i) {
        Rectangle const bounds = Bounds(i);
        bool const selected = i == active;
        bool const hovered = CheckCollisionPointRec(mouse, bounds);

        Color const border = selected ? Color{225, 172, 74, 255}
                                      : (hovered ? Color{120, 124, 130, 255} : Color{70, 74, 80, 255});
        Color const fill = selected ? Color{225, 172, 74, 255}
                                    : (hovered ? Color{40, 43, 48, 255} : Color{30, 32, 36, 255});
        Color const text_color = selected ? Color{18, 20, 23, 255} : Color{190, 194, 199, 255};

        DrawRectangleRec(bounds, fill);
        DrawRectangleLinesEx(bounds, 1.0f, border);

        char const* name = kModeNames[i];
        int const text_width = MeasureText(name, kFontSize);
        DrawText(name, static_cast<int>(bounds.x + (bounds.width - static_cast<float>(text_width)) * 0.5f),
                 static_cast<int>(bounds.y + (bounds.height - static_cast<float>(kFontSize)) * 0.5f),
                 kFontSize, text_color);
    }
}

} // namespace mode_bar

// ------------------------------------------------------------
//  电平表绘制
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

static void DrawGainReductionMeter(float x, float y, float width, float height, float gain_reduction_db) {
    constexpr float kMaxReductionDb = 60.0f;

    DrawText("GAIN REDUCTION", static_cast<int>(x), static_cast<int>(y - 18.0f), 10,
             Color{140, 140, 150, 255});

    DrawRectangle(static_cast<int>(x), static_cast<int>(y), static_cast<int>(width), static_cast<int>(height),
                  Color{25, 25, 30, 255});
    DrawRectangleLines(static_cast<int>(x), static_cast<int>(y), static_cast<int>(width), static_cast<int>(height),
                       Color{60, 60, 70, 255});

    float const norm = std::clamp(-gain_reduction_db / kMaxReductionDb, 0.0f, 1.0f);
    if (norm > 0.0f) {
        int const fill_width = static_cast<int>(norm * (width - 2.0f));
        DrawRectangle(static_cast<int>(x) + 1, static_cast<int>(y) + 1, fill_width,
                      static_cast<int>(height) - 2, Color{225, 172, 74, 255});
    }

    std::string const text = std::format("{:.1f} dB", gain_reduction_db);
    int const text_width = MeasureText(text.c_str(), 12);
    DrawText(text.c_str(), static_cast<int>(x + width - static_cast<float>(text_width) - 8.0f),
             static_cast<int>(y + (height - 12.0f) * 0.5f), 12, Color{235, 235, 238, 255});

    DrawText("0", static_cast<int>(x), static_cast<int>(y + height + 6.0f), 10,
             Color{140, 140, 150, 255});
    std::string const max_text = std::format("-{:.0f} dB", kMaxReductionDb);
    DrawText(max_text.c_str(),
             static_cast<int>(x + width - static_cast<float>(MeasureText(max_text.c_str(), 10))),
             static_cast<int>(y + height + 6.0f), 10, Color{140, 140, 150, 255});
}

} // namespace dynamic_demo

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

    using namespace dynamic_demo;

    DynamicsMode const mode = static_cast<DynamicsMode>(s_mode.load(std::memory_order_relaxed));
    s_processor.SetParameters(
        mode,
        s_threshold_db.load(std::memory_order_relaxed),
        s_ratio.load(std::memory_order_relaxed),
        s_range_db.load(std::memory_order_relaxed),
        s_attack_ms.load(std::memory_order_relaxed),
        s_release_ms.load(std::memory_order_relaxed),
        s_makeup_db.load(std::memory_order_relaxed));

    // miniaudio 的 input/output 是两个独立缓冲区，先拷贝再原地处理
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
    s_gain_reduction_db.store(s_processor.GetGainReductionDb(), std::memory_order_relaxed);
}

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------
int main(void) {
    using namespace dynamic_demo;

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
    InitWindow(kWindowWidth, kWindowHeight, "Dynamics - Compressor");
    SetTargetFPS(60);

    // ── 电平表显示值（UI 线程）──
    float display_input = 0.0f;
    float display_output = 0.0f;
    float display_gain_reduction_db = 0.0f;

    int active_mode = s_mode.load(std::memory_order_relaxed);

    // ── main loop ──
    while (!WindowShouldClose()) {
        Vector2 const mouse = GetMousePosition();

        // ── 模式切换 ──
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            int const hit = mode_bar::HitTest(mouse);
            if (hit >= 0 && hit != active_mode) {
                active_mode = hit;
                s_mode.store(hit, std::memory_order_relaxed);
                SetWindowTitle(std::format("Dynamics - {}", kModeNames[hit]).c_str());
            }
        }
        auto const mode = static_cast<DynamicsMode>(active_mode);

        // ── 电平表衰减 ──
        constexpr float kDecay = 0.86f;
        display_input = std::max(s_input_peak.load(std::memory_order_relaxed), display_input * kDecay);
        display_output = std::max(s_output_peak.load(std::memory_order_relaxed), display_output * kDecay);
        display_gain_reduction_db =
            std::min(s_gain_reduction_db.load(std::memory_order_relaxed), display_gain_reduction_db * kDecay);

        BeginDrawing();
        ClearBackground(Color{18, 20, 23, 255});

        // ── 顶栏 ──
        DrawRectangle(0, 0, kWindowWidth, 64, Color{24, 26, 29, 255});
        DrawRectangle(0, 63, kWindowWidth, 1, Color{48, 52, 57, 255});
        DrawText("DYNAMICS PROCESSOR", 20, 14, 24, Color{232, 233, 235, 255});
        DrawText("GATE / EXPANDER / COMPRESSOR - BLOCK PEAK DETECTION", 21, 43, 10,
                 Color{126, 132, 139, 255});

        // ── 模式按钮 ──
        mode_bar::Draw(active_mode, mouse);

        // ── 旋钮 ──
        DisplayKnobs(mode);

        // ── 电平表 ──
        constexpr float kMeterY = 280.0f;
        constexpr float kMeterWidth = 28.0f;
        constexpr float kMeterHeight = 110.0f;
        DrawVerticalMeter(24.0f, kMeterY, kMeterWidth, kMeterHeight, display_input, "IN",
                          Color{76, 201, 151, 255});
        DrawVerticalMeter(72.0f, kMeterY, kMeterWidth, kMeterHeight, display_output, "OUT",
                          Color{76, 201, 151, 255});
        DrawGainReductionMeter(150.0f, kMeterY + 36.0f, 550.0f, 22.0f, display_gain_reduction_db);

        // ── 底部提示 ──
        DrawText("drag knob: adjust   shift/alt: fine   right click: reset", 20, kWindowHeight - 20, 10,
                 Color{80, 80, 90, 255});

        EndDrawing();
    }

    // ── cleanup ──
    ma_device_stop(&device);
    ma_device_uninit(&device);
    CloseWindow();
    return 0;
}
