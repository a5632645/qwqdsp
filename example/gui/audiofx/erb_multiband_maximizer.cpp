// ------------------------------------------------------------
//  ERB Multiband Maximizer — 12 段多段最大化器
//
//  频段: 20 Hz - 20 kHz 按 ERB 尺度均分 12 段 (相邻高通相减, 完全重构)
//  每段: Drive -> look-ahead limiter -> 段 Ceiling
//  12 段相加后再过一个最终 limiter，保证总输出 ceiling。
//
//  每段最终 drive(dB) = 总 drive(dB) + 20log10(系数)
//  即线性域上的 drive = 系数 * 总 drive。
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

#include "erb_multiband_maximizer_dsp.hpp"

namespace erb_multiband_maximizer_demo {

constexpr int kWindowWidth = 1080;
constexpr int kWindowHeight = 640;
constexpr float kSampleRate = 48000.0f;
constexpr float kMeterFloorDb = -60.0f;
constexpr float kMeterCeilDb = 0.0f;
constexpr float kMaxBandReductionDb = 24.0f;

constexpr Color kBackground{18, 20, 23, 255};
constexpr Color kPanelColor{22, 24, 27, 255};
constexpr Color kPanelLine{38, 41, 46, 255};
constexpr Color kTextColor{232, 233, 235, 255};
constexpr Color kDimColor{126, 132, 139, 255};
constexpr Color kFaintColor{80, 85, 92, 255};

// ------------------------------------------------------------
//  预设：系数作用于各频段的 drive
// ------------------------------------------------------------
struct Preset {
    const char* name;
    const char* description;
    std::array<float, kNumBands> coefficients;
};

static constexpr std::array<Preset, 8> kPresets{
    {
     {"FLAT",
         "same drive on every band",
         {1.00f, 1.00f, 1.00f, 1.00f, 1.00f, 1.00f, 1.00f, 1.00f, 1.00f, 1.00f, 1.00f, 1.00f}},
     {"LOUDNESS",
         "inverse equal-loudness: push lows and highs",
         {1.70f, 1.45f, 1.20f, 1.05f, 0.97f, 0.92f, 0.88f, 0.86f, 0.90f, 0.98f, 1.12f, 1.30f}},
     {"WARMTH",
         "drive the 100-600 Hz body",
         {1.05f, 1.35f, 1.40f, 1.25f, 1.10f, 0.97f, 0.92f, 0.90f, 0.93f, 0.97f, 1.00f, 1.00f}},
     {"PRESENCE",
         "drive 2-6 kHz for forwardness",
         {0.80f, 0.85f, 0.88f, 0.92f, 0.97f, 1.10f, 1.30f, 1.40f, 1.30f, 1.10f, 1.00f, 0.95f}},
     {"DE-MUD",
         "ease 200-600 Hz, keep sub and clarity",
         {1.10f, 0.78f, 0.70f, 0.82f, 1.02f, 1.10f, 1.18f, 1.18f, 1.12f, 1.02f, 1.00f, 1.00f}},
     {"AIR",
         "drive 8 kHz and above",
         {1.00f, 1.00f, 1.00f, 1.00f, 1.00f, 1.00f, 0.98f, 0.98f, 1.03f, 1.12f, 1.30f, 1.50f}},
     {"VOCAL",
         "lift 1-4 kHz, tame sub and sibilance",
         {0.55f, 0.65f, 0.85f, 1.00f, 1.15f, 1.30f, 1.35f, 1.25f, 1.05f, 0.90f, 0.85f, 0.82f}},
     {"MASTERING",
         "gentle smile curve for full mixes",
         {1.12f, 1.06f, 1.00f, 0.97f, 0.94f, 0.91f, 0.89f, 0.90f, 0.94f, 1.00f, 1.06f, 1.12f}},
     }
};

// ------------------------------------------------------------
//  参数（UI 线程写入，音频线程读取）
// ------------------------------------------------------------
static std::atomic<float> s_total_drive_db{6.0f};
static std::atomic<float> s_ceiling_db{-0.1f};
static std::atomic<float> s_contrast{0.0f};
static std::atomic<float> s_release_ms{100.0f};
static std::atomic<float> s_lookahead_ms{5.0f};
static std::atomic<bool> s_bypass{false};
static std::array<std::atomic<float>, kNumBands> s_coefficients;

// ------------------------------------------------------------
//  电平表（音频线程写入，UI 线程读取）
// ------------------------------------------------------------
static std::atomic<float> s_input_peak{0.0f};
static std::atomic<float> s_output_peak{0.0f};
static std::array<std::atomic<float>, kNumBands> s_band_gain_reduction_db;
static std::atomic<float> s_final_gain_reduction_db{0.0f};

static ErbMultibandMaximizer s_processor;

// ------------------------------------------------------------
//  控件与标签
// ------------------------------------------------------------
static int s_preset_index = 0;

static std::array<std::string, kNumBands> s_band_labels;
static std::array<float, kNumBands> s_band_centers;

static Knob s_total_drive_knob;
static Knob s_ceiling_knob;
static Knob s_contrast_knob;
static Knob s_release_knob;
static Knob s_lookahead_knob;

constexpr float kPanelX = 10.0f;
constexpr float kPanelY = 206.0f;
constexpr float kPanelW = static_cast<float>(kWindowWidth) - 20.0f;
constexpr float kPanelH = 366.0f;

// 左侧: 系数-频率曲线
constexpr float kMaxCoefficient = 4.0f;

constexpr float kChartX = 176.0f;
constexpr float kChartY = 232.0f;
constexpr float kChartW = 374.0f;
constexpr float kChartH = 322.0f;
constexpr float kPlotLeft = 208.0f;
constexpr float kPlotRight = 542.0f;
constexpr float kPlotTop = 240.0f;
constexpr float kPlotBottom = 526.0f;

// 右侧: 各段 GR
constexpr float kGrAreaX = 574.0f;
constexpr float kGrAreaW = static_cast<float>(kWindowWidth) - 20.0f - kGrAreaX;
constexpr float kGrColumnW = kGrAreaW / static_cast<float>(kNumBands);
constexpr float kGrMeterW = 24.0f;
constexpr float kGrMeterY = 250.0f;
constexpr float kGrMeterH = 280.0f;

static std::string FormatShortHz(float hz) {
    if (hz < 1000.0f) {
        return std::format("{:.0f}", hz);
    }
    float const k = hz / 1000.0f;
    if (k < 10.0f) {
        std::string text = std::format("{:.1f}", k);
        if (text.ends_with(".0")) {
            text.resize(text.size() - 2);
        }
        return text + "k";
    }
    return std::format("{:.0f}k", k);
}

// contrast: 以 1 为轴心的幂曲线, 让大于 1 的系数更大, 小于 1 的系数更小
static float ShapeCoefficient(float coefficient, float contrast) noexcept {
    if (contrast <= 0.0f) {
        return coefficient;
    }
    float const exponent = std::exp2(contrast);
    return std::clamp(std::pow(std::max(coefficient, 0.0f), exponent), 0.0f, kMaxCoefficient);
}

static float UnshapeCoefficient(float displayed, float contrast) noexcept {
    if (contrast <= 0.0f) {
        return displayed;
    }
    float const exponent = std::exp2(contrast);
    return std::clamp(std::pow(std::max(displayed, 0.0f), 1.0f / exponent), 0.0f, kMaxCoefficient);
}

static void BuildBandLabels() {
    auto const& edges = s_processor.GetBandEdges();
    for (size_t band = 0; band < kNumBands; ++band) {
        float const center = std::sqrt(edges[band] * edges[band + 1]);
        s_band_centers[band] = center;
        s_band_labels[band] = FormatShortHz(center);
    }
}

static void SetupKnob(Knob& knob, Rectangle bound, const char* title, float min, float max, float step, float value,
                      int sensitivity, std::function<std::string(float)> format, std::function<void(float)> on_change) {
    knob.set_bound(bound)
        .set_title(title)
        .set_range(min, max, step, value)
        .set_value(value)
        .set_sensitivity(sensitivity)
        .set_name_font_size(11)
        .set_number_font_size(10)
        .set_fore_color(kTextColor)
        .set_bg_color(Color{24, 26, 29, 255});
    knob.value_to_text_function = std::move(format);
    knob.on_value_change = std::move(on_change);
}

static void ApplyPreset(size_t index) {
    index %= kPresets.size();
    s_preset_index = static_cast<int>(index);
    for (size_t band = 0; band < kNumBands; ++band) {
        s_coefficients[band].store(kPresets[index].coefficients[band], std::memory_order_relaxed);
    }
}

static void CreateKnobs() {
    SetupKnob(
        s_release_knob, Rectangle{20.0f, 76.0f, 96.0f, 112.0f}, "Release", 10.0f, 1000.0f, 5.0f, 100.0f, 2,
        [](float v) { return std::format("{:.0f} ms", v); },
        [](float v) { s_release_ms.store(v, std::memory_order_relaxed); });

    SetupKnob(
        s_lookahead_knob, Rectangle{122.0f, 76.0f, 96.0f, 112.0f}, "Look-ahead", 0.5f, 10.0f, 0.1f, 5.0f, 2,
        [](float v) { return std::format("{:.1f} ms", v); },
        [](float v) { s_lookahead_ms.store(v, std::memory_order_relaxed); });

    SetupKnob(
        s_total_drive_knob, Rectangle{20.0f, 222.0f, 96.0f, 112.0f}, "Drive", 0.0f, 24.0f, 0.5f, 6.0f, 2,
        [](float v) { return std::format("{:+.1f} dB", v); },
        [](float v) { s_total_drive_db.store(v, std::memory_order_relaxed); });

    SetupKnob(
        s_ceiling_knob, Rectangle{16.0f, 420.0f, 70.0f, 112.0f}, "All Ceil", -24.0f, 0.0f, 0.1f, -0.1f, 1,
        [](float v) { return std::format("{:+.1f} dB", v); },
        [](float v) { s_ceiling_db.store(v, std::memory_order_relaxed); });

    SetupKnob(
        s_contrast_knob, Rectangle{90.0f, 420.0f, 70.0f, 112.0f}, "Contrast", 0.0f, 1.0f, 0.01f, 0.0f, 2,
        [](float v) { return std::format("{:.2f}", v); },
        [](float v) { s_contrast.store(v, std::memory_order_relaxed); });
}

// ------------------------------------------------------------
//  绘制
// ------------------------------------------------------------
static float AmplitudeToNorm(float amplitude) noexcept {
    float const db = 20.0f * std::log10(std::max(amplitude, 1.0e-6f));
    return std::clamp((db - kMeterFloorDb) / (kMeterCeilDb - kMeterFloorDb), 0.0f, 1.0f);
}

static void DrawVerticalMeter(float x, float y, float width, float height, float value, const char* label,
                              Color fill_color) {
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
             static_cast<int>(y + height + 6.0f), 10, kDimColor);
}

static void DrawGainReductionMeter(float x, float y, float width, float height, float gain_reduction_db,
                                   const char* label, Color color) {
    DrawText(label, static_cast<int>(x), static_cast<int>(y - 18.0f), 10, kDimColor);

    DrawRectangle(static_cast<int>(x), static_cast<int>(y), static_cast<int>(width), static_cast<int>(height),
                  Color{25, 25, 30, 255});
    DrawRectangleLines(static_cast<int>(x), static_cast<int>(y), static_cast<int>(width), static_cast<int>(height),
                       Color{60, 60, 70, 255});

    float const norm = std::clamp(-gain_reduction_db / kMaxBandReductionDb, 0.0f, 1.0f);
    if (norm > 0.0f) {
        int const fill_width = static_cast<int>(norm * (width - 2.0f));
        DrawRectangle(static_cast<int>(x) + 1, static_cast<int>(y) + 1, fill_width, static_cast<int>(height) - 2,
                      color);
    }

    std::string const text = std::format("{:.1f} dB", gain_reduction_db);
    int const text_width = MeasureText(text.c_str(), 11);
    DrawText(text.c_str(), static_cast<int>(x + width - static_cast<float>(text_width) - 8.0f),
             static_cast<int>(y + (height - 11.0f) * 0.5f), 11, kTextColor);
}

static Color BandColor(size_t band) {
    float const t = static_cast<float>(band) / static_cast<float>(kNumBands - 1);
    return ColorFromHSV(170.0f - 125.0f * t, 0.55f, 0.95f);
}

static void DrawPresetSelector(Rectangle bound) {
    DrawRectangleRec(bound, Color{24, 26, 29, 255});
    DrawRectangleLinesEx(bound, 1.0f, Color{60, 64, 70, 255});

    DrawText("CURVE", static_cast<int>(bound.x + 10.0f), static_cast<int>(bound.y + 7.0f), 10, kDimColor);

    auto const& preset = kPresets[static_cast<size_t>(s_preset_index)];
    int const name_width = MeasureText(preset.name, 14);
    DrawText(preset.name, static_cast<int>(bound.x + (bound.width - static_cast<float>(name_width)) * 0.5f),
             static_cast<int>(bound.y + 26.0f), 14, kTextColor);

    DrawText("<", static_cast<int>(bound.x + 8.0f), static_cast<int>(bound.y + bound.height * 0.5f - 10.0f), 18,
             kDimColor);
    DrawText(">", static_cast<int>(bound.x + bound.width - 18.0f),
             static_cast<int>(bound.y + bound.height * 0.5f - 10.0f), 18, kDimColor);

    std::string const counter = std::format("{}/{}", s_preset_index + 1, kPresets.size());
    int const counter_width = MeasureText(counter.c_str(), 9);
    DrawText(counter.c_str(), static_cast<int>(bound.x + bound.width - static_cast<float>(counter_width) - 10.0f),
             static_cast<int>(bound.y + bound.height - 14.0f), 9, kFaintColor);

    auto const mouse = GetMousePosition();
    if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
        Rectangle const left_button{bound.x, bound.y, 26.0f, bound.height};
        Rectangle const right_button{bound.x + bound.width - 26.0f, bound.y, 26.0f, bound.height};
        size_t const count = kPresets.size();
        if (CheckCollisionPointRec(mouse, left_button)) {
            ApplyPreset((static_cast<size_t>(s_preset_index) + count - 1) % count);
        }
        else if (CheckCollisionPointRec(mouse, right_button)) {
            ApplyPreset((static_cast<size_t>(s_preset_index) + 1) % count);
        }
    }
}

// ------------------------------------------------------------
//  系数-频率曲线（左侧）
// ------------------------------------------------------------
static int s_drag_band = -1;
static float s_drag_start_y = 0.0f;
static float s_drag_start_value = 1.0f;

static float ChartFreqToX(float hz) noexcept {
    static float const log_low = std::log10(kErbLowHz);
    static float const log_high = std::log10(kErbHighHz);
    float const t = (std::log10(hz) - log_low) / (log_high - log_low);
    return kPlotLeft + t * (kPlotRight - kPlotLeft);
}

static float ChartValueToY(float value) noexcept {
    return kPlotBottom - value / kMaxCoefficient * (kPlotBottom - kPlotTop);
}

static int NearestChartBand(float x) noexcept {
    int nearest = 0;
    float best = 1.0e9f;
    for (size_t band = 0; band < kNumBands; ++band) {
        float const dx = std::abs(ChartFreqToX(s_band_centers[band]) - x);
        if (dx < best) {
            best = dx;
            nearest = static_cast<int>(band);
        }
    }
    return nearest;
}

static void HandleCoefficientChartInput() {
    auto const mouse = GetMousePosition();
    Rectangle const plot{kPlotLeft, kPlotTop, kPlotRight - kPlotLeft, kPlotBottom - kPlotTop};

    float const contrast = s_contrast.load(std::memory_order_relaxed);

    if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON) && CheckCollisionPointRec(mouse, plot)) {
        s_drag_band = NearestChartBand(mouse.x);
        s_drag_start_y = mouse.y;
        s_drag_start_value = ShapeCoefficient(
            s_coefficients[static_cast<size_t>(s_drag_band)].load(std::memory_order_relaxed), contrast);
    }

    if (s_drag_band >= 0) {
        if (IsMouseButtonDown(MOUSE_LEFT_BUTTON)) {
            float const value_range = kMaxCoefficient / (kPlotBottom - kPlotTop);
            float const displayed =
                std::clamp(s_drag_start_value + (s_drag_start_y - mouse.y) * value_range, 0.0f, kMaxCoefficient);
            s_coefficients[static_cast<size_t>(s_drag_band)].store(UnshapeCoefficient(displayed, contrast),
                                                                   std::memory_order_relaxed);
        }
        else {
            s_drag_band = -1;
        }
    }
}

static void DrawCoefficientChart() {
    auto const mouse = GetMousePosition();
    Rectangle const plot{kPlotLeft, kPlotTop, kPlotRight - kPlotLeft, kPlotBottom - kPlotTop};

    float const contrast = s_contrast.load(std::memory_order_relaxed);

    int hover_band = -1;
    if (CheckCollisionPointRec(mouse, plot)) {
        int const nearest = NearestChartBand(mouse.x);
        float const node_y = ChartValueToY(
            ShapeCoefficient(s_coefficients[static_cast<size_t>(nearest)].load(std::memory_order_relaxed), contrast));
        if (std::abs(mouse.y - node_y) <= 14.0f) {
            hover_band = nearest;
        }
    }

    DrawText("BAND COEFFICIENT", static_cast<int>(kChartX), static_cast<int>(kChartY - 20.0f), 10, kDimColor);

    DrawRectangle(static_cast<int>(kChartX), static_cast<int>(kChartY), static_cast<int>(kChartW),
                  static_cast<int>(kChartH), Color{19, 21, 24, 255});
    DrawRectangleLinesEx(Rectangle{kChartX, kChartY, kChartW, kChartH}, 1.0f, Color{32, 35, 40, 255});
    DrawRectangleRec(plot, Color{22, 24, 27, 255});

    // 纵向网格 + 系数刻度
    for (int i = 0; i <= 8; ++i) {
        float const value = 0.5f * static_cast<float>(i);
        float const y = ChartValueToY(value);
        Color const line = (i == 2) ? Color{50, 55, 61, 255} : Color{32, 35, 40, 255};
        DrawLine(static_cast<int>(kPlotLeft), static_cast<int>(y), static_cast<int>(kPlotRight), static_cast<int>(y),
                 line);

        if (i % 2 == 0) {
            std::string const label = std::format("{:.0f}", value);
            int const label_width = MeasureText(label.c_str(), 8);
            DrawText(label.c_str(), static_cast<int>(kPlotLeft - static_cast<float>(label_width) - 6.0f),
                     static_cast<int>(y - 4.0f), 8, kFaintColor);
        }
    }

    // 频段边界 + 频率刻度
    auto const& edges = s_processor.GetBandEdges();
    for (size_t edge = 1; edge < kNumBands; ++edge) {
        int const x = static_cast<int>(ChartFreqToX(edges[edge]));
        DrawLine(x, static_cast<int>(kPlotTop), x, static_cast<int>(kPlotBottom), Color{29, 32, 36, 255});
    }
    for (float hz : {20.0f, 100.0f, 1000.0f, 10000.0f, 20000.0f}) {
        int const x = static_cast<int>(ChartFreqToX(hz));
        DrawLine(x, static_cast<int>(kPlotBottom), x, static_cast<int>(kPlotBottom) + 3, Color{60, 64, 70, 255});

        std::string const label = FormatShortHz(hz);
        int const label_width = MeasureText(label.c_str(), 8);
        int const text_x =
            std::clamp(x - label_width / 2, static_cast<int>(kPlotLeft), static_cast<int>(kPlotRight) - label_width);
        DrawText(label.c_str(), text_x, static_cast<int>(kPlotBottom) + 6, 8, kFaintColor);
    }

    // 曲线
    std::array<Vector2, kNumBands> points;
    std::array<float, kNumBands> values;
    for (size_t band = 0; band < kNumBands; ++band) {
        values[band] = ShapeCoefficient(s_coefficients[band].load(std::memory_order_relaxed), contrast);
        points[band] = Vector2{ChartFreqToX(s_band_centers[band]), ChartValueToY(values[band])};
    }

    Color const fill_color{76, 201, 151, 36};
    for (size_t band = 0; band + 1 < kNumBands; ++band) {
        Vector2 const bottom_left{points[band].x, kPlotBottom};
        Vector2 const bottom_right{points[band + 1].x, kPlotBottom};
        DrawTriangle(bottom_left, points[band], points[band + 1], fill_color);
        DrawTriangle(bottom_left, points[band + 1], bottom_right, fill_color);
    }
    for (size_t band = 0; band + 1 < kNumBands; ++band) {
        DrawLineEx(points[band], points[band + 1], 1.8f, Color{76, 201, 151, 255});
    }

    // 节点
    for (size_t band = 0; band < kNumBands; ++band) {
        bool const active = static_cast<int>(band) == s_drag_band;
        bool const hovered = static_cast<int>(band) == hover_band;
        DrawCircleV(points[band], active ? 5.0f : (hovered ? 4.5f : 3.5f), BandColor(band));
        if (active) {
            DrawCircleLinesV(points[band], 7.0f, kTextColor);
        }
    }

    // 当前节点数值
    int const label_band = s_drag_band >= 0 ? s_drag_band : hover_band;
    if (label_band >= 0) {
        std::string const text = std::format("x{:.2f}", values[static_cast<size_t>(label_band)]);
        int const text_width = MeasureText(text.c_str(), 10);
        Vector2 const point = points[static_cast<size_t>(label_band)];
        float const text_x = std::clamp(point.x + 8.0f, kPlotLeft, kPlotRight - static_cast<float>(text_width) - 2.0f);
        float const text_y = std::clamp(point.y - 18.0f, kPlotTop + 2.0f, kPlotBottom - 14.0f);
        DrawText(text.c_str(), static_cast<int>(text_x), static_cast<int>(text_y), 10, kTextColor);
    }
}

static void DrawBypassButton(Rectangle bound) {
    bool const bypass = s_bypass.load(std::memory_order_relaxed);
    auto const mouse = GetMousePosition();
    if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON) && CheckCollisionPointRec(mouse, bound)) {
        s_bypass.store(!bypass, std::memory_order_relaxed);
    }

    Color const accent{225, 172, 74, 255};
    Color const background = bypass ? accent : Color{24, 26, 29, 255};
    Color const foreground = bypass ? Color{20, 22, 26, 255} : kDimColor;
    DrawRectangleRec(bound, background);
    DrawRectangleLinesEx(bound, 1.0f, bypass ? accent : Color{60, 64, 70, 255});

    constexpr char kLabel[] = "BYPASS";
    int const label_width = MeasureText(kLabel, 12);
    DrawText(kLabel, static_cast<int>(bound.x + (bound.width - static_cast<float>(label_width)) * 0.5f),
             static_cast<int>(bound.y + (bound.height - 12.0f) * 0.5f), 12, foreground);
}

static void DrawGrGrid() {
    for (int i = 0; i <= 4; ++i) {
        float const db = kMaxBandReductionDb * 0.25f * static_cast<float>(i);
        float const y = kGrMeterY + (db / kMaxBandReductionDb) * kGrMeterH;
        DrawLine(static_cast<int>(kGrAreaX), static_cast<int>(y), static_cast<int>(kGrAreaX + kGrAreaW),
                 static_cast<int>(y), Color{32, 35, 40, 255});

        std::string const label = std::format("{:.0f}", db);
        int const label_width = MeasureText(label.c_str(), 8);
        DrawText(label.c_str(), static_cast<int>(kGrAreaX - static_cast<float>(label_width) - 6.0f),
                 static_cast<int>(y - 4.0f), 8, kFaintColor);
    }
}

static void DrawBandGrMeter(size_t band, float gain_reduction_db) {
    float const center_x = kGrAreaX + static_cast<float>(band) * kGrColumnW + kGrColumnW * 0.5f;
    float const x = center_x - kGrMeterW * 0.5f;
    Color const color = BandColor(band);

    DrawRectangle(static_cast<int>(x), static_cast<int>(kGrMeterY), static_cast<int>(kGrMeterW),
                  static_cast<int>(kGrMeterH), Color{24, 26, 29, 255});
    DrawRectangleLines(static_cast<int>(x), static_cast<int>(kGrMeterY), static_cast<int>(kGrMeterW),
                       static_cast<int>(kGrMeterH), Color{44, 47, 53, 255});

    float const norm = std::clamp(-gain_reduction_db / kMaxBandReductionDb, 0.0f, 1.0f);
    int const fill_height = static_cast<int>(norm * (kGrMeterH - 2.0f));
    if (fill_height > 0) {
        DrawRectangle(static_cast<int>(x) + 1, static_cast<int>(kGrMeterY) + 1, static_cast<int>(kGrMeterW) - 2,
                      fill_height, color);
    }

    std::string const text = std::format("{:.1f}", gain_reduction_db);
    int const text_width = MeasureText(text.c_str(), 9);
    bool const over_fill =
        static_cast<float>(kGrMeterY + 1.0f + static_cast<float>(fill_height)) > (kGrMeterY + kGrMeterH - 12.0f);
    DrawText(text.c_str(), static_cast<int>(x + (kGrMeterW - static_cast<float>(text_width)) * 0.5f),
             static_cast<int>(kGrMeterY + kGrMeterH - 12.0f), 9,
             over_fill ? Color{20, 22, 26, 255} : Color{200, 203, 208, 255});

    int const label_width = MeasureText(s_band_labels[band].c_str(), 9);
    DrawText(s_band_labels[band].c_str(), static_cast<int>(center_x - static_cast<float>(label_width) * 0.5f),
             static_cast<int>(kGrMeterY + kGrMeterH + 8.0f), 9, kTextColor);
}

} // namespace erb_multiband_maximizer_demo

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

    using namespace erb_multiband_maximizer_demo;

    // 旁通切换时清空限幅器/分频器状态，避免恢复处理时残留 look-ahead 数据
    bool const bypass = s_bypass.load(std::memory_order_relaxed);
    static bool last_bypass = false;
    if (bypass != last_bypass) {
        last_bypass = bypass;
        s_processor.Reset();
    }

    if (bypass) {
        std::copy_n(input, static_cast<size_t>(frameCount) * 2, output);
    }
    else {
        float const total_drive_db = s_total_drive_db.load(std::memory_order_relaxed);
        float const contrast = s_contrast.load(std::memory_order_relaxed);
        std::array<float, kNumBands> drives_db;
        for (size_t band = 0; band < kNumBands; ++band) {
            // 系数在 0 附近会变成 -inf dB，用 -60 dB 作下限
            float const base = std::max(s_coefficients[band].load(std::memory_order_relaxed), 1.0e-3f);
            float const coefficient = ShapeCoefficient(base, contrast);
            drives_db[band] = total_drive_db + qwqdsp::convert::Gain2Db<-60.0f>(coefficient);
        }

        s_processor.SetParameters(drives_db, s_ceiling_db.load(std::memory_order_relaxed),
                                  s_release_ms.load(std::memory_order_relaxed),
                                  s_lookahead_ms.load(std::memory_order_relaxed));

        std::copy_n(input, static_cast<size_t>(frameCount) * 2, output);
        s_processor.Process(output, static_cast<size_t>(frameCount));
    }

    float input_peak = 0.0f;
    float output_peak = 0.0f;
    for (size_t i = 0; i < static_cast<size_t>(frameCount) * 2; ++i) {
        input_peak = std::max(input_peak, std::abs(input[i]));
        output_peak = std::max(output_peak, std::abs(output[i]));
    }

    s_input_peak.store(input_peak, std::memory_order_relaxed);
    s_output_peak.store(output_peak, std::memory_order_relaxed);
    for (size_t band = 0; band < kNumBands; ++band) {
        float const gain_reduction_db = bypass ? 0.0f : s_processor.GetBandGainReductionDb(band);
        s_band_gain_reduction_db[band].store(gain_reduction_db, std::memory_order_relaxed);
    }
    s_final_gain_reduction_db.store(bypass ? 0.0f : s_processor.GetFinalGainReductionDb(), std::memory_order_relaxed);
}

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------
int main(void) {
    using namespace erb_multiband_maximizer_demo;

    s_processor.Init(kSampleRate);
    BuildBandLabels();
    CreateKnobs();
    ApplyPreset(0);

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
    InitWindow(kWindowWidth, kWindowHeight, "ERB Multiband Maximizer");
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
            display_band_gain_reduction[band] = std::min(s_band_gain_reduction_db[band].load(std::memory_order_relaxed),
                                                         display_band_gain_reduction[band] * kDecay);
        }
        display_final_gain_reduction =
            std::min(s_final_gain_reduction_db.load(std::memory_order_relaxed), display_final_gain_reduction * kDecay);

        BeginDrawing();
        ClearBackground(kBackground);

        // ── 顶栏 ──
        DrawRectangle(0, 0, kWindowWidth, 64, Color{24, 26, 29, 255});
        DrawRectangle(0, 63, kWindowWidth, 1, Color{48, 52, 57, 255});
        DrawText("ERB MULTIBAND MAXIMIZER", 20, 14, 24, kTextColor);
        DrawText("12 ERB-SPACED BANDS -> 12 LIMITERS -> SUM -> FINAL LIMITER", 21, 43, 10, kDimColor);

        // ── 全局参数与主电平表 ──
        s_release_knob.display();
        s_lookahead_knob.display();
        DrawBypassButton(Rectangle{228.0f, 100.0f, 100.0f, 32.0f});

        DrawText("MASTER", 350, 76, 9, kDimColor);
        DrawVerticalMeter(350.0f, 90.0f, 24.0f, 88.0f, display_input, "IN", Color{76, 201, 151, 255});
        DrawVerticalMeter(386.0f, 90.0f, 24.0f, 88.0f, display_output, "OUT", Color{76, 201, 151, 255});

        DrawText("BANDS", 450, 88, 9, kDimColor);
        DrawText("12 x ERB (20 Hz - 20 kHz)", 450, 104, 12, Color{220, 223, 228, 255});
        DrawText("11 HP4 DIFFERENCES + FINAL LIMITER", 450, 124, 9, Color{110, 115, 122, 255});
        DrawGainReductionMeter(450.0f, 162.0f, 300.0f, 18.0f, display_final_gain_reduction, "FINAL GR",
                               Color{235, 235, 238, 255});

        DrawText("PRESET", 800, 88, 9, kDimColor);
        DrawText(kPresets[static_cast<size_t>(s_preset_index)].description, 800, 104, 12, kTextColor);
        DrawText("band drive = coefficient x total drive", 800, 126, 10, kFaintColor);

        // ── 频谱段面板 ──
        DrawRectangle(static_cast<int>(kPanelX), static_cast<int>(kPanelY), static_cast<int>(kPanelW),
                      static_cast<int>(kPanelH), kPanelColor);
        DrawRectangleLinesEx(Rectangle{kPanelX, kPanelY, kPanelW, kPanelH}, 1.0f, kPanelLine);

        // 左侧：总 drive / 曲线预设 / 所有 ceiling
        DrawRectangle(14, static_cast<int>(kPanelY) + 4, 150, static_cast<int>(kPanelH) - 8, Color{19, 21, 24, 255});
        DrawRectangleLinesEx(Rectangle{14.0f, kPanelY + 4.0f, 150.0f, kPanelH - 8.0f}, 1.0f, Color{32, 35, 40, 255});
        s_total_drive_knob.display();
        DrawPresetSelector(Rectangle{14.0f, 344.0f, 150.0f, 66.0f});
        s_ceiling_knob.display();
        s_contrast_knob.display();

        // 左侧：系数-频率曲线
        HandleCoefficientChartInput();
        DrawCoefficientChart();

        // 右侧：各 limiter 的 GR
        DrawText("GAIN REDUCTION (dB)", static_cast<int>(kGrAreaX), static_cast<int>(kChartY - 20.0f), 10, kDimColor);
        DrawGrGrid();
        for (size_t band = 0; band < kNumBands; ++band) {
            DrawBandGrMeter(band, display_band_gain_reduction[band]);
        }

        DrawText("drag knob: adjust   drag curve node: band coefficient   shift/alt: fine   right click: reset", 20,
                 kWindowHeight - 20, 10, kFaintColor);

        EndDrawing();
    }

    ma_device_uninit(&device);
    CloseWindow();
}
