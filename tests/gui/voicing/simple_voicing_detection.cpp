#include <cstdio>

#include "miniaudio.h"
#include "raylib.h"
#include "slider.hpp"

#include "subband_voicing_detector.hpp"

// ------------------------------------------------------------
//  常量
// ------------------------------------------------------------

static constexpr float kSampleRate = 48000.0f;
static constexpr int kWindowWidth = 800;
static constexpr int kWindowHeight = 480;

// ---- 概率条画布 ----
static constexpr int kProbX = 30;
static constexpr int kProbY = 20;
static constexpr int kProbW = 740;
static constexpr int kProbH = 24;

// ---- 帧分析参数 ----
static constexpr size_t kFrameSize = 2048; // ~42.7ms

// ------------------------------------------------------------
//  全局状态
// ------------------------------------------------------------
static float g_voicing_prob = 0.0f;
static float g_energy_ratio = 0.0f;
static float g_threshold = 0.5f;
static float g_unvoiced_gain_db = -6.0f;
static float g_voiced_gain_db = 0.0f;

// ---- 检测器 ----
static SubbandVoicingDetector g_detector(kSampleRate);

/// 将检测器参数同步为 FixedVoicingDetector 默认值
static void InitDetectorDefaults() noexcept {
    g_detector.setLpFreq(500.0f);
    g_detector.setHpFreq(5000.0f);
    g_detector.setDelta(0.0001f);
    g_detector.setCenter(0.0f);
    g_detector.setWidth(0.1f);
    g_detector.setRmsTau(0.005f);
}

// ------------------------------------------------------------
//  miniaudio 回调 — 捕获 → 检测 → 增益 → 输出
// ------------------------------------------------------------

extern "C" void MaDuplexCallback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;

    auto* src = static_cast<const float*>(pInput);
    auto* dst = static_cast<float*>(pOutput);
    if (src == nullptr || dst == nullptr)
        return;

    static size_t s_frame_count = 0;

    for (ma_uint32 i = 0; i < frameCount; ++i) {
        float sample = src[i];

        g_detector.processSample(sample);
        ++s_frame_count;

        if (s_frame_count == kFrameSize) {
            float energy_ratio, rms_linear;
            float prob = g_detector.frameResult(energy_ratio, rms_linear);

            g_voicing_prob = prob;
            g_energy_ratio = energy_ratio;

            s_frame_count = 0;
        }

        // ----- 输出: 当前增益 × 采样 -----
        float prob = g_detector.lastProbability();
        float gain_db = (prob >= g_threshold) ? g_voiced_gain_db : g_unvoiced_gain_db;
        dst[i] = sample * std::pow(10.0f, gain_db / 20.0f);
    }
}

// ------------------------------------------------------------
//  绘图
// ------------------------------------------------------------

/// 绘制概率条
static void DrawProbabilityBar(float prob, float threshold) {
    // 背景
    DrawRectangle(kProbX, kProbY, kProbW, kProbH, {30, 30, 30, 255});

    // 概率填充
    int fill_w = static_cast<int>(prob * kProbW);
    Color fill_color;
    if (prob >= threshold) {
        // 浊音 — 绿色渐变
        float t = std::clamp((prob - threshold) / (1.0f - threshold), 0.0f, 1.0f);
        fill_color = {static_cast<unsigned char>(60 - 40 * t), static_cast<unsigned char>(180 + 75 * t),
                      static_cast<unsigned char>(60 - 40 * t), 255};
    }
    else {
        // 清音 — 蓝色渐变
        float t = std::clamp(prob / threshold, 0.0f, 1.0f);
        fill_color = {static_cast<unsigned char>(60 + 40 * t), static_cast<unsigned char>(60 + 40 * t),
                      static_cast<unsigned char>(180 + 75 * t), 255};
    }
    if (fill_w > 0)
        DrawRectangle(kProbX, kProbY, fill_w, kProbH, fill_color);

    // 阈值线
    int thresh_x = kProbX + static_cast<int>(threshold * kProbW);
    DrawLine(thresh_x, kProbY, thresh_x, kProbY + kProbH, {255, 220, 60, 255});

    // 边框
    DrawRectangleLines(kProbX, kProbY, kProbW, kProbH, {80, 80, 80, 255});

    // 标签
    char label[64];
    snprintf(label, sizeof(label), "Voicing: %.1f%%", prob * 100.0f);
    DrawText(label, kProbX + 4, kProbY + 4, 14, WHITE);
}

// ---- 映射曲线画布 ----
static constexpr int kCurveX = 30;
static constexpr int kCurveY = 400;
static constexpr int kCurveW = 740;
static constexpr int kCurveH = 60;

/// 绘制 ratio → prob 映射曲线
static void DrawMappingCurve(float center, float width, float current_ratio, float threshold) {
    // 背景
    DrawRectangle(kCurveX, kCurveY, kCurveW, kCurveH, {20, 20, 20, 255});

    // 网格: 概率 0.0, 0.5, 1.0
    for (int row = 0; row <= 2; ++row) {
        int y = kCurveY + kCurveH - static_cast<int>(static_cast<float>(row) * 0.5f * kCurveH);
        DrawLine(kCurveX, y, kCurveX + kCurveW, y, {50, 50, 50, 255});
    }

    // 阈值线 (概率域)
    int thresh_y = kCurveY + kCurveH - static_cast<int>(threshold * kCurveH);
    DrawLine(kCurveX, thresh_y, kCurveX + kCurveW, thresh_y, {80, 50, 20, 255});

    // 采样曲线: log10(ratio) ∈ [-3, 3]
    static constexpr int kNumSamples = 200;
    static constexpr float kLogMin = -3.0f;
    static constexpr float kLogMax = 3.0f;
    Vector2 points[kNumSamples];

    for (int i = 0; i < kNumSamples; ++i) {
        float t = static_cast<float>(i) / (kNumSamples - 1);
        float log_r = kLogMin + t * (kLogMax - kLogMin);
        float ratio = std::pow(10.0f, log_r);
        float p = RatioToProbability(ratio, center, width);

        float x = kCurveX + t * kCurveW;
        float y = kCurveY + kCurveH - p * kCurveH;
        points[i] = {x, std::clamp(y, static_cast<float>(kCurveY), static_cast<float>(kCurveY + kCurveH - 1))};
    }
    DrawLineStrip(points, kNumSamples, {200, 200, 100, 255});

    // 当前 ratio 位置标记
    float cur_log_r = std::log10(current_ratio + 1e-10f);
    float cur_t = std::clamp((cur_log_r - kLogMin) / (kLogMax - kLogMin), 0.0f, 1.0f);
    int marker_x = kCurveX + static_cast<int>(cur_t * kCurveW);
    float cur_p = RatioToProbability(current_ratio, center, width);
    int marker_y = kCurveY + kCurveH - static_cast<int>(cur_p * kCurveH);

    // 垂直线
    DrawLine(marker_x, kCurveY, marker_x, kCurveY + kCurveH, {255, 100, 100, 180});
    // 标记点
    DrawCircle(marker_x, marker_y, 4, {255, 60, 60, 255});

    // 标签
    char buf[64];
    snprintf(buf, sizeof(buf), "ratio=%.2f  prob=%.2f", current_ratio, cur_p);
    DrawText(buf, kCurveX + 4, kCurveY + 2, 10, {180, 180, 180, 255});
    DrawText("log10(ratio)", kCurveX + kCurveW - 80, kCurveY + kCurveH - 12, 9, {120, 120, 120, 255});
    DrawText("0", kCurveX - 14, kCurveY + kCurveH - 12, 9, {120, 120, 120, 255});
    DrawText("1", kCurveX - 14, kCurveY + 2, 9, {120, 120, 120, 255});
}

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------

int main(void) {
    InitDetectorDefaults();

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, "V/UV Detection - Sub-band Energy Ratio  |  miniaudio + raylib");
    SetTargetFPS(60);

    // ----- Knobs -----

    // ----- Row 1: 阈值 + 子带缩放 -----
    Knob threshold_knob;
    threshold_knob.set_title("Threshold");
    threshold_knob.set_range(0.0f, 1.0f, 0.01f, 0.5f);
    threshold_knob.set_bound(30, 120, 200, 50);
    threshold_knob.value_to_text_function = [](float v) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "%.2f", v);
        return buf;
    };
    threshold_knob.on_value_change = [](float v) { g_threshold = v; };

    Knob lp_scale_knob;
    lp_scale_knob.set_title("LP Scale");
    lp_scale_knob.set_range(0.0f, 5.0f, 0.05f, 1.0f);
    lp_scale_knob.set_bound(260, 120, 200, 50);
    lp_scale_knob.value_to_text_function = [](float v) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "%.2f", v);
        return buf;
    };
    lp_scale_knob.on_value_change = [](float v) { g_detector.setLpScale(v); };

    Knob hp_scale_knob;
    hp_scale_knob.set_title("HP Scale");
    hp_scale_knob.set_range(0.0f, 5.0f, 0.05f, 1.0f);
    hp_scale_knob.set_bound(490, 120, 200, 50);
    hp_scale_knob.value_to_text_function = [](float v) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "%.2f", v);
        return buf;
    };
    hp_scale_knob.on_value_change = [](float v) { g_detector.setHpScale(v); };

    // ----- Row 2: 滤波器频率 + Delta -----
    Knob lp_freq_knob;
    lp_freq_knob.set_title("LP Freq");
    lp_freq_knob.set_range(100.0f, 3000.0f, 10.0f, 500.0f);
    lp_freq_knob.set_bound(30, 190, 200, 50);
    lp_freq_knob.value_to_text_function = [](float v) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "%.0f Hz", v);
        return buf;
    };
    lp_freq_knob.on_value_change = [](float v) { g_detector.setLpFreq(v); };

    Knob hp_freq_knob;
    hp_freq_knob.set_title("HP Freq");
    hp_freq_knob.set_range(500.0f, 8000.0f, 10.0f, 5000.0f);
    hp_freq_knob.set_bound(260, 190, 200, 50);
    hp_freq_knob.value_to_text_function = [](float v) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "%.0f Hz", v);
        return buf;
    };
    hp_freq_knob.on_value_change = [](float v) { g_detector.setHpFreq(v); };

    Knob delta_knob;
    delta_knob.set_title("Delta");
    delta_knob.set_range(0.0f, 0.01f, 0.0001f, 0.0001f);
    delta_knob.set_bound(490, 190, 200, 50);
    delta_knob.value_to_text_function = [](float v) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "%.4f", v);
        return buf;
    };
    delta_knob.on_value_change = [](float v) { g_detector.setDelta(v); };

    // ----- Row 3: 概率映射 + RMS 时间 -----
    Knob center_knob;
    center_knob.set_title("Center");
    center_knob.set_range(-2.0f, 4.0f, 0.05f, 0.0f);
    center_knob.set_bound(30, 260, 200, 50);
    center_knob.value_to_text_function = [](float v) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "%.2f", v);
        return buf;
    };
    center_knob.on_value_change = [](float v) { g_detector.setCenter(v); };

    Knob width_knob;
    width_knob.set_title("Width");
    width_knob.set_range(0.1f, 5.0f, 0.05f, 0.1f);
    width_knob.set_bound(260, 260, 200, 50);
    width_knob.value_to_text_function = [](float v) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "%.2f", v);
        return buf;
    };
    width_knob.on_value_change = [](float v) { g_detector.setWidth(v); };

    Knob rms_tau_knob;
    rms_tau_knob.set_title("RMS Tau");
    rms_tau_knob.set_range(0.005f, 0.2f, 0.005f, 0.005f);
    rms_tau_knob.set_bound(490, 260, 200, 50);
    rms_tau_knob.value_to_text_function = [](float v) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "%.0f ms", v * 1000.0f);
        return buf;
    };
    rms_tau_knob.on_value_change = [](float v) { g_detector.setRmsTau(v); };

    // ----- Row 4: 增益 -----
    Knob unvoiced_gain_knob;
    unvoiced_gain_knob.set_title("Unvoiced Gain");
    unvoiced_gain_knob.set_range(-60.0f, 12.0f, 0.5f, -6.0f);
    unvoiced_gain_knob.set_bound(30, 330, 200, 50);
    unvoiced_gain_knob.value_to_text_function = [](float v) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "%.1f dB", v);
        return buf;
    };
    unvoiced_gain_knob.on_value_change = [](float v) { g_unvoiced_gain_db = v; };

    Knob voiced_gain_knob;
    voiced_gain_knob.set_title("Voiced Gain");
    voiced_gain_knob.set_range(-60.0f, 12.0f, 0.5f, 0.0f);
    voiced_gain_knob.set_bound(260, 330, 200, 50);
    voiced_gain_knob.value_to_text_function = [](float v) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "%.1f dB", v);
        return buf;
    };
    voiced_gain_knob.on_value_change = [](float v) { g_voiced_gain_db = v; };

    // ----- miniaudio 双工 (捕获+播放) -----
    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.playback.format = ma_format_f32;
    config.playback.channels = 1;
    config.sampleRate = static_cast<ma_uint32>(kSampleRate);
    config.dataCallback = MaDuplexCallback;
    config.pUserData = nullptr;
    config.periodSizeInMilliseconds = 10;

    ma_device device;
    ma_result result = ma_device_init(nullptr, &config, &device);
    if (result == MA_SUCCESS) {
        ma_device_start(&device);
    }
    else {
        TraceLog(LOG_WARNING, "miniaudio 双工设备初始化失败，以静默模式运行");
    }

    // ----- 主循环 (分析/判决/增益全在音频回调中完成) -----
    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(BLACK);

        // 概率条
        float threshold = g_threshold;
        float prob = g_voicing_prob;
        DrawProbabilityBar(prob, threshold);

        // Knobs
        threshold_knob.display();
        lp_scale_knob.display();
        hp_scale_knob.display();
        lp_freq_knob.display();
        hp_freq_knob.display();
        delta_knob.display();
        center_knob.display();
        width_knob.display();
        rms_tau_knob.display();
        unvoiced_gain_knob.display();
        voiced_gain_knob.display();

        // V/UV 状态大文字
        bool is_voiced = prob >= threshold;
        const char* status_str = is_voiced ? "VOICED" : "UNVOICED";
        Color status_color = is_voiced ? GREEN : SKYBLUE;
        DrawText(status_str, 30, 70, 36, status_color);

        // 映射曲线
        float center = g_detector.center();
        float width = g_detector.width();
        DrawMappingCurve(center, width, g_energy_ratio, threshold);

        DrawFPS(kWindowWidth - 80, 10);
        EndDrawing();
    }

    // ----- 清理 -----
    if (result == MA_SUCCESS) {
        ma_device_stop(&device);
        ma_device_uninit(&device);
    }
    CloseWindow();
    return 0;
}
