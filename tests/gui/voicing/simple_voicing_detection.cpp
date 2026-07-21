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

// ------------------------------------------------------------
//  全局状态
// ------------------------------------------------------------
static float g_voicing_prob = 0.0f;
static float g_threshold = 0.5f;
static float g_unvoiced_gain_db = 10.0f;
static float g_voiced_gain_db = -10.0f;

// ---- 检测器 ----
static SubbandVoicingDetector g_detector(kSampleRate);

/// 将检测器参数同步为 FixedVoicingDetector 默认值
static void InitDetectorDefaults() noexcept {
    g_detector.setLpFreq(2000.0f);
    g_detector.setHpFreq(5000.0f);
    g_detector.setHpScale(1.0f);
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

    for (ma_uint32 i = 0; i < frameCount; ++i) {
        float sample = src[i];

        g_detector.processSample(sample);

        g_voicing_prob = g_detector.simpleProbability();

        // ----- 输出: 当前增益 × 采样 -----
        float gain_db = (g_voicing_prob >= g_threshold) ? g_voiced_gain_db : g_unvoiced_gain_db;
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

    // ----- Row 2: 滤波器频率 -----
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

    // ----- Row 3: RMS 时间 -----
    Knob rms_tau_knob;
    rms_tau_knob.set_title("RMS Tau");
    rms_tau_knob.set_range(0.001f, 0.2f, 0.001f, 0.005f);
    rms_tau_knob.set_bound(260, 260, 200, 50);
    rms_tau_knob.value_to_text_function = [](float v) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "%.0f ms", v * 1000.0f);
        return buf;
    };
    rms_tau_knob.on_value_change = [](float v) { g_detector.setRmsTau(v); };

    // ----- Row 4: 增益 -----
    Knob unvoiced_gain_knob;
    unvoiced_gain_knob.set_title("Unvoiced Gain");
    unvoiced_gain_knob.set_range(-60.0f, 12.0f, 0.5f, 10.0f);
    unvoiced_gain_knob.set_bound(30, 330, 200, 50);
    unvoiced_gain_knob.value_to_text_function = [](float v) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "%.1f dB", v);
        return buf;
    };
    unvoiced_gain_knob.on_value_change = [](float v) { g_unvoiced_gain_db = v; };

    Knob voiced_gain_knob;
    voiced_gain_knob.set_title("Voiced Gain");
    voiced_gain_knob.set_range(-60.0f, 12.0f, 0.5f, -10.0f);
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
        hp_scale_knob.display();
        lp_freq_knob.display();
        hp_freq_knob.display();
        rms_tau_knob.display();
        unvoiced_gain_knob.display();
        voiced_gain_knob.display();

        // V/UV 状态大文字
        bool is_voiced = prob >= threshold;
        const char* status_str = is_voiced ? "VOICED" : "UNVOICED";
        Color status_color = is_voiced ? GREEN : SKYBLUE;
        DrawText(status_str, 30, 70, 36, status_color);

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
