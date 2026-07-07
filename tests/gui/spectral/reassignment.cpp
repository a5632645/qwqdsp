#include <algorithm>
#include <cmath>
#include <cstdio>
#include <span>

#include "miniaudio.h"
#include "raylib.h"

#include "reassignment/freq_reassignment_frame.hpp"
#include "reassignment/scrolling_image.hpp"
#include "reassignment/spectrogram_column.hpp"
#include "reassignment/spectrogram_frame.hpp"
#include "reassignment/tf_derivative_reassignment_frame.hpp"
#include "reassignment/tf_derivative_reassignment_frame_peak_filter.hpp"
#include "reassignment/tf_reassignment_frame.hpp"
#include "reassignment/time_reassignment_frame.hpp"

#include "reassignment/magma_colormap.hpp"
#include "reassignment/viridis_colormap.hpp"

// ── 窗口与 UI 常量 ──
static constexpr int kWindowWidth = 800;
static constexpr int kWindowHeight = 400;

// ── 画布布局 ──
static constexpr int kCanvasX = 60;
static constexpr int kCanvasY = 30;
static constexpr int kCanvasW = kWindowWidth - kCanvasX - 20;
static constexpr int kCanvasH = kWindowHeight - kCanvasY - 60;

// ── 频谱图参数 ──
static constexpr float kDbFloor = -60.0f;
static constexpr float kFreqMin = 20.0f;
static constexpr float kFreqMax = 20000.0f;

// ── 频率刻度标签 ──
static constexpr float kFreqTickLabels[] = {20, 200, 2000, 20000};
static const Color kGridColor = {50, 50, 50, 255};
static const Color kTextColor = {180, 180, 180, 255};
static const Color kBgColor = {20, 20, 20, 255};

static constexpr int kSampleRate = 48000;
static constexpr int kFftSize = 4096;
static constexpr int kHopSize = kFftSize / 4;

// 滚动时长: 全屏可见时间范围
static constexpr float kScrollSeconds = 3.0f;
static constexpr int kImageWidth = static_cast<int>(kScrollSeconds * kSampleRate / kHopSize);

// ----------------------------------------
// 全局状态
// ----------------------------------------

// #define USE_FREQ_REASSIGNMENT
static TimeReassignmentFrame<ViridisColormap> frame_;
static SpectrogramColumn column_;
static ScrollingImage image_;

// ----------------------------------------
// miniaudio 回调
// ----------------------------------------

extern "C" void MaCaptureCallback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;
    float const* src = reinterpret_cast<float const*>(pInput);
    float* dst = reinterpret_cast<float*>(pOutput);
    std::copy_n(src, frameCount, dst);

    column_.ProcessAudio({src, frameCount}, frame_, [&](std::span<const Color> col) { image_.PushColumn(col); });
}

// ----------------------------------------
// draw
// ----------------------------------------

static void DrawBackground() {
    DrawRectangle(kCanvasX, kCanvasY, kCanvasW, kCanvasH, kBgColor);
}

static void DrawSpectrogram() {
    image_.Draw(kCanvasX, kCanvasY, kCanvasW, kCanvasH);
}

static void DrawFreqGrid() {
    const float logMin = std::log10(kFreqMin);
    const float logMax = std::log10(kFreqMax);

    // 水平频率网格线 + 左侧标签
    for (float freq : kFreqTickLabels) {
        float norm = (std::log10(freq) - logMin) / (logMax - logMin);
        int y = kCanvasY + kCanvasH - static_cast<int>(norm * kCanvasH);
        DrawLine(kCanvasX, y, kCanvasX + kCanvasW, y, kGridColor);

        char label[16];
        if (freq >= 1000.0f)
            snprintf(label, sizeof(label), "%.0fk", freq / 1000.0f);
        else
            snprintf(label, sizeof(label), "%.0f", freq);
        int tw = MeasureText(label, 10);
        DrawText(label, kCanvasX - tw - 6, y - 5, 10, kTextColor);
    }
}

// ----------------------------------------
// main
// ----------------------------------------

int main(void) {
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, "Spectrogram - miniaudio + qwqdsp + raylib");
    SetTargetFPS(60);

    // ── miniaudio 全双工 ──
    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.playback.format = ma_format_f32;
    config.playback.channels = 1;
    config.sampleRate = static_cast<ma_uint32>(kSampleRate);
    config.dataCallback = MaCaptureCallback;
    config.pUserData = nullptr;
    config.periodSizeInMilliseconds = 10;

    ma_device device;
    ma_result result = ma_device_init(nullptr, &config, &device);
    if (result == MA_SUCCESS)
        ma_device_start(&device);
    else
        TraceLog(LOG_WARNING, "miniaudio 捕获设备初始化失败，以静默模式运行");

    // ── 初始化 ──
    column_.Init(kCanvasH, kSampleRate, kFftSize, kHopSize);
#ifdef USE_FREQ_REASSIGNMENT
    frame_.Init(kSampleRate, kFftSize, 1, kCanvasH, kFreqMin, kFreqMax, kDbFloor);
#else
    frame_.Init(kSampleRate, kFftSize, kHopSize, 4, kCanvasH, kFreqMin, kFreqMax, kDbFloor);
#endif
    image_.Init(kImageWidth, kCanvasH);

    // ── 主循环 ──
    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(BLACK);

        DrawBackground();
        DrawSpectrogram();
        DrawFreqGrid();

        DrawFPS(kWindowWidth - 80, 10);
        EndDrawing();
    }

    // ── 清理 ──
    if (result == MA_SUCCESS) {
        ma_device_stop(&device);
        ma_device_uninit(&device);
    }
    image_.Unload();
    CloseWindow();
    return 0;
}
