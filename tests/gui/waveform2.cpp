#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdio>
#include <span>
#include <vector>

#include "miniaudio.h"
#include "raylib.h"
#include "slider.hpp"

// ── 常量 ──
static constexpr float kSampleRate = 48000.0f;
static constexpr int kWindowWidth = 640;
static constexpr int kWindowHeight = 320;
static constexpr int kCanvasX = 30;
static constexpr int kCanvasY = 20;
static constexpr int kCanvasW = 580;
static constexpr int kCanvasH = 220;
static constexpr int kMaxTimeMs = 100;

// 最大时基采样数 + 查找范围
static constexpr size_t kMaxSamples = static_cast<size_t>(kSampleRate * kMaxTimeMs / 1000.0f);
static constexpr size_t kSearchRange = static_cast<size_t>(kSampleRate * 50 / 1000.0f); // 零交叉查找范围 (±150 采样)
static constexpr size_t kBufSize = kMaxSamples + kSearchRange;

// ── 全局状态 (音频线程 ↔ 主线程) ──
static std::vector<float> g_circ_buf(kBufSize);
static std::atomic<size_t> g_circ_pos{0};
static std::atomic<float> g_timebase_ms{20.0f};

// ============================================================
//  miniaudio 回调 — 仅写入环形缓冲，不处理触发
// ============================================================
extern "C" void MaCaptureCallback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;
    (void)pOutput;

    auto* src = static_cast<const float*>(pInput);
    if (src == nullptr)
        return;

    for (ma_uint32 i = 0; i < frameCount; ++i) {
        size_t pos = g_circ_pos.load(std::memory_order_relaxed);
        g_circ_buf[pos] = src[i];
        g_circ_pos.store((pos + 1) % kBufSize, std::memory_order_release);
    }
}

// ============================================================
//  绘制
// ============================================================

/// 绘制深色背景 + 网格 (振幅标度)
static void DrawGrid() {
    DrawRectangle(kCanvasX, kCanvasY, kCanvasW, kCanvasH, {20, 20, 20, 255});

    static constexpr Color kGridColor{50, 50, 50, 255};
    static constexpr float kAmpTicks[] = {-1.0f, -0.5f, 0.0f, 0.5f, 1.0f};
    for (float amp : kAmpTicks) {
        int y = kCanvasY + kCanvasH - static_cast<int>((amp + 1.0f) * 0.5f * kCanvasH);
        DrawLine(kCanvasX, y, kCanvasX + kCanvasW, y, kGridColor);
    }

    // 0 轴线用稍亮的颜色突出
    int y0 = kCanvasY + kCanvasH / 2;
    DrawLine(kCanvasX, y0, kCanvasX + kCanvasW, y0, {80, 80, 80, 255});
}

/// 绘制波形 (LineStrip)
static void DrawWaveform(std::span<const float> buf, size_t len) {
    if (len < 2)
        return;

    static constexpr Color kWaveColor{120, 200, 255, 255};
    std::vector<Vector2> points(len);
    float inv_len_m1 = 1.0f / static_cast<float>(len - 1);
    for (size_t i = 0; i < len; ++i) {
        float x = static_cast<float>(kCanvasX) + (static_cast<float>(i) * inv_len_m1) * kCanvasW;
        float y = kCanvasY + kCanvasH - ((buf[i] + 1.0f) * 0.5f * kCanvasH);
        points[i] = {x, std::clamp(y, static_cast<float>(kCanvasY), static_cast<float>(kCanvasY + kCanvasH - 1))};
    }
    DrawLineStrip(points.data(), static_cast<int>(len), kWaveColor);
}

// ============================================================
//  main
// ============================================================

int main(void) {
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, "Waveform2 - Zero-crossing trigger + raylib");
    SetTargetFPS(60);

    // ── Knob: 时基 ──
    Knob timebase_knob;
    timebase_knob.set_title("Timebase");
    timebase_knob.set_range(1.0f, static_cast<float>(kMaxTimeMs), 1.0f, 20.0f);
    timebase_knob.set_bound(30, 255, 280, 50);
    timebase_knob.value_to_text_function = [](float v) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "%.0f ms", v);
        return buf;
    };
    timebase_knob.on_value_change = [](float v) { g_timebase_ms.store(v, std::memory_order_relaxed); };

    // ── miniaudio 回环捕获 ──
    ma_device_config config = ma_device_config_init(ma_device_type_loopback);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.sampleRate = static_cast<ma_uint32>(kSampleRate);
    config.dataCallback = MaCaptureCallback;
    config.pUserData = nullptr;
    config.periodSizeInMilliseconds = 10;

    ma_device device;
    ma_result result = ma_device_init(nullptr, &config, &device);
    if (result == MA_SUCCESS) {
        ma_device_start(&device);
    }
    else {
        TraceLog(LOG_WARNING, "miniaudio 捕获设备初始化失败，以静默模式运行");
    }

    // 主线程本地显示缓冲
    std::vector<float> local_buf(kBufSize);
    size_t local_len = static_cast<size_t>(kSampleRate * 20.0f / 1000.0f);

    // ── 主循环 ──
    while (!WindowShouldClose()) {
        // 更新 local_len (时基可能已被 Knob 改变)
        local_len = static_cast<size_t>(kSampleRate * g_timebase_ms.load(std::memory_order_relaxed) / 1000.0f);
        local_len = std::min(local_len, kMaxSamples);
        size_t half = local_len / 2;

        // 从环形缓冲中搜索正向过零点，居中显示
        size_t circ_pos = g_circ_pos.load(std::memory_order_acquire);
        if (circ_pos > 0) {
            // 查找中心：时基/2 之前的位置
            size_t search_center = (circ_pos - half + kBufSize) % kBufSize;

            // 在查找范围内搜索正向过零点
            int found_offset = 0; // 0 = 未找到，居中在 search_center
            int search_half = static_cast<int>(kSearchRange) / 2;
            // 不要搜索到 circ_pos 附近（可能正在被写入）
            int max_back = std::min(search_half, static_cast<int>(half) - 2);
            if (max_back > 0) {
                for (int off = -max_back; off <= max_back; ++off) {
                    size_t idx = (search_center + static_cast<size_t>(off) + kBufSize) % kBufSize;
                    size_t prev = (idx + kBufSize - 1) % kBufSize;
                    if (g_circ_buf[prev] < 0.0f && g_circ_buf[idx] >= 0.0f) {
                        found_offset = off;
                        break;
                    }
                }
            }

            // 以找到的过零点（或 search_center）为中心，拷贝 display_len 个采样
            size_t center = (search_center + static_cast<size_t>(found_offset) + kBufSize) % kBufSize;
            size_t start = (center + kBufSize - half) % kBufSize;
            for (size_t j = 0; j < local_len; ++j) {
                local_buf[j] = g_circ_buf[(start + j) % kBufSize];
            }
        }

        BeginDrawing();
        ClearBackground(BLACK);

        DrawGrid();
        if (local_len > 0) {
            DrawWaveform(local_buf, local_len);
        }

        // 状态文字
        char info[64];
        snprintf(info, sizeof(info), "Zero-crossing trigger  |  Timebase: %.0f ms",
                 g_timebase_ms.load(std::memory_order_relaxed));
        DrawText(info, kCanvasX, kCanvasY + kCanvasH + 4, 10, {180, 180, 180, 255});

        timebase_knob.display();

        DrawFPS(kWindowWidth - 80, 10);
        EndDrawing();
    }

    // ── 清理 ──
    if (result == MA_SUCCESS) {
        ma_device_stop(&device);
        ma_device_uninit(&device);
    }
    CloseWindow();
    return 0;
}
