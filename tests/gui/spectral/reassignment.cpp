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
#include "reassignment/tf_derivative_reassignment_frame_conv.hpp"
#include "reassignment/tf_phase_vocoder_reassignment_frame.hpp"
#include "reassignment/tf_phase_vocoder_reassignment_frame_conv.hpp"
#include "reassignment/tf_reassignment_frame.hpp"
#include "reassignment/time_reassignment_frame.hpp"

#include "reassignment/jet_colormap.hpp"
#include "reassignment/magma_colormap.hpp"
#include "reassignment/viridis_colormap.hpp"

// ── 音频参数 ──
static constexpr int kSampleRate = 48000;
static constexpr int kFftSize = 4096;
static constexpr int kHopSize = kFftSize / 16;

// ── zeroPad ──
// 简单类（标准谱图/时间重分配）：较大 zeroPad 提升频率分辨率
static constexpr int kZeroPadSimple = 4;
// 全 TF 类（子列缓冲已提供时间精度）
static constexpr int kZeroPadFull = 1;

// ── 频谱图显示 ──
static constexpr float kDbFloor = -72.0f;
static constexpr float kFreqMin = 20.0f;
static constexpr float kFreqMax = 20000.0f;

// ── 滚动图像 ──
static constexpr float kScrollSeconds = 3.0f;
static constexpr int kImageWidth = static_cast<int>(kScrollSeconds * kSampleRate / kHopSize);

// ── UI 布局 ──
static constexpr int kWindowWidth = 800;
static constexpr int kWindowHeight = 400;
static constexpr int kCanvasX = 60;
static constexpr int kCanvasY = 30;
static constexpr int kCanvasW = kWindowWidth - kCanvasX - 20;
static constexpr int kCanvasH = kWindowHeight - kCanvasY - 60;

// ── UI 颜色与标签 ──
static constexpr float kFreqTickLabels[] = {20, 200, 2000, 20000};
static const Color kGridColor = {50, 50, 50, 255};
static const Color kTextColor = {180, 180, 180, 255};
static const Color kBgColor = {20, 20, 20, 255};

// ----------------------------------------
// 全局状态
// ----------------------------------------

enum class FrameType : int {
    kSpectrogram,
    kFreqReassignment,
    kTimeReassignment,
    kTfReassignment,
    kTfPhaseVocoder,
    kTfPhaseVocoderPeak,
    kTfDerivative,
    kTfDerivativePeak,
    kTfPhaseVocoderConv,
    kTfDerivativeConv,
    kCount
};

static FrameType g_frame_type = FrameType::kTfPhaseVocoderPeak;

static constexpr const char* kFrameNames[] = {
    "Spectrogram",
    "Frequency Reassignment",
    "Time Reassignment",
    "Time-Frequency Reassignment",
    "Phase Vocoder Reassignment",
    "Phase Vocoder + Peak Filter",
    "Derivative Reassignment",
    "Derivative + Peak Filter",
    "Phase Vocoder + Convergence",
    "Derivative + Convergence",
};
static_assert(std::size(kFrameNames) == static_cast<int>(FrameType::kCount));

static SpectrogramFrame<MagmaColormap> f_sp;
static FreqReassignmentFrame<MagmaColormap> f_freq;
static TimeReassignmentFrame<MagmaColormap> f_time;
static TfReassignmentFrame<MagmaColormap> f_tf;
static TfPhaseVocoderReassignmentFrame<MagmaColormap, false> f_pv;
static TfPhaseVocoderReassignmentFrame<MagmaColormap, true> f_pv_pk;
static TfDerivativeReassignmentFrame<MagmaColormap, false> f_deriv;
static TfDerivativeReassignmentFrame<MagmaColormap, true> f_deriv_pk;
static TfPhaseVocoderReassignmentFrameConv<MagmaColormap> f_pv_conv;
static TfDerivativeReassignmentFrameConv<MagmaColormap> f_deriv_conv;

static SpectrogramColumn column_;
static ScrollingImage image_;

// ----------------------------------------
// miniaudio 回调
// ----------------------------------------

extern "C" void MaCaptureCallback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;
    (void)pOutput;
    float const* src = reinterpret_cast<float const*>(pInput);
    auto push = [&](std::span<const Color> col) { image_.PushColumn(col); };

    switch (g_frame_type) {
        case FrameType::kSpectrogram:
            column_.ProcessAudio({src, frameCount}, f_sp, push);
            break;
        case FrameType::kFreqReassignment:
            column_.ProcessAudio({src, frameCount}, f_freq, push);
            break;
        case FrameType::kTimeReassignment:
            column_.ProcessAudio({src, frameCount}, f_time, push);
            break;
        case FrameType::kTfReassignment:
            column_.ProcessAudio({src, frameCount}, f_tf, push);
            break;
        case FrameType::kTfPhaseVocoder:
            column_.ProcessAudio({src, frameCount}, f_pv, push);
            break;
        case FrameType::kTfPhaseVocoderPeak:
            column_.ProcessAudio({src, frameCount}, f_pv_pk, push);
            break;
        case FrameType::kTfDerivative:
            column_.ProcessAudio({src, frameCount}, f_deriv, push);
            break;
        case FrameType::kTfDerivativePeak:
            column_.ProcessAudio({src, frameCount}, f_deriv_pk, push);
            break;
        case FrameType::kTfPhaseVocoderConv:
            column_.ProcessAudio({src, frameCount}, f_pv_conv, push);
            break;
        case FrameType::kTfDerivativeConv:
            column_.ProcessAudio({src, frameCount}, f_deriv_conv, push);
            break;
        default:
            break;
    }
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
    if (result == MA_SUCCESS)
        ma_device_start(&device);
    else
        TraceLog(LOG_WARNING, "miniaudio 捕获设备初始化失败，以静默模式运行");

    // ── 初始化 ──
    column_.Init(kCanvasH, kSampleRate, kFftSize, kHopSize);
    f_sp.Init(kSampleRate, kFftSize, kZeroPadSimple, kCanvasH, kFreqMin, kFreqMax, kDbFloor);
    f_freq.Init(kSampleRate, kFftSize, kZeroPadFull, kCanvasH, kFreqMin, kFreqMax, kDbFloor);
    f_time.Init(kSampleRate, kFftSize, kHopSize, kZeroPadSimple, kCanvasH, kFreqMin, kFreqMax, kDbFloor);
    f_tf.Init(kSampleRate, kFftSize, kHopSize, kZeroPadFull, kCanvasH, kFreqMin, kFreqMax, kDbFloor);
    f_pv.Init(kSampleRate, kFftSize, kHopSize, kZeroPadFull, kCanvasH, kFreqMin, kFreqMax, kDbFloor);
    f_pv_pk.Init(kSampleRate, kFftSize, kHopSize, kZeroPadFull, kCanvasH, kFreqMin, kFreqMax, kDbFloor);
    f_deriv.Init(kSampleRate, kFftSize, kHopSize, kZeroPadFull, kCanvasH, kFreqMin, kFreqMax, kDbFloor);
    f_deriv_pk.Init(kSampleRate, kFftSize, kHopSize, kZeroPadFull, kCanvasH, kFreqMin, kFreqMax, kDbFloor);
    f_pv_conv.Init(kSampleRate, kFftSize, kHopSize, kZeroPadFull, kCanvasH, kFreqMin, kFreqMax, kDbFloor);
    f_deriv_conv.Init(kSampleRate, kFftSize, kHopSize, kZeroPadFull, kCanvasH, kFreqMin, kFreqMax, kDbFloor);
    image_.Init(kImageWidth, kCanvasH);

    // ── 主循环 ──
    while (!WindowShouldClose()) {
        // ── 键盘切换 ──
        int ch = GetCharPressed();
        if (ch >= '1' && ch <= '9') {
            auto t = static_cast<FrameType>(ch - '1');
            if (t >= FrameType::kSpectrogram && t < FrameType::kCount)
                g_frame_type = t;
        }
        else if (ch == '0') {
            g_frame_type = FrameType::kTfDerivativeConv;
        }

        BeginDrawing();
        ClearBackground(BLACK);

        DrawBackground();
        DrawSpectrogram();
        DrawFreqGrid();

        // ── 参数信息 ──
        {
            char buf[96];
            constexpr int kOverlapPct = 100 - kHopSize * 100 / kFftSize;
            const int active_zero_pad =
                (g_frame_type == FrameType::kSpectrogram || g_frame_type == FrameType::kTimeReassignment)
                    ? kZeroPadSimple
                    : kZeroPadFull;
            snprintf(buf, sizeof(buf), "FFT: %d | Hop: %d (%d%%) | ZeroPad: %d | Noise Floor: %.0f dB", kFftSize,
                     kHopSize, kOverlapPct, active_zero_pad, kDbFloor);
            DrawText(buf, kCanvasX + 4, 10, 10, kTextColor);
        }

        // ── 算法切换面板 ──
        {
            constexpr Color kActiveColor = {255, 220, 60, 255};
            constexpr Color kInactiveColor = {120, 120, 120, 255};
            constexpr int kRowH = 14;
            constexpr int kCols = 4;
            constexpr int kCellW = kCanvasW / kCols;
            const int active = static_cast<int>(g_frame_type);

            // 0 1 2 3
            // 4 5 6
            // 7 8 9
            constexpr int kCountPerRow[] = {4, 3, 3};
            int idx = 0;
            for (int row = 0; row < 3; ++row) {
                for (int col = 0; col < kCountPerRow[row]; ++col, ++idx) {
                    char buf[48];
                    int key = (idx < 9) ? ('1' + idx) : '0';
                    snprintf(buf, sizeof(buf), "%c:%s", key, kFrameNames[idx]);
                    int x = kCanvasX + col * kCellW;
                    int y = kCanvasY + kCanvasH + 4 + row * kRowH;
                    DrawText(buf, x, y, 10, (idx == active) ? kActiveColor : kInactiveColor);
                }
            }
        }
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
