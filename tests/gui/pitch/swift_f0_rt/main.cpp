#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <deque>
#include <span>
#include <vector>

#include "miniaudio.h"
#include "raylib.h"

#include <qwqdsp/colormap/colormap.hpp>

#include "../../spectral/reassignment/colormap_adapter.hpp"
#include "../../spectral/reassignment/nc_reassignment_frame.hpp"
#include "../../spectral/reassignment/scrolling_image.hpp"
#include "../../spectral/reassignment/spectrogram_column.hpp"

#include "decimator.hpp"
#include "pitch_frontend.hpp"

// ------------------------------------------------------------
//  常量
// ------------------------------------------------------------

static constexpr int kSampleRate = 48000;

// ---- 背景: NC 方法（无时间重分配）----
static constexpr int kFftSize = 4096;
static constexpr int kHopSize = kFftSize / 16;
static constexpr int kNcZeroPad = 2;
static constexpr float kDbFloor = -72.0f;
static constexpr float kFreqMin = 20.0f;
static constexpr float kFreqMax = 20000.0f;

// ---- Swift F0: 48k → 16k ----
static constexpr int kDecimateFactor = kSampleRate / qwqdsp_swift_f0::kSampleRate; // 3
static constexpr float kDecimateAttenDb = 80.0f;
static constexpr int kDecimateTapsPerPhase = 32;

// ---- 滚动图像 ----
static constexpr float kScrollSeconds = 3.0f;
static constexpr int kImageWidth = static_cast<int>(kScrollSeconds * kSampleRate / kHopSize);

// ---- UI 布局 ----
static constexpr int kWindowWidth = 900;
static constexpr int kWindowHeight = 480;
static constexpr int kCanvasX = 60;
static constexpr int kCanvasY = 30;
static constexpr int kCanvasW = kWindowWidth - kCanvasX - 20;
static constexpr int kCanvasH = kWindowHeight - kCanvasY - 90;

// ---- UI 颜色 ----
static constexpr float kFreqTickLabels[] = {50, 200, 500, 2000, 20000};
static constexpr Color kGridColor = {50, 50, 50, 255};
static constexpr Color kTextColor = {180, 180, 180, 255};
static constexpr Color kBgColor = {20, 20, 20, 255};
// 背景 Resonator 色图为 黑→蓝→青→绿→黄→红，基频线用品红以避开整个色域
static constexpr Color kPitchColor = {255, 40, 255, 255};
static constexpr Color kPitchOutlineColor = {0, 0, 0, 255};

using ColorMap = ColormapAdapter<qwqdsp_colormap::Resonator>;

// ------------------------------------------------------------
//  全局状态
// ------------------------------------------------------------

static SpectrogramColumn g_column;
static NcReassignmentFrame<ColorMap> g_nc;
static ScrollingImage g_image;

static swift_f0_rt::Decimator g_decimator;
static swift_f0_rt::LogMagStft g_stft;
static swift_f0_rt::FrameRing g_frame_ring;
static swift_f0_rt::ResultRing g_result_ring;
static swift_f0_rt::SwiftF0Worker g_worker;

/// 背景瀑布图已推入的列总数（音频线程写，主线程读）
static std::atomic<uint64_t> g_total_columns{0};
/// 16kHz 侧已产出帧数（音频线程写）
static std::atomic<uint64_t> g_frame_count{0};

/// 显示用置信度阈值（主线程写，绘制读）
static std::atomic<float> g_conf_threshold{0.5f};

/// GUI 线程持有的近期基频结果（按时间递增）
struct PitchPoint {
    int64_t center_sample{};
    float pitch_hz{};
    float confidence{};
};
static std::deque<PitchPoint> g_pitch_hist;

/// 音频线程工作缓冲（避免回调内分配）
static std::vector<float> g_decimated;

// ------------------------------------------------------------
//  音频回调
// ------------------------------------------------------------

extern "C" void MaCaptureCallback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;
    (void)pOutput;
    float const* src = reinterpret_cast<float const*>(pInput);
    if (src == nullptr) {
        return;
    }

    // ---- 背景: NC 方法 ----
    auto push = [&](std::span<const Color> col) {
        g_image.PushColumn(col);
        g_total_columns.fetch_add(1, std::memory_order_relaxed);
    };
    g_column.ProcessAudio({src, frameCount}, g_nc, push);

    // ---- Swift F0 前端: 48k → 16k → STFT → log-mag 帧 ----
    g_decimator.Process({src, frameCount}, g_decimated);
    if (!g_decimated.empty()) {
        g_stft.Process(g_decimated, [&](std::span<const float> lm, int64_t center_16k) {
            swift_f0_rt::Frame f;
            std::copy(lm.begin(), lm.end(), f.lm.begin());
            // 抽取输出第 k 个样本对应输入样本 k*factor - 群延迟，换算回 48k 时间轴
            f.center_sample = center_16k * kDecimateFactor - g_decimator.DelaySamples();
            f.index = g_frame_count.fetch_add(1, std::memory_order_relaxed);
            g_frame_ring.Push(f);
        });
    }
}

// ------------------------------------------------------------
//  绘制
// ------------------------------------------------------------

/// 频率 → 画布 y（与背景 NC 帧相同的对数映射）
static int FreqToY(float hz) noexcept {
    float const log_min = std::log10(kFreqMin);
    float const log_max = std::log10(kFreqMax);
    float const norm = (std::log10(std::clamp(hz, kFreqMin, kFreqMax)) - log_min) / (log_max - log_min);
    return kCanvasY + kCanvasH - static_cast<int>(norm * static_cast<float>(kCanvasH));
}

static void DrawFreqGrid() {
    float const log_min = std::log10(kFreqMin);
    float const log_max = std::log10(kFreqMax);

    for (float freq : kFreqTickLabels) {
        float const norm = (std::log10(freq) - log_min) / (log_max - log_min);
        int const y = kCanvasY + kCanvasH - static_cast<int>(norm * static_cast<float>(kCanvasH));
        DrawLine(kCanvasX, y, kCanvasX + kCanvasW, y, kGridColor);

        char label[16];
        if (freq >= 1000.0f) {
            snprintf(label, sizeof(label), "%.0fk", freq / 1000.0f);
        }
        else {
            snprintf(label, sizeof(label), "%.0f", freq);
        }
        int const tw = MeasureText(label, 10);
        DrawText(label, kCanvasX - tw - 6, y - 5, 10, kTextColor);
    }
}

/**
 * @brief 叠加基频轨迹
 *
 * 背景第 c 列（从 0 计数）在输入累计 fft_size + c*hop 个样本时产出，
 * 其分析窗中心位于 c*hop + fft_size/2。基频帧已知分析窗中心的绝对样本位置，
 * 因此可直接换算成背景列坐标，与瀑布图严格对齐。
 */
static void DrawPitchTrack() {
    uint64_t const total_columns = g_total_columns.load(std::memory_order_relaxed);
    if (total_columns < static_cast<uint64_t>(kImageWidth)) {
        return;
    }
    float const threshold = g_conf_threshold.load(std::memory_order_relaxed);
    float const left_column = static_cast<float>(total_columns - static_cast<uint64_t>(kImageWidth));

    // 两遍绘制：先暗色粗描边，再亮色细主线，保证在任意背景亮度上可辨
    for (int pass = 0; pass < 2; ++pass) {
        Color const color = (pass == 0) ? kPitchOutlineColor : kPitchColor;
        float const width = (pass == 0) ? 3.0f : 1.5f;

        bool has_prev = false;
        Vector2 prev{};
        for (auto const& p : g_pitch_hist) {
            float const col = (static_cast<float>(p.center_sample) - 0.5f * static_cast<float>(kFftSize))
                            / static_cast<float>(kHopSize);
            float const xf = (col - left_column) * static_cast<float>(kCanvasW) / static_cast<float>(kImageWidth);
            if (xf < 0.0f || xf >= static_cast<float>(kCanvasW) || p.confidence < threshold) {
                has_prev = false;
                continue;
            }
            Vector2 const cur{kCanvasX + xf, static_cast<float>(FreqToY(p.pitch_hz))};
            if (has_prev) {
                DrawLineEx(prev, cur, width, color);
            }
            prev = cur;
            has_prev = true;
        }
    }
}

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------

int main(void) {
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, "Swift F0 RT Pitch - NC background | miniaudio + qwqdsp + raylib");
    SetTargetFPS(60);

    // ---- miniaudio 回环捕获 ----
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

    // ---- 初始化背景 ----
    g_column.Init(kCanvasH, kSampleRate, kFftSize, kHopSize);
    g_nc.Init(kSampleRate, kFftSize, kNcZeroPad, kCanvasH, kFreqMin, kFreqMax, kDbFloor);
    g_image.Init(kImageWidth, kCanvasH);

    // ---- 初始化 Swift F0 前端 ----
    g_decimated.reserve(4096); // 预留容量，避免音频回调内分配
    g_decimator.Init(kDecimateFactor, kDecimateAttenDb, kDecimateTapsPerPhase);
    g_stft.Init();
    g_worker.Start(&g_frame_ring, &g_result_ring);

    // ---- 主循环 ----
    while (!WindowShouldClose()) {
        // ---- 收取推理结果 ----
        {
            swift_f0_rt::PitchResult r;
            while (g_result_ring.Pop(r)) {
                g_pitch_hist.push_back(PitchPoint{r.center_sample, r.pitch_hz, r.confidence});
            }
            // 丢弃已滚出瀑布图左侧的点
            uint64_t const total_columns = g_total_columns.load(std::memory_order_relaxed);
            if (total_columns > static_cast<uint64_t>(kImageWidth)) {
                // 最左可见列的分析窗中心样本位置
                int64_t const left_sample = static_cast<int64_t>(total_columns - static_cast<uint64_t>(kImageWidth))
                                              * kHopSize
                                          + kFftSize / 2;
                while (!g_pitch_hist.empty() && g_pitch_hist.front().center_sample < left_sample) {
                    g_pitch_hist.pop_front();
                }
            }
        }

        BeginDrawing();
        ClearBackground(BLACK);

        // ---- 背景 + 基频叠加 ----
        DrawRectangle(kCanvasX, kCanvasY, kCanvasW, kCanvasH, kBgColor);
        g_image.Draw(kCanvasX, kCanvasY, kCanvasW, kCanvasH);
        DrawPitchTrack();
        DrawFreqGrid();

        // ---- 信息 ----
        {
            char buf[160];
            constexpr int kOverlapPct = 100 - kHopSize * 100 / kFftSize;
            snprintf(buf, sizeof(buf), "Background: NC (no time reassignment) | FFT %d Hop %d (%d%%) ZeroPad %d | %.0f dB",
                     kFftSize, kHopSize, kOverlapPct, kNcZeroPad, kDbFloor);
            DrawText(buf, kCanvasX + 4, 10, 10, kTextColor);
        }
        {
            char buf[192];
            snprintf(buf, sizeof(buf),
                     "SwiftF0: 48k->16k x%d | n_fft %d hop %d | block %d | infer %.1f ms | queue %d | frames %llu",
                     kDecimateFactor, qwqdsp_swift_f0::kNFFT, qwqdsp_swift_f0::kHopLength, swift_f0_rt::kBlock,
                     g_worker.LastInferenceMs(), g_frame_ring.Depth(),
                     static_cast<unsigned long long>(g_worker.OutputFrame()));
            DrawText(buf, kCanvasX + 4, kCanvasY + kCanvasH + 8, 10, kTextColor);
        }
        {
            float const threshold = g_conf_threshold.load(std::memory_order_relaxed);
            char buf[96];
            snprintf(buf, sizeof(buf), "Confidence threshold: %.2f   (A / D to adjust)", threshold);
            DrawText(buf, kCanvasX + 4, kCanvasY + kCanvasH + 24, 10, kTextColor);
        }
        {
            // 最新一帧的基频读数
            float hz = 0.0f, conf = 0.0f;
            bool valid = false;
            if (!g_pitch_hist.empty()) {
                auto const& p = g_pitch_hist.back();
                hz = p.pitch_hz;
                conf = p.confidence;
                valid = true;
            }
            char buf[96];
            if (valid) {
                snprintf(buf, sizeof(buf), "Latest: %.1f Hz   conf %.2f", hz, conf);
            }
            else {
                snprintf(buf, sizeof(buf), "Latest: --");
            }
            DrawText(buf, kCanvasX + 4, kCanvasY + kCanvasH + 40, 12,
                     (valid && conf >= g_conf_threshold.load(std::memory_order_relaxed)) ? kPitchColor : kTextColor);
        }
        DrawFPS(kWindowWidth - 80, 10);
        EndDrawing();

        // ---- 阈值调节 ----
        if (IsKeyDown(KEY_A)) {
            float v = g_conf_threshold.load(std::memory_order_relaxed) - 0.01f;
            g_conf_threshold.store(std::clamp(v, 0.0f, 1.0f), std::memory_order_relaxed);
        }
        if (IsKeyDown(KEY_D)) {
            float v = g_conf_threshold.load(std::memory_order_relaxed) + 0.01f;
            g_conf_threshold.store(std::clamp(v, 0.0f, 1.0f), std::memory_order_relaxed);
        }
    }

    // ---- 清理 ----
    if (result == MA_SUCCESS) {
        ma_device_stop(&device);
        ma_device_uninit(&device);
    }
    g_worker.Stop();
    g_image.Unload();
    CloseWindow();
    return 0;
}
