#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <format>
#include <mutex>
#include <numbers>
#include <span>
#include <string>
#include <vector>

#include <qwqdsp/colormap/colormap.hpp>

#include "miniaudio.h"
#include "raylib.h"

#include "gui_test.hpp"
#include "reassignment/colormap_adapter.hpp"
#include "reassignment/spectrogram_column.hpp"
#include "reassignment/windowless_nc_frame.hpp"

// ── 分析器参数 ──
static constexpr int kSampleRate = 48000;
static constexpr int kFftSize = 4096;
static constexpr int kHopSize = kFftSize / 16; // 256 样本 ≈ 5.33 ms 一列
static constexpr int kOctaveDivisions = 12;    // 每倍频程 12 个滤波器
static constexpr int kBandCount = 120;         // 10 个倍频程：20.87 Hz – 20.16 kHz
static constexpr float kBandRefHz = 1000.0f;   // 1/12 倍频程网格的参考中心频率
static constexpr int kBandRefIndex = 67;       // 1000 Hz 在网格中的下标
static constexpr float kHalfBand = 1.0f / (2.0f * kOctaveDivisions); // 半带宽（倍频程）
static constexpr float kDbFloor = -84.0f;      // 帧内部底噪（无信号行）

// ── 显示参数 ──
static constexpr float kDbTop = 0.0f;
static constexpr float kDbBottom = -72.0f;
static constexpr float kDbGridStep = 6.0f;
static constexpr float kPeakDecayDbPerSec = 18.0f;

// ── 平滑参数（按列更新，列率 = kSampleRate / kHopSize）──
static constexpr float kAttackMs = 10.0f;
static constexpr float kReleaseMs = 220.0f;

// ── 窗口与布局 ──
static constexpr int kWindowWidth = 980;
static constexpr int kWindowHeight = 480;
static constexpr int kCanvasX = 58;
static constexpr int kCanvasY = 56;
static constexpr int kCanvasW = 858;
static constexpr int kCanvasH = 380;

// ── 颜色 ──
static const Color kBgColor = {8, 8, 10, 255};
static const Color kPanelColor = {14, 14, 18, 255};
static const Color kPanelBorder = {58, 58, 70, 255};
static const Color kGridColor = {38, 38, 48, 255};
static const Color kGridStrongColor = {56, 56, 70, 255};
static const Color kTextColor = {196, 196, 208, 255};
static const Color kPeakColor = {255, 236, 180, 255};

// ── 频率刻度 ──
static constexpr float kFreqTicks[] = {20.0f, 30.0f, 40.0f, 60.0f, 100.0f, 200.0f,  300.0f, 500.0f,
                                       1000.0f, 2000.0f, 3000.0f, 4000.0f, 6000.0f, 10000.0f, 20000.0f};

// ── 截图验证开关（验证时改 true；交付前为 false）──
static constexpr bool kGuiTest = false;
static constexpr char const* kShotPath = "temp/nc_octave_bar/shots/bar.png";

// ----------------------------------------
// 1/12 倍频程网格
// ----------------------------------------

/**
 * @brief 生成 1/12 倍频程中心频率表（升序，1000 Hz 为参考中心）
 *
 * 下标 67 对应 1000 Hz，下标 0/119 对应 20.87 Hz / 20.16 kHz。
 */
static std::array<float, kBandCount> MakeBandCenters() noexcept {
    std::array<float, kBandCount> centers{};
    for (int i = 0; i < kBandCount; ++i) {
        centers[static_cast<size_t>(i)] =
            kBandRefHz * std::exp2(static_cast<float>(i - kBandRefIndex) / static_cast<float>(kOctaveDivisions));
    }
    return centers;
}

static const std::array<float, kBandCount> g_band_centers = MakeBandCenters();
static const float g_freq_min = g_band_centers.front() * std::exp2(-kHalfBand);
static const float g_freq_max = g_band_centers.back() * std::exp2(kHalfBand);

// ----------------------------------------
// 全局状态
// ----------------------------------------

using ColorMap = ColormapAdapter<qwqdsp_colormap::Magma>;

static WindowlessNcFrame<ColorMap> g_analyzer;
static SpectrogramColumn g_column;

static std::mutex g_level_mutex;                            // 保护下面两行（音频线程写，GUI 读）
static std::array<float, kBandCount> g_latest_db{};         // 已做攻击/释放平滑的逐带电平
static std::array<float, kBandCount> g_smoothed_db{};       // 音频线程的平滑状态

static std::array<float, kBandCount> g_peak_db{};           // GUI 线程的峰值保持
static float g_attack_coef = 0.0f;
static float g_release_coef = 0.0f;

static ma_device g_device;
static bool g_device_ready = false;

// ----------------------------------------
// 坐标映射
// ----------------------------------------

/// @brief 频率(Hz) → 面板 x（对数轴）
static float FreqToX(float freq) noexcept {
    const float log_min = std::log10(g_freq_min);
    const float log_max = std::log10(g_freq_max);
    const float t = (std::log10(freq) - log_min) / (log_max - log_min);
    return static_cast<float>(kCanvasX) + t * static_cast<float>(kCanvasW);
}

/// @brief 电平(dB) → 面板 y（0 dB 在顶，kDbBottom 在底）
static int DbToY(float db) noexcept {
    const float t = (std::clamp(db, kDbBottom, kDbTop) - kDbBottom) / (kDbTop - kDbBottom);
    return kCanvasY + kCanvasH - static_cast<int>(std::lround(t * static_cast<float>(kCanvasH)));
}

// ----------------------------------------
// 分析链
// ----------------------------------------

/**
 * @brief 每产出一列：把帧的 dB 数据做攻击/释放平滑后发布给 GUI
 *
 * 颜色列不使用；柱状显示直接读 GetColumnDb()，避免从颜色反查电平。
 */
static void PushColumn(std::span<const Color> /*colors*/) noexcept {
    const std::span<const float> db = g_analyzer.GetColumnDb();
    std::lock_guard<std::mutex> lock(g_level_mutex);
    for (int i = 0; i < kBandCount; ++i) {
        const float target = db[static_cast<size_t>(i)];
        const float prev = g_smoothed_db[static_cast<size_t>(i)];
        const float alpha = (target > prev) ? g_attack_coef : g_release_coef;
        g_smoothed_db[static_cast<size_t>(i)] = prev + alpha * (target - prev);
    }
    g_latest_db = g_smoothed_db;
}

static void Dsp_Init() noexcept {
    g_column.Init(kBandCount, kSampleRate, kFftSize, kHopSize);
    g_analyzer.InitBands(kSampleRate, kFftSize, kHopSize, g_band_centers, kDbFloor, 1.0f);

    const float column_dt = static_cast<float>(kHopSize) / static_cast<float>(kSampleRate);
    g_attack_coef = 1.0f - std::exp(-column_dt / (kAttackMs * 0.001f));
    g_release_coef = 1.0f - std::exp(-column_dt / (kReleaseMs * 0.001f));

    g_smoothed_db.fill(kDbFloor);
    g_latest_db.fill(kDbFloor);
    g_peak_db.fill(kDbBottom);
}

extern "C" void MaCaptureCallback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;
    (void)pOutput;
    float const* src = reinterpret_cast<float const*>(pInput);
    g_column.ProcessAudio({src, frameCount}, g_analyzer, PushColumn);
}

/**
 * @brief 生成确定性测试信号并跑完整分析链（截图验证用）
 *
 * 六个 1/12 倍频程中心频率上的纯音给可断言的柱高，白噪声给自然的频段填充
 * （频段宽度随频率上升 → 每倍频程约 +3 dB）。不初始化音频设备，结果可复现。
 */
static void FeedTestSignal() noexcept {
    struct Tone {
        float freq;
        float amp;
    };
    constexpr Tone kTones[] = {
        {31.25f, 0.22f}, {125.0f, 0.18f}, {500.0f, 0.12f},
        {2000.0f, 0.09f}, {8000.0f, 0.07f}, {16000.0f, 0.10f},
    };
    constexpr float kNoiseRms = 0.02f;
    constexpr int kSeconds = 2;
    constexpr int kBlock = 512;

    std::vector<float> block(kBlock);
    std::uint32_t rng = 0x9E3779B9u;
    const int total = kSeconds * kSampleRate;
    for (int start = 0; start < total; start += kBlock) {
        const int count = std::min(kBlock, total - start);
        for (int i = 0; i < count; ++i) {
            const double t = static_cast<double>(start + i) / static_cast<double>(kSampleRate);
            double s = 0.0;
            for (const Tone& tone : kTones) {
                s += static_cast<double>(tone.amp)
                   * std::sin(2.0 * std::numbers::pi_v<double> * static_cast<double>(tone.freq) * t
                              + static_cast<double>(tone.freq) * 0.001);
            }
            rng = rng * 1664525u + 1013904223u;
            const float u = static_cast<float>(rng >> 8) / static_cast<float>(1u << 24) * 2.0f - 1.0f;
            s += static_cast<double>(u) * static_cast<double>(kNoiseRms) * std::numbers::sqrt3_v<float>;
            block[static_cast<size_t>(i)] = static_cast<float>(s);
        }
        g_column.ProcessAudio(std::span<const float>(block.data(), static_cast<size_t>(count)), g_analyzer,
                              PushColumn);
    }
}

// ----------------------------------------
// 绘制
// ----------------------------------------

static Color MixColor(Color a, Color b, float t) noexcept {
    t = std::clamp(t, 0.0f, 1.0f);
    auto mix = [t](unsigned char x, unsigned char y) {
        return static_cast<unsigned char>(std::lround(static_cast<float>(x) + (static_cast<float>(y) - x) * t));
    };
    return Color{mix(a.r, b.r), mix(a.g, b.g), mix(a.b, b.b), 255};
}

/// @brief 柱色：低电平深蓝紫 → 品红 → 高电平亮红
static Color BarColor(float db) noexcept {
    const float t = std::clamp((db - kDbBottom) / (kDbTop - kDbBottom), 0.0f, 1.0f);
    if (t < 0.5f)
        return MixColor(Color{56, 44, 138, 255}, Color{196, 62, 198, 255}, t * 2.0f);
    return MixColor(Color{196, 62, 198, 255}, Color{255, 72, 76, 255}, (t - 0.5f) * 2.0f);
}

/// @brief 面板底、网格与刻度
static void DrawPanel() noexcept {
    DrawRectangle(kCanvasX, kCanvasY, kCanvasW, kCanvasH, kPanelColor);

    // 水平 dB 网格
    for (float db = kDbTop; db >= kDbBottom - 0.5f; db -= kDbGridStep) {
        const int y = DbToY(db);
        const Color color = (db == kDbTop) ? kGridStrongColor : kGridColor;
        DrawLine(kCanvasX, y, kCanvasX + kCanvasW, y, color);
    }

    // 垂直频率网格
    const float log_min = std::log10(g_freq_min);
    const float log_max = std::log10(g_freq_max);
    for (float freq : kFreqTicks) {
        if (freq < g_freq_min || freq > g_freq_max)
            continue;
        const float norm = (std::log10(freq) - log_min) / (log_max - log_min);
        const int x = kCanvasX + static_cast<int>(std::lround(norm * static_cast<float>(kCanvasW)));
        DrawLine(x, kCanvasY, x, kCanvasY + kCanvasH, kGridColor);
    }

    DrawRectangleLinesEx(Rectangle{static_cast<float>(kCanvasX), static_cast<float>(kCanvasY),
                                   static_cast<float>(kCanvasW), static_cast<float>(kCanvasH)},
                         1.0f, kPanelBorder);
}

/// @brief 左右两侧的 dB 刻度标签
static void DrawDbLabels() noexcept {
    for (float db = kDbTop; db >= kDbBottom - 0.5f; db -= kDbGridStep) {
        const int y = DbToY(db);
        const std::string label = std::format("{:.0f}", db);
        const int width = MeasureText(label.c_str(), 10);
        DrawText(label.c_str(), kCanvasX - width - 8, y - 5, 10, kTextColor);
        DrawText(label.c_str(), kCanvasX + kCanvasW + 8, y - 5, 10, kTextColor);
    }
}

/// @brief 频率刻度标签
static void DrawFreqLabels() noexcept {
    const int y = kCanvasY + kCanvasH + 6;
    for (float freq : kFreqTicks) {
        if (freq < g_freq_min || freq > g_freq_max)
            continue;
        const float norm = (std::log10(freq) - std::log10(g_freq_min)) / (std::log10(g_freq_max) - std::log10(g_freq_min));
        const int x = kCanvasX + static_cast<int>(std::lround(norm * static_cast<float>(kCanvasW)));
        const std::string label = (freq >= 1000.0f) ? std::format("{:.0f}k", freq / 1000.0f) : std::format("{:.0f}", freq);
        const int width = MeasureText(label.c_str(), 10);
        DrawText(label.c_str(), x - width / 2, y, 10, kTextColor);
    }
}

/// @brief 柱状频谱 + 峰值保持
static void DrawBars(const std::array<float, kBandCount>& db) noexcept {
    for (int i = 0; i < kBandCount; ++i) {
        const float fc = g_band_centers[static_cast<size_t>(i)];
        const int x0 = static_cast<int>(std::lround(FreqToX(fc * std::exp2(-kHalfBand))));
        const int x1 = static_cast<int>(std::lround(FreqToX(fc * std::exp2(kHalfBand))));
        const int width = std::max(1, x1 - x0 - 1);

        const float value = std::clamp(db[static_cast<size_t>(i)], kDbBottom, kDbTop);
        const int y = DbToY(value);
        const int height = std::clamp(kCanvasY + kCanvasH - y, 1, kCanvasH);

        const Color color = BarColor(value);
        DrawRectangle(x0, y, width, height, color);
        DrawRectangle(x0, y, width, 1, MixColor(color, Color{255, 255, 255, 255}, 0.35f));

        const float peak = std::clamp(g_peak_db[static_cast<size_t>(i)], kDbBottom, kDbTop);
        const int peak_y = DbToY(peak);
        if (peak_y < y - 1)
            DrawRectangle(x0, peak_y - 1, width, 2, kPeakColor);
    }
}

/// @brief 标题
static void DrawTitle() noexcept {
    DrawText("NC 1/12-Octave Spectrum", 16, 10, 22, kTextColor);
}

static void DrawAnalyzer() noexcept {
    std::array<float, kBandCount> db{};
    {
        std::lock_guard<std::mutex> lock(g_level_mutex);
        db = g_latest_db;
    }

    const float dt = std::min(GetFrameTime(), 0.1f);
    for (int i = 0; i < kBandCount; ++i) {
        float& peak = g_peak_db[static_cast<size_t>(i)];
        peak = std::max(peak - kPeakDecayDbPerSec * dt, db[static_cast<size_t>(i)]);
    }

    DrawPanel();
    DrawBars(db);
    DrawDbLabels();
    DrawFreqLabels();
    DrawTitle();
}

// ----------------------------------------
// 音频设备
// ----------------------------------------

static void StartCapture() noexcept {
    ma_device_config config = ma_device_config_init(ma_device_type_loopback);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.sampleRate = static_cast<ma_uint32>(kSampleRate);
    config.dataCallback = MaCaptureCallback;
    config.pUserData = nullptr;
    config.periodSizeInMilliseconds = 10;

    if (ma_device_init(nullptr, &config, &g_device) == MA_SUCCESS) {
        ma_device_start(&g_device);
        g_device_ready = true;
    }
    else {
        TraceLog(LOG_WARNING, "miniaudio capture device init failed, running silent");
    }
}

static void StopCapture() noexcept {
    if (g_device_ready) {
        ma_device_stop(&g_device);
        ma_device_uninit(&g_device);
        g_device_ready = false;
    }
}

// ----------------------------------------
// main
// ----------------------------------------

int main(void) {
    qwqdsp_gui_test::ShotMode<kGuiTest> shot(kShotPath);
    shot.ApplyWindowFlags();

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, "NC 1/12-Octave Spectrum");
    SetTargetFPS(60);

    Dsp_Init();

    if (shot.IsEnabled()) {
        FeedTestSignal(); // 跳过设备，喂确定性信号
    }
    else {
        StartCapture();
    }

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(kBgColor);
        DrawAnalyzer();
        EndDrawing();

        if (shot.Tick())
            break;
    }

    StopCapture();
    CloseWindow();
    return 0;
}
