#include <array>
#include <vector>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <numbers>

#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/window/blackman_harris.hpp>
#include <qwqdsp/window/helper.hpp>

#include "miniaudio.h"
#include "raylib.h"

// ════════════════════════════════════════════════════════════
//  常量
// ════════════════════════════════════════════════════════════

// ── 窗口 & 画布 ──
static constexpr int kWindowWidth = 1280;
static constexpr int kWindowHeight = 720;
static constexpr int kCanvasX = 50;
static constexpr int kCanvasY = 30;
static constexpr int kCanvasW = kWindowWidth - kCanvasX - 20;
static constexpr int kCanvasH = kWindowHeight - kCanvasY - 60;

// ── DSP ──
static constexpr int kSampleRate = 48000;
static constexpr int kHop = 64;

// ── 多分辨率 ──
// 分频点 (Hz, 升序排列，介于 kFreqMin ~ kFreqMax 之间)
static constexpr std::array<float, 3> kCrossovers = {1500.0f, 3000.0f, 6000.0f};
static constexpr int kNumResolutions = kCrossovers.size() + 1;
// 自动生成 FFT 尺寸 — 索引 0 最大 FFT（最低频段）
static constexpr std::array<int, kNumResolutions> kFftSizes = [] {
    std::array<int, kNumResolutions> arr{};
    int v = 512;
    for (int i = kNumResolutions - 1; i >= 0; --i) {
        arr[i] = v;
        v *= 2;
    }
    return arr;
}();
// ── 显示 ──
static constexpr float kFreqMin = 20.0f;
static constexpr float kFreqMax = 20000.0f;
static constexpr float kDbFloor = -80.0f;

// ── 平滑 ──
static constexpr float kAttackMs = 1.0f;
static constexpr float kReleaseMs = 150.0f;

// ── 颜色 ──
static const Color kGridColor = {50, 50, 50, 255};
static const Color kTextColor = {180, 180, 180, 255};
static const Color kCurveColor = {255, 255, 255, 255};
static const Color kBgColor = {20, 20, 20, 255};
// 参考曲线颜色（每个 FFT 一条）
static const Color kRefColors[] = {
    {120, 200, 255, 160}, // 4096 — 青半透明
    {255, 180,  80, 160}, // 2048 — 橙半透明
    { 80, 220, 120, 160}, // 1024 — 绿半透明
    {240, 120, 200, 160}, // 512  — 粉半透明
};

// ════════════════════════════════════════════════════════════
//  显示查找表
// ════════════════════════════════════════════════════════════

struct DisplayLutEntry {
    int channel_idx;
    float bin_idx;
};
static std::vector<DisplayLutEntry> g_display_lut;

// ════════════════════════════════════════════════════════════
//  FFT 通道
// ════════════════════════════════════════════════════════════

struct FftChannel {
    std::vector<float> buf;          // 环形输入 buffer
    std::vector<float> window;       // 窗系数（保持不变）
    std::vector<float> win_buf;      // 加窗后工作区
    std::vector<float> gain;         // FFT 幅度输出
    qwqdsp_spectral::RealFFT fft;
    std::vector<float> mag_dB;      // 平滑后的幅度 dB
    std::vector<float> smoothed_dB; // 平滑状态（每个 bin 独立）
};
static std::array<FftChannel, kNumResolutions> g_channels;

// ── 绘制快照（音频线程发布，GUI 线程只读） ──
static std::array<std::vector<float>, 2> g_curve_dB;
static std::atomic<int> g_curve_read_idx{0};

// ── 平滑系数 ──
static float g_attack_alpha = 0.0f;
static float g_release_alpha = 0.0f;

static std::array<float, kHop> g_hop_buf{};
static int g_hop_count = 0;

static float SampleMagDb(const std::vector<float>& mag_dB, float bin_idx) noexcept {
    int num_bins = static_cast<int>(mag_dB.size());
    int i0 = std::clamp(static_cast<int>(bin_idx), 0, num_bins - 1);
    int i1 = std::min(i0 + 1, num_bins - 1);
    float frac = bin_idx - static_cast<float>(i0);
    return mag_dB[i0] + frac * (mag_dB[i1] - mag_dB[i0]);
}

static void PublishDisplayCurve() {
    int write_idx = 1 - g_curve_read_idx.load(std::memory_order_relaxed);
    auto& curve = g_curve_dB[write_idx];
    for (int px = 0; px < kCanvasW; ++px) {
        const auto& lut = g_display_lut[px];
        curve[px] = SampleMagDb(g_channels[lut.channel_idx].mag_dB, lut.bin_idx);
    }
    g_curve_read_idx.store(write_idx, std::memory_order_release);
}

// 参考曲线每 FFT 每像素的 bin 索引查找表
static std::array<std::vector<float>, kNumResolutions> g_ref_bin_lut;

// ════════════════════════════════════════════════════════════
//  ReferenceView — 独立参考分析器
// ════════════════════════════════════════════════════════════

class ReferenceView {
public:
    void Init() {
        for (int i = 0; i < kNumResolutions; ++i) {
            int n = kFftSizes[i];
            auto& ch = channels_[i];
            ch.buf.assign(n, 0.0f);
            ch.window.resize(n);
            qwqdsp_window::BlackmanHarris::Window(ch.window, true);
            qwqdsp_window::Helper::Normalize(ch.window);
            ch.win_buf.resize(n);
            ch.gain.assign(n / 2 + 1, 0.0f);
            ch.fft.Init(n);
            ch.smoothed_dB.assign(n / 2 + 1, kDbFloor);
        }
    }

    void Process(float const* src) {
        for (int i = 0; i < kNumResolutions; ++i) {
            auto& ch = channels_[i];
            int n = static_cast<int>(ch.buf.size());

            std::memmove(ch.buf.data(), ch.buf.data() + kHop, (n - kHop) * sizeof(float));
            std::memcpy(ch.buf.data() + n - kHop, src, kHop * sizeof(float));

            for (int j = 0; j < n; ++j)
                ch.win_buf[j] = ch.buf[j] * ch.window[j];

            ch.fft.FFTGainPhase(ch.win_buf, ch.gain);
            for (int j = 0; j < n / 2 + 1; ++j) {
                float raw_dB = 20.0f * std::log10(ch.gain[j] + 1e-12f);
                float prev = ch.smoothed_dB[j];
                float alpha = (raw_dB > prev) ? g_attack_alpha : g_release_alpha;
                ch.smoothed_dB[j] = prev + alpha * (raw_dB - prev);
            }
        }
    }

    void Draw() const;

private:
    struct Channel {
        std::vector<float> buf;
        std::vector<float> window;
        std::vector<float> win_buf;
        std::vector<float> gain;
        qwqdsp_spectral::RealFFT fft;
        std::vector<float> smoothed_dB;
    };
    std::array<Channel, kNumResolutions> channels_;
};

static ReferenceView g_ref_view;

void ReferenceView::Draw() const {
    for (int i = 0; i < kNumResolutions; ++i) {
        const auto& ch = channels_[i];
        int num_bins = static_cast<int>(ch.smoothed_dB.size());

        std::vector<Vector2> pts(kCanvasW);
        for (int px = 0; px < kCanvasW; ++px) {
            float bin_idx = g_ref_bin_lut[i][px];
            int i0 = std::clamp(static_cast<int>(bin_idx), 0, num_bins - 1);
            int i1 = std::min(i0 + 1, num_bins - 1);
            float frac = bin_idx - static_cast<float>(i0);
            float db_val = ch.smoothed_dB[i0] + frac * (ch.smoothed_dB[i1] - ch.smoothed_dB[i0]);
            float normY = std::clamp((db_val - kDbFloor) / (-kDbFloor), 0.0f, 1.0f);
            pts[px] = {static_cast<float>(kCanvasX + px),
                       kCanvasY + kCanvasH - normY * kCanvasH};
        }
        DrawLineStrip(pts.data(), kCanvasW, kRefColors[i]);
    }
}

// ════════════════════════════════════════════════════════════
//  初始化
// ════════════════════════════════════════════════════════════

static void Dsp_Init() {
    // ── 平滑系数 ──
    float T = static_cast<float>(kHop) / static_cast<float>(kSampleRate);
    g_attack_alpha = 1.0f - std::exp(-T / (kAttackMs * 0.001f));
    g_release_alpha = 1.0f - std::exp(-T / (kReleaseMs * 0.001f));

    // ── FFT 通道 ──
    for (int i = 0; i < kNumResolutions; ++i) {
        int n = kFftSizes[i];
        auto& ch = g_channels[i];
        ch.buf.assign(n, 0.0f);
        ch.window.resize(n);
        qwqdsp_window::BlackmanHarris::Window(ch.window, true);
        qwqdsp_window::Helper::Normalize(ch.window);
        ch.win_buf.resize(n);
        ch.gain.assign(n / 2 + 1, 0.0f);
        ch.fft.Init(n);
        ch.mag_dB.assign(n / 2 + 1, kDbFloor);
        ch.smoothed_dB.assign(n / 2 + 1, kDbFloor);
    }

    for (auto& curve : g_curve_dB) {
        curve.assign(kCanvasW, kDbFloor);
    }

    // ── 显示查找表 ──
    const float kLogFreqMin = std::log10(kFreqMin);
    const float kLogFreqMax = std::log10(kFreqMax);
    const float kLogFreqSpan = kLogFreqMax - kLogFreqMin;

    g_display_lut.resize(kCanvasW);
    for (int px = 0; px < kCanvasW; ++px) {
        float t = static_cast<float>(px) / static_cast<float>(kCanvasW - 1);
        float logFreq = kLogFreqMin + t * kLogFreqSpan;
        float freq = std::pow(10.0f, logFreq);

        // 根据分频点确定所在频段
        int seg = 0;
        for (int i = 0; i < kNumResolutions - 1; ++i) {
            if (freq >= kCrossovers[i])
                seg = i + 1;
        }

        int fftSize = kFftSizes[seg];
        int numBins = fftSize / 2 + 1;
        float binIdx = freq / static_cast<float>(kSampleRate) * static_cast<float>(fftSize);
        binIdx = std::clamp(binIdx, 0.0f, static_cast<float>(numBins - 1));

        g_display_lut[px] = {seg, binIdx};
    }

    // ── 参考曲线查找表（每个 FFT 全频段） ──
    for (int i = 0; i < kNumResolutions; ++i) {
        int fftSize = kFftSizes[i];
        int numBins = fftSize / 2 + 1;
        g_ref_bin_lut[i].resize(kCanvasW);
        for (int px = 0; px < kCanvasW; ++px) {
            float t = static_cast<float>(px) / static_cast<float>(kCanvasW - 1);
            float logFreq = kLogFreqMin + t * kLogFreqSpan;
            float freq = std::pow(10.0f, logFreq);
            float binIdx = freq / static_cast<float>(kSampleRate) * static_cast<float>(fftSize);
            g_ref_bin_lut[i][px] = std::clamp(binIdx, 0.0f, static_cast<float>(numBins - 1));
        }
    }

    // ── 参考分析器 ──
    g_ref_view.Init();

    PublishDisplayCurve();
}

// ════════════════════════════════════════════════════════════
//  DSP 处理
// ════════════════════════════════════════════════════════════

static void Dsp_Process(float const* src) {
    for (int i = 0; i < kNumResolutions; ++i) {
        auto& ch = g_channels[i];
        int n = static_cast<int>(ch.buf.size());

        // 左移 kHop，尾部填入新采样
        std::memmove(ch.buf.data(), ch.buf.data() + kHop, (n - kHop) * sizeof(float));
        std::memcpy(ch.buf.data() + n - kHop, src, kHop * sizeof(float));

        // 加窗
        for (int j = 0; j < n; ++j) {
            ch.win_buf[j] = ch.buf[j] * ch.window[j];
        }

        // FFT → 幅度 dB → 指数平滑
        ch.fft.FFTGainPhase(ch.win_buf, ch.gain);
        for (int j = 0; j < n / 2 + 1; ++j) {
            float raw_dB = 20.0f * std::log10(ch.gain[j] + 1e-12f);
            float prev = ch.smoothed_dB[j];
            float alpha = (raw_dB > prev) ? g_attack_alpha : g_release_alpha;
            ch.smoothed_dB[j] = prev + alpha * (raw_dB - prev);
            ch.mag_dB[j] = ch.smoothed_dB[j];
        }
    }
    g_ref_view.Process(src);
    PublishDisplayCurve();
}

// ════════════════════════════════════════════════════════════
//  音频回调
// ════════════════════════════════════════════════════════════

extern "C" void MaCaptureCallback(ma_device* pDevice, void* pOutput,
                                  const void* pInput, ma_uint32 frameCount) {
    float const* src = reinterpret_cast<float const*>(pInput);
    float* dst = reinterpret_cast<float*>(pOutput);
    int remaining = static_cast<int>(frameCount);
    if (src == nullptr) {
        if (dst != nullptr) {
            std::fill_n(dst, remaining, 0.0f);
        }
        return;
    }
    if (dst != nullptr) {
        std::copy_n(src, remaining, dst);
    }

    while (remaining > 0) {
        int need = kHop - g_hop_count;
        int can_do = std::min(need, remaining);
        std::copy_n(src, can_do, g_hop_buf.begin() + g_hop_count);
        g_hop_count += can_do;
        remaining -= can_do;
        src += can_do;

        if (g_hop_count == kHop) {
            Dsp_Process(g_hop_buf.data());
            g_hop_count = 0;
        }
    }
}

// ════════════════════════════════════════════════════════════
//  绘制
// ════════════════════════════════════════════════════════════

static void DrawBackground() {
    DrawRectangle(kCanvasX, kCanvasY, kCanvasW, kCanvasH, kBgColor);
}

static void DrawYAxisGrid() {
    static constexpr float kYTicks[] = {0.0f, -20.0f, -40.0f, -60.0f, -80.0f};
    for (float db : kYTicks) {
        float norm = (db - kDbFloor) / (-kDbFloor);
        int y = kCanvasY + kCanvasH - static_cast<int>(norm * kCanvasH);
        DrawLine(kCanvasX, y, kCanvasX + kCanvasW, y, kGridColor);
        char label[16];
        snprintf(label, sizeof(label), "%.0f", db);
        int tw = MeasureText(label, 10);
        DrawText(label, kCanvasX - tw - 6, y - 5, 10, kTextColor);
    }
}

static void DrawXAxisGrid() {
    const float kLogFreqMin = std::log10(kFreqMin);
    const float kLogFreqMax = std::log10(kFreqMax);
    const float kLogFreqSpan = kLogFreqMax - kLogFreqMin;

    // decade 子网格
    static constexpr float kDecadeStart[] = {20, 200, 2000};
    for (float start : kDecadeStart) {
        float step = start;
        for (float freq = start; freq <= start * 10.0f + 0.5f; freq += step) {
            float norm = (std::log10(freq) - kLogFreqMin) / kLogFreqSpan;
            float sx = kCanvasX + norm * static_cast<float>(kCanvasW);
            DrawLine(static_cast<int>(sx), kCanvasY,
                     static_cast<int>(sx), kCanvasY + kCanvasH, kGridColor);
        }
    }

    // 标签
    static constexpr float kXTickLabels[] = {20, 200, 2000, 20000};
    for (float freq : kXTickLabels) {
        float norm = (std::log10(freq) - kLogFreqMin) / kLogFreqSpan;
        float sx = kCanvasX + norm * static_cast<float>(kCanvasW);
        char label[16];
        if (freq >= 1000.0f)
            snprintf(label, sizeof(label), "%.0fk", freq / 1000.0f);
        else
            snprintf(label, sizeof(label), "%.0f", freq);
        int tw = MeasureText(label, 10);
        DrawText(label, static_cast<int>(sx) - tw / 2,
                 kCanvasY + kCanvasH + 4, 10, kTextColor);
    }
}

static void DrawCurve() {
    std::vector<Vector2> pts(kCanvasW);
    int curve_idx = g_curve_read_idx.load(std::memory_order_acquire);
    const auto& curve = g_curve_dB[curve_idx];
    for (int px = 0; px < kCanvasW; ++px) {
        float dbVal = curve[px];
        float normY = (dbVal - kDbFloor) / (-kDbFloor);
        normY = std::clamp(normY, 0.0f, 1.0f);
        float y = kCanvasY + kCanvasH - normY * kCanvasH;
        pts[px] = {static_cast<float>(kCanvasX + px), y};
    }
    if (kCanvasW > 1)
        DrawLineStrip(pts.data(), kCanvasW, kCurveColor);
}

// ════════════════════════════════════════════════════════════
//  main
// ════════════════════════════════════════════════════════════

int main(void) {
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight,
               "Multi-Resolution Spectrum - miniaudio + qwqdsp + raylib");
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

    // ── DSP 初始化 ──
    Dsp_Init();

    ma_device device;
    ma_result result = ma_device_init(nullptr, &config, &device);
    if (result == MA_SUCCESS) {
        ma_device_start(&device);
    } else {
        TraceLog(LOG_WARNING, "miniaudio 捕获设备初始化失败，以静默模式运行");
    }

    // ── 主循环 ──
    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(BLACK);

        DrawBackground();
        DrawYAxisGrid();
        // g_ref_view.Draw();
        DrawCurve();
        DrawXAxisGrid();

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
