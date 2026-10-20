/// @brief 用时间折叠把长窗频谱搬进固定大小的 FFT, 并按对数频率绘制 min/max 频谱带

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstring>
#include <format>
#include <numbers>
#include <span>

#include <qwqdsp/spectral/complex_fft_adv.hpp>
#include <qwqdsp/window/blackman.hpp>
#include <qwqdsp/window/helper.hpp>

#include "gui_test.hpp"

#include "miniaudio.h"
#include "raylib.h"

/// 截图验证时改为 true
static constexpr bool kGuiTest = false;
/// 截图输出路径（相对工作区根）
static constexpr char const* kShotPath = "temp/portnoff_spectrum/shots/helper.png";

// ------------------------------------------------------------
// 分析参数
// ------------------------------------------------------------

static constexpr int kSampleRate = 48000;

/// FFT 大小固定; 折叠因子决定实际分析窗长与恢复的谱线密度
static constexpr int kFftSize = 1024;
static constexpr int kFoldFactor = 8;
static constexpr int kWindowSize = kFftSize * kFoldFactor;
/// 等效 FFT 长度由 FFT 大小与折叠因子决定
static constexpr int kSpecBins = kWindowSize / 2 + 1;

static constexpr int kHopSize = kFftSize / 2;
static constexpr float kBinHz = static_cast<float>(kSampleRate) / static_cast<float>(kWindowSize);

// ------------------------------------------------------------
// 显示参数
// ------------------------------------------------------------

static constexpr int kWidth = 960;
static constexpr int kHeight = 540;
static constexpr int kCanvasX = 64;
static constexpr int kCanvasY = 40;
static constexpr int kCanvasW = kWidth - kCanvasX - 24;
static constexpr int kCanvasH = kHeight - kCanvasY - 64;

static constexpr float kFreqMin = 20.0f;
static constexpr float kFreqMax = 20000.0f;
static constexpr float kDbFloor = -100.0f;
static constexpr float kDbTop = 10.0f;
static constexpr float kDbRawFloor = -160.0f;
static constexpr float kMinGain = 1.0e-8f;

static constexpr float kAttackMs = 1.0f;
static constexpr float kReleaseMs = 100.0f;

// 低频显示平滑: 半宽(以 bin 计)在 kSmoothFullHz 以下全开, 到 kSmoothZeroHz 关闭
static constexpr float kSmoothHalfBins = 2.0f;
static constexpr float kSmoothFullHz = 500.0f;
static constexpr float kSmoothZeroHz = 2500.0f;

static const Color kBgColor = {16, 18, 22, 255};
static const Color kCanvasColor = {22, 26, 32, 255};
static const Color kGridColor = {48, 54, 64, 255};
static const Color kTextColor = {170, 178, 190, 255};
static const Color kCurveColor = {110, 210, 255, 255};
static const Color kMinColor = {60, 120, 160, 220};
static const Color kBandColor = {40, 110, 160, 120};
static const Color kCenterGridColor = {80, 92, 108, 255};

// ------------------------------------------------------------
// 分析状态
// ------------------------------------------------------------

static std::array<float, kWindowSize> g_window;
static std::array<float, kWindowSize> g_history;

/// 按下标恢复谐波时每个折叠段的相位 e^(-j2π·Δ·l)
static std::array<std::array<std::complex<float>, kFoldFactor>, kFoldFactor> g_seg_phase;
/// 频偏在折叠块内的逐点相位 e^(-j2π·Δ·i/N)
static std::array<std::array<std::complex<float>, kFftSize>, kFoldFactor> g_sample_phase;

static std::array<std::complex<float>, kFftSize> g_fold;
static std::array<std::complex<float>, kFftSize> g_fft_out;
static qwqdsp_spectral::ComplexFftAdv g_fft;

/// 音频线程拥有; GUI 线程只读发布缓冲
static std::array<float, kSpecBins> g_smoothed_db;
static std::array<std::array<float, kSpecBins>, 2> g_publish_db;
static std::atomic<int> g_publish_index = 0;

static float g_attack_alpha = 0.0f;
static float g_release_alpha = 0.0f;

static std::array<float, kHopSize> g_hop_buf;
static int g_hop_count = 0;

// ------------------------------------------------------------
// 显示查找表
// ------------------------------------------------------------

static std::array<int, kCanvasW> g_bin_lo;
static std::array<int, kCanvasW> g_bin_hi;
static std::array<float, kCanvasW> g_bin_pos;
static std::array<float, kCanvasW> g_mix;
static std::array<float, kSpecBins> g_slope;
static std::array<float, kSpecBins> g_render_db;
static std::array<float, kSpecBins + 1> g_db_prefix;
static std::array<float, kCanvasW> g_curve_min_y;
static std::array<float, kCanvasW> g_curve_max_y;

static float DbToY(float db) {
    // 低于 dB floor 的部分不夹到画布底边, 交给绘制裁剪区处理
    float const norm = (db - kDbFloor) / (kDbTop - kDbFloor);
    return static_cast<float>(kCanvasY + kCanvasH) - norm * static_cast<float>(kCanvasH);
}

static void Dsp_Init() {
    qwqdsp_window::Blackman::Window(g_window, true);
    qwqdsp_window::Helper::Normalize(g_window);
    g_history.fill(0.0f);
    g_smoothed_db.fill(kDbRawFloor);
    for (auto& buf : g_publish_db) {
        buf.fill(kDbRawFloor);
    }

    constexpr float twopi = std::numbers::pi_v<float> * 2.0f;
    for (int l = 0; l < kFoldFactor; ++l) {
        float delta = static_cast<float>(l) / static_cast<float>(kFoldFactor);
        for (int seg = 0; seg < kFoldFactor; ++seg) {
            float phase = -twopi * delta * static_cast<float>(seg);
            g_seg_phase[l][seg] = {std::cos(phase), std::sin(phase)};
        }
        for (int i = 0; i < kFftSize; ++i) {
            float phase = -twopi * delta * static_cast<float>(i) / static_cast<float>(kFftSize);
            g_sample_phase[l][i] = {std::cos(phase), std::sin(phase)};
        }
    }

    g_fft.Init(kFftSize);

    float const dt = static_cast<float>(kHopSize) / static_cast<float>(kSampleRate);
    g_attack_alpha = 1.0f - std::exp(-dt / (kAttackMs * 0.001f));
    g_release_alpha = 1.0f - std::exp(-dt / (kReleaseMs * 0.001f));
}

static void InitDisplay() {
    float const log_min = std::log(kFreqMin);
    float const log_span = std::log(kFreqMax) - log_min;
    // 预计算每列覆盖的 bin 范围, 供 min/max 采样
    for (int px = 0; px < kCanvasW; ++px) {
        float const f0 = std::exp(log_min + log_span * static_cast<float>(px) / static_cast<float>(kCanvasW));
        float const f1 = std::exp(log_min + log_span * static_cast<float>(px + 1) / static_cast<float>(kCanvasW));

        int lo = static_cast<int>(std::ceil(f0 / kBinHz - 0.5f));
        int hi = static_cast<int>(std::floor(f1 / kBinHz + 0.5f));
        lo = std::clamp(lo, 1, kSpecBins - 1);
        hi = std::clamp(hi, 1, kSpecBins - 1);
        if (hi < lo) {
            int const nearest = static_cast<int>(std::lround(0.5f * (f0 + f1) / kBinHz));
            lo = hi = std::clamp(nearest, 1, kSpecBins - 1);
        }
        g_bin_lo[px] = lo;
        g_bin_hi[px] = hi;
        float const f_mid = std::exp(log_min + log_span * (static_cast<float>(px) + 0.5f)
                                     / static_cast<float>(kCanvasW));
        g_bin_pos[px] = f_mid / kBinHz;
        // 用连续的小数 bin 跨度做混合权重, 避免整数 bin 计数造成的跳变
        float const span_bins = (f1 - f0) / kBinHz;
        g_mix[px] = std::clamp(span_bins - 1.0f, 0.0f, 1.0f);
    }
    g_curve_min_y.fill(DbToY(kDbFloor));
    g_curve_max_y.fill(DbToY(kDbFloor));
}

// ------------------------------------------------------------
// 折叠分析
// ------------------------------------------------------------

static void ProcessHop(float const* input) {
    std::memmove(g_history.data(), g_history.data() + kHopSize,
                 (kWindowSize - kHopSize) * sizeof(float));
    std::copy_n(input, kHopSize, g_history.data() + kWindowSize - kHopSize);

    // 时间折叠: y[i] = Σ_l x[lN+i]·w[lN+i]; 其 N 点 DFT 等于长窗信号在 2πk/N 处的采样。
    // 再对每个折叠段乘 e^(-j2π·Δ·l)、对块内样本乘 e^(-j2π·Δ·i/N), 把采样网格平移 Δ,
    // 于是 Δ=0..(kFoldFactor-1)/kFoldFactor 交织出长窗的全部谱线。
    for (int l = 0; l < kFoldFactor; ++l) {
        for (int i = 0; i < kFftSize; ++i) {
            float sum_re = 0.0f;
            float sum_im = 0.0f;
            for (int seg = 0; seg < kFoldFactor; ++seg) {
                int idx = seg * kFftSize + i;
                float value = g_history[idx] * g_window[idx];
                auto const phase = g_seg_phase[l][seg];
                sum_re += value * phase.real();
                sum_im += value * phase.imag();
            }
            auto const phase = g_sample_phase[l][i];
            g_fold[i] = {
                sum_re * phase.real() - sum_im * phase.imag(),
                sum_re * phase.imag() + sum_im * phase.real(),
            };
        }

        g_fft.FFT(std::span<std::complex<float> const>{g_fold},
                  std::span<std::complex<float>>{g_fft_out});
        for (int k = 0; k <= kFftSize / 2; ++k) {
            int bin = k * kFoldFactor + l;
            if (bin >= kSpecBins) {
                break;
            }
            float const gain = std::abs(g_fft_out[k]);
            float const db = 20.0f * std::log10(std::max(gain, kMinGain));
            float const prev = g_smoothed_db[bin];
            float const alpha = (db > prev) ? g_attack_alpha : g_release_alpha;
            g_smoothed_db[bin] = prev + alpha * (db - prev);
        }
    }

    int write_index = 1 - g_publish_index.load(std::memory_order_relaxed);
    g_publish_db[write_index] = g_smoothed_db;
    g_publish_index.store(write_index, std::memory_order_release);
}

extern "C" void MaCaptureCallback(ma_device* pDevice, void* pOutput,
                                  const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;
    (void)pOutput;
    float const* src = static_cast<float const*>(pInput);
    if (src == nullptr) {
        return;
    }

    int remaining = static_cast<int>(frameCount);
    while (remaining > 0) {
        int need = kHopSize - g_hop_count;
        int count = std::min(need, remaining);
        std::copy_n(src, count, g_hop_buf.data() + g_hop_count);
        g_hop_count += count;
        remaining -= count;
        src += count;

        if (g_hop_count == kHopSize) {
            ProcessHop(g_hop_buf.data());
            g_hop_count = 0;
        }
    }
}

static void FeedTestSignal() {
    std::array<float, kHopSize> hop{};
    int sample_index = 0;
    for (int h = 0; h < 80; ++h) {
        for (int i = 0; i < kHopSize; ++i, ++sample_index) {
            float const t = static_cast<float>(sample_index) / static_cast<float>(kSampleRate);
            hop[i] = std::sin(2.0f * std::numbers::pi_v<float> * 440.0f * t)
                   + 0.5f * std::sin(2.0f * std::numbers::pi_v<float> * 3000.0f * t)
                   + 0.25f * std::sin(2.0f * std::numbers::pi_v<float> * 10000.0f * t);
        }
        ProcessHop(hop.data());
    }
}

// ------------------------------------------------------------
// 曲线
// ------------------------------------------------------------

/// 每列覆盖的 bin 范围内取 min/max, 竖线段显示该列的动态范围。
/// 低频 bin 在像素上很稀疏, 用单调三次 Hermite 在 bin 间插值后再与 min/max 混合。
static void BuildCurve(std::span<float const> spec_db) {
    g_db_prefix[0] = 0.0f;
    for (int bin = 0; bin < kSpecBins; ++bin) {
        g_db_prefix[bin + 1] = g_db_prefix[bin] + spec_db[bin];
    }
    g_render_db[0] = spec_db[0];
    for (int bin = 1; bin < kSpecBins; ++bin) {
        float const f = static_cast<float>(bin) * kBinHz;
        float const fade = std::clamp((kSmoothZeroHz - f) / (kSmoothZeroHz - kSmoothFullHz), 0.0f, 1.0f);
        if (fade <= 0.0f) {
            g_render_db[bin] = spec_db[bin];
            continue;
        }
        int const half = std::max(1, static_cast<int>(std::lround(kSmoothHalfBins * fade)));
        int const lo = std::max(1, bin - half);
        int const hi = std::min(kSpecBins - 1, bin + half);
        g_render_db[bin] = (g_db_prefix[hi + 1] - g_db_prefix[lo]) / static_cast<float>(hi - lo + 1);
    }

    g_slope.fill(0.0f);
    for (int bin = 1; bin + 1 < kSpecBins; ++bin) {
        float const d0 = g_render_db[bin] - g_render_db[bin - 1];
        float const d1 = g_render_db[bin + 1] - g_render_db[bin];
        g_slope[bin] = (d0 * d1 > 0.0f) ? 2.0f * d0 * d1 / (d0 + d1) : 0.0f;
    }

    for (int px = 0; px < kCanvasW; ++px) {
        int const lo = g_bin_lo[px];
        int const hi = g_bin_hi[px];
        float min_db = g_render_db[lo];
        float max_db = g_render_db[lo];
        for (int bin = lo + 1; bin <= hi; ++bin) {
            min_db = std::min(min_db, g_render_db[bin]);
            max_db = std::max(max_db, g_render_db[bin]);
        }

        float const pos = g_bin_pos[px];
        int const bin0 = std::clamp(static_cast<int>(pos), 1, kSpecBins - 2);
        float const t = std::clamp(pos - static_cast<float>(bin0), 0.0f, 1.0f);
        float const t2 = t * t;
        float const t3 = t2 * t;
        float const h00 = 2.0f * t3 - 3.0f * t2 + 1.0f;
        float const h10 = t3 - 2.0f * t2 + t;
        float const h01 = -2.0f * t3 + 3.0f * t2;
        float const h11 = t3 - t2;
        float const smooth_db = h00 * g_render_db[bin0] + h10 * g_slope[bin0]
                              + h01 * g_render_db[bin0 + 1] + h11 * g_slope[bin0 + 1];

        float const mix = g_mix[px];
        min_db = smooth_db + mix * (min_db - smooth_db);
        max_db = smooth_db + mix * (max_db - smooth_db);

        g_curve_min_y[px] = DbToY(min_db);
        g_curve_max_y[px] = DbToY(max_db);
    }
}

static float FreqToX(float freq) {
    float t = (std::log(freq) - std::log(kFreqMin)) / (std::log(kFreqMax) - std::log(kFreqMin));
    return static_cast<float>(kCanvasX) + t * static_cast<float>(kCanvasW);
}

// ------------------------------------------------------------
// 绘制
// ------------------------------------------------------------

static void DrawGrid() {
    DrawRectangle(kCanvasX, kCanvasY, kCanvasW, kCanvasH, kCanvasColor);

    static constexpr float kYTicks[] = {60.0f, 40.0f, 20.0f, 0.0f, -20.0f, -40.0f, -60.0f};
    for (float db : kYTicks) {
        float y = DbToY(db);
        Color const color = (db == 0.0f) ? kCenterGridColor : kGridColor;
        DrawLineV({static_cast<float>(kCanvasX), y},
                  {static_cast<float>(kCanvasX + kCanvasW), y}, color);
        auto label = std::format("{:.0f}", db);
        DrawText(label.c_str(), kCanvasX - MeasureText(label.c_str(), 12) - 8,
                 static_cast<int>(y) - 6, 12, kTextColor);
    }
    DrawText("dB", 18, kCanvasY - 18, 12, kTextColor);

    // 每个 decade 十条对数网格线
    static constexpr float kDecadeStart[] = {20.0f, 200.0f, 2000.0f};
    for (float start : kDecadeStart) {
        float step = start;
        for (float freq = start; freq <= start * 10.0f + 0.5f; freq += step) {
            float x = FreqToX(freq);
            DrawLineV({x, static_cast<float>(kCanvasY)},
                      {x, static_cast<float>(kCanvasY + kCanvasH)}, kGridColor);
        }
    }

    static constexpr float kXTicks[] = {20.0f, 200.0f, 2000.0f, 20000.0f};
    for (float freq : kXTicks) {
        float x = FreqToX(freq);
        auto label = std::format("{:.0f}", freq >= 1000.0f ? freq / 1000.0f : freq);
        if (freq >= 1000.0f) {
            label += "k";
        }
        DrawText(label.c_str(), static_cast<int>(x) - MeasureText(label.c_str(), 12) / 2,
                 kCanvasY + kCanvasH + 6, 12, kTextColor);
    }
    DrawText("Hz", 6, kCanvasY + kCanvasH + 6, 12, kTextColor);
}

static void DrawSpectrum() {
    // 竖线段显示每列的 min/max 动态范围, 再叠 min/max 包络
    static std::array<Vector2, kCanvasW> max_points;
    static std::array<Vector2, kCanvasW> min_points;
    float const x0 = static_cast<float>(kCanvasX);
    for (int px = 0; px < kCanvasW; ++px) {
        float const x = x0 + static_cast<float>(px);
        float const y_min = g_curve_min_y[px];
        float const y_max = g_curve_max_y[px];
        DrawLineV({x, y_min}, {x, y_max}, kBandColor);
        min_points[px] = {x, y_min};
        max_points[px] = {x, y_max};
    }

    DrawLineStrip(min_points.data(), kCanvasW, kMinColor);
    DrawLineStrip(max_points.data(), kCanvasW, kCurveColor);
}

int main(void) {
    qwqdsp_gui_test::ShotMode<kGuiTest> shot(kShotPath);
    shot.ApplyWindowFlags();

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    auto title = std::format("Portnoff Spectrum - {}-point spectrum from {}-point FFT",
                             kWindowSize, kFftSize);
    InitWindow(kWidth, kHeight, title.c_str());
    SetTargetFPS(60);

    Dsp_Init();
    InitDisplay();

    ma_device device{};
    auto start_audio = [&device]() -> ma_result {
        ma_device_config config = ma_device_config_init(ma_device_type_loopback);
        config.capture.format = ma_format_f32;
        config.capture.channels = 1;
        config.sampleRate = static_cast<ma_uint32>(kSampleRate);
        config.dataCallback = MaCaptureCallback;
        config.pUserData = nullptr;
        config.periodSizeInMilliseconds = 10;

        ma_result result = ma_device_init(nullptr, &config, &device);
        if (result == MA_SUCCESS) {
            ma_device_start(&device);
        }
        else {
            TraceLog(LOG_WARNING, "loopback capture init failed: %s", ma_result_description(result));
        }
        return result;
    };

    ma_result result = MA_DEVICE_NOT_STARTED;
    if (shot.IsEnabled()) {
        FeedTestSignal();
    }
    else {
        result = start_audio();
    }

    while (!WindowShouldClose()) {
        int read_index = g_publish_index.load(std::memory_order_acquire);
        BuildCurve(g_publish_db[read_index]);

        BeginDrawing();
        ClearBackground(kBgColor);
        DrawText("Portnoff fold spectrum", kCanvasX, 12, 18, kTextColor);
        auto info = std::format("FFT {} x{} -> {}  |  {:.1f} Hz/bin  |  min-max",
                                kFftSize, kFoldFactor, kWindowSize, kBinHz);
        DrawText(info.c_str(), kCanvasX + 210, 16, 12, kTextColor);

        DrawGrid();
        BeginScissorMode(kCanvasX, kCanvasY, kCanvasW, kCanvasH);
        DrawSpectrum();
        EndScissorMode();
        DrawFPS(kWidth - 90, 12);
        EndDrawing();

        if (shot.Tick()) {
            break;
        }
    }

    if (result == MA_SUCCESS) {
        ma_device_stop(&device);
        ma_device_uninit(&device);
    }
    CloseWindow();
    return 0;
}
