// ────────────────────────────────────────────────────────────
//  polyphase_filter_bank_view.cpp
//  实时多相分析滤波器组 (DFT 调制) + 对数频率插值曲线可视化
//
//  信号链:
//    miniaudio 双工 (loopback) → 帧累加器
//    → 多相分析滤波器组 (M=512, L=16384, Blackman 窗) ← 音频线程
//    → 子带幅度 (dB) → attack/release 平滑 → s_db_latest[257]
//                                                      ↓
//                        主线程读取 → 对数频率插值 → 曲线渲染
// ────────────────────────────────────────────────────────────

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <numbers>
#include <span>

#include "miniaudio.h"
#include "raylib.h"

#include <qwqdsp/filter/window_fir.hpp>
#include <qwqdsp/window/blackman.hpp>
#include <qwqdsp/window/hann.hpp>
#include <qwqdsp/window/helper.hpp>

// ════════════════════════════════════════════════════════════
//  用户可调常量
// ════════════════════════════════════════════════════════════

// ── DSP ──
static constexpr float kSampleRate = 48000.0f;
static constexpr int kNumSubbands = 512;                           // M
static constexpr int kPolyphaseLen = 32;                           // P — 每相位抽头数
static constexpr int kPrototypeLen = kNumSubbands * kPolyphaseLen; // L = M × P

// ── 平滑 (一阶 IIR, attack/release) ──
static constexpr float kAttackMs = 1.0f;
static constexpr float kReleaseMs = 150.0f;

// ── 窗口 & 画布 ──
static constexpr int kWindowWidth = 640;
static constexpr int kWindowHeight = 320;
static constexpr int kMarginLeft = 65;
static constexpr int kMarginRight = 20;
static constexpr int kMarginTop = 45;
static constexpr int kMarginBottom = 55;
static constexpr int kCanvasWidth = kWindowWidth - kMarginLeft - kMarginRight;
static constexpr int kCanvasHeight = kWindowHeight - kMarginTop - kMarginBottom;
static constexpr float kDbFloor = -80.0f;
static constexpr float kDbRange = 80.0f; // Y 轴动态范围

// ── 频率显示范围 & 网格 ──
static constexpr float kDisplayFreqMin = 20.0f;
static constexpr float kDisplayFreqMax = 20000.0f;
static constexpr float kXTickLabels[] = {20, 200, 2000, 20000};
static const Color kGridColor = {50, 50, 50, 255};
static const Color kTextColor = {180, 180, 180, 255};
static const Color kCurveColor = {120, 200, 255, 255};

// ════════════════════════════════════════════════════════════
//  音频线程 → 主线程 共享 dB 缓冲区
// ════════════════════════════════════════════════════════════
static constexpr int kDisplayBands = kNumSubbands / 2 + 1;
static std::array<float, kDisplayBands> s_db_latest;

// ── 回调耗时统计 (音频线程写入, 主线程读取) ──
static float s_cb_us_last = 0.0f; // 最近一次 (µs)
static float s_cb_us_avg = 0.0f;  // 指数平滑平均
static float s_cb_us_max = 0.0f;  // 观测峰值 (主线程每秒重置)

// ════════════════════════════════════════════════════════════
//  多相分析滤波器组
// ════════════════════════════════════════════════════════════
class PolyphaseAnalysisFB {
    // 多相系数 & 延迟线: [phase * kPolyphaseLen + tap]
    alignas(32) std::array<float, kNumSubbands * kPolyphaseLen> coeffs_;
    alignas(32) std::array<float, kNumSubbands * kPolyphaseLen> delay_;
    // FFT 工作缓冲区
    std::array<float, kNumSubbands> fft_real_;
    std::array<float, kNumSubbands> fft_imag_;
    // 平滑后的子带 dB 值 (bin 0 ~ M/2)
    std::array<float, kNumSubbands / 2 + 1> smoothed_db_;
    // 一阶平滑系数
    float attack_alpha_;
    float release_alpha_;
public:
    // ── 初始化: 设计原型滤波器 → 多相分解 → 计算平滑系数 ──
    void Init() {
        // 1. 原型低通 FIR (截止 π/M)
        std::array<float, kPrototypeLen> ir;
        qwqdsp_filter::WindowFIR::Lowpass(ir, std::numbers::pi_v<float> / static_cast<float>(kNumSubbands));
        qwqdsp_window::Blackman::ApplyWindow(ir, false); // 对称窗 → FIR 设计
        qwqdsp_filter::WindowFIR::Normalize(ir);         // DC 增益 = 1

        // 2. 多相分解: coeffs_[k][j] = ir[k + j*M]
        //    k = 相位索引 (0..M-1), j = 子滤波器抽头 (0..P-1)
        for (int k = 0; k < kNumSubbands; ++k) {
            for (int j = 0; j < kPolyphaseLen; ++j) {
                int src = k + j * kNumSubbands;
                coeffs_[k * kPolyphaseLen + j] = (src < kPrototypeLen) ? ir[src] : 0.0f;
            }
        }

        // 3. 重置状态
        Reset();

        // 4. 平滑时间常数 → α
        float T = static_cast<float>(kNumSubbands) / kSampleRate; // 帧周期
        attack_alpha_ = 1.0f - std::exp(-T / (kAttackMs * 0.001f));
        release_alpha_ = 1.0f - std::exp(-T / (kReleaseMs * 0.001f));
    }

    void Reset() {
        std::fill(delay_.begin(), delay_.end(), 0.0f);
        std::fill(smoothed_db_.begin(), smoothed_db_.end(), kDbFloor);
    }

    // ── 处理一帧 (M 个输入样本) ──
    void Process(const float* input) {
        // Step 1 — 多相滤波 (转置 FIR)
        for (int k = 0; k < kNumSubbands; ++k) {
            float* d = &delay_[k * kPolyphaseLen];
            const float* c = &coeffs_[k * kPolyphaseLen];

            // 反向分配: branch k 接收 input[M-1-k]
            // (commutator 模型 — 最新样本进 branch 0)
            float x = input[kNumSubbands - 1 - k];

            // 转置直接 II 型 FIR
            float u = c[0] * x + d[0];
            for (int j = 0; j < kPolyphaseLen - 1; ++j) {
                d[j] = d[j + 1] + c[j + 1] * x;
            }
            d[kPolyphaseLen - 1] = c[kPolyphaseLen - 1] * x;

            fft_real_[k] = u;
            fft_imag_[k] = 0.0f;
        }

        // Step 2 — M 点复数 FFT (基 2 DIT)
        FFT();

        // Step 3 — 子带幅度 → dB → 一阶平滑
        constexpr int kDisplayBands = kNumSubbands / 2 + 1; // = 9
        for (int k = 0; k < kDisplayBands; ++k) {
            float mag_sq = fft_real_[k] * fft_real_[k] + fft_imag_[k] * fft_imag_[k];
            float dB = 10.0f * std::log10(mag_sq + 1e-12f);

            float prev = smoothed_db_[k];
            float alpha = (dB > prev) ? attack_alpha_ : release_alpha_;
            smoothed_db_[k] = prev + alpha * (dB - prev);
        }
    }

    // ── 存取器 ──
    std::span<const float> GetDB() const noexcept {
        return smoothed_db_;
    }
    int NumBands() const noexcept {
        return smoothed_db_.size();
    }
private:
    // ── 基 2 按时间抽取 (DIT) FFT, 原地计算 ──
    //    要求 kNumSubbands 为 2 的幂
    void FFT() noexcept {
        constexpr int n = kNumSubbands;

        // 位反转重排
        for (int i = 1, j = 0; i < n; ++i) {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
            if (i < j) {
                std::swap(fft_real_[i], fft_real_[j]);
                std::swap(fft_imag_[i], fft_imag_[j]);
            }
        }

        // 蝶形运算
        for (int len = 2; len <= n; len <<= 1) {
            float ang = 2.0f * std::numbers::pi_v<float> / static_cast<float>(len);
            float w_real = std::cos(ang);
            float w_imag = -std::sin(ang); // exp(-j·ang)

            for (int i = 0; i < n; i += len) {
                float wr = 1.0f, wi = 0.0f;
                int half = len >> 1;
                for (int j = 0; j < half; ++j) {
                    int even = i + j;
                    int odd = i + j + half;

                    float tr = wr * fft_real_[odd] - wi * fft_imag_[odd];
                    float ti = wr * fft_imag_[odd] + wi * fft_real_[odd];

                    fft_real_[odd] = fft_real_[even] - tr;
                    fft_imag_[odd] = fft_imag_[even] - ti;
                    fft_real_[even] = fft_real_[even] + tr;
                    fft_imag_[even] = fft_imag_[even] + ti;

                    float nwr = wr * w_real - wi * w_imag;
                    float nwi = wr * w_imag + wi * w_real;
                    wr = nwr;
                    wi = nwi;
                }
            }
        }
    }
};

// ── 全局实例 ──
static PolyphaseAnalysisFB s_pfb;

// ── 帧累加器 (音频线程局部变量) ──
static std::array<float, kNumSubbands> s_frame_acc{};
static int s_frame_count = 0;

// ════════════════════════════════════════════════════════════
//  miniaudio 回调 (audio 线程) — 采样累加 → 滤波器组 → 共享 dB
// ════════════════════════════════════════════════════════════
extern "C" void MaCaptureCallback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
    (void)pOutput;
    auto t0 = std::chrono::steady_clock::now();

    float const* src = static_cast<float const*>(pInput);

    // 逐采样累加 → 满一帧后运行滤波器组
    for (ma_uint32 i = 0; i < frameCount; ++i) {
        s_frame_acc[s_frame_count++] = src[i];
        if (s_frame_count == kNumSubbands) {
            s_frame_count = 0;
            s_pfb.Process(s_frame_acc.data());
            auto db = s_pfb.GetDB();
            std::copy_n(db.begin(), db.size(), s_db_latest.begin());
        }
    }

    // ── 耗时统计 ──
    auto t1 = std::chrono::steady_clock::now();
    float us = std::chrono::duration<float, std::micro>(t1 - t0).count();
    s_cb_us_last = us;
    s_cb_us_avg += 0.05f * (us - s_cb_us_avg);
    if (us > s_cb_us_max)
        s_cb_us_max = us;
}

// ════════════════════════════════════════════════════════════
//  可视化 (raylib)
// ════════════════════════════════════════════════════════════

static void DrawAxes() {
    constexpr float kBinFreq = kSampleRate / kNumSubbands;
    float logFMin = std::log10(kDisplayFreqMin);
    float logFRange = std::log10(kDisplayFreqMax) - logFMin;

    // ── Y 轴网格 & 刻度 (dB) ──
    for (float db = 0.0f; db >= kDbFloor; db -= 20.0f) {
        float yy = kMarginTop + static_cast<float>(kCanvasHeight) * (1.0f - (db - kDbFloor) / kDbRange);
        int y = static_cast<int>(yy);
        DrawLine(kMarginLeft, y, kMarginLeft + kCanvasWidth, y, kGridColor);
        DrawText(TextFormat("%.0f", db), 10, y - 8, 12, kTextColor);
    }

    // ── X 轴对数频率网格竖线 ──
    constexpr float kDecadeStart[] = {20, 200, 2000};
    for (float start : kDecadeStart) {
        float step = start;
        for (float freq = start; freq <= start * 10.0f && freq <= kDisplayFreqMax; freq += step) {
            if (freq < kDisplayFreqMin)
                continue;
            float norm = (std::log10(freq) - logFMin) / logFRange;
            int x = kMarginLeft + static_cast<int>(norm * kCanvasWidth);
            DrawLine(x, kMarginTop, x, kMarginTop + kCanvasHeight, kGridColor);
        }
    }

    // ── X 轴频率标签 ──
    for (float freq : kXTickLabels) {
        if (freq < kDisplayFreqMin || freq > kDisplayFreqMax)
            continue;
        float norm = (std::log10(freq) - logFMin) / logFRange;
        int x = kMarginLeft + static_cast<int>(norm * kCanvasWidth);
        const char* lbl = (freq >= 1000.0f) ? TextFormat("%.0fk", freq / 1000.0f) : TextFormat("%.0f", freq);
        int tw = MeasureText(lbl, 10);
        DrawText(lbl, x - tw / 2, kMarginTop + kCanvasHeight + 4, 10, kTextColor);
    }

    // ── 坐标轴 ──
    DrawLine(kMarginLeft, kMarginTop, kMarginLeft, kMarginTop + kCanvasHeight, WHITE);
    DrawLine(kMarginLeft, kMarginTop + kCanvasHeight, kMarginLeft + kCanvasWidth, kMarginTop + kCanvasHeight, WHITE);
}

static void DrawSubbandCurve(std::span<const float> dB) {
    constexpr float kBinFreq = kSampleRate / kNumSubbands;
    float logFMin = std::log10(kDisplayFreqMin);
    float logFRange = std::log10(kDisplayFreqMax) - logFMin;

    // 固定大小点数组 (避免每帧分配)
    static std::array<Vector2, kCanvasWidth> pts;

    for (int px = 0; px < kCanvasWidth; ++px) {
        float normX = static_cast<float>(px) / (kCanvasWidth - 1.0f);
        float freq = std::pow(10.0f, logFMin + normX * logFRange);
        float binF = freq / kBinFreq;

        int seg = static_cast<int>(binF); // floor
        seg = std::clamp(seg, 0, static_cast<int>(dB.size()) - 2);
        int base = std::clamp(seg - 1, 0, static_cast<int>(dB.size()) - 3);

        // 三点 Lagrange 抛物线插值 (通过 base, base+1, base+2)
        float d = binF - static_cast<float>(base + 1);
        float dbVal =
            dB[base] * d * (d - 1.0f) * 0.5f + dB[base + 1] * (1.0f - d * d) + dB[base + 2] * d * (d + 1.0f) * 0.5f;
        dbVal = std::clamp(dbVal, kDbFloor, 0.0f);

        float normY = (dbVal - kDbFloor) / kDbRange;
        float y = kMarginTop + kCanvasHeight - normY * kCanvasHeight;

        pts[px].x = static_cast<float>(kMarginLeft + px);
        pts[px].y = y;
    }

    DrawLineStrip(pts.data(), kCanvasWidth, kCurveColor);
}

// ════════════════════════════════════════════════════════════
//  main
// ════════════════════════════════════════════════════════════
int main(void) {
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, "Polyphase Analysis Filter Bank — qwqdsp + raylib");
    SetTargetFPS(60);

    // ── 初始化滤波器组 & 共享缓冲区 ──
    s_pfb.Init();
    s_db_latest.fill(kDbFloor);

    // ── miniaudio 回环捕获 (loopback) ──
    ma_device_config config = ma_device_config_init(ma_device_type_loopback);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.sampleRate = static_cast<ma_uint32>(kSampleRate);
    config.dataCallback = MaCaptureCallback;
    config.pUserData = nullptr;
    config.periodSizeInMilliseconds = 10;

    ma_device device;
    ma_result result = ma_device_init(nullptr, &config, &device);
    bool audio_ok = (result == MA_SUCCESS);
    if (audio_ok) {
        ma_device_start(&device);
    }
    else {
        TraceLog(LOG_WARNING, "miniaudio init failed — running in silent mode");
    }

    // ── 主循环 ──
    while (!WindowShouldClose()) {
        // 从共享缓冲区取最新 dB 快照 (栈拷贝, ~1 µs)
        auto local_db = s_db_latest;

        // ── 渲染 ──
        BeginDrawing();
        ClearBackground({18, 18, 18, 255});

        DrawAxes();
        DrawSubbandCurve(local_db);

        // 峰值信息
        float peak_db = -999.0f;
        int peak_bin = 0;
        for (int i = 0; i < static_cast<int>(local_db.size()); ++i) {
            if (local_db[i] > peak_db) {
                peak_db = local_db[i];
                peak_bin = i;
            }
        }
        float peak_freq = static_cast<float>(peak_bin) * kSampleRate / static_cast<float>(kNumSubbands);

        // 每秒重置峰值
        static int s_cb_reset = 0;
        if (++s_cb_reset >= 60) {
            s_cb_reset = 0;
            s_cb_us_max = 0.0f;
        }

        DrawText(TextFormat("M=%d  L=%d  |  "
                            "Peak: bin %d = %.1f Hz @ %.1f dB  |  "
                            "CB: %.0f us (avg %.0f, max %.0f)  |  FPS: %d",
                            kNumSubbands, kPrototypeLen, peak_bin, peak_freq, peak_db, s_cb_us_last, s_cb_us_avg,
                            s_cb_us_max, GetFPS()),
                 10, 10, 14, WHITE);

        EndDrawing();
    }

    // ── 清理 ──
    if (audio_ok) {
        ma_device_stop(&device);
        ma_device_uninit(&device);
    }
    CloseWindow();
    return 0;
}
