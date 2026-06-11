// ────────────────────────────────────────────────────────────
//  polyphase_filter_bank_view2.cpp
//  实时余弦调制 (CMFB / Pseudo-QMF) 分析滤波器组 + 对数插值曲线
//
//  信号链:
//    miniaudio → 帧累加器
//    → CMFB (M=512, 2M 多相 + DCT-IV) ← 音频线程
//    → |Xₖ| → dB → 平滑 → s_db_latest[512]
//                                      ↓
//                    主线程 → 对数插值 → 曲线渲染
//
//  参考实现: CosineModulatedFilterBank (AI gen, 已整合至此)
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
#include <vector>

#include "../../playing/miniaudio.h"
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
static constexpr int kNumSubbands = 256;                             // M
static constexpr int kPolyphaseLen = 32;                             // P₁ — 原每相位抽头数
static constexpr int kPrototypeLen = kNumSubbands * kPolyphaseLen;   // L = M × P₁ = 16384
static constexpr int kPolyphase2M = 2 * kNumSubbands;                // 2M = 1024
static constexpr int kPolyphase2MLen = kPrototypeLen / kPolyphase2M; // P₂ = L / 2M = 16

// ── 平滑 (一阶 IIR, attack/release) ──
static constexpr float kAttackMs = 1.0f;
static constexpr float kReleaseMs = 150.0f;

// ── 窗口 & 画布 ──
static constexpr int kWindowWidth = 1280;
static constexpr int kWindowHeight = 720;
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
static constexpr int kDisplayBands = kNumSubbands; // M 个实数子带
static std::array<float, kDisplayBands> s_db_latest;

// ── 回调耗时统计 (音频线程写入, 主线程读取) ──
static float s_cb_us_last = 0.0f; // 最近一次 (µs)
static float s_cb_us_avg = 0.0f;  // 指数平滑平均
static float s_cb_us_max = 0.0f;  // 观测峰值 (主线程每秒重置)

// ════════════════════════════════════════════════════════════
//  CMFB 分析滤波器组 (余弦调制 / Pseudo-QMF)
//  2M 路多相 + (-1)ᵐ 符号交替 + 余弦对称合并 + DCT-IV
// ════════════════════════════════════════════════════════════
class CosineModulatedFB {
    // 2M 多相系数: [phase × kPolyphase2MLen + tap]
    alignas(32) std::array<float, kPolyphase2M * kPolyphase2MLen> coeffs_;
    // 共享环形缓冲区 (原型滤波器全长)
    alignas(32) std::array<float, kPrototypeLen> ring_buf_{};
    int wp_ = 0; // 环形缓冲区写指针
    // 多相滤波输出 [2M]
    std::array<float, kPolyphase2M> out_;
    // 余弦调制矩阵 (堆分配, M × 2M)
    std::vector<float> dct_mat_;
    // 平滑后的 dB 值 (M 个实数子带)
    std::array<float, kNumSubbands> smoothed_db_;
    float attack_alpha_;
    float release_alpha_;
public:
    // ── 初始化 ──
    void Init() {
        // 1. 原型低通 FIR (截止 π/(2M), 子带带宽减半)
        std::array<float, kPrototypeLen> ir;
        qwqdsp_filter::WindowFIR::Lowpass(ir, std::numbers::pi_v<float> / static_cast<float>(2 * kNumSubbands));
        qwqdsp_window::Hann::ApplyWindow(ir, false);
        qwqdsp_filter::WindowFIR::Normalize(ir);

        // 2. 2M 多相分解: coeffs[l·P + m] = ir[2M·m + l] × (-1)^m
        //    l = 0..2M-1 (相位), m = 0..P₂-1 (子滤波器抽头)
        for (int l = 0; l < kPolyphase2M; ++l) {
            for (int m = 0; m < kPolyphase2MLen; ++m) {
                int src = kPolyphase2M * m + l;
                float sign = (m % 2 == 0) ? 1.0f : -1.0f; // (-1)^m
                coeffs_[l * kPolyphase2MLen + m] = (src < kPrototypeLen) ? ir[src] * sign : 0.0f;
            }
        }

        // 3. 预计算余弦调制矩阵 [M × 2M]
        //    C[k][i] = cos(π/M·(k+½)·(i − (L-1)/2))
        constexpr float kCenter = static_cast<float>(kPrototypeLen - 1) * 0.5f;
        dct_mat_.resize(static_cast<size_t>(kNumSubbands) * kPolyphase2M);
        for (int k = 0; k < kNumSubbands; ++k) {
            float phase_k =
                (static_cast<float>(k) + 0.5f) * std::numbers::pi_v<float> / static_cast<float>(kNumSubbands);
            for (int i = 0; i < kPolyphase2M; ++i) {
                float angle = phase_k * (static_cast<float>(i) - kCenter);
                dct_mat_[k * kPolyphase2M + i] = 2 * std::cos(angle);
            }
        }

        Reset();

        // 4. 平滑系数
        float T = static_cast<float>(kNumSubbands) / kSampleRate;
        attack_alpha_ = 1.0f - std::exp(-T / (kAttackMs * 0.001f));
        release_alpha_ = 1.0f - std::exp(-T / (kReleaseMs * 0.001f));
    }

    void Reset() {
        std::fill(ring_buf_.begin(), ring_buf_.end(), 0.0f);
        std::fill(smoothed_db_.begin(), smoothed_db_.end(), kDbFloor);
        wp_ = 0;
    }

    // ── 处理一帧 (M 个输入样本) ──
    void Process(const float* input) {
        // Step 1 — 写入 M 个新采样到环形缓冲区
        for (int i = 0; i < kNumSubbands; ++i)
            ring_buf_[(wp_ + i) % kPrototypeLen] = input[i];

        // Step 2 — M 对直接型 FIR, 从共享环形缓冲区读取
        //   分支 l:   最新采样在 (wp+M-1-l)     (当前帧)
        //   分支 l+M: 最新采样在 (wp+M-1-l-M)   (上一帧同位置)
        //   后续抽头: 步进 -2M
        for (int l = 0; l < kNumSubbands; ++l) {
            // ── 分支 l (当前帧) ──
            {
                const float* c = &coeffs_[l * kPolyphase2MLen];
                float sum = 0.0f;
                int idx = (wp_ + kNumSubbands - 1 - l + kPrototypeLen) % kPrototypeLen;
                for (int m = 0; m < kPolyphase2MLen; ++m) {
                    sum += c[m] * ring_buf_[idx];
                    idx = (idx - kPolyphase2M + kPrototypeLen) % kPrototypeLen;
                }
                out_[l] = sum;
            }

            // ── 分支 l+M (上一帧同位置) ──
            {
                int l2 = l + kNumSubbands;
                const float* c = &coeffs_[l2 * kPolyphase2MLen];
                float sum = 0.0f;
                int idx = (wp_ + kNumSubbands - 1 - l - kNumSubbands + kPrototypeLen) % kPrototypeLen;
                for (int m = 0; m < kPolyphase2MLen; ++m) {
                    sum += c[m] * ring_buf_[idx];
                    idx = (idx - kPolyphase2M + kPrototypeLen) % kPrototypeLen;
                }
                out_[l2] = sum;
            }
        }

        // ── 推进写指针 ──
        wp_ = (wp_ + kNumSubbands) % kPrototypeLen;

        // Step 2 — 直接余弦调制: X[k] = Σ out[i] · cos(π/M·(k+½)·(i−(L-1)/2))
        for (int k = 0; k < kNumSubbands; ++k) {
            const float* row = &dct_mat_[k * kPolyphase2M];
            float sum = 0.0f;
            for (int i = 0; i < kPolyphase2M; ++i) {
                sum += row[i] * out_[i];
            }

            // Step 4 — 幅度 → dB → 一阶平滑
            // 修正：幅度转 dB 需要乘 20，或者 log10(sum * sum)
            float dB = 20.0f * std::log10(std::abs(sum) + 1e-12f);
            float prev = smoothed_db_[k];
            float alpha = (dB > prev) ? attack_alpha_ : release_alpha_;
            smoothed_db_[k] = prev + alpha * (dB - prev);
        }
    }

    std::span<const float> GetDB() const noexcept {
        return smoothed_db_;
    }
    int NumBands() const noexcept {
        return kNumSubbands;
    }
};

// ── 全局实例 ──
static CosineModulatedFB s_pfb;

// ── 帧累加器 (音频线程局部变量) ──
static std::array<float, kNumSubbands> s_frame_acc{};
static int s_frame_count = 0;

// ════════════════════════════════════════════════════════════
//  miniaudio 回调 (audio 线程) — 采样累加 → 滤波器组 → 共享 dB
// ════════════════════════════════════════════════════════════
extern "C" void MaCaptureCallback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
    auto t0 = std::chrono::steady_clock::now();

    float const* src = static_cast<float const*>(pInput);
    float* dst = static_cast<float*>(pOutput);

    // 直通 (loopback)
    std::copy_n(src, frameCount, dst);

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
    // CMFB: 子带间距 = fs / (2M), 中心频率 = (k+0.5) × fs/(2M)
    constexpr float kCMFMBinSpacing = kSampleRate / (2.0f * kNumSubbands);
    float logFMin = std::log10(kDisplayFreqMin);
    float logFRange = std::log10(kDisplayFreqMax) - logFMin;

    // 固定大小点数组 (避免每帧分配)
    static std::array<Vector2, kCanvasWidth> pts;

    for (int px = 0; px < kCanvasWidth; ++px) {
        float normX = static_cast<float>(px) / (kCanvasWidth - 1.0f);
        float freq = std::pow(10.0f, logFMin + normX * logFRange);
        float binF = freq / kCMFMBinSpacing - 0.5f;

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
    InitWindow(kWindowWidth, kWindowHeight, "CMFB (Cosine Modulated) Filter Bank — qwqdsp + raylib");
    SetTargetFPS(60);

    // ── 初始化滤波器组 & 共享缓冲区 ──
    s_pfb.Init();
    s_db_latest.fill(kDbFloor);

    // ── miniaudio 全双工 (loopback) ──
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
