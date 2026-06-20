#include <atomic>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <vector>

#include "raylib.h"
#include "miniaudio.h"

#include "qwqdsp/spectral/ipp_real_fft.hpp"

// 频谱驱动 (81 bin)
#include "spectrum_driver.h"

static constexpr int kWindowWidth = 1100;
static constexpr int kWindowHeight = 640;
static constexpr float kSampleRate = 44100;

enum class SpectrumViewMode {
    Bins,
    Mel,
};

// ── 频谱分析全局状态 ──
static SpectrumContext gSpectrumCtx;
static constexpr uint16_t kSptLen = SPECTRUM_LEN;
static constexpr uint8_t kBinNum = SPECTRUM_USER_BIN_NUM;
static constexpr uint8_t kMelBandNum = SPECTRUM_MEL_FILTER_NUM;

static int16_t gPcmBuffer[kSptLen]; // int16 PCM 帧
static float gFloatRing[kSptLen];   // float 环形暂存
static size_t gFloatPos = 0;        // 当前写入位置

static uint8_t gLevels[kBinNum];          // 量化电平
static uint8_t gMelLevels[kMelBandNum];   // Mel量化电平
static std::atomic<bool> gNewData{false}; // 渲染线程可消费标志
static SpectrumViewMode gViewMode = SpectrumViewMode::Bins;

// ── 渲染颜色 ──
static constexpr Color kBarLow  = { 0,   255, 0,   255 };   // 绿
static constexpr Color kBarMid  = { 255, 255, 0,   255 };   // 黄
static constexpr Color kBarHigh = { 255, 0,   0,   255 };   // 红
static constexpr Color kBarBg   = { 20,   20,  20,  255 };

// ── rfft_api: 使用 IPP RealFFT (浮点 CCS) 包装为 int32 Ooura 兼容格式 ──
//    输入 buf[0..n-1] = int32 时域数据 (经汉明窗加权)
//    IPP 输出 CCS 格式: [DC_re, DC_im(=0), bin1_re, bin1_im, ..., Ny_re, Ny_im(=0)]
//    前 n 个值恰好与 Ooura 打包布局兼容 (spectrum_process 会置零 buf[0,1]):
//        buf[0]   = DC real
//        buf[1]   = DC imag (=0, spectrum_process 会置零)
//        buf[2k]  = bin_k real
//        buf[2k+1] = bin_k imag   (k = 1 .. n/2-1)
//    返回 1 成功, 0 失败
int rfft_api(int32_t *buf, int n, int inverse) {
    static qwqdsp_spectral::IppRealFFT s_fft;
    static int s_prev_n = 0;

    if (inverse != 1) return 0;
    if (s_prev_n != n) {
        s_fft.Init(n);
        s_prev_n = n;
    }

    // 1) int32 → float 时域
    std::vector<float> time(n);
    for (int i = 0; i < n; ++i)
        time[i] = static_cast<float>(buf[i]);

    // 2) IPP 浮点 RealFFT → CCS 格式 (输出 n+2 个元素)
    std::vector<float> ccs(n + 2);
    s_fft.FFT(time.data(), ccs.data());

    // 3) CCS 前 n 个值 → int32 buf
    //    CCS[0]=DC_re, CCS[1]=0(DC_im), CCS[2..n-1]=bin1..bin(n/2-1) 实/虚交错
    //    spectrum_process 会将 buf[0,1] 置零, 然后从 idx=0 开始按实虚对遍历
    for (int i = 0; i < n; ++i)
        buf[i] = static_cast<int32_t>(ccs[i]);

    return 1;
}

// ── float → int16 (饱和) ──
static inline int16_t FloatToInt16(float v) {
    if (v >  1.0f) return  32767;
    if (v < -1.0f) return -32768;
    return static_cast<int16_t>(v * 32767.0f);
}

static inline uint8_t DbfsToLevel(float dbfs) {
    if (dbfs <= FFT_DBFS_FLOOR) return 0;
    if (dbfs >= 0.0f) return SPECTRUM_MAX_LEVEL;
    float level = (dbfs - FFT_DBFS_FLOOR) * SPECTRUM_MAX_LEVEL / (0.0f - FFT_DBFS_FLOOR);
    return static_cast<uint8_t>(level + 0.5f);
}

// ============================================================
//  miniaudio 数据回调 — 累积采样 → 满帧后执行 FFT
// ============================================================
extern "C" void MaCaptureCallback(ma_device* pDevice, void* pOutput,
                                  const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;

    // 全双工：pInput 来自麦克风，pOutput 送扬声器
    auto* input  = static_cast<const float*>(pInput);
    auto* output = static_cast<float*>(pOutput);

    if (pInput == nullptr || pOutput == nullptr) return;

    // 直通：麦克风输入 → 扬声器输出 (监听)
    std::memcpy(output, input, frameCount * sizeof(float));

    // 累积 float 采样，满一帧后执行频谱分析
    for (ma_uint32 i = 0; i < frameCount; ++i) {
        gFloatRing[gFloatPos++] = input[i];
        if (gFloatPos >= kSptLen) {
            gFloatPos = 0;

            // float → int16
            for (uint16_t j = 0; j < kSptLen; ++j) {
                gPcmBuffer[j] = FloatToInt16(gFloatRing[j]);
            }

            // 频谱处理 + 电平量化
            if (spectrum_process(&gSpectrumCtx, gPcmBuffer)) {
                // 将量化电平拷贝到渲染缓冲区 (release 语义保证数据可见)
                std::memcpy(gLevels, gSpectrumCtx.level_cur, kBinNum);
#if SPECTRUM_ENABLE_DEBUG_MEL
                spectrum_update_mel_dbfs(&gSpectrumCtx);
                for (uint8_t band = 0; band < kMelBandNum; ++band) {
                    gMelLevels[band] = DbfsToLevel(gSpectrumCtx.mel_dbfs[band]);
                }
#endif
                gNewData.store(true, std::memory_order_release);
            }
        }
    }
}

// ============================================================
//  绘制 81-bin 频谱图
// ============================================================
static void DrawSpectrumBars() {
    const int   canvasX    = 50;
    const int   canvasY    = 30;
    const int   canvasW    = kWindowWidth - canvasX - 20;
    const int   canvasH    = kWindowHeight - canvasY - 60;
    const float barSpacing = 2.0f;
    const float barTotalW  = static_cast<float>(canvasW) / static_cast<float>(kBinNum);
    const float barW       = barTotalW - barSpacing;

    // ── 背景 ──
    DrawRectangle(canvasX, canvasY, canvasW, canvasH, kBarBg);

    // ── 水平参考线 ──
    constexpr int kRefTicks[] = { 0, 23, 46, 70, 93, SPECTRUM_MAX_LEVEL };
    for (int ref : kRefTicks) {
        int y = canvasY + canvasH - (ref * canvasH / SPECTRUM_MAX_LEVEL);
        DrawLine(canvasX, y, canvasX + canvasW, y, Color{ 50, 50, 50, 255 });
    }

    // ── 消费新数据 (acquire) ──
    (void)gNewData.exchange(false, std::memory_order_acquire);

    // ── 绘制每个 bin ──
    for (uint8_t i = 0; i < kBinNum; ++i) {
        float x   = canvasX + i * barTotalW;
        int   lev = gLevels[i];
        int   barH = (lev * canvasH) / SPECTRUM_MAX_LEVEL;
        int   y    = canvasY + canvasH - barH;

        // 颜色: 低→绿, 中→黄, 高→红
        Color color;
        if      (lev < 42)  color = kBarLow;
        else if (lev < 84)  color = kBarMid;
        else                color = kBarHigh;

        DrawRectangle(static_cast<int>(x), y,
                      static_cast<int>(barW), barH, color);
    }

    // ── 频率标签 (底部) ──
    for (uint8_t i = 0; i < kBinNum; ++i) {
        float cx = canvasX + i * barTotalW + barTotalW * 0.5f;
        char  label[16];
        const int freqHz = (static_cast<int>(i) + 1) * SPECTRUM_SAMPLE_RATE / SPECTRUM_LEN;
        if (freqHz >= 1000)
            snprintf(label, sizeof(label), "%.1fk", freqHz / 1000.0f);
        else
            snprintf(label, sizeof(label), "%d", freqHz);
        int tw = MeasureText(label, 10);
        DrawText(label, static_cast<int>(cx) - tw / 2,
                 canvasY + canvasH + 4, 10, GRAY);
    }

    // ── Y 轴刻度 ──
    for (int ref : kRefTicks) {
        int y   = canvasY + canvasH - (ref * canvasH / SPECTRUM_MAX_LEVEL);
        float db = FFT_DBFS_FLOOR + (ref * (0.0f - FFT_DBFS_FLOOR) / SPECTRUM_MAX_LEVEL);
        char lbl[16];
        snprintf(lbl, sizeof(lbl), "%.0f dB", db);
        DrawText(lbl, 2, y - 5, 10, GRAY);
    }
}

// ============================================================
//  绘制 32-band Mel滤波器组输出
// ============================================================
static void DrawMelBands() {
    const int   canvasX    = 50;
    const int   canvasY    = 30;
    const int   canvasW    = kWindowWidth - canvasX - 20;
    const int   canvasH    = kWindowHeight - canvasY - 60;
    const float barSpacing = 4.0f;
    const float barTotalW  = static_cast<float>(canvasW) / static_cast<float>(kMelBandNum);
    const float barW       = barTotalW - barSpacing;

    DrawRectangle(canvasX, canvasY, canvasW, canvasH, kBarBg);

    constexpr int kRefTicks[] = { 0, 23, 46, 70, 93, SPECTRUM_MAX_LEVEL };
    for (int ref : kRefTicks) {
        int y = canvasY + canvasH - (ref * canvasH / SPECTRUM_MAX_LEVEL);
        DrawLine(canvasX, y, canvasX + canvasW, y, Color{ 50, 50, 50, 255 });
    }

    (void)gNewData.exchange(false, std::memory_order_acquire);

    for (uint8_t i = 0; i < kMelBandNum; ++i) {
        float x = canvasX + i * barTotalW;
        int lev = gMelLevels[i];
        int barH = (lev * canvasH) / SPECTRUM_MAX_LEVEL;
        int y = canvasY + canvasH - barH;

        Color color;
        if      (lev < 42)  color = kBarLow;
        else if (lev < 84)  color = kBarMid;
        else                color = kBarHigh;

        DrawRectangle(static_cast<int>(x), y,
                      static_cast<int>(barW), barH, color);
    }

    for (uint8_t i = 0; i < kMelBandNum; ++i) {
        float cx = canvasX + i * barTotalW + barTotalW * 0.5f;
        char label[16];
        const int freqHz = static_cast<int>(gSpectrumCtx.mel_center_hz[i] + 0.5f);
        if (freqHz >= 1000)
            snprintf(label, sizeof(label), "%.1fk", freqHz / 1000.0f);
        else
            snprintf(label, sizeof(label), "%d", freqHz);
        int tw = MeasureText(label, 10);
        DrawText(label, static_cast<int>(cx) - tw / 2,
                 canvasY + canvasH + 4, 10, GRAY);
    }

    for (int ref : kRefTicks) {
        int y = canvasY + canvasH - (ref * canvasH / SPECTRUM_MAX_LEVEL);
        float db = FFT_DBFS_FLOOR + (ref * (0.0f - FFT_DBFS_FLOOR) / SPECTRUM_MAX_LEVEL);
        char lbl[16];
        snprintf(lbl, sizeof(lbl), "%.0f dB", db);
        DrawText(lbl, 2, y - 5, 10, GRAY);
    }
}

// ============================================================
//  主函数
// ============================================================
int main(void) {
    // ── 初始化频谱驱动 ──
    spectrum_init(&gSpectrumCtx);

    // ── raylib 窗口 ──
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, "Spectrogram - miniaudio + qwqdsp + raylib");
    SetTargetFPS(60);

    // ── miniaudio 全双工 (捕获 + 播放) ──
    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format   = ma_format_f32;
    config.capture.channels = 1;
    config.playback.format  = ma_format_f32;
    config.playback.channels = 1;
    config.sampleRate       = static_cast<ma_uint32>(kSampleRate);
    config.dataCallback     = MaCaptureCallback;
    config.pUserData        = nullptr;
    config.periodSizeInMilliseconds = 10;

    ma_device device;
    ma_result result = ma_device_init(nullptr, &config, &device);
    if (result == MA_SUCCESS) {
        ma_device_start(&device);
    }
    else {
        TraceLog(LOG_WARNING, "miniaudio 捕获设备初始化失败，以静默模式运行");
    }

    // ── 主循环 ──
    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_SPACE)) {
            gViewMode = (gViewMode == SpectrumViewMode::Bins) ? SpectrumViewMode::Mel : SpectrumViewMode::Bins;
        }

        BeginDrawing();
        ClearBackground(BLACK);
        if (gViewMode == SpectrumViewMode::Mel) {
            DrawMelBands();
        }
        else {
            DrawSpectrumBars();
        }
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
