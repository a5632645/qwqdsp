#include <algorithm>
#include <cmath>
#include <cstdio>
#include <format>
#include <span>
#include <string>
#include <utility>

#include "miniaudio.h"
#include "raylib.h"

#include "reassignment/colormap_adapter.hpp"
#include "reassignment/freq_reassignment_frame.hpp"
#include "reassignment/scrolling_image.hpp"
#include "reassignment/spectrogram_column.hpp"
#include "reassignment/spectrogram_frame.hpp"
#include "reassignment/tf_derivative_reassignment_frame.hpp"
#include "reassignment/tf_derivative_reassignment_frame_conv.hpp"
#include "reassignment/tf_phase_vocoder_reassignment_frame.hpp"
#include "reassignment/tf_phase_vocoder_reassignment_frame_conv.hpp"
#include "reassignment/nc_reassignment_frame.hpp"
#include "reassignment/nc_time_reassignment_frame.hpp"
#include "reassignment/windowless_nc_frame.hpp"
#include "reassignment/windowless_nc_reassign_frame.hpp"
#include "reassignment/tf_reassignment_frame.hpp"
#include "reassignment/time_reassignment_frame.hpp"
#include <qwqdsp/colormap/colormap.hpp>

// ── 音频参数 ──
static constexpr int kSampleRate = 48000;
static constexpr int kFftSize = 4096;
static constexpr int kHopSize = kFftSize / 16;

// ── zeroPad ──
// 简单类（标准谱图/时间重分配）：较大 zeroPad 提升频率分辨率
static constexpr int kZeroPadSimple = 2;
// 全 TF 类（子列缓冲已提供时间精度）
static constexpr int kZeroPadFull = 1;
// ── NC 方法 ──
static constexpr int kNcZeroPad = 2;
// ── Windowless NC：NC bin 带宽缩放系数（>1 带宽更宽/分辨率更低，<1 相反）──
static constexpr float kNcBandwidthScale = 1.0f;
// ── Windowless NC + TF 重分配：窗长下限（周期数）。低频 NC 窗长被 0.075 s 上限夹到
//    20–30 Hz 只剩 1.5–2.25 个周期，实测（labs/nc_reassign）加到 4 个周期最优 ──
static constexpr float kNcMinPeriodsFloor = 4.0f;
// ── Windowless NC + TF 重分配：重分配截止频率。高频 NC bin 自身的时频分辨率已经很强，
//    搬到瞬时频率/群延迟反而更差（lab 分带线上能量：不重分配在 1–12 kHz 已 = 0.998，
//    与全程重分配持平；20–100 Hz 只有 0.418）。截止以上退化为 bin 中心 + 窗中心，
//    并与无重分配的 WindowlessNcFrame 逐像素一致 ──
static constexpr float kNcReassignMaxHz = 1000.0f;

// ── 频谱图显示 ──
static constexpr float kDbFloor = -72.0f;
static constexpr float kWindowLessNcDbFloor = -72.0f;
static constexpr float kFreqMin = 20.0f;
static constexpr float kFreqMax = 20000.0f;

// ── 滚动图像 ──
static constexpr float kScrollSeconds = 3.0f;
static constexpr int kImageWidth = static_cast<int>(kScrollSeconds * kSampleRate / kHopSize);

// ── UI 布局 ──
// 顶部 3 行参数文字（右上角留出 FPS）；画布；底部 4×4 算法切换面板。
static constexpr int kWindowWidth = 800;
static constexpr int kWindowHeight = 470;
static constexpr int kCanvasX = 60;
static constexpr int kCanvasY = 44;
static constexpr int kCanvasW = kWindowWidth - kCanvasX - 20;
static constexpr int kSelectorCols = 4;
static constexpr int kSelectorRows = 4;
static constexpr int kSelectorRowH = 24;
static constexpr int kSelectorGap = 6;        // 画布与面板之间
static constexpr int kBottomMargin = 14;      // 面板下方留白
static constexpr int kCanvasH =
    kWindowHeight - kCanvasY - kSelectorGap - kSelectorRows * kSelectorRowH - kBottomMargin;
// FPS 计数器（raylib 的 DrawFPS 用 20 号字，比参数文字高一倍，故单占右上角）
static constexpr int kFpsX = kWindowWidth - 80;
static constexpr int kFpsY = 4;

// ── UI 颜色与标签 ──
static constexpr float kFreqTickLabels[] = {20, 200, 2000, 20000};
static const Color kGridColor = {50, 50, 50, 255};
static const Color kTextColor = {180, 180, 180, 255};
static const Color kBgColor = {20, 20, 20, 255};

using ColorMap = ColormapAdapter<qwqdsp_colormap::Magma>;

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
    kNcMethod,
    kNcTimeMethod,
    kWindowlessNc,
    kWindowlessNcTf,
    kWindowlessNcFreq,
    kWindowlessNcTime,
    kCount
};

static FrameType g_frame_type = FrameType::kWindowlessNc;

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
    "NC Method",
    "NC Time",
    "Windowless NC",
    "Windowless NC TF",
    "Windowless NC Freq",
    "Windowless NC Time",
};
static_assert(std::size(kFrameNames) == static_cast<int>(FrameType::kCount));

static SpectrogramFrame<ColorMap> f_sp;
static FreqReassignmentFrame<ColorMap> f_freq;
static TimeReassignmentFrame<ColorMap> f_time;
static TfReassignmentFrame<ColorMap> f_tf;
static TfPhaseVocoderReassignmentFrame<ColorMap, false> f_pv;
static TfPhaseVocoderReassignmentFrame<ColorMap, true> f_pv_pk;
static TfDerivativeReassignmentFrame<ColorMap, false> f_deriv;
static TfDerivativeReassignmentFrame<ColorMap, true> f_deriv_pk;
static TfPhaseVocoderReassignmentFrameConv<ColorMap> f_pv_conv;
static TfDerivativeReassignmentFrameConv<ColorMap> f_deriv_conv;
static NcReassignmentFrame<ColorMap> f_nc;
static NcTimeReassignmentFrame<ColorMap> f_nc_time;
static WindowlessNcFrame<ColorMap> f_windowless;
static WindowlessNcReassignFrame<ColorMap> f_windowless_tf;
static WindowlessNcReassignFrame<ColorMap, NcReassignMode::kFreq> f_windowless_nc_freq;
static WindowlessNcReassignFrame<ColorMap, NcReassignMode::kTime> f_windowless_nc_time;

static SpectrogramColumn column_;
static ScrollingImage image_;

// ----------------------------------------
// 各算法真实使用的参数
// ----------------------------------------

/**
 * @brief 一个算法实际使用的参数
 *
 * 与 main() 里各帧 Init 的实参一一对应（改 Init 时同步改这张表）。
 */
struct FrameParams {
    int zero_pad;          // 零填充倍数；0 = 该算法不做 FFT
    bool rect_window;      // true = 矩形窗（NC 类直接对原始帧做 FFT，不加窗）
    bool windowless;       // true = 无窗递归滑动 DFT（不用 FFT/窗/零填充）
    int sub_columns;       // 时间轴子列数 = fftSize / hopSize；0 = 不做时间重分配
    float min_periods;     // 窗长周期数下限（仅无窗 NC 法）；0 = 无下限
    float reassign_max_hz; // 重分配截止频率(Hz，仅无窗 NC 重分配法)；0 = 不重分配
    float db_floor;        // 颜色映射的幅度下限
    const char* estimator; // 频率/时间估计算子
};

static constexpr FrameParams kFrameParams[] = {
    // kSpectrogram
    {kZeroPadSimple, false, false, 0, 0.0f, 0.0f, kDbFloor,
     "freq: nearest bin (low) + max-hold (high)"},
    // kFreqReassignment
    {kZeroPadFull, false, false, 0, 0.0f, 0.0f, kDbFloor,
     "freq: cross-spectrum phase diff of x[n-1]"},
    // kTimeReassignment
    {kZeroPadSimple, false, false, kFftSize / kHopSize, 0.0f, 0.0f, kDbFloor,
     "time: group delay of roll(X_h,1)"},
    // kTfReassignment
    {kZeroPadFull, false, false, kFftSize / kHopSize, 0.0f, 0.0f, kDbFloor,
     "freq: cross-spectrum phase diff | time: group delay of roll(X_h,1)"},
    // kTfPhaseVocoder
    {kZeroPadFull, false, false, kFftSize / kHopSize, 0.0f, 0.0f, kDbFloor,
     "freq: phase vocoder dphi | time: time-weighted window"},
    // kTfPhaseVocoderPeak
    {kZeroPadFull, false, false, kFftSize / kHopSize, 0.0f, 0.0f, kDbFloor,
     "freq: phase vocoder dphi | time: time-weighted window | + Loris sidelobe reject"},
    // kTfDerivative
    {kZeroPadFull, false, false, kFftSize / kHopSize, 0.0f, 0.0f, kDbFloor,
     "freq: derivative window | time: time-weighted window"},
    // kTfDerivativePeak
    {kZeroPadFull, false, false, kFftSize / kHopSize, 0.0f, 0.0f, kDbFloor,
     "freq: derivative window | time: time-weighted window | + Loris sidelobe reject"},
    // kTfPhaseVocoderConv
    {kZeroPadFull, false, false, kFftSize / kHopSize, 0.0f, 0.0f, kDbFloor,
     "freq: phase vocoder dphi | time: time-weighted window | + conv weight (drop 0.3-0.7)"},
    // kTfDerivativeConv
    {kZeroPadFull, false, false, kFftSize / kHopSize, 0.0f, 0.0f, kDbFloor,
     "freq: derivative window | time: time-weighted window | + conv weight (drop 0.3-0.7)"},
    // kNcMethod
    {kNcZeroPad, true, false, 0, 0.0f, 0.0f, kWindowLessNcDbFloor,
     "NC: negative correlation of the bin pair"},
    // kNcTimeMethod
    {kNcZeroPad, true, false, kFftSize / kHopSize, 0.0f, 0.0f, kWindowLessNcDbFloor,
     "NC: negative correlation of the bin pair | time: group delay of roll(X,1)"},
    // kWindowlessNc
    {0, false, true, 0, kNcMinPeriodsFloor, 0.0f, kWindowLessNcDbFloor,
     "NC: sliding DFT per bin | time axis: window center (no reassignment)"},
    // kWindowlessNcTf
    {0, false, true, 0, kNcMinPeriodsFloor, kNcReassignMaxHz, kWindowLessNcDbFloor,
     "+ freq/time reassignment"},
    // kWindowlessNcFreq
    {0, false, true, 0, kNcMinPeriodsFloor, kNcReassignMaxHz, kWindowLessNcDbFloor,
     "+ freq reassignment"},
    // kWindowlessNcTime
    {0, false, true, 0, kNcMinPeriodsFloor, kNcReassignMaxHz, kWindowLessNcDbFloor,
     "+ time reassignment"},
};
static_assert(std::size(kFrameParams) == static_cast<int>(FrameType::kCount));
static_assert(static_cast<int>(FrameType::kCount) == kSelectorCols * kSelectorRows);

/**
 * @brief 无窗 NC 各 bin 实际使用的窗长范围（样本）
 *
 * @return {最小窗长, 最大窗长}；非无窗方法返回 {0, 0}
 */
static std::pair<int, int> WindowlessWindowRange() noexcept {
    switch (g_frame_type) {
        case FrameType::kWindowlessNc:
            return {f_windowless.MinWindowSamples(), f_windowless.MaxWindowSamples()};
        case FrameType::kWindowlessNcTf:
            return {f_windowless_tf.MinWindowSamples(), f_windowless_tf.MaxWindowSamples()};
        case FrameType::kWindowlessNcFreq:
            return {f_windowless_nc_freq.MinWindowSamples(), f_windowless_nc_freq.MaxWindowSamples()};
        case FrameType::kWindowlessNcTime:
            return {f_windowless_nc_time.MinWindowSamples(), f_windowless_nc_time.MaxWindowSamples()};
        default:
            return {0, 0};
    }
}

// ----------------------------------------
// 算法切换面板
// ----------------------------------------

/**
 * @brief 第 index 个算法按钮在切换面板里的矩形
 *
 * @param index 算法下标（0 .. FrameType::kCount-1）
 */
static Rectangle SelectorItemRect(int index) noexcept {
    constexpr int kCellW = kCanvasW / kSelectorCols;
    const int col = index % kSelectorCols;
    const int row = index / kSelectorCols;
    return {
        static_cast<float>(kCanvasX + col * kCellW),
        static_cast<float>(kCanvasY + kCanvasH + kSelectorGap + row * kSelectorRowH),
        static_cast<float>(kCellW) - 2.0f,
        static_cast<float>(kSelectorRowH) - 2.0f,
    };
}

/**
 * @brief 画算法切换面板，并处理左键点击切换
 *
 * 交互与 audiofx/filters.cpp 的 drawSelector 一致：选中项反白（黑字）、未选中项
 * 描边、悬停提亮，左键点击即选中。算法有 16 个，单排放不下完整名字，故排成
 * 4 列 × 4 行。
 */
static void DrawAlgorithmSelector() noexcept {
    constexpr Color kSelectedFore = {255, 220, 60, 255};
    constexpr Color kHoverFore = {235, 235, 235, 255};
    constexpr Color kIdleFore = {120, 120, 120, 255};
    constexpr int kFontSize = 11;
    constexpr int kTextPadX = 6;
    constexpr int kTextPadY = 7;

    const Vector2 mouse = GetMousePosition();
    const bool clicked = IsMouseButtonPressed(MOUSE_LEFT_BUTTON);

    for (int i = 0; i < static_cast<int>(FrameType::kCount); ++i) {
        const Rectangle cell = SelectorItemRect(i);
        const bool hovered = CheckCollisionPointRec(mouse, cell);
        if (hovered && clicked)
            g_frame_type = static_cast<FrameType>(i);

        const bool active = (i == static_cast<int>(g_frame_type));
        const Color fore = active ? kSelectedFore : (hovered ? kHoverFore : kIdleFore);
        if (active)
            DrawRectangleRec(cell, fore);
        else
            DrawRectangleLinesEx(cell, 1.0f, fore);
        DrawText(kFrameNames[i], static_cast<int>(cell.x) + kTextPadX, static_cast<int>(cell.y) + kTextPadY,
                 kFontSize, active ? BLACK : fore);
    }
}

/**
 * @brief 画当前算法的参数信息（3 行）
 *
 * 第 1 行：变换/窗/步进/底噪；第 2 行：时间轴或 NC bin 几何；第 3 行：估计算子。
 * 取值全部来自 kFrameParams 与各帧自身（无窗 NC 的窗长范围由帧给出），
 * 避免显示与实际算法脱节。
 */
static void DrawFrameParams() {
    constexpr int kFontSize = 10;
    constexpr int kTextX = kCanvasX + 4;
    constexpr int kLineY0 = 4;
    constexpr int kLineY1 = 16;
    constexpr int kLineY2 = 28;

    const FrameParams& params = kFrameParams[static_cast<int>(g_frame_type)];

    std::string line1;
    if (params.windowless) {
        line1 = std::format(
            "FFT: none (recursive sliding DFT) | window: none | column every {} samples | floor: {:.0f} dB",
            kHopSize, params.db_floor);
    }
    else {
        const int fft_len = kFftSize * params.zero_pad;
        const double bin_hz = static_cast<double>(kSampleRate) / static_cast<double>(fft_len);
        const double overlap_pct =
            100.0 * (1.0 - static_cast<double>(kHopSize) / static_cast<double>(kFftSize));
        line1 = std::format(
            "FFT: {} x{} = {} (bin {:.2f} Hz) | window: {} | hop: {} ({:.1f}% overlap) | floor: {:.0f} dB",
            kFftSize, params.zero_pad, fft_len, bin_hz, params.rect_window ? "rect" : "BH3-3T",
            kHopSize, overlap_pct, params.db_floor);
    }

    std::string line2;
    if (params.windowless) {
        const auto [n_min, n_max] = WindowlessWindowRange();
        const double ms_per_sample = 1000.0 / static_cast<double>(kSampleRate);
        line2 = std::format("NC: {} bins (1 per pixel) | bandwidth x{:.1f} | N: {}-{} samples ({:.1f}-{:.1f} ms)",
                            kCanvasH, kNcBandwidthScale, n_min, n_max,
                            static_cast<double>(n_min) * ms_per_sample,
                            static_cast<double>(n_max) * ms_per_sample);
        if (params.min_periods > 0.0f)
            line2 += std::format(" | N >= {:.0f} periods", params.min_periods);
        if (params.reassign_max_hz > 0.0f)
            line2 += std::format(" | reassign < {:.0f} Hz", params.reassign_max_hz);
    }
    else {
        // NC 类的「bin 对」间距 = 零填充倍数（见 NcReassignmentFrame / NcTimeReassignmentFrame）
        if (params.rect_window)
            line2 = std::format("NC bin pair: rect FFT (k, k+{}) | ", params.zero_pad);
        if (params.sub_columns > 0)
            line2 += std::format("time axis: {} sub-columns (fftSize/hop = {}/{})", params.sub_columns, kFftSize,
                                 kHopSize);
        else
            line2 += "time axis: none (one column per hop)";
    }

    DrawText(line1.c_str(), kTextX, kLineY0, kFontSize, kTextColor);
    DrawText(line2.c_str(), kTextX, kLineY1, kFontSize, kTextColor);
    DrawText(params.estimator, kTextX, kLineY2, kFontSize, kTextColor);
}

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
        case FrameType::kNcMethod:
            column_.ProcessAudio({src, frameCount}, f_nc, push);
            break;
        case FrameType::kNcTimeMethod:
            column_.ProcessAudio({src, frameCount}, f_nc_time, push);
            break;
        case FrameType::kWindowlessNc:
            column_.ProcessAudio({src, frameCount}, f_windowless, push);
            break;
        case FrameType::kWindowlessNcTf:
            column_.ProcessAudio({src, frameCount}, f_windowless_tf, push);
            break;
        case FrameType::kWindowlessNcFreq:
            column_.ProcessAudio({src, frameCount}, f_windowless_nc_freq, push);
            break;
        case FrameType::kWindowlessNcTime:
            column_.ProcessAudio({src, frameCount}, f_windowless_nc_time, push);
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
    f_nc.Init(kSampleRate, kFftSize, kNcZeroPad, kCanvasH, kFreqMin, kFreqMax, kWindowLessNcDbFloor);
    f_nc_time.Init(kSampleRate, kFftSize, kHopSize, kNcZeroPad, kCanvasH, kFreqMin, kFreqMax, kWindowLessNcDbFloor);
    f_windowless.Init(kSampleRate, kFftSize, kHopSize, kNcZeroPad, kCanvasH, kFreqMin, kFreqMax, kWindowLessNcDbFloor,
                      kNcBandwidthScale);
    f_windowless_tf.Init(kSampleRate, kFftSize, kHopSize, kNcZeroPad, kCanvasH, kFreqMin, kFreqMax,
                         kWindowLessNcDbFloor, kNcBandwidthScale, kNcMinPeriodsFloor, kNcReassignMaxHz);
    f_windowless_nc_freq.Init(kSampleRate, kFftSize, kHopSize, kNcZeroPad, kCanvasH, kFreqMin, kFreqMax,
                              kWindowLessNcDbFloor, kNcBandwidthScale, kNcMinPeriodsFloor, kNcReassignMaxHz);
    f_windowless_nc_time.Init(kSampleRate, kFftSize, kHopSize, kNcZeroPad, kCanvasH, kFreqMin, kFreqMax,
                              kWindowLessNcDbFloor, kNcBandwidthScale, kNcMinPeriodsFloor, kNcReassignMaxHz);
    image_.Init(kImageWidth, kCanvasH);

    // ── 主循环 ──
    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(BLACK);

        DrawBackground();
        DrawSpectrogram();
        DrawFreqGrid();

        // ── 参数信息（按当前算法取真实参数）──
        DrawFrameParams();

        // ── 算法切换面板（左键点击切换）──
        DrawAlgorithmSelector();

        DrawFPS(kFpsX, kFpsY);
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
