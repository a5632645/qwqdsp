// ============================================================
// ADAA 振荡器 / waveshaper —— 实时试听 + 波形/频谱对照
// ============================================================
// 两种模式（顶部按钮或 Tab 切换）：
//   Oscillator ：抗混叠振荡器（锯齿 / 三角 / 梯形方波），相位增量 = f0/fs
//   Waveshaper ：正弦波 → 静态非线性（硬削波 / 波折叠 / tanh / 多项式）
// 每种模式都能在 8 档算法间切换（数字键 1..8），实时听差别、看频谱：
//   trivial / AA-FIR-1(=DPW-2) / AA-FIR-2(=DPW-3) / AA-IIR-1 / AA-IIR-2
//   / OVS-FIR-8（FIR 多相解析过采样）/ OVS-IIR-8（IIR 全通和多相解析过采样）/ polyBLEP
// 最后两档只在 Waveshaper 模式可用（解析过采样是「信号 → 非线性」结构）；polyBLEP 只对振荡器
// 有意义。频谱面板上灰色是 trivial 的谱，彩色是当前方法的谱 —— 混叠分量直接看得出来。
//
// 「多项式」shape 走三次多项式 f(x)=x+d·x²+d·x³：其余算法用它的实数退化，两个解析过采样档
// 则把 H(z)=z+z²+z³ 作用在复数解析信号上取实部（Vicanek 归一化 Re(H(d·z))/d，见 adaa_dsp.hpp）。
// d 由 drive 旋钮单独映射：d = drive/20（drive 默认 6 → d = 0.3；drive 拉满 20 → d = 1，即
// labs/adaa_iir 的参照点）。
//
// DSP 在本目录的 `adaa_dsp.hpp`（与 labs/adaa_iir 的 Python 参考实现逐点比对过）。
// polyBLEP 用库组件 `qwqdsp_oscillator::PolyBlep`。
//
// 音频输出用 raylib 的 AudioStream（与 synth/polyblep.cpp 同路子，故目标加 NO_MINIAUDIO）。

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <span>

#include "raylib.h"
#include "slider.hpp"

#include "qwqdsp/convert.hpp"
#include "qwqdsp/oscillator/polyblep.hpp"
#include "qwqdsp/spectral/real_fft_adv.hpp"
#include "qwqdsp/window/blackman_harris.hpp"

#include "adaa_dsp.hpp"

namespace {

// ------------------------------------------------------------
// 常量与布局
// ------------------------------------------------------------

constexpr int kWidth = 1100;
constexpr int kHeight = 630;
constexpr float kFs = 48000.0f;

constexpr int kModeRowY = 0;
constexpr int kModeRowH = 26;

constexpr int kLeftW = 250;          ///< 左侧控件列宽度
constexpr int kPanelX = kLeftW + 10;
constexpr int kPanelW = kWidth - kPanelX - 10;

constexpr int kScopeY = 34;
constexpr int kScopeH = 290;
constexpr int kSpecY = kScopeY + kScopeH + 8;
constexpr int kSpecH = kHeight - kSpecY - 26;

constexpr int kMethodBtnY = 34;
constexpr int kMethodBtnH = 24;
constexpr int kNumMethods = 8;       ///< 5 个 ADAA 档 + 2 个解析过采样档 + polyBLEP

// 方法索引（0..4 = adaa::Method；解析过采样与 polyBLEP 不在 adaa::Method 里）
constexpr int kMethodFirOvs = 5;     ///< FIR 多相解析过采样（L=8）
constexpr int kMethodIirOvs = 6;     ///< IIR 全通和多相解析过采样（L=8）
constexpr int kMethodPolyBlep = 7;   ///< 库组件 polyBLEP（仅振荡器）
constexpr int kMethodLastAdaa = 4;   ///< adaa::Method::AaIIR2

constexpr int kWaveBtnY = kMethodBtnY + kNumMethods * (kMethodBtnH + 2) + 8;
constexpr int kWaveBtnH = 24;

constexpr int kKnobY = kWaveBtnY + 4 * (kWaveBtnH + 2) + 10;   ///< 4 个 shape 按钮

constexpr int kScopeCap = 4096;      ///< 示波器环形缓冲
constexpr int kFftSize = 4096;
constexpr int kNumBins = kFftSize / 2 + 1;
constexpr float kSpecFloorDb = -100.0f;
constexpr float kSpecTopDb = 10.0f;

constexpr float kSpecMinHz = 20.0f;
constexpr float kSpecMaxHz = 20000.0f;

enum class Mode { Oscillator = 0, Waveshaper, NumModes };
enum class OscWave { Saw = 0, Triangle, Square, NumWaves };
enum class ShapeKind { HardClip = 0, Wavefold, Tanh, Poly, NumShapes };
constexpr int kNumPwlShapes = 3;     ///< 前 3 个 shape 是分段线性；Poly 是三次多项式

constexpr const char* kModeNames[] = {"Oscillator", "Waveshaper"};
constexpr const char* kMethodLabels[] = {
    "trivial", "AA-FIR-1 (DPW-2)", "AA-FIR-2 (DPW-3)",
    "AA-IIR-1 (butter2)", "AA-IIR-2 (ellip-10)",
    "OVS-FIR-8 (analytic)", "OVS-IIR-8 (analytic)", "polyBLEP*",
};
constexpr const char* kOscWaveLabels[] = {"saw", "triangle", "square(trap)"};
constexpr const char* kShapeLabels[] = {
    "hard clip", "wavefold", "tanh(0.3)", "poly z+z^2+z^3 (analytic)",
};

/// 「多项式」shape 的旋钮 → d 映射：高次谐波 ∝ d^(n−1)，d 大时波形幅度很快超出示波器量程，
/// 故取 drive/20（drive 默认 6 → d = 0.3；drive 拉满 20 → d = 1，即 Python 参照点）。
constexpr float kPolyDriveScale = 1.0f / 20.0f;

[[nodiscard]] float PolyDrive(float drive) noexcept { return drive * kPolyDriveScale; }

// ------------------------------------------------------------
// 显示用的无锁环形缓冲（音频线程写，UI 线程读）
// ------------------------------------------------------------

template <size_t N>
class ScopeBuffer {
public:
    void Push(float v) noexcept {
        buf_[w_.load(std::memory_order_relaxed) % N] = v;
        w_.fetch_add(1, std::memory_order_relaxed);
    }

    /// 把最近 ``n`` 个样本拷进 ``dst``（n ≤ N）。
    void CopyLatest(float* dst, size_t n) const noexcept {
        const size_t w = w_.load(std::memory_order_relaxed);
        const size_t start = w - n;
        for (size_t i = 0; i < n; ++i) {
            dst[i] = buf_[(start + i) % N];
        }
    }

    [[nodiscard]] uint64_t Written() const noexcept { return w_.load(std::memory_order_relaxed); }

private:
    std::array<float, N> buf_{};
    std::atomic<uint64_t> w_{0};
};

// ------------------------------------------------------------
// 全局参数（UI 线程写、音频线程读；演示里用放松的原子即可）
// ------------------------------------------------------------

std::atomic<Mode> g_mode{Mode::Oscillator};
std::atomic<int> g_method{static_cast<int>(adaa::Method::AaIIR2)};
std::atomic<int> g_wave{static_cast<int>(OscWave::Saw)};
std::atomic<int> g_shape{static_cast<int>(ShapeKind::HardClip)};
std::atomic<float> g_freq{440.0f};
std::atomic<float> g_drive{6.0f};      ///< waveshaper 的输入增益
std::atomic<float> g_level{0.5f};
std::atomic<float> g_reset_req{0.0f};  ///< 换算法时请求复位

// ------------------------------------------------------------
// DSP 对象（只在音频线程里访问）
// ------------------------------------------------------------

adaa::Oscillator g_osc;
adaa::Shaper g_shaper;
adaa::PeriodicWave g_waves[static_cast<int>(OscWave::NumWaves)];
adaa::Pwl g_shapes[kNumPwlShapes];
adaa::Cubic g_cubic;             ///< 「多项式」shape（每 buffer 由 drive 更新）
adaa::AnalyticOvsFir g_ovs_fir;  ///< FIR 多相解析过采样
adaa::AnalyticOvsIir g_ovs_iir;  ///< IIR 全通和多相解析过采样
adaa::Filter g_butter = adaa::MakeButter2();
adaa::Filter g_iir2 = adaa::MakeElliptic10(0.49, 80.0, 0.1);   // 10 阶椭圆，rs=80dB
qwqdsp_oscillator::PolyBlep<qwqdsp_oscillator::blep_coeff::BlackmanNutallApprox> g_polyblep;

// 显示缓冲：当前方法 / trivial 参考
ScopeBuffer<kScopeCap> g_scope_cur;
ScopeBuffer<kScopeCap> g_scope_naive;

// 频率/相位状态
double g_sine_phase{};
int g_active_method = -1;    ///< 音频线程记录当前已配置的方法（换算法才重配 DSP）

// ------------------------------------------------------------
// 音频线程
// ------------------------------------------------------------

void PushScope(float cur, float naive) noexcept {
    g_scope_cur.Push(cur);
    g_scope_naive.Push(naive);
}

void AudioOutCallback(void* buffer, unsigned int frames) noexcept {
    struct T { float l; float r; };
    std::span<T> out{static_cast<T*>(buffer), frames};

    const Mode mode = g_mode.load(std::memory_order_relaxed);
    const int method = g_method.load(std::memory_order_relaxed);
    const float f0 = g_freq.load(std::memory_order_relaxed);
    const float drive = g_drive.load(std::memory_order_relaxed);
    const float level = g_level.load(std::memory_order_relaxed);

    if (g_reset_req.exchange(0.0f, std::memory_order_relaxed) != 0.0f) {
        g_osc.Reset();
        g_shaper.Reset();
        g_ovs_fir.Reset();
        g_ovs_iir.Reset();
        g_active_method = -1;                  // 换模式/换算法后必须重配另一个对象
    }

    const auto& filter =
        (method == static_cast<int>(adaa::Method::AaIIR2)) ? g_iir2 : g_butter;

    if (mode == Mode::Oscillator) {
        g_osc.SetWave(&g_waves[static_cast<size_t>(g_wave.load(std::memory_order_relaxed))]);
        // ⚠ 只有换算法时才重配：SetMethod → SetFilter 会复位 ŷ 状态，每个 buffer 都调会
        //   在块边界留下瞬态（相位不动、状态被清零），听起来是 ~fps 的嗡声。
        if (method != g_active_method) {
            if (method != kMethodPolyBlep) {
                // 解析过采样档对振荡器没意义 → 回退到 trivial（仍走 Oscillator 以推进相位）
                const auto m = (method <= kMethodLastAdaa)
                    ? static_cast<adaa::Method>(method) : adaa::Method::Trivial;
                g_osc.SetMethod(m, filter);
            }
            g_active_method = method;
        }
        g_osc.SetDelta(static_cast<float>(f0) / kFs);
        g_polyblep.SetFreq(f0, kFs);
        const int wave = g_wave.load(std::memory_order_relaxed);

        for (auto& s : out) {
            float y = 0.0f;
            if (method == kMethodPolyBlep) {
                switch (static_cast<OscWave>(wave)) {
                    case OscWave::Saw: y = g_polyblep.Sawtooth(); break;
                    case OscWave::Triangle: y = g_polyblep.Triangle(); break;
                    default: y = g_polyblep.Sqaure(); break;
                }
            }
            else {
                y = g_osc.Process();
            }
            const float naive = g_waves[static_cast<size_t>(wave)].Value(g_osc.Phase());
            PushScope(y, naive);
            const float o = y * level;
            s.l = o;
            s.r = o;
        }
    }
    else {
        const int shape = g_shape.load(std::memory_order_relaxed);
        const bool poly = (shape == static_cast<int>(ShapeKind::Poly));
        const bool analytic = (method == kMethodFirOvs || method == kMethodIirOvs);
        const float d = PolyDrive(drive);
        g_cubic.b1 = 1.0f;
        g_cubic.b2 = d;
        g_cubic.b3 = d;
        if (poly) {
            g_shaper.SetCubic(&g_cubic);
        }
        else {
            g_shaper.SetShape(&g_shapes[static_cast<size_t>(shape)]);
        }
        const adaa::Pwl* const pwl = poly ? nullptr : &g_shapes[static_cast<size_t>(shape)];
        if (analytic) {
            if (poly) {
                g_ovs_fir.SetPolyDrive(d);
                g_ovs_iir.SetPolyDrive(d);
            }
            else {
                g_ovs_fir.SetPwlShape(pwl);
                g_ovs_iir.SetPwlShape(pwl);
            }
        }
        if (method != g_active_method) {
            if (method <= kMethodLastAdaa) {
                g_shaper.SetMethod(static_cast<adaa::Method>(method), filter);
            }
            g_active_method = method;
        }
        const double dp = static_cast<double>(f0) / kFs;

        for (auto& s : out) {
            const float sin_v = static_cast<float>(
                std::sin(2.0 * 3.14159265358979323846 * g_sine_phase));
            g_sine_phase += dp;
            if (g_sine_phase > 1.0) { g_sine_phase -= 1.0; }
            // 「多项式」shape 的自变量取单位正弦（drive 只经 d 进入多项式）；其余 shape 沿用
            // drive·sin 的既有做法。
            const float xin = poly ? sin_v : drive * sin_v;

            float y = 0.0f;
            if (analytic) {
                y = (method == kMethodFirOvs) ? g_ovs_fir.Process(xin)
                                              : g_ovs_iir.Process(xin);
            }
            else if (method == kMethodPolyBlep) {
                // polyBLEP 档对 waveshaper 没有意义 → 回退到 trivial
                y = poly ? g_cubic.Value(xin) : pwl->Value(xin);
            }
            else {
                y = g_shaper.Process(xin);
            }
            const float naive = poly ? g_cubic.Value(xin) : pwl->Value(xin);
            PushScope(y, naive);
            const float o = y * level;
            s.l = o;
            s.r = o;
        }
    }
}

// ------------------------------------------------------------
// 面板辅助
// ------------------------------------------------------------

float FreqToX(float hz) {
    const float t = (std::log10(hz) - std::log10(kSpecMinHz))
                  / (std::log10(kSpecMaxHz) - std::log10(kSpecMinHz));
    return static_cast<float>(kPanelX) + t * static_cast<float>(kPanelW);
}

float DbToY(float db) {
    const float t = (db - kSpecFloorDb) / (kSpecTopDb - kSpecFloorDb);
    return static_cast<float>(kSpecY + kSpecH) - t * static_cast<float>(kSpecH);
}

/// 画一个按钮，返回是否被点击。
bool DrawButton(Rectangle r, const char* text, bool active, bool enabled = true) {
    const bool hovered = CheckCollisionPointRec(GetMousePosition(), r);
    const bool clicked = enabled && hovered && IsMouseButtonPressed(MOUSE_LEFT_BUTTON);
    const Color bg = active ? RAYWHITE : (enabled ? Color{40, 40, 40, 255} : Color{24, 24, 24, 255});
    const Color fg = active ? BLACK : (enabled ? RAYWHITE : Color{110, 110, 110, 255});
    DrawRectangleRec(r, bg);
    if (!active) {
        DrawRectangleLinesEx(r, 1.0f, enabled ? Color{120, 120, 120, 255} : Color{70, 70, 70, 255});
    }
    DrawText(text, static_cast<int>(r.x) + 6, static_cast<int>(r.y) + 5, 12, fg);
    return clicked;
}

// ------------------------------------------------------------
// 示波器
// ------------------------------------------------------------

std::array<float, kScopeCap> g_scope_tmp{};

void DrawScope(float f0) {
    DrawRectangleLines(kPanelX, kScopeY, kPanelW, kScopeH, Color{70, 70, 70, 255});

    // 自适应窗口：大约画 4 个周期（至少 64 点、至多 kScopeCap）
    const float period = kFs / std::max(1.0f, f0);
    const size_t n = static_cast<size_t>(
        std::clamp(period * 4.0f, 64.0f, static_cast<float>(kScopeCap)));
    if (g_scope_cur.Written() < n) { return; }

    g_scope_cur.CopyLatest(g_scope_tmp.data(), n);

    // 中轴线
    const float mid_y = static_cast<float>(kScopeY) + kScopeH * 0.5f;
    DrawLineV({static_cast<float>(kPanelX), mid_y},
              {static_cast<float>(kPanelX + kPanelW), mid_y}, Color{50, 50, 50, 255});

    const float sx = static_cast<float>(kPanelW) / static_cast<float>(n - 1);
    const float scale = kScopeH * 0.42f;
    Vector2 prev{};
    for (size_t i = 0; i < n; ++i) {
        const Vector2 p{static_cast<float>(kPanelX) + static_cast<float>(i) * sx,
                        mid_y - std::clamp(g_scope_tmp[i], -1.5f, 1.5f) * scale};
        if (i > 0) {
            DrawLineV(prev, p, Color{120, 200, 255, 255});
        }
        // 逐样本点：削波这类「近垂直跳变」用折线画会在某些列留缝，点补上
        DrawPixelV(p, Color{160, 220, 255, 255});
        prev = p;
    }
    DrawText(std::format("scope: last {:.0f} periods ({:.2f} ms)", 4.0,
                         static_cast<double>(n) / kFs * 1000.0).c_str(),
             kPanelX + 8, kScopeY + 6, 12, Color{170, 170, 170, 255});
}

// ------------------------------------------------------------
// 频谱
// ------------------------------------------------------------

qwqdsp_spectral::RealFftAdv g_fft;
std::array<float, kFftSize> g_fft_in{};
std::array<float, kNumBins> g_gain_cur{};
std::array<float, kNumBins> g_gain_naive{};
float g_spec_peak_db = 0.0f;

void UpdateSpectrum() {
    if (g_scope_cur.Written() < kFftSize) { return; }
    g_scope_cur.CopyLatest(g_fft_in.data(), kFftSize);
    qwqdsp_window::BlackmanHarris::ApplyWindow(g_fft_in, true);
    g_fft.FFTGainPhase(g_fft_in, g_gain_cur);
    g_scope_naive.CopyLatest(g_fft_in.data(), kFftSize);
    qwqdsp_window::BlackmanHarris::ApplyWindow(g_fft_in, true);
    g_fft.FFTGainPhase(g_fft_in, g_gain_naive);

    // 用 trivial 的峰值当 0 dB 参考（两种方法共用同一参考，便于直接比较）
    float peak = 1e-12f;
    for (const float v : g_gain_naive) { peak = std::max(peak, v); }
    g_spec_peak_db = 20.0f * std::log10(peak);
}

void DrawSpectrumCurve(const std::array<float, kNumBins>& gain, Color color) {
    const float bin_hz = kFs / static_cast<float>(kFftSize);
    Vector2 prev{};
    bool has_prev = false;
    for (int plot_x = 0; plot_x < kPanelW; ++plot_x) {
        const float t = static_cast<float>(plot_x) / static_cast<float>(kPanelW - 1);
        const float hz = std::pow(10.0f, std::log10(kSpecMinHz)
                                  + t * (std::log10(kSpecMaxHz) - std::log10(kSpecMinHz)));
        const int bin = std::clamp(static_cast<int>(hz / bin_hz), 0, kNumBins - 1);
        const float db = 20.0f * std::log10(std::max(gain[bin], 1e-12f)) - g_spec_peak_db;
        const Vector2 p{static_cast<float>(kPanelX + plot_x), DbToY(std::clamp(db, kSpecFloorDb, kSpecTopDb))};
        if (has_prev) { DrawLineV(prev, p, color); }
        prev = p;
        has_prev = true;
    }
}

void DrawSpectrum() {
    DrawRectangleLines(kPanelX, kSpecY, kPanelW, kSpecH, Color{70, 70, 70, 255});

    // 网格：20/100/1k/10k Hz 与 0/-20/.../-100 dB
    for (const float hz : {100.0f, 1000.0f, 10000.0f}) {
        const float x = FreqToX(hz);
        DrawLineV({x, static_cast<float>(kSpecY)}, {x, static_cast<float>(kSpecY + kSpecH)},
                  Color{45, 45, 45, 255});
        DrawText(std::format("{}", static_cast<int>(hz)).c_str(),
                 static_cast<int>(x) + 2, kSpecY + kSpecH - 14, 10, Color{140, 140, 140, 255});
    }
    for (int db = 0; db >= -100; db -= 20) {
        const float y = DbToY(static_cast<float>(db));
        DrawLineV({static_cast<float>(kPanelX), y},
                  {static_cast<float>(kPanelX + kPanelW), y}, Color{45, 45, 45, 255});
        DrawText(std::format("{}", db).c_str(), kPanelX + 3, static_cast<int>(y) - 11, 10,
                 Color{140, 140, 140, 255});
    }

    DrawSpectrumCurve(g_gain_naive, Color{110, 110, 110, 255});
    DrawSpectrumCurve(g_gain_cur, Color{255, 190, 90, 255});

    DrawText("spectrum: 0 dB = trivial peak   gray = trivial   orange = current method",
             kPanelX + 8, kSpecY + 5, 12, Color{170, 170, 170, 255});
}

} // namespace

// ------------------------------------------------------------

int main() {
    // DSP 初始化
    g_waves[static_cast<int>(OscWave::Saw)] = adaa::MakeSaw();
    g_waves[static_cast<int>(OscWave::Triangle)] = adaa::MakeTriangle();
    g_waves[static_cast<int>(OscWave::Square)] = adaa::MakeSquare(128.0f);
    g_shapes[static_cast<int>(ShapeKind::HardClip)] = adaa::MakeHardClip(1.0f);
    g_shapes[static_cast<int>(ShapeKind::Wavefold)] = adaa::MakeWavefold(0.7f);
    g_shapes[static_cast<int>(ShapeKind::Tanh)] = adaa::MakeTanhFit(1.0f, 0.3f);
    g_fft.Init(kFftSize);

    InitWindow(kWidth, kHeight, "ADAA oscillator / waveshaper");
    SetTargetFPS(60);

    InitAudioDevice();
    const bool audio_ok = IsAudioDeviceReady();
    AudioStream stream{};
    if (audio_ok) {
        SetAudioStreamBufferSizeDefault(512);
        stream = LoadAudioStream(static_cast<unsigned int>(kFs), 32, 2);
        SetAudioStreamCallback(stream, AudioOutCallback);
        PlayAudioStream(stream);
    }
    else {
        TraceLog(LOG_WARNING, "音频设备不可用，仅显示（截图/无头模式下正常）");
    }

    // ---- 控件 ----
    Knob pitch_knob;
    pitch_knob.on_value_change = [](float v) {
        g_freq.store(qwqdsp::convert::Pitch2Freq(v));
        g_reset_req.store(1.0f);
    };
    pitch_knob.set_bound(0, kKnobY, 80, 70);
    pitch_knob.set_range(0.0f, 127.0f, 0.1f, 69.0f);   // 69 = A4 = 440 Hz
    pitch_knob.set_bg_color(BLACK);
    pitch_knob.set_fore_color(RAYWHITE);
    pitch_knob.set_title("pitch");
    pitch_knob.value_to_text_function = [](float v) {
        return std::format("{:.1f}\n{:.0f} Hz", static_cast<double>(v),
                           static_cast<double>(qwqdsp::convert::Pitch2Freq(v)));
    };

    Knob drive_knob;
    drive_knob.on_value_change = [](float v) { g_drive.store(v); };
    drive_knob.set_bound(85, kKnobY, 80, 70);
    drive_knob.set_range(0.5f, 20.0f, 0.1f, 6.0f);
    drive_knob.set_bg_color(BLACK);
    drive_knob.set_fore_color(RAYWHITE);
    drive_knob.set_title("drive");
    drive_knob.value_to_text_function = [](float v) {
        if (g_shape.load() == static_cast<int>(ShapeKind::Poly)) {
            return std::format("{:.1f}x\nd={:.2f}", static_cast<double>(v),
                               static_cast<double>(PolyDrive(v)));
        }
        return std::format("{:.1f}x", static_cast<double>(v));
    };

    Knob level_knob;
    level_knob.on_value_change = [](float v) { g_level.store(v); };
    level_knob.set_bound(170, kKnobY, 80, 70);
    level_knob.set_range(0.0f, 1.0f, 0.01f, 0.5f);
    level_knob.set_bg_color(BLACK);
    level_knob.set_fore_color(RAYWHITE);
    level_knob.set_title("level");
    level_knob.value_to_text_function = [](float v) {
        return std::format("{:.2f}", static_cast<double>(v));
    };

    while (!WindowShouldClose()) {
        const Mode mode = g_mode.load();
        const int method = g_method.load();

        // 快捷键
        if (IsKeyPressed(KEY_TAB)) {
            g_mode.store(mode == Mode::Oscillator ? Mode::Waveshaper : Mode::Oscillator);
            g_reset_req.store(1.0f);
        }
        for (int i = 0; i < kNumMethods; ++i) {
            const bool is_poly = (i == kMethodPolyBlep);
            const bool is_analytic = (i == kMethodFirOvs || i == kMethodIirOvs);
            if (IsKeyPressed(KEY_ONE + i)
                && !(mode == Mode::Waveshaper && is_poly)
                && !(mode == Mode::Oscillator && is_analytic)) {
                g_method.store(i);
                g_reset_req.store(1.0f);
            }
        }
        // 参数变了 → 音频线程侧复位（避免状态不匹配的爆音）
        if (static_cast<Mode>(g_mode.load()) != mode) { g_reset_req.store(1.0f); }

        BeginDrawing();
        {
            ClearBackground(BLACK);

            // ---- 顶部模式按钮 ----
            const float half = static_cast<float>(kWidth) * 0.5f;
            if (DrawButton({0, kModeRowY, half - 1, kModeRowH}, kModeNames[0],
                           mode == Mode::Oscillator)) {
                g_mode.store(Mode::Oscillator);
                g_reset_req.store(1.0f);
            }
            if (DrawButton({half + 1, kModeRowY, half - 1, kModeRowH}, kModeNames[1],
                           mode == Mode::Waveshaper)) {
                g_mode.store(Mode::Waveshaper);
                g_reset_req.store(1.0f);
            }

            // ---- 方法按钮 ----
            for (int i = 0; i < kNumMethods; ++i) {
                const bool is_poly = (i == kMethodPolyBlep);
                const bool is_analytic = (i == kMethodFirOvs || i == kMethodIirOvs);
                const bool enabled = !(is_poly && mode == Mode::Waveshaper)
                                  && !(is_analytic && mode == Mode::Oscillator);
                const Rectangle r{0, static_cast<float>(kMethodBtnY + i * (kMethodBtnH + 2)),
                                  static_cast<float>(kLeftW), static_cast<float>(kMethodBtnH)};
                if (DrawButton(r, kMethodLabels[i], method == i, enabled)) {
                    g_method.store(i);
                    g_reset_req.store(1.0f);
                }
            }

            // ---- 波形 / 非线性按钮 ----
            const char* const* labels = (mode == Mode::Oscillator) ? kOscWaveLabels : kShapeLabels;
            const int count = (mode == Mode::Oscillator) ? static_cast<int>(OscWave::NumWaves)
                                                         : static_cast<int>(ShapeKind::NumShapes);
            const int selected = (mode == Mode::Oscillator) ? g_wave.load() : g_shape.load();
            for (int i = 0; i < count; ++i) {
                const Rectangle r{0, static_cast<float>(kWaveBtnY + i * (kWaveBtnH + 2)),
                                  static_cast<float>(kLeftW), static_cast<float>(kWaveBtnH)};
                if (DrawButton(r, labels[i], selected == i)) {
                    if (mode == Mode::Oscillator) { g_wave.store(i); }
                    else { g_shape.store(i); }
                }
            }

            // ---- 旋钮 ----
            pitch_knob.display();
            drive_knob.SetEnable(mode == Mode::Waveshaper);
            drive_knob.display();
            level_knob.display();

            // ---- 面板 ----
            DrawScope(g_freq.load());
            UpdateSpectrum();
            DrawSpectrum();

            // ---- 角落信息 ----
            const float f0 = g_freq.load();
            const char* note = nullptr;
            if (method == kMethodPolyBlep) {
                note = (mode == Mode::Oscillator)
                    ? "polyBLEP"
                    : "polyBLEP n/a for waveshaper (falls back to trivial)";
            }
            else if (g_shape.load() == static_cast<int>(ShapeKind::Poly)
                     && (method == kMethodFirOvs || method == kMethodIirOvs)) {
                note = "analytic OVS: complex cubic shaper, d = drive/20";
            }
            else {
                note = kMethodLabels[method];
            }
            DrawText(std::format("{}   fs={:.0f} Hz   f0={:.1f} Hz   [Tab] mode   [1..8] method",
                                 note, static_cast<double>(kFs), static_cast<double>(f0)).c_str(),
                     kPanelX + 8, kHeight - 15, 12, Color{170, 170, 170, 255});
            DrawFPS(kWidth - 80, kModeRowH + 6);
            if (!audio_ok) {
                DrawText("audio device unavailable", kPanelX + 8, kScopeY + kScopeH - 18, 12, RED);
            }
        }
        EndDrawing();
    }

    if (audio_ok) {
        if (IsAudioStreamPlaying(stream)) { StopAudioStream(stream); }   // shot 模式没建流
        if (IsAudioStreamValid(stream)) { UnloadAudioStream(stream); }
        CloseAudioDevice();
    }
    CloseWindow();
    return 0;
}
