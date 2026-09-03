#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <numbers>
#include <vector>

#include "raylib.h"
#include "miniaudio.h"

#include "qwqdsp/convert.hpp"
#include "qwqdsp/colormap/resonator.hpp"
#include "qwqdsp/filter/rbj.hpp"
#include "qwqdsp/filter/biquad.hpp"
#include "qwqdsp/oscillator/vic_sine_osc.hpp"
#include "simd.hpp"
#include <immintrin.h>

// ============================================================
//  常量 & 运行时参数
// ============================================================
// 每个谐振器的级联 biquad 节数 (可在 CMake 命令行添加 -DRESONATOR_BQ_COUNT=8 覆盖)
#ifndef RESONATOR_BQ_COUNT
#define RESONATOR_BQ_COUNT 3
#endif
static constexpr int    kResonatorBQ = RESONATOR_BQ_COUNT;
static constexpr int    kMaxBQ       = 8;  // 栈数组上限
static constexpr int    kWindowWidth      = 1100;
static constexpr int    kWindowHeight     = 640;
static constexpr float  kSampleRate       = 48000.0f;
static constexpr int    kSimdGroup        = 8;            // Float256 宽度
static constexpr int    kMaxResonators    = 1024;          // 最大谐振器数
static constexpr int    kMaxGroups        = kMaxResonators / kSimdGroup;
static constexpr int    kSpecCols         = 800;
static constexpr int    kSpecPadX         = 70;
static constexpr int    kSpecPadY         = 20;
static constexpr int    kSpecPadBottom    = 40;
static constexpr float  kMinDisplayDb     = -80.0f;
static constexpr float  kDisplayMaxFreq   = 16000.0f;
static constexpr float  kFreqMin          = 20.0f;

// 运行时可调参数 (主线程写，音频线程读)
static std::atomic<int>    s_param_num_res{256};       // 默认 256
static std::atomic<float>  s_param_min_bw{40.0f};
static std::atomic<size_t> s_param_hop_size{128};      // 默认 128
static std::atomic<int>    s_param_freq_dist{0};       // 0=log, 1=linear, 2=mel
static std::atomic<bool>   s_param_changed{false};
static std::atomic<bool>   s_pending_reinit{false};

// 实际生效的值 (音频线程使用)
static int    s_cur_num_groups  = 256 / 8;
static std::atomic<size_t> s_cur_hop_size{128};

// ============================================================
//  频率分段 & 谐振器状态 (仅在音频回调线程中访问)
// ============================================================
struct Band {
    float start_hz;
    float center_hz;
    float end_hz;
    float bw_hz;    // = end - start
};

static Band                 s_bands[kMaxResonators];
static size_t               s_samples_since_hop = 0;

// ============================================================
//  SIMD 谐振器 — Float256 同时处理 8 个谐振器
// ============================================================
struct Resonator8 {
    // VicSineOsc 系数 & 状态 (8 lanes)
    simd::Float256 k1_, k2_;
    simd::Float256 u_, v_;

    // 级联 biquad (实部 & 虚部)
    struct BQ8 {
        simd::Float256 b0_, b1_, b2_, a1_, a2_;
        simd::Float256 s1_, s2_;
    };
    BQ8 re_filt_[kResonatorBQ];
    BQ8 im_filt_[kResonatorBQ];

    // 输出 & 上一帧
    simd::Float256 re_, im_;
    simd::Float256 last_re_, last_im_;

    void Init(const Band bands[8], float fs) noexcept {
        alignas(32) float k1a[8], k2a[8];
        // 栈数组用 kMaxBQ 固定大小
        alignas(32) float b0a[kMaxBQ][8], b1a[kMaxBQ][8], b2a[kMaxBQ][8], a1a[kMaxBQ][8], a2a[kMaxBQ][8];

        constexpr int nbq = kResonatorBQ;

        // 计算 Butterworth Q 值 (2*nbq 阶)
        float bq_q[kMaxBQ];
        for (int j = 0; j < nbq; ++j) {
            float theta = (2.0f * j + 1.0f) * std::numbers::pi_v<float> / (4.0f * nbq);
            bq_q[j] = 1.0f / (2.0f * std::cos(theta));
            bq_q[j] = 0.5f;
        }

        for (int i = 0; i < 8; ++i) {
            float w = bands[i].center_hz / fs * (std::numbers::pi_v<float> * 2.0f);
            float t  = std::tan(w * 0.5f);
            k1a[i] = t;
            k2a[i] = 2.0f * t / (1.0f + t * t);

            float cutoff_w = (bands[i].bw_hz * 0.5f) / fs * (std::numbers::pi_v<float> * 2.0f);
            qwqdsp_filter::RBJ rbj;
            for (int j = 0; j < nbq; ++j) {
                rbj.Lowpass(cutoff_w, bq_q[j]);
                b0a[j][i] = rbj.b0;   b1a[j][i] = rbj.b1;   b2a[j][i] = rbj.b2;
                a1a[j][i] = rbj.a1;   a2a[j][i] = rbj.a2;
            }
        }

        k1_ = simd::Float256{k1a[0], k1a[1], k1a[2], k1a[3], k1a[4], k1a[5], k1a[6], k1a[7]};
        k2_ = simd::Float256{k2a[0], k2a[1], k2a[2], k2a[3], k2a[4], k2a[5], k2a[6], k2a[7]};
        u_ = simd::Float256{1,1,1,1,1,1,1,1};
        v_ = simd::Float256{0,0,0,0,0,0,0,0};

        for (int j = 0; j < nbq; ++j) {
            re_filt_[j].b0_ = simd::Float256{b0a[j][0], b0a[j][1], b0a[j][2], b0a[j][3], b0a[j][4], b0a[j][5], b0a[j][6], b0a[j][7]};
            re_filt_[j].b1_ = simd::Float256{b1a[j][0], b1a[j][1], b1a[j][2], b1a[j][3], b1a[j][4], b1a[j][5], b1a[j][6], b1a[j][7]};
            re_filt_[j].b2_ = simd::Float256{b2a[j][0], b2a[j][1], b2a[j][2], b2a[j][3], b2a[j][4], b2a[j][5], b2a[j][6], b2a[j][7]};
            re_filt_[j].a1_ = simd::Float256{a1a[j][0], a1a[j][1], a1a[j][2], a1a[j][3], a1a[j][4], a1a[j][5], a1a[j][6], a1a[j][7]};
            re_filt_[j].a2_ = simd::Float256{a2a[j][0], a2a[j][1], a2a[j][2], a2a[j][3], a2a[j][4], a2a[j][5], a2a[j][6], a2a[j][7]};
            re_filt_[j].s1_ = simd::Float256{0,0,0,0,0,0,0,0};
            re_filt_[j].s2_ = simd::Float256{0,0,0,0,0,0,0,0};
            im_filt_[j] = re_filt_[j];
        }

        re_ = im_ = simd::Float256{0,0,0,0,0,0,0,0};
        last_re_ = last_im_ = simd::Float256{0,0,0,0,0,0,0,0};
    }

    void Tick(float x) noexcept {
        simd::Float256 xv{x,x,x,x,x,x,x,x};

        // 8 × VicSineOsc
        simd::Float256 w = u_ - k1_ * v_;
        v_ = v_ + k2_ * w;
        u_ = w - k1_ * v_;

        // 8 × 复下混频: x * (u + j*v)
        simd::Float256 re = xv * u_;
        simd::Float256 im = xv * v_;

        // 级联 biquad (TDF2)
        for (int j = 0; j < kResonatorBQ; ++j) {
            auto& f = re_filt_[j];
            simd::Float256 out = re * f.b0_ + f.s1_;
            f.s1_ = re * f.b1_ - out * f.a1_ + f.s2_;
            f.s2_ = re * f.b2_ - out * f.a2_;
            re = out;
        }
        for (int j = 0; j < kResonatorBQ; ++j) {
            auto& f = im_filt_[j];
            simd::Float256 out = im * f.b0_ + f.s1_;
            f.s1_ = im * f.b1_ - out * f.a1_ + f.s2_;
            f.s2_ = im * f.b2_ - out * f.a2_;
            im = out;
        }

        last_re_ = re_;  last_im_ = im_;
        re_ = re;        im_ = im;
    }
};

static Resonator8 s_resonator_groups[kMaxGroups];

// 提取 8 个标量增益 & 重分配频率 (hop 时调用，标量开销可忽略)
static void ExtractGroup(int gi, float gains[8], float freq_offsets[8]) noexcept {
    auto& g = s_resonator_groups[gi];
    for (int i = 0; i < 8; ++i) {
        float re = g.re_[i], im = g.im_[i];
        float lre = g.last_re_[i], lim = g.last_im_[i];
        gains[i] = std::sqrt(re * re + im * im);
        std::complex<float> common{re, im};
        std::complex<float> freq{lre, lim};
        float a = std::arg(freq * std::conj(common));
        freq_offsets[i] = std::fmod(a / (std::numbers::pi_v<float> * 2.0f), 1.0f);
    }
}

// ============================================================
//  频谱图缓冲区 (音频线程写列，主线程读列)
// ============================================================
static int                 s_spec_rows = 0;
static std::vector<float>  s_spec_buf;             // [cols * rows]
static std::atomic<size_t> s_spec_col{0};           // 下一个要写入的列 (循环)

// 指数平滑状态 (每显示行一个，仅回调线程访问)
static std::vector<float>  s_smoothed;
static constexpr float     kSmoothAlpha = 0.3f;     // EMA 系数

// 对数频率坐标预计算
static float s_log_f_min  = 0.0f;
static float s_log_f_range = 0.0f;

static std::atomic<bool>   g_capturing{false};

// ── 频率 → 归一化坐标 (根据当前分布类型) ──
static float FreqToNorm(float freq_hz) noexcept {
    int dist = s_param_freq_dist.load(std::memory_order_relaxed);
    if (dist == 1) {  // linear
        return (freq_hz - kFreqMin) / (kDisplayMaxFreq - kFreqMin);
    } else if (dist == 2) {  // mel
        auto mel = [](float f) { return 1127.0f * std::log(1.0f + f / 700.0f); };
        return (mel(freq_hz) - mel(kFreqMin)) / (mel(kDisplayMaxFreq) - mel(kFreqMin));
    }
    // default: log
    return (std::log10(freq_hz) - s_log_f_min) / s_log_f_range;
}

// ── 重初始化谐振器组 (主线程调用) ──
static void ReinitializeResonators(int num_res, float min_bw, int dist) noexcept {
    // 须为 8 的倍数
    num_res = (num_res / 8) * 8;
    if (num_res < 8) num_res = 8;
    if (num_res > kMaxResonators) num_res = kMaxResonators;
    int num_groups = num_res / 8;

    // 更新频率坐标基准
    s_log_f_min   = std::log10(kFreqMin);
    s_log_f_range = std::log10(kDisplayMaxFreq) - s_log_f_min;

    // 按分布类型计算频段
    if (dist == 1) {  // linear
        float step = (kDisplayMaxFreq - kFreqMin) / num_res;
        for (int i = 0; i < num_res; ++i) {
            float lo = kFreqMin + i * step;
            float hi = kFreqMin + (i + 1) * step;
            s_bands[i].start_hz   = lo;
            s_bands[i].end_hz     = hi;
            s_bands[i].center_hz  = (lo + hi) * 0.5f;
            s_bands[i].bw_hz      = std::max(hi - lo, min_bw);
        }
    } else if (dist == 2) {  // mel
        auto mel = [](float f) { return 1127.0f * std::log(1.0f + f / 700.0f); };
        auto inv_mel = [](float m) { return 700.0f * (std::exp(m / 1127.0f) - 1.0f); };
        float m_min = mel(kFreqMin);
        float m_max = mel(kDisplayMaxFreq);
        float step  = (m_max - m_min) / num_res;
        for (int i = 0; i < num_res; ++i) {
            float lo = inv_mel(m_min + i * step);
            float hi = inv_mel(m_min + (i + 1) * step);
            s_bands[i].start_hz   = lo;
            s_bands[i].end_hz     = hi;
            s_bands[i].center_hz  = (lo + hi) * 0.5f;
            s_bands[i].bw_hz      = std::max(hi - lo, min_bw);
        }
    } else {  // log (default)
        double log_l = std::log(static_cast<double>(kFreqMin));
        double log_h = std::log(static_cast<double>(kDisplayMaxFreq));
        double step  = (log_h - log_l) / num_res;
        for (int i = 0; i < num_res; ++i) {
            float lo = static_cast<float>(std::exp(log_l + i * step));
            float hi = static_cast<float>(std::exp(log_l + (i + 1) * step));
            s_bands[i].start_hz   = lo;
            s_bands[i].end_hz     = hi;
            s_bands[i].center_hz  = (lo + hi) * 0.5f;
            s_bands[i].bw_hz      = std::max(hi - lo, min_bw);
        }
    }

    // 重新初始化谐振器组
    for (int g = 0; g < num_groups; ++g)
        s_resonator_groups[g].Init(s_bands + g * 8, kSampleRate);

    // 提交生效值
    s_cur_num_groups = num_groups;
    s_param_changed.store(false, std::memory_order_release);
}

static float Arg(std::complex<float> z) noexcept {
    float a = std::arg(z);
    a /= 2.0f * std::numbers::pi_v<float>;
    return std::fmod(a, 1.0f);
}

// ============================================================
//  miniaudio 数据回调 — 每采样 Tick 所有谐振器
// ============================================================
extern "C" void MaCaptureCallback(ma_device* pDevice, void* pOutput,
                                  const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;
    (void)pOutput;

    // loopback：pInput 来自系统音频输出
    auto* input  = static_cast<const float*>(pInput);

    if (pInput == nullptr) return;

    // ── 检查是否需要重初始化 (主线程设标志，音频线程执行) ──
    if (s_pending_reinit.exchange(false, std::memory_order_acquire)) {
        ReinitializeResonators(s_param_num_res.load(std::memory_order_relaxed),
                               s_param_min_bw.load(std::memory_order_relaxed),
                               s_param_freq_dist.load(std::memory_order_relaxed));
        s_cur_hop_size.store(s_param_hop_size.load(std::memory_order_relaxed),
                             std::memory_order_relaxed);
        // 清空平滑器 & 频谱缓存
        std::fill(s_smoothed.begin(), s_smoothed.end(), kMinDisplayDb);
        std::fill(s_spec_buf.begin(), s_spec_buf.end(), kMinDisplayDb);
        s_spec_col.store(0, std::memory_order_relaxed);
    }

    // ── 每采样 Tick 所有谐振器 ──
    int num_groups = s_cur_num_groups;
    for (ma_uint32 i = 0; i < frameCount; ++i) {
        float x = input[i];
        for (int g = 0; g < num_groups; ++g)
            s_resonator_groups[g].Tick(x);
    }

    s_samples_since_hop += frameCount;

    // ── 滚动降采样: 每 N 跳写一列, N=1 为正常速度 ──
    static int s_hop_counter = 0;
    static constexpr int kScrollDecimate = 1;

    // ── 每 hop_size 采样产出一列 ──
    size_t hop_size = s_cur_hop_size.load(std::memory_order_relaxed);
    while (s_samples_since_hop >= hop_size) {
        s_samples_since_hop -= hop_size;
        if (++s_hop_counter % kScrollDecimate != 0) continue;

        int rows = s_spec_rows;
        float col_data[1024];

        // ── 收集所有谐振器的 re, im ──
        alignas(32) float all_re[kMaxResonators], all_im[kMaxResonators];
        for (int g = 0; g < num_groups; ++g) {
            auto& grp = s_resonator_groups[g];
            for (int i = 0; i < 8; ++i) {
                all_re[g * 8 + i] = grp.re_[i];
                all_im[g * 8 + i] = grp.im_[i];
            }
        }
        // 延迟一阶：shifted_re[i] = all_re[i-1], shifted_re[0] = 0
        alignas(32) float shifted_re[kMaxResonators], shifted_im[kMaxResonators];
        shifted_re[0] = 1.0f;
        shifted_im[0] = 0.0f;
        for (int i = 1; i < kMaxResonators; ++i) {
            shifted_re[i] = all_re[i - 1];
            shifted_im[i] = all_im[i - 1];
        }

        // ── 1) 用重分配频率将各谐振器能量精确分派到显示行 ──
        std::fill(col_data, col_data + rows, kMinDisplayDb);
        alignas(32) float gains[8], freq_offsets[8];
        for (int g = 0; g < num_groups; ++g) {
            ExtractGroup(g, gains, freq_offsets);
            int base = g * 8;
            for (int i = 0; i < 8; ++i) {
                float inst_hz = s_bands[base + i].center_hz
                              + freq_offsets[i] * kSampleRate;
                if (inst_hz < kFreqMin) inst_hz = kFreqMin;
                if (inst_hz > kDisplayMaxFreq) continue;

                auto time = std::complex{all_re[base + i], all_im[base + i]} * std::complex{shifted_re[base + i], -shifted_im[base + i]};
                // 0~1
                auto timeT = 0.5f - Arg(time);
                // 时域重分配权值: |timeT| 越接近 0.5 能量越集中
                float time_weight = 1.0f - 2 * std::abs(timeT - 0.5f);

                float norm = FreqToNorm(inst_hz);
                int r = std::clamp(static_cast<int>(norm * rows), 0, rows - 1);
                float raw = qwqdsp::convert::Gain2Db<kMinDisplayDb>(gains[i]);
                if (raw > col_data[r]) col_data[r] = raw;
            }
        }

        // ── 2) 指数平滑 ──
        // for (int r = 0; r < rows; ++r) {
        //     s_smoothed[r] += kSmoothAlpha * (col_data[r] - s_smoothed[r]);
        //     col_data[r] = s_smoothed[r];
        // }

        // 写入频谱图环形缓冲区
        size_t idx = s_spec_col.load(std::memory_order_relaxed);
        std::copy_n(col_data, rows, s_spec_buf.begin() + idx * rows);
        s_spec_col.store((idx + 1) % kSpecCols, std::memory_order_release);
    }
}

// ============================================================
//  频谱图颜色映射  (Black → Blue → Cyan → Green → Yellow → Red)
//  已提取为 qwqdsp/colormap/resonator.hpp 中的 ResistorMap 表
// ============================================================
static Color MapDbToColor(float db) {
    float t = std::clamp((db - kMinDisplayDb) / (-kMinDisplayDb), 0.0f, 1.0f);
    auto rgb = qwqdsp_colormap::Resonator::Map(t);
    return Color{rgb[0], rgb[1], rgb[2], 255};
}

// ============================================================
//  主函数
// ============================================================
int main(void) {
    // ── raylib 窗口 ──
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, "Spectrogram - miniaudio + qwqdsp + raylib");
    SetTargetFPS(60);

    // ── 频谱图尺寸 ──
    const int spec_rows = kWindowHeight - kSpecPadY - kSpecPadBottom;
    const int spec_cols = kSpecCols;
    const int spec_x    = kSpecPadX;
    const int spec_y    = kSpecPadY;
    const int spec_w    = kWindowWidth - kSpecPadX;
    const int spec_h    = spec_rows;
    s_spec_rows = spec_rows;

    // ── 初始计算频段 & 谐振器 ──
    ReinitializeResonators(256, 40.0f, 0);
    s_cur_hop_size.store(128, std::memory_order_relaxed);
    s_param_hop_size.store(128, std::memory_order_relaxed);

    // 频谱图共享缓冲区 & 平滑器
    s_spec_buf.assign(spec_cols * spec_rows, kMinDisplayDb);
    s_smoothed.assign(spec_rows, kMinDisplayDb);

    // ── 显示用纹理 ──
    std::vector<Color> pixels(spec_cols * spec_rows, BLACK);
    Texture2D texture{};
    {
        Image img;
        img.data     = pixels.data();
        img.width    = spec_cols;
        img.height   = spec_rows;
        img.mipmaps  = 1;
        img.format   = PIXELFORMAT_UNCOMPRESSED_R8G8B8A8;
        texture = LoadTextureFromImage(img);
    }

    // ── miniaudio 回环捕获 (loopback) ──
    ma_device_config config = ma_device_config_init(ma_device_type_loopback);
    config.capture.format   = ma_format_f32;
    config.capture.channels = 1;
    config.sampleRate       = (ma_uint32)kSampleRate;
    config.dataCallback     = MaCaptureCallback;
    config.pUserData        = nullptr;
    config.periodSizeInMilliseconds = 10;

    ma_device device;
    ma_result result = ma_device_init(nullptr, &config, &device);
    if (result == MA_SUCCESS) {
        ma_device_start(&device);
        g_capturing.store(true, std::memory_order_relaxed);
    } else {
        TraceLog(LOG_WARNING, "miniaudio 捕获设备初始化失败，以静默模式运行");
    }

    // ── 主循环 (仅渲染，FFT 在音频线程中完成) ──
    while (!WindowShouldClose()) {
        // ── 读取音频线程产出的频谱列，更新像素缓冲区 ──
        {
            size_t base = s_spec_col.load(std::memory_order_acquire);

            for (int c = 0; c < spec_cols; ++c) {
                // 环形缓冲中第 c 列（从左到右）：base → 最左（最旧），base-1 → 最右（最新）
                size_t src = (base + c) % spec_cols;
                float* col = s_spec_buf.data() + src * spec_rows;

                for (int r = 0; r < spec_rows; ++r)
                    // 翻转 Y：r=0 (最低频) → 纹理底部，r=max (最高频) → 纹理顶部
                    pixels[c + (spec_rows - 1 - r) * spec_cols] = MapDbToColor(col[r]);
            }

            UpdateTexture(texture, pixels.data());
        }

        // ── 绘制 ──
        BeginDrawing();
        ClearBackground(Color{10, 10, 12, 255});

        DrawTexture(texture, spec_x, spec_y, WHITE);

        DrawRectangleLines(spec_x - 1, spec_y - 1,
                           spec_w + 2, spec_h + 2, Color{60, 60, 70, 255});

        // ── 频率轴标签 ──
        {
            constexpr float freq_ticks[] = { 31, 62, 125, 250, 500, 1000, 2000, 4000, 8000, 16000 };

            for (auto f : freq_ticks) {
                float norm = FreqToNorm(f);
                int y = spec_y + spec_rows - 1 - (int)(norm * spec_rows);
                y = std::clamp(y, spec_y, spec_y + spec_h - 1);

                DrawLine(spec_x - 4, y, spec_x, y, Color{80, 80, 90, 255});
                const char* label = (f >= 1000) ? TextFormat("%dk", f / 1000)
                                                : TextFormat("%d", (int)f);
                DrawText(label, spec_x - MeasureText(label, 10) - 8, y - 5, 10, LIGHTGRAY);
            }
        }

        // ── 颜色指示条 (dBFS) ──
        {
            int bar_x = spec_x + kSpecCols + 16;
            int bar_y = spec_y;
            int bar_w = 28;
            int bar_h = spec_h;

            // 标题
            DrawText("dBFS", bar_x, bar_y - 14, 10, LIGHTGRAY);

            // 渐变条
            for (int r = 0; r < bar_h; ++r) {
                float db = (float)r / (bar_h - 1) * kMinDisplayDb;  // top=0 → bottom=-80
                Color c = MapDbToColor(db);
                DrawRectangle(bar_x, bar_y + r, bar_w, 1, c);
            }
            DrawRectangleLines(bar_x, bar_y, bar_w, bar_h, Color{60, 60, 70, 255});

            // dB 刻度标签
            constexpr float kTickDb[] = { 0, -20, -40, -60, -80 };
            for (auto db : kTickDb) {
                float t = db / kMinDisplayDb;   // 0 → 1
                int y = bar_y + (int)(t * bar_h);
                const char* lbl = (db == 0) ? "0" : TextFormat("%d", (int)db);
                DrawText(lbl, bar_x + bar_w + 4, y - 5, 10, LIGHTGRAY);
                DrawLine(bar_x, y, bar_x - 3, y, Color{80, 80, 90, 255});
            }
        }

        // ── 时间轴标签 ──
        {
            int label_y = spec_y + spec_h + 6;
            DrawText("< older", spec_x, label_y, 10, DARKGRAY);
            DrawText("newer >", spec_x + spec_w - 40, label_y, 10, DARKGRAY);
        }

        DrawText("Spectrogram", spec_x, 4, 14, LIGHTGRAY);

        {
            static constexpr const char* kDistNames[] = { "Log", "Linear", "Mel" };
            int d = s_param_freq_dist.load(std::memory_order_relaxed);
            size_t hop = s_cur_hop_size.load();
            DrawText(TextFormat("[1,2,3]%s  [+/-]BW:%.0f  [<>]Resonators:%d  [[]]Hop:%zu  [SPACE] Puase  [ESC] Quit  FPS:%d",
                     kDistNames[d], s_param_min_bw.load(), s_cur_num_groups * 8, hop, GetFPS()),
                     spec_x + 150, 5, 10, WHITE);
        }

        if (g_capturing.load(std::memory_order_relaxed))
            DrawText("● REC", spec_x + spec_w - 80, 4, 12, GREEN);
        else
            DrawText("■ STOP", spec_x + spec_w - 60, 4, 12, RED);

        EndDrawing();

        // ── 参数热键 ──
        bool need_reinit = false;
        if (IsKeyPressed(KEY_SPACE)) {
            if (g_capturing.load()) {
                if (result == MA_SUCCESS) ma_device_stop(&device);
                g_capturing.store(false, std::memory_order_relaxed);
            } else {
                if (result == MA_SUCCESS) ma_device_start(&device);
                g_capturing.store(true, std::memory_order_relaxed);
            }
        }
        if (IsKeyPressed(KEY_ONE))        { s_param_freq_dist.store(0); need_reinit = true; }
        if (IsKeyPressed(KEY_TWO))        { s_param_freq_dist.store(1); need_reinit = true; }
        if (IsKeyPressed(KEY_THREE))      { s_param_freq_dist.store(2); need_reinit = true; }
        if (IsKeyPressed(KEY_EQUAL)) {  // +
            float bw = s_param_min_bw.load();
            bw = std::min(bw + 5.0f, 200.0f);
            s_param_min_bw.store(bw);
            need_reinit = true;
        }
        if (IsKeyPressed(KEY_MINUS)) {  // -
            float bw = s_param_min_bw.load();
            bw = std::max(bw - 5.0f, 5.0f);
            s_param_min_bw.store(bw);
            need_reinit = true;
        }
        if (IsKeyPressed(KEY_COMMA)) {  // < 减少谐振器
            int n = s_param_num_res.load();
            n = std::max(n - 32, 32);
            s_param_num_res.store(n);
            need_reinit = true;
        }
        if (IsKeyPressed(KEY_PERIOD)) {  // > 增加谐振器
            int n = s_param_num_res.load();
            n = std::min(n + 32, kMaxResonators);
            s_param_num_res.store(n);
            need_reinit = true;
        }

        if (need_reinit) {
            s_pending_reinit.store(true, std::memory_order_release);
        }
        // [`[`/`]` hop 变化直接写入原子，无需重初始化
        if (IsKeyPressed(KEY_LEFT_BRACKET)) {
            size_t h = s_param_hop_size.load();
            h = std::max(h - 32, size_t(32));
            s_param_hop_size.store(h);
            s_cur_hop_size.store(h, std::memory_order_relaxed);
        }
        if (IsKeyPressed(KEY_RIGHT_BRACKET)) {
            size_t h = s_param_hop_size.load();
            h = std::min(h + 32, size_t(2048));
            s_param_hop_size.store(h);
            s_cur_hop_size.store(h, std::memory_order_relaxed);
        }
    }

    // ── 清理 ──
    if (result == MA_SUCCESS) {
        ma_device_stop(&device);
        ma_device_uninit(&device);
    }
    UnloadTexture(texture);
    CloseWindow();
    return 0;
}
