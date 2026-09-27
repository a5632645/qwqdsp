#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <numbers>

#include "miniaudio.h"
#include "raylib.h"
#include "slider.hpp"

#include <qwqdsp/filter/one_pole_tpt.hpp>
#include <qwqdsp/fx/delay_line.hpp>

static constexpr int kWindowWidth = 480;
static constexpr int kWindowHeight = 240;
static constexpr float kSampleRate = 48000.0f;
static constexpr float kMaxDelayMs = 1000.0f;

using DelayLine = qwqdsp_fx::DelayLine<qwqdsp_fx::DelayLineInterp::None>;

// ------------------------------------------------------------
//  DelayTap — 带线性交叉淡化的延迟读取
// ------------------------------------------------------------
class DelayTap {
public:
    void SetTarget(float samples) noexcept {
        pending_ = samples;
    }

    float Read(DelayLine& dl) noexcept {
        float a = dl.GetAfterPush(static_cast<int>(prev_));
        float b = dl.GetAfterPush(static_cast<int>(target_));
        float out = a + (b - a) * gradient_;

        gradient_ += kDt;
        if (gradient_ >= 1.0f) {
            gradient_ -= 1.0f;
            prev_ = target_;
            target_ = pending_;
        }
        return out;
    }
private:
    float prev_ = 0.0f;
    float target_ = 0.0f;
    float pending_ = 0.0f;
    float gradient_ = 1.0f;
    static constexpr float kDt = 1.0f / 2048.0f;
};
// ------------------------------------------------------------
//  Echo — 立体声回声效果器
// ------------------------------------------------------------
class Echo {
public:
    void Init(float sample_rate) {
        sr_ = sample_rate;
        delay_l_.Init(kMaxDelayMs, sr_);
        delay_r_.Init(kMaxDelayMs, sr_);
    }

    void Process(const float* input, float* output_l, float* output_r, size_t frame_count) noexcept {
        // 仅参数变更时重新加载、更新滤波器和 Tap
        if (params_dirty_.exchange(false, std::memory_order_relaxed)) {
            float predelay = predelay_ms.load(std::memory_order_relaxed);
            float predelay_rr = predelay_r_ratio.load(std::memory_order_relaxed);
            float echo_time = echo_time_ms.load(std::memory_order_relaxed);
            float echo_time_rr = echo_time_r_ratio.load(std::memory_order_relaxed);
            float repeat_val = repeat.load(std::memory_order_relaxed);
            float cutoff_val = cutoff.load(std::memory_order_relaxed);
            float hpf_cutoff_val = hpf_cutoff.load(std::memory_order_relaxed);

            float predelay_l = predelay;
            float predelay_r = predelay * (1.0f + predelay_rr);
            float echo_time_l = echo_time;
            float echo_time_r = echo_time * (1.0f + echo_time_rr);

            repeat_r_ = 0.0f;
            if (echo_time_l > 0.0f && echo_time_r > 0.0f) [[likely]] {
                repeat_r_ = std::pow(repeat_val, echo_time_r / echo_time_l);
            }

            float w = 2.0f * std::numbers::pi_v<float> * cutoff_val / sr_;
            float wh = 2.0f * std::numbers::pi_v<float> * hpf_cutoff_val / sr_;
            lpf_l_.MakeLowpass(w);
            lpf_r_.MakeLowpass(w);
            hpf_l_.MakeLowpass(wh);
            hpf_r_.MakeLowpass(wh);

            tap_pd_l_.SetTarget(predelay_l * sr_ / 1000.0f);
            tap_pd_r_.SetTarget(predelay_r * sr_ / 1000.0f);
            tap_echo_l_.SetTarget(echo_time_l * sr_ / 1000.0f);
            tap_echo_r_.SetTarget(echo_time_r * sr_ / 1000.0f);
        }

        float p_in = 0.0f, p_out_l = 0.0f, p_out_r = 0.0f;
        float repeat_val = repeat.load(std::memory_order_relaxed);

        for (size_t i = 0; i < frame_count; ++i) {
            float in = input[i];
            p_in = std::max(p_in, std::abs(in));

            delay_l_.Push(in + fb_l_);
            output_l[i] = tap_pd_l_.Read(delay_l_);

            float el = tap_echo_l_.Read(delay_l_);
            el = hpf_l_.TickHighpass(el);
            el = lpf_l_.TickLowpass(el);
            fb_l_ = el * repeat_val;

            delay_r_.Push(in + fb_r_);
            output_r[i] = tap_pd_r_.Read(delay_r_);

            float er = tap_echo_r_.Read(delay_r_);
            er = hpf_r_.TickHighpass(er);
            er = lpf_r_.TickLowpass(er);
            fb_r_ = er * repeat_r_;

            p_out_l = std::max(p_out_l, std::abs(output_l[i]));
            p_out_r = std::max(p_out_r, std::abs(output_r[i]));
        }

        peak_in.store(p_in, std::memory_order_relaxed);
        peak_out_l.store(p_out_l, std::memory_order_relaxed);
        peak_out_r.store(p_out_r, std::memory_order_relaxed);
    }

    std::atomic<float> predelay_ms{20.0f};
    std::atomic<float> predelay_r_ratio{0.25f};
    std::atomic<float> echo_time_ms{200.0f};
    std::atomic<float> echo_time_r_ratio{0.25f};
    std::atomic<float> repeat{0.5f};
    std::atomic<float> cutoff{8000.0f};
    std::atomic<float> hpf_cutoff{100.0f};

    std::atomic<float> peak_in{0.0f};
    std::atomic<float> peak_out_l{0.0f};
    std::atomic<float> peak_out_r{0.0f};

    // 参数变更标记
    std::atomic<bool> params_dirty_{true};
private:
    float sr_ = 48000.0f;
    float repeat_r_ = 0.0f;
    DelayLine delay_l_;
    DelayLine delay_r_;
    DelayTap tap_pd_l_;
    DelayTap tap_pd_r_;
    DelayTap tap_echo_l_;
    DelayTap tap_echo_r_;
    qwqdsp_filter::OnePoleTPT lpf_l_;
    qwqdsp_filter::OnePoleTPT lpf_r_;
    qwqdsp_filter::OnePoleTPT hpf_l_;
    qwqdsp_filter::OnePoleTPT hpf_r_;
    float fb_l_ = 0.0f;
    float fb_r_ = 0.0f;
};

static Echo s_echo;

// ------------------------------------------------------------
//  miniaudio 回调
// ------------------------------------------------------------
extern "C" void MaCallback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;

    auto* input = static_cast<const float*>(pInput);
    auto* output = static_cast<float*>(pOutput);

    if (input == nullptr || output == nullptr)
        return;

    constexpr size_t kChunk = 512;
    float buf_l[kChunk];
    float buf_r[kChunk];
    size_t remaining = static_cast<size_t>(frameCount);
    size_t offset = 0;

    while (remaining > 0) {
        size_t nf = std::min(remaining, kChunk);
        s_echo.Process(input + offset, buf_l, buf_r, nf);

        for (size_t i = 0; i < nf; ++i) {
            output[2 * (offset + i)]     = buf_l[i];
            output[2 * (offset + i) + 1] = buf_r[i];
        }

        offset += nf;
        remaining -= nf;
    }
}

// ------------------------------------------------------------
//  工具：创建旋钮
// ------------------------------------------------------------
static Knob MakeKnob(int& kx, int ky, int w, int h, const char* title, float vmin, float vmax, float vstep,
                     float vdef) {
    Knob knob;
    knob.set_bound(kx, ky, w, h)
        .set_title(title)
        .set_range(vmin, vmax, vstep, vdef)
        .set_value(vdef)
        .set_name_font_size(9)
        .set_number_font_size(9)
        .set_fore_color(Color{200, 200, 210, 255})
        .set_bg_color(Color{30, 30, 35, 255});
    kx += w + 14;
    return knob;
}

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------
int main(void) {
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, "Echo Effect — Stereo");
    SetTargetFPS(60);

    // ── 初始化 Echo ──
    s_echo.Init(kSampleRate);

    // ── miniaudio full-duplex（单声道入 → 立体声出）──
    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.playback.format = ma_format_f32;
    config.playback.channels = 2;
    config.sampleRate = (ma_uint32)kSampleRate;
    config.dataCallback = MaCallback;
    config.pUserData = nullptr;
    config.periodSizeInMilliseconds = 10;

    ma_device device;
    ma_result result = ma_device_init(nullptr, &config, &device);
    bool audio_ok = (result == MA_SUCCESS);
    if (audio_ok) {
        ma_device_start(&device);
    }

    // ── 旋钮布局 ──
    constexpr int kKnobW = 58, kKnobH = 60;
    constexpr int kGap4 = 18;
    int kx, ky;

    // 第 1 行：延迟参数
    kx = (kWindowWidth - (4 * kKnobW + 3 * kGap4)) / 2;
    ky = 42;

    Knob knob_pd = MakeKnob(kx, ky, kKnobW, kKnobH, "Predelay", 0, 250, 1, 20);
    Knob knob_pd_ro = MakeKnob(kx, ky, kKnobW, kKnobH, "R Pre %", 0.0f, 50.0f, 1.0f, 25.0f);
    Knob knob_et = MakeKnob(kx, ky, kKnobW, kKnobH, "EchoTime", 0, 500, 5, 200);
    Knob knob_et_ro = MakeKnob(kx, ky, kKnobW, kKnobH, "R Echo %", -50.0f, 50.0f, 1.0f, 25.0f);

    // 第 2 行：反馈与滤波
    kx = (kWindowWidth - (3 * kKnobW + 2 * 20)) / 2;
    ky = 42 + kKnobH + 24;

    Knob knob_rpt = MakeKnob(kx, ky, kKnobW, kKnobH, "Repeat", 0.0f, 0.9f, 0.01f, 0.5f);
    Knob knob_cut = MakeKnob(kx, ky, kKnobW, kKnobH, "Cutoff", 200.0f, 20000.0f, 100.0f, 8000.0f);
    Knob knob_hpf = MakeKnob(kx, ky, kKnobW, kKnobH, "Cutoff HPF", 20.0f, 500.0f, 10.0f, 100.0f);

    // ── 旋钮回调绑定 ──
    auto bind_knob = [](Knob& knob, std::atomic<float>& target, auto fmt, float scale = 1.0f) {
        knob.value_to_text_function = fmt;
        knob.on_value_change = [&target, scale](float v) {
            target.store(v / scale, std::memory_order_relaxed);
            s_echo.params_dirty_.store(true, std::memory_order_relaxed);
        };
    };

    bind_knob(knob_pd, s_echo.predelay_ms, [](float v) { return TextFormat("%.0f ms", v); });
    bind_knob(knob_pd_ro, s_echo.predelay_r_ratio, [](float v) { return TextFormat("+%.0f%%", v); }, 100.0f);
    bind_knob(knob_et, s_echo.echo_time_ms, [](float v) { return TextFormat("%.0f ms", v); });
    bind_knob(knob_et_ro, s_echo.echo_time_r_ratio, [](float v) { return TextFormat("%+.0f%%", v); }, 100.0f);
    bind_knob(knob_rpt, s_echo.repeat, [](float v) { return TextFormat("%.0f%%", v * 100.0f); });
    bind_knob(knob_cut, s_echo.cutoff, [](float v) { return TextFormat("%.0f Hz", v); });
    bind_knob(knob_hpf, s_echo.hpf_cutoff, [](float v) { return TextFormat("%.0f Hz", v); });

    // 强制同步初始值
    s_echo.predelay_ms.store(knob_pd.get_value(), std::memory_order_relaxed);
    s_echo.predelay_r_ratio.store(knob_pd_ro.get_value() / 100.0f, std::memory_order_relaxed);
    s_echo.echo_time_ms.store(knob_et.get_value(), std::memory_order_relaxed);
    s_echo.echo_time_r_ratio.store(knob_et_ro.get_value() / 100.0f, std::memory_order_relaxed);
    s_echo.repeat.store(knob_rpt.get_value(), std::memory_order_relaxed);
    s_echo.cutoff.store(knob_cut.get_value(), std::memory_order_relaxed);
    s_echo.hpf_cutoff.store(knob_hpf.get_value(), std::memory_order_relaxed);

    // 触发首次参数加载
    s_echo.params_dirty_.store(true, std::memory_order_relaxed);

    // ── 峰值电平显示状态（带 hold 衰减）──
    float disp_peak_in = 0.0f, disp_peak_out_l = 0.0f, disp_peak_out_r = 0.0f;

    // ── main loop ──
    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(Color{10, 10, 12, 255});

        // ── 旋钮 ──
        knob_pd.display();
        knob_pd_ro.display();
        knob_et.display();
        knob_et_ro.display();
        knob_rpt.display();
        knob_cut.display();
        knob_hpf.display();

        // ── 峰值电平表 ──
        {
            constexpr float kDecay = 0.88f;
            constexpr float kBarW = 18.0f;
            constexpr float kBarX1 = 14.0f;
            constexpr float kBarX2 = 414.0f;
            constexpr float kBarX3 = 434.0f;
            constexpr float kBarY = 14.0f;
            constexpr float kBarH = 212.0f;

            // 更新显示峰值（瞬间上升，缓慢衰减）
            auto update_peak = [](float& disp, std::atomic<float>& src) {
                float raw = src.load(std::memory_order_relaxed);
                disp = (raw > disp) ? raw : disp * kDecay;
            };
            update_peak(disp_peak_in, s_echo.peak_in);
            update_peak(disp_peak_out_l, s_echo.peak_out_l);
            update_peak(disp_peak_out_r, s_echo.peak_out_r);

            // 将线性幅度转为 dB 归一化 [0,1]（-60dB ~ 0dB）
            auto to_norm = [](float amp) -> float {
                float db = 20.0f * std::log10(std::max(amp, 1e-10f));
                return std::clamp((db + 60.0f) / 60.0f, 0.0f, 1.0f);
            };

            // 根据归一化值取颜色（绿→黄→红）
            auto bar_color = [](float norm) -> Color {
                if (norm < 0.667f) {
                    float t = norm / 0.667f;
                    return Color{(unsigned char)(t * 255.0f), 255, 0, 255};
                }
                else {
                    float t = (norm - 0.667f) / 0.333f;
                    return Color{255, (unsigned char)((1.0f - t) * 255.0f), 0, 255};
                }
            };

            // 绘制单个电平条
            auto draw_bar = [](float x, float w, float y, float h, float norm, const char* label, Color col) {
                // 背景
                DrawRectangle(static_cast<int>(x), static_cast<int>(y), static_cast<int>(w), static_cast<int>(h),
                              Color{25, 25, 30, 255});
                // 边框
                DrawRectangleLines(static_cast<int>(x), static_cast<int>(y), static_cast<int>(w), static_cast<int>(h),
                                   Color{60, 60, 70, 255});
                // 填充（从底部向上）
                if (norm > 0.001f) {
                    int fill_h = static_cast<int>(norm * h);
                    DrawRectangle(static_cast<int>(x) + 1, static_cast<int>(y + h - fill_h), static_cast<int>(w) - 2,
                                  fill_h, col);
                }
                // 标签
                DrawText(label, static_cast<int>(x), static_cast<int>(y + h + 4), 9, Color{140, 140, 150, 255});
            };

            float n_in = to_norm(disp_peak_in);
            float n_outl = to_norm(disp_peak_out_l);
            float n_outr = to_norm(disp_peak_out_r);

            draw_bar(kBarX1, kBarW, kBarY, kBarH, n_in, "IN", bar_color(n_in));
            draw_bar(kBarX2, kBarW, kBarY, kBarH, n_outl, "L", bar_color(n_outl));
            draw_bar(kBarX3, kBarW, kBarY, kBarH, n_outr, "R", bar_color(n_outr));
        }

        EndDrawing();
    }

    // ── cleanup ──
    if (audio_ok) {
        ma_device_stop(&device);
        ma_device_uninit(&device);
    }
    CloseWindow();
    return 0;
}
