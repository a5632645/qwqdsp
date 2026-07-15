#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <numbers>

#include "miniaudio.h"
#include "raylib.h"
#include "slider.hpp"

#include <qwqdsp/fx/plate_reverb.hpp>

static constexpr int kWindowWidth = 620;
static constexpr int kWindowHeight = 280;
static constexpr float kSampleRate = 48000.0f;

// ------------------------------------------------------------
//  PlateReverbWrapper — 带原子参数的板混响包装
// ------------------------------------------------------------
class PlateReverbWrapper {
public:
    void Init(float sample_rate) {
        sr_ = sample_rate;
        reverb_.Init(sample_rate);
    }

    void Process(const float* input, float* output_l, float* output_r, size_t frame_count) noexcept {
        if (params_dirty_.exchange(false, std::memory_order_relaxed)) {
            float m = mix_.load(std::memory_order_relaxed);
            float pd = predelay_.load(std::memory_order_relaxed);
            float lp = lowpass_.load(std::memory_order_relaxed);
            float d = decay_.load(std::memory_order_relaxed);
            float sz = size_.load(std::memory_order_relaxed);
            float damp = damping_.load(std::memory_order_relaxed);

            reverb_.SetMix(m);
            reverb_.SetPredelay(pd / 1000.0f);   // ms -> s
            reverb_.SetLowpass(lp);
            reverb_.SetDecay(d);
            reverb_.SetSize(sz);
            reverb_.SetDamping(damp);
        }

        float p_in_l = 0.0f, p_in_r = 0.0f, p_out_l = 0.0f, p_out_r = 0.0f;

        for (size_t i = 0; i < frame_count; ++i) {
            float in_l = input[2 * i];
            float in_r = input[2 * i + 1];
            p_in_l = std::max(p_in_l, std::abs(in_l));
            p_in_r = std::max(p_in_r, std::abs(in_r));

            float l, r;
            reverb_.Process(in_l, in_r, &l, &r);
            output_l[i] = l;
            output_r[i] = r;

            p_out_l = std::max(p_out_l, std::abs(l));
            p_out_r = std::max(p_out_r, std::abs(r));
        }

        peak_in_l.store(p_in_l, std::memory_order_relaxed);
        peak_in_r.store(p_in_r, std::memory_order_relaxed);
        peak_out_l.store(p_out_l, std::memory_order_relaxed);
        peak_out_r.store(p_out_r, std::memory_order_relaxed);
    }

    std::atomic<float> mix_{0.5f};
    std::atomic<float> predelay_{20.0f};   // ms
    std::atomic<float> lowpass_{20000.0f}; // Hz
    std::atomic<float> decay_{2000.0f}; // ms
    std::atomic<float> size_{1.0f};
    std::atomic<float> damping_{5000.0f};  // Hz

    std::atomic<float> peak_in_l{0.0f};
    std::atomic<float> peak_in_r{0.0f};
    std::atomic<float> peak_out_l{0.0f};
    std::atomic<float> peak_out_r{0.0f};

    std::atomic<bool> params_dirty_{true};
private:
    float sr_ = 48000.0f;
    qwqdsp_fx::PlateReverb reverb_;
};

static PlateReverbWrapper s_reverb;

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
        s_reverb.Process(input + offset, buf_l, buf_r, nf);

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
    InitWindow(kWindowWidth, kWindowHeight, "Plate Reverb — PlateReverb");
    SetTargetFPS(60);

    // ── 初始化 PlateReverb ──
    s_reverb.Init(kSampleRate);

    // ── miniaudio full-duplex（立体声入 → 立体声出）──
    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 2;
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

    // 第 1 行：混响主体参数
    kx = (kWindowWidth - (6 * kKnobW + 5 * kGap4)) / 2;
    ky = 42;

    Knob knob_mix   = MakeKnob(kx, ky, kKnobW, kKnobH, "Mix",      0.0f,    1.0f,     0.01f, 0.5f);
    Knob knob_pd    = MakeKnob(kx, ky, kKnobW, kKnobH, "Predelay", 0.0f,  100.0f,     1.0f,  20.0f);
    Knob knob_decay = MakeKnob(kx, ky, kKnobW, kKnobH, "Decay",  100.0f, 10000.0f, 50.0f, 2000.0f);
    Knob knob_size  = MakeKnob(kx, ky, kKnobW, kKnobH, "Size",     0.5f,    2.0f,     0.01f, 1.0f);
    Knob knob_lp    = MakeKnob(kx, ky, kKnobW, kKnobH, "Lowpass", 16.0f, 20000.0f, 100.0f,  20000.0f);
    Knob knob_damp  = MakeKnob(kx, ky, kKnobW, kKnobH, "Damping", 16.0f, 20000.0f, 100.0f,  5000.0f);

    // ── 旋钮回调绑定 ──
    auto bind_knob = [](Knob& knob, std::atomic<float>& target, auto fmt, float scale = 1.0f) {
        knob.value_to_text_function = fmt;
        knob.on_value_change = [&target, scale](float v) {
            target.store(v / scale, std::memory_order_relaxed);
            s_reverb.params_dirty_.store(true, std::memory_order_relaxed);
        };
    };

    bind_knob(knob_mix,   s_reverb.mix_,      [](float v) { return TextFormat("%.0f%%", v * 100.0f); });
    bind_knob(knob_pd,    s_reverb.predelay_,  [](float v) { return TextFormat("%.0f ms", v); });
    bind_knob(knob_decay, s_reverb.decay_,     [](float v) { return TextFormat("%.0f ms", v); });
    bind_knob(knob_size,  s_reverb.size_,      [](float v) { return TextFormat("%.2f", v); });
    bind_knob(knob_lp,    s_reverb.lowpass_,   [](float v) { return TextFormat("%.0f Hz", v); });
    bind_knob(knob_damp,  s_reverb.damping_,   [](float v) { return TextFormat("%.0f Hz", v); });

    // 强制同步初始值
    s_reverb.mix_.store(knob_mix.get_value(), std::memory_order_relaxed);
    s_reverb.predelay_.store(knob_pd.get_value(), std::memory_order_relaxed);
    s_reverb.decay_.store(knob_decay.get_value() / 1.0f, std::memory_order_relaxed);
    s_reverb.size_.store(knob_size.get_value(), std::memory_order_relaxed);
    s_reverb.lowpass_.store(knob_lp.get_value(), std::memory_order_relaxed);
    s_reverb.damping_.store(knob_damp.get_value(), std::memory_order_relaxed);

    // 触发首次参数加载
    s_reverb.params_dirty_.store(true, std::memory_order_relaxed);

    // ── 峰值电平显示状态（带 hold 衰减）──
    float disp_peak_in_l = 0.0f, disp_peak_in_r = 0.0f, disp_peak_out_l = 0.0f, disp_peak_out_r = 0.0f;

    // ── main loop ──
    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(Color{10, 10, 12, 255});

        // ── 旋钮 ──
        knob_mix.display();
        knob_pd.display();
        knob_decay.display();
        knob_size.display();
        knob_lp.display();
        knob_damp.display();

        // ── 峰值电平表 ──
        {
            constexpr float kDecay = 0.88f;
            constexpr float kBarW = 18.0f;
            constexpr float kGap = 4.0f;
            constexpr float kBarX1 = 14.0f;
            constexpr float kBarX2 = kBarX1 + kBarW + kGap;
            constexpr float kBarX3 = 574.0f;
            constexpr float kBarX4 = kBarX3 + kBarW + kGap;
            constexpr float kBarY = 14.0f;
            constexpr float kBarH = 252.0f;

            // 更新显示峰值（瞬间上升，缓慢衰减）
            auto update_peak = [](float& disp, std::atomic<float>& src) {
                float raw = src.load(std::memory_order_relaxed);
                disp = (raw > disp) ? raw : disp * kDecay;
            };
            update_peak(disp_peak_in_l, s_reverb.peak_in_l);
            update_peak(disp_peak_in_r, s_reverb.peak_in_r);
            update_peak(disp_peak_out_l, s_reverb.peak_out_l);
            update_peak(disp_peak_out_r, s_reverb.peak_out_r);

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
                DrawRectangle(static_cast<int>(x), static_cast<int>(y), static_cast<int>(w), static_cast<int>(h),
                              Color{25, 25, 30, 255});
                DrawRectangleLines(static_cast<int>(x), static_cast<int>(y), static_cast<int>(w), static_cast<int>(h),
                                   Color{60, 60, 70, 255});
                if (norm > 0.001f) {
                    int fill_h = static_cast<int>(norm * h);
                    DrawRectangle(static_cast<int>(x) + 1, static_cast<int>(y + h - fill_h), static_cast<int>(w) - 2,
                                  fill_h, col);
                }
                DrawText(label, static_cast<int>(x), static_cast<int>(y + h + 4), 9, Color{140, 140, 150, 255});
            };

            float n_in_l  = to_norm(disp_peak_in_l);
            float n_in_r  = to_norm(disp_peak_in_r);
            float n_outl  = to_norm(disp_peak_out_l);
            float n_outr  = to_norm(disp_peak_out_r);

            draw_bar(kBarX1, kBarW, kBarY, kBarH, n_in_l,  "L in",  bar_color(n_in_l));
            draw_bar(kBarX2, kBarW, kBarY, kBarH, n_in_r,  "R in",  bar_color(n_in_r));
            draw_bar(kBarX3, kBarW, kBarY, kBarH, n_outl,  "L out", bar_color(n_outl));
            draw_bar(kBarX4, kBarW, kBarY, kBarH, n_outr,  "R out", bar_color(n_outr));
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
