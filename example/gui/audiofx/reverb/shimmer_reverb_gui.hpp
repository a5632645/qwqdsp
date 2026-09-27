#pragma once

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <numbers>

#include "raylib.h"
#include "slider.hpp"

// 包含同目录下的 shimmer DSP
#include "shimmer_reverb.hpp"

// ------------------------------------------------------------
//  ShimmerReverbGui — Shimmer 混响的完整 GUI 包装器
// ------------------------------------------------------------
struct ShimmerReverbGui {
    // ── 原子参数（UI 线程写入，音频线程读取）──
    std::atomic<float> mix_{0.5f};
    std::atomic<float> predelay_{20.0f};   // ms
    std::atomic<float> lowpass_{20000.0f}; // Hz
    std::atomic<float> decay_{2000.0f};    // ms
    std::atomic<float> size_{1.0f};
    std::atomic<float> damping_{5000.0f};   // Hz
    std::atomic<float> pitch_shift_{12.0f}; // semitones
    std::atomic<bool> params_dirty_{true};

    static constexpr const char* kWindowTitle = "Shimmer Reverb — ShimmerReverb";

    static constexpr int kWindowWidth = 620;
    static constexpr int kWindowHeight = 300;

    // --------------------------------------------------------
    //  Init
    // --------------------------------------------------------
    void Init(float sample_rate) {
        sr_ = sample_rate;
        reverb_.Init(sample_rate);
    }

    // --------------------------------------------------------
    //  Process — 音频线程回调（纯 DSP）
    // --------------------------------------------------------
    void Process(const float* input, float* output_l, float* output_r, size_t frame_count) noexcept {
        if (params_dirty_.exchange(false, std::memory_order_relaxed)) {
            float m = mix_.load(std::memory_order_relaxed);
            float pd = predelay_.load(std::memory_order_relaxed);
            float lp = lowpass_.load(std::memory_order_relaxed);
            float d = decay_.load(std::memory_order_relaxed);
            float sz = size_.load(std::memory_order_relaxed);
            float damp = damping_.load(std::memory_order_relaxed);
            float ps = pitch_shift_.load(std::memory_order_relaxed);

            reverb_.SetMix(m);
            reverb_.SetPredelay(pd / 1000.0f); // ms -> s
            reverb_.SetLowpass(lp);
            reverb_.SetDecay(d);
            reverb_.SetSize(sz);
            reverb_.SetDamping(damp);
            reverb_.SetPitchShift(ps);
        }

        for (size_t i = 0; i < frame_count; ++i) {
            float in_l = input[2 * i];
            float in_r = input[2 * i + 1];
            float l, r;
            reverb_.Process(in_l, in_r, &l, &r);
            output_l[i] = l;
            output_r[i] = r;
        }
    }

    // --------------------------------------------------------
    //  CreateKnobs — 两行布局
    // --------------------------------------------------------
    void CreateKnobs() {
        constexpr int kKnobW = 58, kKnobH = 60;
        constexpr int kGap4 = 18;
        int kx, ky;

        // ── 第 1 行：混响主体参数（6 个）──
        kx = (kWindowWidth - (6 * kKnobW + 5 * kGap4)) / 2;
        ky = 42;

        knobs_[0] = MakeKnob(kx, ky, kKnobW, kKnobH, "Mix", 0.0f, 1.0f, 0.01f, 0.5f);
        knobs_[1] = MakeKnob(kx, ky, kKnobW, kKnobH, "Predelay", 0.0f, 100.0f, 1.0f, 20.0f);
        knobs_[2] = MakeKnob(kx, ky, kKnobW, kKnobH, "Decay", 100.0f, 10000.0f, 50.0f, 2000.0f);
        knobs_[3] = MakeKnob(kx, ky, kKnobW, kKnobH, "Size", 0.05f, 2.0f, 0.01f, 1.0f);
        knobs_[4] = MakeKnob(kx, ky, kKnobW, kKnobH, "Lowpass", 16.0f, 20000.0f, 100.0f, 20000.0f);
        knobs_[5] = MakeKnob(kx, ky, kKnobW, kKnobH, "Damping", 16.0f, 20000.0f, 100.0f, 5000.0f);

        // ── 第 2 行：Pitch Shift（居中）──
        kx = (kWindowWidth - kKnobW) / 2;
        ky = 42 + kKnobH + 22;

        knobs_[6] = MakeKnob(kx, ky, kKnobW, kKnobH, "Pitch", -24.0f, 24.0f, 1.0f, 12.0f);

        // ── 旋钮回调绑定 ──
        std::atomic<float>* param_ptrs[kNumKnobs] = {&mix_,     &predelay_, &decay_,      &size_,
                                                     &lowpass_, &damping_,  &pitch_shift_};

        for (size_t i = 0; i < kNumKnobs; ++i) {
            auto* target = param_ptrs[i];
            auto* dirty = &params_dirty_;
            knobs_[i].on_value_change = [target, dirty](float v) {
                target->store(v, std::memory_order_relaxed);
                dirty->store(true, std::memory_order_relaxed);
            };
        }

        // 设置数值显示格式
        knobs_[0].value_to_text_function = [](float v) { return TextFormat("%.0f%%", v * 100.0f); };
        knobs_[1].value_to_text_function = [](float v) { return TextFormat("%.0f ms", v); };
        knobs_[2].value_to_text_function = [](float v) { return TextFormat("%.0f ms", v); };
        knobs_[3].value_to_text_function = [](float v) { return TextFormat("%.2f", v); };
        knobs_[4].value_to_text_function = [](float v) { return TextFormat("%.0f Hz", v); };
        knobs_[5].value_to_text_function = [](float v) { return TextFormat("%.0f Hz", v); };
        knobs_[6].value_to_text_function = [](float v) { return TextFormat("%+.0f st", v); };

        // 强制同步初始值
        for (size_t i = 0; i < kNumKnobs; ++i)
            param_ptrs[i]->store(knobs_[i].get_value(), std::memory_order_relaxed);

        // 触发首次参数加载
        params_dirty_.store(true, std::memory_order_relaxed);
    }

    // --------------------------------------------------------
    //  DisplayKnobs — 每帧绘制所有旋钮
    // --------------------------------------------------------
    void DisplayKnobs() {
        for (auto& k : knobs_)
            k.display();
    }

    static constexpr size_t kNumKnobs = 7;
    Knob knobs_[kNumKnobs];
private:
    // --------------------------------------------------------
    //  工具：创建旋钮
    // --------------------------------------------------------
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

    float sr_ = 48000.0f;
    ShimmerReverb reverb_;
};
