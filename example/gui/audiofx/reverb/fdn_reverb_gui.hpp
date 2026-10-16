#pragma once

#include <atomic>
#include <cstddef>

#include "raylib.h"
#include "slider.hpp"

#include <qwqdsp/fx/fdn_reverb.hpp>

// ------------------------------------------------------------
//  FdnReverbGui — FDN 混响的 GUI 包装器
//
//  预设（延迟长度 / 反馈矩阵 / 输入输出增益）由模板参数注入，界面上只给衰减与搁架参数。
// ------------------------------------------------------------
template <qwqdsp_fx::CFdnPreset TPreset>
struct FdnReverbGui {
    // ── 原子参数（UI 线程写入，音频线程读取）──
    std::atomic<float> mix_{qwqdsp_fx::FdnReverb<TPreset>::kDefaultMix};
    std::atomic<float> decay_{qwqdsp_fx::FdnReverb<TPreset>::kDefaultDecayMs};         // ms
    std::atomic<float> high_decay_{qwqdsp_fx::FdnReverb<TPreset>::kDefaultHighDecayMs}; // ms
    std::atomic<float> shelf_{qwqdsp_fx::FdnReverb<TPreset>::kDefaultShelf};            // Hz
    std::atomic<bool> params_dirty_{true};

    static constexpr int kWindowWidth = 620;
    static constexpr int kWindowHeight = 280;

    // --------------------------------------------------------
    //  Init
    // --------------------------------------------------------
    void Init(float sample_rate) {
        reverb_.Init(sample_rate);
    }

    // --------------------------------------------------------
    //  Process — 音频线程回调（纯 DSP，不计算峰值）
    // --------------------------------------------------------
    void Process(const float* input, float* output_l, float* output_r, size_t frame_count) noexcept {
        if (params_dirty_.exchange(false, std::memory_order_relaxed)) {
            reverb_.SetMix(mix_.load(std::memory_order_relaxed));
            reverb_.SetDecay(decay_.load(std::memory_order_relaxed));
            reverb_.SetHighDecay(high_decay_.load(std::memory_order_relaxed));
            reverb_.SetShelfFrequency(shelf_.load(std::memory_order_relaxed));
        }

        for (size_t i = 0; i < frame_count; ++i) {
            float l, r;
            reverb_.Process(input[2 * i], input[2 * i + 1], &l, &r);
            output_l[i] = l;
            output_r[i] = r;
        }
    }

    // --------------------------------------------------------
    //  CreateKnobs — 创建所有旋钮并绑定回调
    // --------------------------------------------------------
    void CreateKnobs() {
        using Reverb = qwqdsp_fx::FdnReverb<TPreset>;
        constexpr int kKnobW = 58, kKnobH = 60;
        int kx = (kWindowWidth - (kNumKnobs * kKnobW + (kNumKnobs - 1) * kGap)) / 2;
        constexpr int ky = 42;

        knobs_[0] = MakeKnob(kx, ky, kKnobW, kKnobH, "Mix", 0.0f, 1.0f, 0.01f, Reverb::kDefaultMix);
        knobs_[1] = MakeKnob(kx, ky, kKnobW, kKnobH, "Decay", 200.0f, 10000.0f, 10.0f, Reverb::kDefaultDecayMs);
        knobs_[2] = MakeKnob(kx, ky, kKnobW, kKnobH, "HF Decay", 200.0f, 10000.0f, 10.0f, Reverb::kDefaultHighDecayMs);
        knobs_[3] = MakeKnob(kx, ky, kKnobW, kKnobH, "HF Freq", 500.0f, 20000.0f, 100.0f, Reverb::kDefaultShelf);

        std::atomic<float>* param_ptrs[kNumKnobs] = {&mix_, &decay_, &high_decay_, &shelf_};
        for (size_t i = 0; i < kNumKnobs; ++i) {
            auto* target = param_ptrs[i];
            auto* dirty = &params_dirty_;
            knobs_[i].on_value_change = [target, dirty](float v) {
                target->store(v, std::memory_order_relaxed);
                dirty->store(true, std::memory_order_relaxed);
            };
        }

        knobs_[0].value_to_text_function = [](float v) { return TextFormat("%.0f%%", v * 100.0f); };
        knobs_[1].value_to_text_function = [](float v) { return TextFormat("%.0f ms", v); };
        knobs_[2].value_to_text_function = [](float v) { return TextFormat("%.0f ms", v); };
        knobs_[3].value_to_text_function = [](float v) { return TextFormat("%.0f Hz", v); };

        for (size_t i = 0; i < kNumKnobs; ++i) {
            param_ptrs[i]->store(knobs_[i].get_value(), std::memory_order_relaxed);
        }
        params_dirty_.store(true, std::memory_order_relaxed);
    }

    // --------------------------------------------------------
    //  DisplayKnobs — 每帧绘制所有旋钮
    // --------------------------------------------------------
    void DisplayKnobs() {
        for (auto& k : knobs_) {
            k.display();
        }
    }

    static constexpr size_t kNumKnobs = 4;
    Knob knobs_[kNumKnobs];
private:
    static constexpr int kGap = 18;

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

    qwqdsp_fx::FdnReverb<TPreset> reverb_;
};

/// @brief 三档 FDN 混响的 GUI 包装器别名
using Fdn4SmallGui = FdnReverbGui<qwqdsp_fx::Fdn4PresetSmall>;
using Fdn4MediumGui = FdnReverbGui<qwqdsp_fx::Fdn4PresetMedium>;
using Fdn4LargeGui = FdnReverbGui<qwqdsp_fx::Fdn4PresetLarge>;
