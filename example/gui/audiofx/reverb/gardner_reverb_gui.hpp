#pragma once

#include <atomic>
#include <cstddef>

#include "raylib.h"
#include "slider.hpp"

#include <qwqdsp/fx/gardner_reverb.hpp>

// ------------------------------------------------------------
//  GardnerReverbGui — Gardner 扩散混响的 GUI 包装器
//
//  混响类型由模板参数给出（小房间 / 大厅），参数范围与默认值取自该类型自身的常量。
// ------------------------------------------------------------
template <typename TReverb>
struct GardnerReverbGui {
    // ── 原子参数（UI 线程写入，音频线程读取）──
    std::atomic<float> mix_{TReverb::kDefaultMix};
    std::atomic<float> decay_{TReverb::kDefaultDecayMs};     // ms
    std::atomic<float> damping_{TReverb::kDefaultDamping};   // Hz
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
            reverb_.SetDamping(damping_.load(std::memory_order_relaxed));
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
        constexpr int kKnobW = 58, kKnobH = 60;
        int kx = (kWindowWidth - (kNumKnobs * kKnobW + (kNumKnobs - 1) * kGap)) / 2;
        constexpr int ky = 42;

        knobs_[0] = MakeKnob(kx, ky, kKnobW, kKnobH, "Mix", 0.0f, 1.0f, 0.01f, TReverb::kDefaultMix);
        knobs_[1] = MakeKnob(kx, ky, kKnobW, kKnobH, "Decay", TReverb::kMinDecayMs, TReverb::kMaxDecayMs, 10.0f,
                             TReverb::kDefaultDecayMs);
        knobs_[2] = MakeKnob(kx, ky, kKnobW, kKnobH, "Damping", 500.0f, 20000.0f, 100.0f, TReverb::kDefaultDamping);

        std::atomic<float>* param_ptrs[kNumKnobs] = {&mix_, &decay_, &damping_};
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
        knobs_[2].value_to_text_function = [](float v) { return TextFormat("%.0f Hz", v); };

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

    static constexpr size_t kNumKnobs = 3;
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

    TReverb reverb_;
};

/// @brief Gardner 小房间扩散混响的 GUI 包装器
using GardnerRoomGui = GardnerReverbGui<qwqdsp_fx::GardnerRoomReverb>;

/// @brief Gardner 大厅扩散混响的 GUI 包装器
using GardnerHallGui = GardnerReverbGui<qwqdsp_fx::GardnerHallReverb>;
