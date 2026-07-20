#pragma once

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>

#include "raylib.h"
#include "slider.hpp"

#include "pitch_quantizer_rt.hpp"
#include "../../../nogui/audiofx/pitch_quantize/scale_helper.hpp"

// ------------------------------------------------------------
//  PitchQuantizerRTGui — 实时音高量化器的完整 GUI 包装器
// ------------------------------------------------------------
struct PitchQuantizerRTGui {
    // ── 原子参数（UI 线程写入，音频线程读取）──
    std::atomic<int> root_note_{0};        // 0=C, 1=C#, ..., 11=B
    std::atomic<int> scale_type_{0};       // 0=Major, 1=MinorNatural, ...
    std::atomic<float> mix_{1.0f};
    std::atomic<float> color_{1.0f};
    std::atomic<float> low_cut_hz_{0.0f};
    std::atomic<float> high_cut_hz_{20000.0f};
    std::atomic<bool> params_dirty_{true};

    static constexpr const char* kWindowTitle = "Pitch Quantizer — Real-time";

    static constexpr int kWindowWidth = 480;
    static constexpr int kWindowHeight = 180;

    // 根音名称
    static constexpr const char* kRootNames[12] = {
        "C", "C#", "D", "D#", "E", "F",
        "F#", "G", "G#", "A", "A#", "B"
    };

    // 音阶名称（与 ScaleHelper::Type 顺序一致）
    static constexpr size_t kNumScaleTypes = 14;
    static constexpr const char* kScaleNames[kNumScaleTypes] = {
        // 核心
        "Major", "Minor",
        // 教会调式
        "Dor", "Phr", "Lyd", "Mix", "Loc",
        // 古典小调
        "HMin", "MMin",
        // 五声与特殊
        "MPen", "mPen", "Blue", "WT",
        // 半音阶
        "Chr"
    };

    // --------------------------------------------------------
    //  Init
    // --------------------------------------------------------
    void Init(float sample_rate) {
        constexpr size_t kFftSize = 2048;
        constexpr size_t kHopSize = kFftSize / 4; // 512
        pq_.Init(sample_rate, kHopSize, kFftSize);
        sr_ = sample_rate;
    }

    // --------------------------------------------------------
    //  Process — 音频线程回调（纯 DSP）
    // --------------------------------------------------------
    void Process(const float* input, float* output, size_t frame_count) noexcept {
        // 参数变更时同步到 DSP
        if (params_dirty_.exchange(false, std::memory_order_relaxed)) {
            int root = root_note_.load(std::memory_order_relaxed);
            int st = scale_type_.load(std::memory_order_relaxed);
            pq_.SetKey(root, static_cast<ScaleHelper::Type>(st));
        }

        pq_.mix_ = mix_.load(std::memory_order_relaxed);
        pq_.SetColor(color_.load(std::memory_order_relaxed));
        pq_.SetLowCut(low_cut_hz_.load(std::memory_order_relaxed));
        pq_.SetHighCut(high_cut_hz_.load(std::memory_order_relaxed));

        // 单声道：拷贝输入到输出缓冲，然后 in-place 处理
        std::copy(input, input + frame_count, output);
        pq_.Process(std::span<float>{output, frame_count});
    }

    // --------------------------------------------------------
    //  CreateKnobs — 创建控件
    // --------------------------------------------------------
    void CreateKnobs() {
        // 底部一行：4 个旋钮居中排列
        constexpr int kKnobW = 58, kKnobH = 60;
        constexpr int kGap = 14;
        constexpr int kRowW = 4 * kKnobW + 3 * kGap;
        int kx = (kWindowWidth - kRowW) / 2;
        int ky = 90;

        low_cut_knob_ = MakeKnob(kx, ky, kKnobW, kKnobH, "Low Cut", 20.0f, 20000.0f, 10.0f, 20.0f);
        low_cut_knob_.value_to_text_function = [](float v) { return TextFormat("%.0f Hz", v); };
        low_cut_knob_.on_value_change = [this](float v) {
            low_cut_hz_.store(v, std::memory_order_relaxed);
        };
        low_cut_hz_.store(low_cut_knob_.get_value(), std::memory_order_relaxed);
        kx += kKnobW + kGap;

        high_cut_knob_ = MakeKnob(kx, ky, kKnobW, kKnobH, "High Cut", 20.0f, 20000.0f, 100.0f, 20000.0f);
        high_cut_knob_.value_to_text_function = [](float v) { return TextFormat("%.0f Hz", v); };
        high_cut_knob_.on_value_change = [this](float v) {
            high_cut_hz_.store(v, std::memory_order_relaxed);
        };
        high_cut_hz_.store(high_cut_knob_.get_value(), std::memory_order_relaxed);
        kx += kKnobW + kGap;

        color_knob_ = MakeKnob(kx, ky, kKnobW, kKnobH, "Color", 0.0f, 1.0f, 0.01f, 1.0f);
        color_knob_.value_to_text_function = [](float v) { return TextFormat("%.0f%%", v * 100.0f); };
        color_knob_.on_value_change = [this](float v) {
            color_.store(v, std::memory_order_relaxed);
        };
        color_.store(color_knob_.get_value(), std::memory_order_relaxed);
        kx += kKnobW + kGap;

        mix_knob_ = MakeKnob(kx, ky, kKnobW, kKnobH, "Mix", 0.0f, 1.0f, 0.01f, 1.0f);
        mix_knob_.value_to_text_function = [](float v) { return TextFormat("%.0f%%", v * 100.0f); };
        mix_knob_.on_value_change = [this](float v) {
            mix_.store(v, std::memory_order_relaxed);
        };
        mix_.store(mix_knob_.get_value(), std::memory_order_relaxed);

        // 选择器区域
        selector_root_x_ = 14;
        selector_root_y_ = 14;
        selector_root_w_ = static_cast<float>(kWindowWidth) - 14.0f - 14.0f;

        selector_scale_y_ = selector_root_y_ + 30.0f + 10.0f;
        selector_scale_w_ = selector_root_w_;

        // 强制同步初始值
        params_dirty_.store(true, std::memory_order_relaxed);
    }

    // --------------------------------------------------------
    //  DisplayKnobs — 每帧绘制所有控件
    // --------------------------------------------------------
    void DisplayKnobs() {
        auto const mouse_pos = GetMousePosition();

        // ---- 根音选择器 ----
        DrawText("Root", selector_root_x_, selector_root_y_ - 12, 10, Color{140, 140, 150, 255});
        DrawSelector(
            Rectangle{selector_root_x_, selector_root_y_, selector_root_w_, 30.0f},
            12, kRootNames,
            static_cast<size_t>(root_note_.load(std::memory_order_relaxed)),
            mouse_pos,
            [this](size_t i) {
                root_note_.store(static_cast<int>(i), std::memory_order_relaxed);
                params_dirty_.store(true, std::memory_order_relaxed);
            });

        // ---- 音阶选择器 ----
        DrawText("Scale", selector_root_x_, selector_scale_y_ - 12, 10, Color{140, 140, 150, 255});
        DrawSelector(
            Rectangle{selector_root_x_, selector_scale_y_, selector_scale_w_, 30.0f},
            kNumScaleTypes, kScaleNames,
            static_cast<size_t>(scale_type_.load(std::memory_order_relaxed)),
            mouse_pos,
            [this](size_t i) {
                scale_type_.store(static_cast<int>(i), std::memory_order_relaxed);
                params_dirty_.store(true, std::memory_order_relaxed);
            });

        // ---- 底部一行：Low Cut / High Cut / Color / Mix ----
        low_cut_knob_.display();
        high_cut_knob_.display();
        color_knob_.display();
        mix_knob_.display();
    }

    // --------------------------------------------------------
    //  Reset
    // --------------------------------------------------------
    void Reset() noexcept {
        pq_.Reset();
    }

private:
    // --------------------------------------------------------
    //  工具：创建旋钮（与 reverb 风格一致）
    // --------------------------------------------------------
    static Knob MakeKnob(int kx, int ky, int w, int h, const char* title, float vmin, float vmax, float vstep,
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
        return knob;
    }

    // --------------------------------------------------------
    //  工具：绘制选择器
    // --------------------------------------------------------
    static void DrawSelector(Rectangle bound, size_t num_items, const char* const* names,
                             size_t selected, Vector2 mouse_pos,
                             std::function<void(size_t)> on_click) {
        float const each_width = bound.width / static_cast<float>(num_items);

        for (size_t i = 0; i < num_items; ++i) {
            float const x = bound.x + static_cast<float>(i) * each_width;
            float const w = each_width - 2.0f;

            Rectangle const item{x, bound.y, w, bound.height};
            bool const hovered = CheckCollisionPointRec(mouse_pos, item);

            if (i == selected) {
                DrawRectangleRec(item, WHITE);
                DrawText(names[i], static_cast<int>(x) + 4, static_cast<int>(bound.y) + 8, 12, BLACK);
            }
            else {
                Color const fore = hovered ? Color{180, 180, 190, 255} : Color{100, 100, 110, 255};
                DrawRectangleLinesEx(item, 1.0f, fore);
                DrawText(names[i], static_cast<int>(x) + 4, static_cast<int>(bound.y) + 8, 12, fore);
            }

            if (hovered && IsMouseButtonPressed(MOUSE_LEFT_BUTTON) && on_click != nullptr) {
                on_click(i);
            }
        }
    }

    float sr_ = 48000.0f;
    PitchQuantizerRT pq_;

    Knob low_cut_knob_;
    Knob high_cut_knob_;
    Knob color_knob_;
    Knob mix_knob_;

    // 选择器位置
    float selector_root_x_ = 14.0f;
    float selector_root_y_ = 14.0f;
    float selector_root_w_ = 0.0f;
    float selector_scale_y_ = 0.0f;
    float selector_scale_w_ = 0.0f;
};
