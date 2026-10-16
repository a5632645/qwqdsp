#pragma once

// ------------------------------------------------------------
//  NaiveEchoGui — 最朴素单声道反馈延迟的 GUI 包装器
//
//  模型：y[n] = x[n] + d[n - D]，写入推迟线的是 x[n] + g·d[n - D]。
//  延迟取整数样本、读位置不做插值、参数不做任何平滑。
// ------------------------------------------------------------

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>

#include "raylib.h"
#include "slider.hpp"

#include "adpcm_delay_line.hpp"
#include "quantized_delay_line.hpp"

// 延迟线类型由存储格式给出（Storage::Line），因此压缩格式可以自带块存储实现
template <typename Storage>
struct NaiveEchoGui {
    static constexpr const char* kStorageName = Storage::kName;
    static constexpr const char* kStorageDescription = Storage::kDescription;

    static constexpr int kWindowWidth = 360;
    static constexpr int kWindowHeight = 210;
    static constexpr float kMaxDelayMs = 1000.0f;
    static constexpr size_t kNumKnobs = 2;

    // ── 原子参数（UI 线程写入，音频线程每块读取一次）──
    std::atomic<float> delay_ms_{200.0f};
    std::atomic<float> feedback_{0.5f};

    Knob knobs_[kNumKnobs];

    // --------------------------------------------------------
    //  Init
    // --------------------------------------------------------
    void Init(float sample_rate) {
        sample_rate_ = sample_rate;
        line_.Init(static_cast<size_t>(kMaxDelayMs * sample_rate * 0.001f));
    }

    // --------------------------------------------------------
    //  Process — 音频线程回调
    //
    //  整个 echo 都在该格式的运算域里跑：输入转换、延迟读取、反馈乘、饱和加、干+湿
    //  混音、写回；只有进出 miniaudio 缓冲时做浮点换算。
    // --------------------------------------------------------
    void Process(const float* input, float* output, size_t frame_count) noexcept {
        size_t const delay_samples = SampleDelay(delay_ms_.load(std::memory_order_relaxed));
        // 反馈量也进同一格式
        Value const gain = Storage::FromFloat(feedback_.load(std::memory_order_relaxed));

        for (size_t i = 0; i < frame_count; ++i) {
            Value const in = ConvertIn(input[i]);
            Value const delayed = line_.GetBeforePush(delay_samples);
            line_.Push(Storage::Add(in, ConvertMul(gain, delayed)));
            output[i] = Storage::ToFloat(Storage::Add(in, delayed));
        }
    }

    // --------------------------------------------------------
    //  CreateKnobs — 创建所有旋钮并绑定回调
    // --------------------------------------------------------
    void CreateKnobs() {
        constexpr int kKnobW = 58, kKnobH = 60;
        int kx = (kWindowWidth - (static_cast<int>(kNumKnobs) * kKnobW + static_cast<int>(kNumKnobs - 1) * kGap)) / 2;
        constexpr int ky = 42;

        knobs_[0] = MakeKnob(kx, ky, kKnobW, kKnobH, "Delay", 1.0f, kMaxDelayMs, 5.0f, 200.0f);
        knobs_[1] = MakeKnob(kx, ky, kKnobW, kKnobH, "Feedback", 0.0f, 0.99f, 0.01f, 0.5f);

        std::atomic<float>* param_ptrs[kNumKnobs] = {&delay_ms_, &feedback_};
        for (size_t i = 0; i < kNumKnobs; ++i) {
            knobs_[i].on_value_change = [target = param_ptrs[i]](float v) {
                target->store(v, std::memory_order_relaxed);
            };
        }
        knobs_[0].value_to_text_function = [](float v) { return TextFormat("%.0f ms", v); };
        knobs_[1].value_to_text_function = [](float v) { return TextFormat("%.0f%%", v * 100.0f); };

        for (size_t i = 0; i < kNumKnobs; ++i) {
            param_ptrs[i]->store(knobs_[i].get_value(), std::memory_order_relaxed);
        }
    }

    // --------------------------------------------------------
    //  DisplayKnobs — 每帧绘制所有旋钮
    // --------------------------------------------------------
    void DisplayKnobs() {
        for (auto& k : knobs_) {
            k.display();
        }
    }
private:
    using Line = typename Storage::Line;
    using Value = typename Storage::Value;

    /// @brief 浮点输入 → 该格式（16 位格式带一阶误差反馈）
    Value ConvertIn(float x) noexcept {
        if constexpr (Storage::kUseErrorFeedback) {
            return Storage::FromFloatShaped(x, ef_in_);
        }
        else {
            return Storage::FromFloat(x);
        }
    }

    /// @brief 该格式的定点乘（16 位格式带一阶误差反馈）
    Value ConvertMul(Value a, Value b) noexcept {
        if constexpr (Storage::kUseErrorFeedback) {
            return Storage::MulShaped(a, b, ef_mul_);
        }
        else {
            return Storage::Mul(a, b);
        }
    }

    static constexpr int kGap = 18;

    /// @brief 毫秒延迟换算为样本数，钳制到该存储实现可用的延迟范围
    size_t SampleDelay(float ms) const noexcept {
        size_t const n = static_cast<size_t>(ms * sample_rate_ * 0.001f);
        // 块压缩格式只能读到已经写完的整块，最短延迟因此大于 1 个样本
        return std::clamp<size_t>(n, Line::kMinDelaySamples, line_.Capacity());
    }

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
        kx += w + kGap;
        return knob;
    }

    float sample_rate_ = 48000.0f;
    // 16 位格式的一阶误差反馈余量（1/256 LSB 细网格）
    int64_t ef_in_ = 0;
    int64_t ef_mul_ = 0;
    Line line_;
};
