#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <numbers>

#include "miniaudio.h"
#include "raylib.h"
#include "slider.hpp"

#include "reverb/peak_meter.hpp"
#include "reverb/plate_reverb_gui.hpp"
#include "reverb/shimmer_reverb_gui.hpp"

// ------------------------------------------------------------
//  ReverbGuiConcept — 混响 GUI 包装器必须满足的接口
//
// 任何用作 ReverbGui 的类型都需要提供以下成员：
//   - static constexpr const char* kWindowTitle
//   - static constexpr int         kWindowWidth
//   - static constexpr int         kWindowHeight
//   - void Init(float sample_rate)
//   - void Process(const float* input, float* output_l,
//                  float* output_r, size_t frame_count) noexcept
//   - void CreateKnobs()
//   - void DisplayKnobs()
// ------------------------------------------------------------
template <typename T>
concept ReverbGuiConcept = requires(T& obj, const float* in, float* ol, float* or_, size_t n, float sr) {
    // ── 静态编译期常量 ──
    { T::kWindowTitle } -> std::convertible_to<const char*>;
    { T::kWindowWidth } -> std::convertible_to<int>;
    { T::kWindowHeight } -> std::convertible_to<int>;

    // ── 方法签名 ──
    { obj.Init(sr) } -> std::same_as<void>;
    { obj.Process(in, ol, or_, n) } noexcept -> std::same_as<void>;
    { obj.CreateKnobs() } -> std::same_as<void>;
    { obj.DisplayKnobs() } -> std::same_as<void>;
};

using ReverbGui = ShimmerReverbGui;
static_assert(ReverbGuiConcept<ReverbGui>,
              "ReverbGui 必须满足 ReverbGuiConcept：需要 kWindowTitle/kWindowWidth/kWindowHeight "
              "静态常量，以及 Init/Process/CreateKnobs/DisplayKnobs 方法");

// 全局实例
static ReverbGui s_reverb;
static PeakMeter s_peak_meter;

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

        // 处理前：测量输入峰值
        s_peak_meter.MeasureInput(input + offset, nf);

        // DSP 处理
        s_reverb.Process(input + offset, buf_l, buf_r, nf);

        // 处理后：测量输出峰值
        s_peak_meter.MeasureOutput(buf_l, buf_r, nf);

        for (size_t i = 0; i < nf; ++i) {
            output[2 * (offset + i)] = buf_l[i];
            output[2 * (offset + i) + 1] = buf_r[i];
        }

        offset += nf;
        remaining -= nf;
    }
}

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------
int main(void) {
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(ReverbGui::kWindowWidth, ReverbGui::kWindowHeight, ReverbGui::kWindowTitle);
    SetTargetFPS(60);

    static constexpr float kSampleRate = 48000.0f;

    // ── 初始化混响 ──
    s_reverb.Init(kSampleRate);
    s_reverb.CreateKnobs();

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

    // ── main loop ──
    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(Color{10, 10, 12, 255});

        // ── 旋钮 ──
        s_reverb.DisplayKnobs();

        // ── 峰值电平表（hold 衰减更新 + 绘制）──
        {
            constexpr float kDecay = 0.88f;
            constexpr float kBarW = 18.0f;
            constexpr float kGap = 4.0f;
            constexpr float kGroupGap = 12.0f;
            // 所有电平条放在右侧：4 条 ×18 + 3 个间距(4+12+4) = 92px 偏移
            constexpr float kBarX = static_cast<float>(ReverbGui::kWindowWidth) - 92.0f - 14.0f;
            constexpr float kBarY = 14.0f;
            constexpr float kBarH = 252.0f;

            s_peak_meter.Update(kDecay);
            s_peak_meter.Draw(kBarX, kBarY, kBarH, kBarW, kGap, kGroupGap);
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
