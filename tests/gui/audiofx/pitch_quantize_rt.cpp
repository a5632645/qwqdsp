#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <numbers>

#include "miniaudio.h"
#include "raylib.h"
#include "slider.hpp"

#include "pitch_quantize/pitch_quantizer_rt_gui.hpp"

static constexpr float kSampleRate = 48000.0f;

// 全局实例
static PitchQuantizerRTGui s_pq;

// ------------------------------------------------------------
//  miniaudio 回调（单声道入 → 单声道出）
// ------------------------------------------------------------
extern "C" void MaCallback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;

    auto* input = static_cast<const float*>(pInput);
    auto* output = static_cast<float*>(pOutput);

    if (input == nullptr || output == nullptr)
        return;

    constexpr size_t kChunk = 512;
    float buf[kChunk];
    size_t remaining = static_cast<size_t>(frameCount);
    size_t offset = 0;

    while (remaining > 0) {
        size_t nf = std::min(remaining, kChunk);
        s_pq.Process(input + offset, buf, nf);

        for (size_t i = 0; i < nf; ++i)
            output[offset + i] = buf[i];

        offset += nf;
        remaining -= nf;
    }
}

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------
int main(void) {
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(PitchQuantizerRTGui::kWindowWidth, PitchQuantizerRTGui::kWindowHeight,
               PitchQuantizerRTGui::kWindowTitle);
    SetTargetFPS(60);

    // ──初始化音高量化器 ──
    s_pq.Init(kSampleRate);
    s_pq.CreateKnobs();

    // ── miniaudio full-duplex（单声道入 → 单声道出）──
    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.playback.format = ma_format_f32;
    config.playback.channels = 1;
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

        // ── 控件 ──
        s_pq.DisplayKnobs();

        // ── 底部提示 ──
        DrawText("right click knob: reset", 8,
                 PitchQuantizerRTGui::kWindowHeight - 18, 10, Color{80, 80, 90, 255});

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
