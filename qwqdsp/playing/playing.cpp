#include <algorithm>
#include <array>
#include <cstddef>
#include <span>

#include "raylib.h"
#include "../playing/slider.hpp"

#include "qwqdsp/osciilor/poly_sine_osc.hpp"

static constexpr int kWidth = 500;
static constexpr int kHeight = 400;
static constexpr float kFs = 48000.0f;

static qwqdsp::oscillor::PolySineOsc dsp;

static void AudioInputCallback(void* _buffer, unsigned int frames) {
    struct T {
        float l;
        float r;
    };
    std::span buffer{reinterpret_cast<T*>(_buffer), frames};
    for (auto& s : buffer) {
        dsp.Tick();
        s.l = dsp.Sine();
        s.r = dsp.Cosine();
    }
}

int main(void) {
    InitWindow(kWidth, kHeight, "formant");

    InitAudioDevice();
    SetAudioStreamBufferSizeDefault(512);
    AudioStream stream = LoadAudioStream(48000, 32, 2);
    SetAudioStreamCallback(stream, AudioInputCallback);
    PlayAudioStream(stream);
    
    Rectangle dsf_bound;
    dsf_bound.x = 0;
    dsf_bound.y = 0;
    dsf_bound.width = 50;
    dsf_bound.height = 50;

    Knob w;
    w.on_value_change = [](float v) {
        dsp.SetW(v);
    };
    dsf_bound.y += dsf_bound.height;
    w.set_bound(dsf_bound);
    w.set_range(0.0f, 3.14f, 0.001f, 0.01f);
    w.set_bg_color(BLACK);
    w.set_fore_color(RAYWHITE);
    w.set_title("w");
    
    SetTargetFPS(30);
    while (!WindowShouldClose()) {
        BeginDrawing();
        {
            ClearBackground(BLACK);

            w.display();
        }
        EndDrawing();
    }

    UnloadAudioStream(stream);
    CloseAudioDevice();
    CloseWindow();
}
