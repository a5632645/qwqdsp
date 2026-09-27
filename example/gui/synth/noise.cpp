#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <span>

#include "raylib.h"
#include "slider.hpp"

#include "qwqdsp/oscillator/noise.hpp"

static constexpr int kWidth = 500;
static constexpr int kHeight = 400;

enum NoiseType {
    White = 0,
    Pink,
    PinkHQ,
    Brown,
    Clicks,
    Clicks2,
    Stair,
    NumNoiseTypes
};
static constexpr const char* kNoiseNames[]{
    "white",
    "pink",
    "pink_hq",
    "brown",
    "clicks",
    "clicks2",
    "stair"
};

static NoiseType noise_type = NoiseType::White;
static qwqdsp_oscillator::WhiteNoise white_noise;
static qwqdsp_oscillator::PinkNoise pink_noise;
static qwqdsp_oscillator::PinkNoiseHQ pink_hq_noise;
static qwqdsp_oscillator::BrownNoise brown_noise;
static qwqdsp_oscillator::Clicks clicks;
static qwqdsp_oscillator::Clicks2 clicks2;
static qwqdsp_oscillator::Stair stair;

static void AudioInputCallback(void* _buffer, unsigned int frames) {
    struct T {
        float l;
        float r;
    };
    std::span buffer{reinterpret_cast<T*>(_buffer), frames};

    switch (noise_type) {
    case White:
        for (auto& s : buffer) {
            s.l = white_noise.Next();
            s.r = s.l;
        }
        break;
    case Pink:
        for (auto& s : buffer) {
            s.l = pink_noise.Next();
            s.r = s.l;
        }
        break;
    case PinkHQ:
        for (auto& s : buffer) {
            s.l = pink_hq_noise.Next();
            s.r = s.l;
        }
        break;
    case Brown:
        for (auto& s : buffer) {
            s.l = brown_noise.Next();
            s.r = s.l;
        }
        break;
    case Clicks:
        for (auto& s : buffer) {
            s.l = clicks.Next();
            s.r = s.l;
        }
        break;
    case Clicks2:
        for (auto& s : buffer) {
            s.l = clicks2.Next();
            s.r = s.l;
        }
        break;
    case Stair:
        for (auto& s : buffer) {
            s.l = stair.Next();
            s.r = s.l;
        }
        break;
    default:
        assert(false);
    }
}

int main(void) {
    InitWindow(kWidth, kHeight, "noise oscillator");

    InitAudioDevice();
    SetAudioStreamBufferSizeDefault(512);
    AudioStream stream = LoadAudioStream(48000, 32, 2);
    SetAudioStreamCallback(stream, AudioInputCallback);
    PlayAudioStream(stream);

    Rectangle bound;
    bound.x = 0;
    bound.y = 0;
    bound.width = 50;
    bound.height = 50;

    Knob probability_knob;
    probability_knob.on_value_change = [](float v) {
        clicks.SetProbability(v);
        clicks2.SetProbability(v);
        stair.SetProbability(v);
    };
    bound.y += bound.height;
    probability_knob.set_bound(bound);
    probability_knob.set_range(0.0f, 1.0f, 0.001f, 0.01f);
    probability_knob.set_bg_color(BLACK);
    probability_knob.set_fore_color(RAYWHITE);
    probability_knob.set_title("probability");

    bound.y += bound.height;
    auto type_bound = bound;
    type_bound.width = GetScreenWidth();
    type_bound.height = 30;
    float each_width = type_bound.width / static_cast<float>(NumNoiseTypes);
    size_t num_types = static_cast<size_t>(NumNoiseTypes);

    SetTargetFPS(30);
    while (!WindowShouldClose()) {
        BeginDrawing();
        {
            ClearBackground(BLACK);

            probability_knob.display();

            auto mouse_pos = GetMousePosition();
            for (size_t i = 0; i < num_types; ++i) {
                float x = i * each_width;
                float y = type_bound.y;
                float w = each_width - 2;
                float h = type_bound.height;
                if (CheckCollisionPointRec(mouse_pos, Rectangle{x,y,w,h}) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
                    noise_type = static_cast<NoiseType>(i);
                }

                if (i != noise_type) {
                    DrawRectangleLines(x, y, w, h, WHITE);
                    DrawText(kNoiseNames[i], x, y, 12, WHITE);
                }
                else {
                    DrawRectangle(x, y, w, h, WHITE);
                    DrawText(kNoiseNames[i], x, y, 12, BLACK);
                }
            }
        }
        EndDrawing();
    }

    UnloadAudioStream(stream);
    CloseAudioDevice();
    CloseWindow();
}
