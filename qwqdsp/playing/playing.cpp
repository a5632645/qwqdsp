#include <mutex>
#include <span>
#include <semaphore>
#include <raylib.h>

#include "../playing/slider.hpp"
#include "qwqdsp/osciilor/noise.hpp"
#include "qwqdsp/filter/parallel_allpass.hpp"
#include "qwqdsp/osciilor/vic_sine_osc.hpp"
#include "qwqdsp/filter/iir_hilbert.hpp"
#include "qwqdsp/convert.hpp"

static qwqdsp::oscillor::PinkNoise noise;
static qwqdsp::filter::ParallelAllpass allpass;
static qwqdsp::filter::IIRHilbertDeeper<> hilbert;
static qwqdsp::filter::IIRHilbertDeeper<> hilbert2;
static qwqdsp::oscillor::VicSineOsc quadosc;

static std::binary_semaphore lock{1};

static void AudioInputCallback(void* _buffer, unsigned int frames) {
    lock.acquire();

    struct T {
        float l;
        float r;
    };
    std::span buffer{reinterpret_cast<T*>(_buffer), frames};

    for (auto& x : buffer) {
        float s = noise.Next();
        auto[up, down] = allpass.Tick(s);
        auto cpx_up = hilbert2.Tick(up);

        quadosc.Tick();
        auto cpx_down = hilbert.Tick(down);
        cpx_down *= quadosc.GetCpx();
        cpx_up *= quadosc.GetCpx();
        s = cpx_down.real() - cpx_up.real();

        s *= 0.5f;
        x.l = s;
        x.r = x.l;
    }

    lock.release();
}

int main(void) {
    InitWindow(640, 480, "formant");

    InitAudioDevice();
    SetAudioStreamBufferSizeDefault(512);
    AudioStream stream = LoadAudioStream(48000, 32, 2);
    SetAudioStreamCallback(stream, AudioInputCallback);
    PlayAudioStream(stream);

    Knob cutoff;
    cutoff.on_value_change = [](float v) {
        lock.acquire();
        allpass.BuildButterworth(15, v);
        lock.release();
    };
    cutoff.set_bound(0, 0, 80, 80);
    cutoff.set_range(0.01f, 3.0f, 0.01f, std::numbers::pi_v<float> / 2);
    cutoff.set_bg_color(BLACK);
    cutoff.set_fore_color(RAYWHITE);
    cutoff.set_title("cutoff");

    Knob shift;
    shift.on_value_change = [](float v) {
        lock.acquire();
        quadosc.SetFreq(v, 48000);
        lock.release();
    };
    shift.set_bound(0, 80, 80, 80);
    shift.set_range(-500.0f, 500.0f, 0.1f, 0.0f);
    shift.set_bg_color(BLACK);
    shift.set_fore_color(RAYWHITE);
    shift.set_title("shift");
    
    SetTargetFPS(30);
    while (!WindowShouldClose()) {
        BeginDrawing();
        {
            ClearBackground(BLACK);
            cutoff.display();
            shift.display();
        }
        EndDrawing();
    }

    UnloadAudioStream(stream);
    CloseAudioDevice();
    CloseWindow();
}
