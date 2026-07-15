#include <AudioFile.h>
#include <audio_ops.hpp>
#include <qwqdsp/fx/pitch_shifter.hpp>
#include <qwqdsp/fx/pitch_shifter2.hpp>
#include <work_dir.hpp>
#include <qwqdsp/window/blackman_harris.hpp>

// ------------------------------------------------------------
// demos
// ------------------------------------------------------------

static void PitchShift_Grain() {
    auto wav_path = qwqdsp_support::WormholeWav();
    AudioFile<float> file{wav_path};
    auto& x_vec = file.samples.front();

    qwqdsp_fx::PitchShifter dsp;
    dsp.SetPitchShift(7.0f);
    for (float& x : x_vec) {
        x = dsp.Tick(x);
    }

    qwqdsp_support::AudioOps::Normalize(x_vec);

    file.setNumChannels(1);
    file.save(qwqdsp_support::OutputFile("pitch_shift_grain.wav"));
}

static void PitchShift_PhaseVocoder() {
    auto wav_path = qwqdsp_support::WormholeWav();
    AudioFile<float> file{wav_path};
    auto& x_vec = file.samples.front();
    auto x2_vec = x_vec;

    qwqdsp_fx::PhaseVocoder dsp;
    dsp.pitch_shift = 7.0f;
    dsp.Process(x_vec.data(), x2_vec.data(), x_vec.size());

    qwqdsp_support::AudioOps::Normalize(x_vec);

    file.setNumChannels(1);
    file.save(qwqdsp_support::OutputFile("pitch_shift_pv.wav"));
}

int main() {
    PitchShift_Grain();
    PitchShift_PhaseVocoder();
}
