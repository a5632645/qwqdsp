#include <qwqdsp/fx/uniform_convolution.hpp>
#include <qwqdsp/oscillator/noise.hpp>
#include <work_dir.hpp>
#include <AudioFile.h>
#include <qwqdsp/convert.hpp>
#include <audio_ops.hpp>

int main() {
    AudioFile<float> file{qwqdsp_support::InputFile("drumloop.wav")};

    constexpr float decay_ms = 2000;
    int samples = decay_ms * file.getSampleRate() / 1000.0f;
    std::vector<float> ir;
    ir.resize(samples);

    float g = 1.0f;
    float mul = qwqdsp::convert::Samples2DecayDb(samples, -60.0f);
    qwqdsp_oscillator::Clicks2 noise;
    noise.SetProbability(0.1f);
    for (int i = 0; i < samples; ++i) {
        ir[i] = g * noise.Next();
        g *= mul;
    }

    qwqdsp_fx::UniformConvolution conv;
    conv.Init(1024);
    conv.SetIR(ir);
    conv.Process(file.samples.front());

    qwqdsp_support::AudioOps::Normalize(file.samples.front());

    file.save(qwqdsp_support::OutputFile("conv_noise.wav"));
}
