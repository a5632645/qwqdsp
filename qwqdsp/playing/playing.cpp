#include "qwqdsp/fx/resample_iir.hpp"
#include "qwqdsp/fx/resample.hpp"
#include "qwqdsp/fx/resample_iir_dynamic.hpp"
#include "AudioFile.h"

int main() {
    constexpr auto kPath = R"(C:\Users\Kawai\Music\sweep.wav)";
    AudioFile<float> infile;
    infile.load(kPath);
    auto& input = infile.samples.front();

    AudioFile<float>::AudioBuffer output;
    output.resize(1);
    auto& out = output.front();

    static constexpr float kOutputFs = 44100.0f;

    qwqdsp::fx::Resample resample;
    resample.Init(infile.getSampleRate(), kOutputFs, 180, 1023);
    out = resample.Process(infile.samples.front());

    AudioFile<float> outfile;
    outfile.setAudioBuffer(output);
    outfile.setSampleRate(kOutputFs);
    outfile.setBitDepth(32);
    outfile.save(R"(C:\Users\Kawai\Music\resample.wav)");
}