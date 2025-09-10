#include "qwqdsp/fx/resample.hpp"
#include "AudioFile.h"

int main() {
constexpr auto kPath = R"(C:\Users\Kawai\Music\sweep.wav)";
    AudioFile<float> infile;
    infile.load(kPath);

    constexpr auto target_fs = 44100;
    qwqdsp::fx::Resample resample;
    resample.Init(infile.getSampleRate(), target_fs, 130, 255);

    AudioFile<float>::AudioBuffer output;
    output.resize(1);
    output.front() = resample.Process(infile.samples.front());

    AudioFile<float> outfile;
    outfile.setAudioBuffer(output);
    outfile.setSampleRate(target_fs);
    outfile.setBitDepth(32);
    outfile.save(R"(C:\Users\Kawai\Music\gunge_slice-synth.wav)");
}