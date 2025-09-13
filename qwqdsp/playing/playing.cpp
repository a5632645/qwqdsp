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

    static constexpr float kOutputFs = 48000.0f;

    // qwqdsp::fx::ResampleIIR resample;
    // resample.Init(infile.getSampleRate(), kOutputFs, 100, 127);
    // resample.Init(infile.getSampleRate(), kOutputFs);
    // output.front() = resample.Process(input);

    qwqdsp::fx::ResampleIIRDynamic resample;
    resample.Init(infile.getSampleRate(), 10000.0f);
    // resample.Set(infile.getSampleRate(), 48000.0f);
    resample.SetPitchShift(-24.0f);

    auto it_in = input.begin();
    while (it_in != input.end()) {
        while (!resample.IsReady()) {
            resample.Push(*it_in);
            ++it_in;
            if (it_in == input.end()) {
                break;
            }
        }

        while (resample.IsReady()) {
            out.push_back(resample.Read());
        }
    }

    AudioFile<float> outfile;
    outfile.setAudioBuffer(output);
    outfile.setSampleRate(kOutputFs);
    outfile.setBitDepth(32);
    outfile.save(R"(C:\Users\Kawai\Music\resample.wav)");
}