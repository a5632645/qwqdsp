#include "qwqdsp/fx/resample_iir.hpp"
#include "qwqdsp/fx/resample.hpp"
#include "qwqdsp/fx/resample_iir_dynamic.hpp"
#include "qwqdsp/fx/new_coeffs.h"
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

    // qwqdsp::fx::ResampleIIR<BestCoeffs<float>, 255> resample;
    // resample.Init(infile.getSampleRate(), kOutputFs);
    // out = resample.Process(infile.samples.front());

    qwqdsp::fx::ResampleIIRDynamic<BestCoeffs<float>, 255> resample;
    resample.Init(infile.getSampleRate());
    resample.Set(infile.getSampleRate(), kOutputFs);
    auto it = input.begin();
    while (it != input.end()) {
        while (!resample.IsReady() && it != input.end()) {
            resample.Push(*it);
            ++it;
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