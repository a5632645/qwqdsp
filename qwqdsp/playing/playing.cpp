#include "AudioFile.h"
#include "qwqdsp/adaptive/nlms.hpp"
#include "qwqdsp/adaptive/rls_filter.hpp"

int main() {
    AudioFile<float> female;
    female.load(R"(C:\Users\Kawai\Downloads\Music\samples_x.wav)");
    AudioFile<float> mix;
    mix.load(R"(C:\Users\Kawai\Downloads\Music\samples_d.wav)");

    AudioFile<float>::AudioBuffer buffer;
    buffer.resize(1);
    auto& output = buffer.front();
    auto& x = female.samples.front();
    auto& target = mix.samples.front();
    size_t const size = std::min(x.size(), target.size());

    // qwqdsp::adaptive::NLMS<512> nlms;
    // nlms.Reset();
    // for (size_t i = 0; i < size; ++i) {
    //     float const pred = nlms.Tick(x[i], target[i]);
    //     output.push_back(target[i] - pred);
    // }

    AudioFile<float> male;
    male.setAudioBuffer(buffer);
    male.setSampleRate(female.getSampleRate());
    male.setBitDepth(32);
    male.save(R"(C:\Users\Kawai\Music\male.wav)");
}