#include "AudioFile.h"
#include "qwqdsp/adaptive/burg_lp.hpp"
#include <numbers>
#include <complex>

int main() {
    AudioFile<float> speech;
    speech.load(R"(C:\Users\Kawai\Music\speech.wav)");
    auto& x = speech.samples.front();

    qwqdsp::adaptive::BurgLP burg;
    burg.Init(x.size());

    float a[4000];
    burg.Process(x, a);

    speech.save(R"(C:\Users\Kawai\Music\speech-whiten.wav)");
}