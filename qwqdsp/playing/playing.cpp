#include "AudioFile.h"
#include "qwqdsp/osciilor/vic_sine_osc.hpp"
#include "qwqdsp/filter/fir_hilbert.hpp"
#include "qwqdsp/filter/fir_hilbert_coeffs.h"

static constexpr auto kInputFile 
= R"(C:\Users\Kawai\Music\speech.wav)";
static constexpr auto kOutputFile 
= R"(C:\Users\Kawai\Music\gunge_slice-shift.wav)";
static constexpr auto kShift = 150; //hz

static void FreqShifter() {
    qwqdsp::filter::FirHilbert<qwqdsp::filter::fircoeff::Hilbert<float>> hilbert;
    qwqdsp::oscillor::VicSineOsc osc_;
    
    AudioFile<float> file;
    if (file.load(kInputFile)) {
        osc_.Reset(0);
        osc_.SetFreq(kShift, file.getSampleRate());

        auto& io = file.samples.front();
        for (auto& s : io) {
            auto analyze_signal = hilbert.Tick(s);
            osc_.Tick();
            auto quad = osc_.GetCpx();
            analyze_signal *= quad;
            s = analyze_signal.real();
        }

        file.setNumChannels(1);
        file.save(kOutputFile);
    }
}

int main() {
    FreqShifter();
}