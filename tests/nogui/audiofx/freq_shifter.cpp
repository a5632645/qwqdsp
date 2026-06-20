#include "AudioFile.h"
#include "qwqdsp/filter/iir_cpx_hilbert.hpp"
#include "qwqdsp/filter/iir_hilbert.hpp"
#include "qwqdsp/oscillator/vic_sine_osc.hpp"
#include <work_dir.hpp>

static constexpr auto kShift = -150; // hz

static void FreqShifter() {
    qwqdsp_filter::IIRHilbertDeeper<> hilbert;
    qwqdsp_oscillator::VicSineOsc osc_;

    AudioFile<float> file;
    if (file.load(qwqdsp_support::WormholeWav())) {
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
        file.save(qwqdsp_support::OutputFile("freq_shift.wav"));
    }
}

static void FreqShifterAntialaising() {
    qwqdsp_filter::IIRHilbertDeeper<> hilbert;
    qwqdsp_filter::IIRHilbertDeeperCpx<> antialaising_filter;
    qwqdsp_oscillator::VicSineOsc osc_;

    AudioFile<float> file;
    if (file.load(qwqdsp_support::WormholeWav())) {
        osc_.Reset(0);
        osc_.SetFreq(kShift, file.getSampleRate());

        auto& io = file.samples.front();
        for (auto& s : io) {
            auto analyze_signal = hilbert.Tick(s);
            osc_.Tick();
            auto quad = osc_.GetCpx();
            analyze_signal *= quad;
            // 移除负频率
            analyze_signal = antialaising_filter.Tick(analyze_signal);
            s = analyze_signal.real();
        }

        file.setNumChannels(1);
        file.save(qwqdsp_support::OutputFile("freq_shift2.wav"));
    }
}

int main() {
    FreqShifter();
    FreqShifterAntialaising();
}
