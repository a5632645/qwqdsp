#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/window/hann.hpp>
#include <random>
#include <qwqdsp/convert.hpp>

// @ref https://apulsoft.ch/blog/allpass-fir-approximation/
int main() {
    qwqdsp_spectral::RealFFT fft;
    fft.Init(1024);

    constexpr int bins = fft.NumBins(1024);

    float ir[1024];
    float gains[bins];
    float phases[bins];

    float ripple_db = 0.1f;
    float max_phase_change = std::acos(2 * std::pow(10.0f, -ripple_db / 20.0f) - 1.0f);

    float bin_phase = 0.0f;
    std::uniform_real_distribution<float> dist{-max_phase_change, max_phase_change};
    std::default_random_engine rng{std::random_device{}()};
    for (int i = 0; i < bins; ++i) {
        gains[i] = 1.0f;
        phases[i] = bin_phase;
        bin_phase += dist(rng);
    }

    fft.IFFTGainPhase(ir, gains, phases);
    fft.TimeDomainShift(ir);
    qwqdsp_window::Hann::ApplyWindow(ir, false);

    fft.Init(2048);
    float pad[2048]{};
    std::copy_n(ir, 1024, pad);

    constexpr int pad_bins = fft.NumBins(2048);
    float pad_g[pad_bins];
    fft.FFTGainPhase(pad, pad_g);

    for (int i = 0; i < pad_bins; ++i) {
        pad_g[i] = qwqdsp::convert::Gain2Db<-100.0f>(pad_g[i]);
    }
}
