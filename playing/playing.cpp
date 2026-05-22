#include <qwqdsp/spectral/ipp_real_fft.hpp>
#include <qwqdsp/window/hann.hpp>
#include <qwqdsp/window/helper.hpp>
#include <qwqdsp/oscillator/vic_sine_osc.hpp>

int main() {
    constexpr int fft_size = 2048;
    constexpr int fft_bins = fft_size / 2 + 1;

    float test[fft_size];
    qwqdsp_oscillator::VicSineOsc osc;
    osc.SetFreq(100.0f, fft_size);
    osc.Reset();
    for (int i = 0; i < fft_size; ++i) {
        test[i] = osc.Tick();
    }

    float window[fft_size];
    qwqdsp_window::Hann::Window(window, true);
    // qwqdsp_window::Helper::Normalize(window);

    for (int i = 0; i < fft_size; ++i) {
        test[i] *= window[i];
        test[i] *= 32768;
    }

    qwqdsp_spectral::IppRealFFT fft;
    fft.Init(fft_size);
    float reim[fft_size+2];
    fft.FFT(test, reim);

    float psd[fft_bins];
    for (int i = 0; i < fft_size+2; i+=2) {
        float re = reim[i];
        float im = reim[i + 1];
        psd[i / 2] = (re * re + im * im);
    }
}
