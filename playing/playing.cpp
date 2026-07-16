#include <qwqdsp/oscillator/vic_sine_osc.hpp>
#include <qwqdsp/spectral/real_fft_adv.hpp>
#include <qwqdsp/window/blackman.hpp>


int main() {
    qwqdsp_spectral::RealFftAdv fft0;
    fft0.Init(256);
    qwqdsp_spectral::RealFftAdv fft1;
    fft1.Init(512);
    qwqdsp_spectral::RealFftAdv fft2;
    fft2.Init(1024);

    float buffer[1024];
    qwqdsp_oscillator::VicSineOsc osc;
    osc.SetFreq(20.5f, 1024.0f);
    for (auto& x : buffer) {
        x = osc.Tick();
    }

    qwqdsp_window::Blackman::ApplyWindow(buffer, true);

    float g0[513]{};
    {
        float xin[1024]{};
        for (int i = 0; i < 1024; ++i) {
            xin[i % 256] += buffer[i];
        }
        fft0.FFTGainPhase({xin, 256}, {g0, 129});
    }

    float g1[513]{};
    {
        float xin[1024]{};
        for (int i = 0; i < 1024; ++i) {
            xin[i % 512] += buffer[i];
        }
        fft1.FFTGainPhase({xin, 512}, {g1, 257});
    }

    float g2[513]{};
    {
        float xin[1024]{};
        for (int i = 0; i < 1024; ++i) {
            xin[i % 1024] += buffer[i];
        }
        fft2.FFTGainPhase({xin, 1024}, {g2, 513});
    }
}
