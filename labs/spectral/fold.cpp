#include <qwqdsp/convert.hpp>
#include <qwqdsp/oscillator/vic_sine_osc.hpp>
#include <qwqdsp/spectral/real_fft_adv.hpp>
#include <qwqdsp/window/window.hpp>
#include <array>

int main() {
    float x_vec[1024];

    qwqdsp_oscillator::VicSineOsc osc;
    osc.SetFreq(5.5, 1024);
    for (auto& x : x_vec) {
        x = osc.Tick();
    }

    float win_vec[1024];
    qwqdsp_window::Blackman::Window(win_vec, true);
    qwqdsp_window::Helper::Normalize(win_vec);
    for (int i = 0; i < 1024; ++i) {
        x_vec[i] *= win_vec[i];
    }

    {
        // ----- fold version -----
        qwqdsp_spectral::RealFftAdv fft;
        fft.Init(256);

        float fold_vec[256];
        for (int i = 0; i < 256; ++i) {
            fold_vec[i] = x_vec[i] + x_vec[i + 256] + x_vec[i + 512] + x_vec[i + 768];
        }

        float fold_g_vec[129];
        fft.FFTGainPhase(fold_vec, fold_g_vec);

        for (auto& g : fold_g_vec) {
            g = qwqdsp::convert::Gain2Db<-100.0f>(g);
        }

        float equal_len_fold_g_vec[513];
        for (int i = 0; i < 513; ++i) {
            equal_len_fold_g_vec[i] = fold_g_vec[i / 4];
        }

        // ----- full version -----
        float full_g_vec[513];
        fft.Init(1024);
        fft.FFTGainPhase(x_vec, full_g_vec);
        for (auto& g : full_g_vec) {
            g = qwqdsp::convert::Gain2Db<-100.0f>(g);
        }

        // 这看起来是下采样的版本
        std::ignore = full_g_vec;
        std::ignore = equal_len_fold_g_vec;
    }
}
