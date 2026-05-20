#include <qwqdsp/filter/cheby_hilbert.hpp>
#include <qwqdsp/spectral/complex_fft.hpp>
#include <qwqdsp/convert.hpp>

int main() {
    qwqdsp_filter::ChebyHilbert dsp;
    dsp.Reset();

    std::complex<float> y[4096];
    for (int i = 0; i < 4096; ++i) {
        y[i] = dsp.Tick(i == 0 ? 1.0f : 0.0f);
    }

    qwqdsp_spectral::ComplexFFT fft;
    fft.Init(4096);

    float gains[4096];
    fft.FFTGainPhase(y, gains);
    for (int i = 0; i < 4096; ++i) {
        gains[i] = qwqdsp::convert::Gain2Db<-150.0f>(gains[i]);
    }
}
