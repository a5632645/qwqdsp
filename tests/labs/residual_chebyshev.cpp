#include <cmath>
#include <complex>
#include <numbers>
#include <qwqdsp/convert.hpp>
#include <qwqdsp/filter/biquad.hpp>
#include <qwqdsp/spectral/real_fft.hpp>

int main() {
    constexpr int N = 4;
    constexpr int kIrLen = 2048;
    constexpr double ripple = 6.0;
    double cutoff = qwqdsp::convert::DigitalFreq2AnalogBilinear(0.5 * 0.5 * 0.4, 0.5) * std::numbers::pi_v<float> * 2;
    double eps = std::sqrt(std::pow(10.0, ripple / 10.0) - 1.0);
    double A = 1.0 / static_cast<double>(2 * N) * std::asinh(1.0 / eps);
    double k_re = std::sinh(A);
    double k_im = std::cosh(A);
    bool highpass = false;

    // s域
    double k = 1;
    std::complex<double> half_spoles[N];
    for (int i = 0; i < N; ++i) {
        double phi = (2.0 * static_cast<double>(i + 1) - 1.0) * std::numbers::pi_v<double> / static_cast<double>(4 * N);
        auto pole = std::complex{k_re * -std::sin(phi), k_im * std::cos(phi)};
        if (!highpass) {
            half_spoles[i] = pole * cutoff;
            k *= std::norm(pole);
            k *= cutoff;
            k *= cutoff;
        }
        else {
            half_spoles[i] = cutoff / pole;
        }
    }
    k /= std::sqrt(1.0 + eps * eps);

    // 双线性变换
    std::complex<double> half_zpoles[N];
    for (int i = 0; i < N; ++i) {
        half_zpoles[i] = (1.0 + half_spoles[i]) / (1.0 - half_spoles[i]);
        k /= std::real((1.0 - half_spoles[i]) * (1.0 - std::conj(half_spoles[i])));
    }

    // 部分分式分解
    std::complex<double> residual[N];
    for (int i = 0; i < N; ++i) {
        auto zpole = half_zpoles[i];

        std::complex<double> up = 1.0;
        std::complex<double> tmp_up = zpole + 1.0;
        if (highpass) {
            tmp_up = zpole - 1.0;
        }
        for (int j = 0; j < N; ++j) {
            up *= tmp_up;
            up *= tmp_up;
        }

        std::complex<double> down = 1.0;
        for (int j = 0; j < N; ++j) {
            if (i == j) {
                down *= (zpole - std::conj(half_zpoles[j]));
            }
            else {
                down *= (zpole - half_zpoles[j]);
                down *= (zpole - std::conj(half_zpoles[j]));
            }
        }

        residual[i] = up / down;
    }

    qwqdsp_filter::Biquad cascade[N];
    float ir[kIrLen]{1.0f};
    for (int i = 0; i < N; ++i) {
        cascade[i].Set(1, highpass ? -2 : 2, 1, -2 * half_zpoles[i].real(), std::norm(half_zpoles[i]));
    }
    for (int i = 0; i < kIrLen; ++i) {
        for (int j = 0; j < N; ++j) {
            ir[i] = cascade[j].Tick(ir[i]);
        }
        ir[i] *= k;
    }

    qwqdsp_filter::Biquad parallel[N];
    float ir2[kIrLen]{1.0f};
    for (int i = 0; i < N; ++i) {
        parallel[i].Set(0, 2 * residual[i].real(), -2 * std::real(residual[i] * std::conj(half_zpoles[i])),
                        -2 * half_zpoles[i].real(), std::norm(half_zpoles[i]));
    }
    for (int i = 0; i < kIrLen; ++i) {
        float x = ir2[i];
        float y = x;
        for (int j = 0; j < N; ++j) {
            y += parallel[j].Tick(x);
        }
        ir2[i] = y * k;
    }

    // fft
    float g1[kIrLen / 2 + 1];
    float g2[kIrLen / 2 + 1];
    qwqdsp_spectral::RealFFT fft;
    fft.Init(kIrLen);
    fft.FFTGainPhase(ir, g1);
    fft.FFTGainPhase(ir2, g2);

    for (int i = 0; i < kIrLen / 2 + 1; ++i) {
        g1[i] = qwqdsp::convert::Gain2Db<-200.0f>(g1[i]);
        g2[i] = qwqdsp::convert::Gain2Db<-200.0f>(g2[i]);
    }
}
