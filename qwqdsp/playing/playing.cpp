#include <numbers>
#include <complex>
#include "AudioFile.h"
#include "qwqdsp/adaptive/burg_lp.hpp"
#include "qwqdsp/convert.hpp"
#include "qwqdsp/spectral/real_fft.hpp"
#include "qwqdsp/window/helper.hpp"

float LPCBurg(const float x[]/*in wave data*/, size_t N, float a[]/*out LPC*/, size_t M)
{
    std::vector<float> b1 (N, 0);
    std::vector<float> b2 (N, 0);
    std::vector<float> aa (M, 0);

    // (3)

    float gain = 0.0;
    for(size_t j = 0; j < N; j++)
    {
        gain += x[j] * x[j];
    }

    gain /= N;
    if(gain <= 0)
    {
        return 0.0;    // warning empty
    }

    // (9)

    b1[0] = x[0];
    b2[N - 2] = x[N - 1];
    for(size_t j = 1; j < N - 1; j++)
    {
        b1[j] = b2[j - 1] = x[j];
    }

    for(size_t i = 0; ; ++i)
    {
        // (7)

        float num = 0.0, denum = 0.0;
        for(size_t j = 0; j < N - i - 1; j++)
        {
            num += b1[j] * b2[j];
            denum += b1[j] * b1[j] + b2[j] * b2[j];
        }

        if(denum <= 0)
        {
            return 0.0;    // warning ill-conditioned
        }

        a[i] = 2.0 * num / denum;

        // (10)

        gain *= 1.0 - a[i] * a[i];

        // (5)

        for(size_t j = 0; j+1 <= i; j++)
        {
            a[j] = aa[j] - a[i] * aa[i - j - 1];
        }

        if(i == M-1) break;

        // (8)  Watch out: i -> i+1

        for(size_t j = 0; j <= i; j++)
        {
            aa[j] = a[j];
        }

        for(size_t j = 0; j <= N - i - 2; j++)
        {
            b1[j] -= aa[i] * b2[j];
            b2[j] = b2[j + 1] - aa[i] * b1[j + 1];
        }
    }

    return gain;
}

        // std::copy(x.begin() + 1, x.end(), eb_.begin());
        // for (size_t kidx = 0; auto& k : latticek) {
        //     ++kidx;

        //     float up{};
        //     float down{};
        //     for (size_t i = 0; i < x.size() - kidx; ++i) {
        //         up += x[i] * eb_[i];
        //         down += x[i] * x[i];
        //         down += eb_[i] * eb_[i];
        //     }
        //     k = -2.0f * up / down;

        //     for (size_t i = 0; i < x.size() - kidx - 1; ++i) {
        //         float const upgo = x[i] + eb_[i] * k;
        //         float const downgo = eb_[i + 1] + x[i + 1] * k;
        //         x[i] = upgo;
        //         eb_[i] = downgo;
        //     }
        // }

int main() {
    AudioFile<float> speech;
    speech.load(R"(C:\Users\Kawai\Music\speech-slice.wav)");
    auto& x = speech.samples.front();
    auto xback = x;

    qwqdsp::adaptive::BurgLP burg;
    burg.Init(x.size());

    constexpr size_t numpoles = 35;

    float k[numpoles];
    burg.Process(x, k);

    std::array<float, numpoles + 1> upgoing{1};
    std::array<float, numpoles + 1> downgoing{1};
    for (size_t kidx = 0; kidx < numpoles; ++kidx) {
        for (size_t i = kidx + 1; i != 0; --i) {
            downgoing[i] = downgoing[i - 1];
        }
        downgoing[0] = 0;

        for (size_t i = 0; i < kidx + 2; ++i) {
            float up = upgoing[i] + k[kidx] * downgoing[i];
            float down = downgoing[i] + k[kidx] * upgoing[i];
            upgoing[i] = up;
            downgoing[i] = down;
        }
    }

    std::array<float, numpoles> tfs;
    for (size_t i = 0; i < numpoles; ++i) {
        for(size_t j = 0; j+1 <= i; j++)
        {
            k[j] = tfs[j] + k[i] * tfs[i - j - 1];
        }

        for(size_t j = 0; j <= i; j++)
        {
            tfs[j] = k[j];
        }
    }

    float mel_begin = 0;
    float mel_end = std::min<float>(20000, speech.getSampleRate() / 2);
    float gain[1024];
    for (size_t i = 0; i < 1024; ++i) {
        float freq = std::lerp(mel_begin, mel_end, i / 1024.0f);
        float omega = freq * std::numbers::pi_v<float> * 2 / speech.getSampleRate();

        auto z_responce = std::complex{1.0f, 0.0f};
        for (int i = 0; i < numpoles; ++i) {
            auto z = std::polar(1.0f, -omega * (i + 1));
            z_responce += upgoing[i + 1] * z;
        }
        z_responce = 1.0f / z_responce;

        gain[i] = std::abs(z_responce);
    }

    float as2[numpoles];
    LPCBurg(xback.data(), xback.size(), as2, numpoles);
    float gain2[1024];
    for (size_t i = 0; i < 1024; ++i) {
        float freq = std::lerp(mel_begin, mel_end, i / 1024.0f);
        float omega = freq * std::numbers::pi_v<float> * 2 / speech.getSampleRate();

        auto z_responce = std::complex{1.0f, 0.0f};
        for (int i = 0; i < numpoles; ++i) {
            auto z = std::polar(1.0f, -omega * (i + 1));
            z_responce -= as2[i] * z;
        }
        z_responce = 1.0f / z_responce;

        gain2[i] = std::abs(z_responce);
    }
}