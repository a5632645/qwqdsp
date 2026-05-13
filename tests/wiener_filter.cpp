#include <qwqdsp/oscillator/noise.hpp>
#include <qwqdsp/spectral/ipp_real_fft.hpp>
#include <qwqdsp/fx/uniform_convolution.hpp>
#include "../playing/AudioFile.h"

constexpr float noise_gain = 1e-3f;

static std::vector<float> PowerSpectrum(const std::vector<float>& x) {
    const size_t n = x.size();
    qwqdsp_spectral::IppRealFFT fft;
    fft.Init(n);
    std::vector<float> fft_out(n + 2);
    fft.FFT(x.data(), fft_out.data());

    std::vector<float> power(n / 2 + 1);
    power[0] = fft_out[0] * fft_out[0];
    for (size_t i = 1; i < n / 2; ++i) {
        float re = fft_out[2 * i];
        float im = fft_out[2 * i + 1];
        power[i] = re * re + im * im;
    }
    power[n / 2] = fft_out[1] * fft_out[1];

    // float max = *std::max_element(power.begin(), power.end());
    // for (auto& x : power) {
    //     x /= max;
    // }

    return power;
}

static std::vector<float> PadTo2(const std::vector<float>& x) {
    int n = x.size();
    int size = 1;
    while (size < n) {
        size *= 2;
    }

    std::vector<float> r;
    r.resize(size);
    std::copy(x.begin(), x.end(), r.begin());
    return r;
}

static std::vector<float> WienerFilter(const std::vector<float>& x, const std::vector<float>& psd_x, const std::vector<float>& psd_noise) {
    std::vector<float> h_filter; h_filter.reserve(psd_x.size());
    for (size_t i = 0; i < psd_x.size(); ++i) {
        float a = psd_x[i];
        float b = psd_noise[i];
        float h = a / (a + b);
        if (!std::isnormal(h)) {
            h = 0;
        }
        h_filter.push_back(h);
    }

    auto pad_x = PadTo2(x);
    qwqdsp_spectral::IppRealFFT fft;
    fft.Init(pad_x.size());
    std::vector<float> pad_x_fft; pad_x_fft.resize(pad_x.size() + 2);
    fft.FFT(pad_x.data(), pad_x_fft.data());

    for (int i = 0; i < pad_x_fft.size(); i += 2) {
        float g = h_filter[i / 2];
        pad_x_fft[i] *= g;
        pad_x_fft[i + 1] *= g;
    }

    fft.IFFT(pad_x_fft.data(), pad_x.data());
    return pad_x;
}

static void Test1() {
    AudioFile<float> file;
    file.load(R"(C:\Users\Kawai\Music\wormhole.wav)");
    auto& x = file.samples.front();

    std::vector<float> noise; noise.reserve(x.size());
    std::vector<float> mix; mix.reserve(x.size());
    qwqdsp_oscillator::WhiteNoise noise_gen;
    for (size_t i = 0; i < x.size(); ++i) {
        float n = noise_gen.Next() * noise_gain;
        noise.push_back(n);
        mix.push_back(n + x[i]);
    }

    auto pad_x = PadTo2(x);
    auto pad_noise = PadTo2(noise);
    auto psd_x = PowerSpectrum(pad_x);
    auto psd_noise = PowerSpectrum(pad_noise);
    auto out = WienerFilter(mix, psd_x, psd_noise);
    
    file.setNumSamplesPerChannel(x.size());
    file.setNumChannels(1);
    file.samples.front() = x;
    file.save("1/x.wav");

    file.setNumSamplesPerChannel(noise.size());
    file.setNumChannels(1);
    file.samples.front() = noise;
    file.save("1/noise.wav");

    file.setNumSamplesPerChannel(out.size());
    file.setNumChannels(1);
    file.samples.front() = out;
    file.save("1/wiener_out.wav");

    file.setNumSamplesPerChannel(mix.size());
    file.setNumChannels(1);
    file.samples.front() = mix;
    file.save("1/wiener_in.wav");
}

static std::vector<float> WienerFilter2(const std::vector<float>& r, const std::vector<float>& psd_r, const std::vector<float>& psd_noise) {
    std::vector<float> h_filter; h_filter.reserve(psd_r.size());
    for (size_t i = 0; i < psd_r.size(); ++i) {
        float a = psd_r[i];
        float b = psd_noise[i];
        float h = (std::sqrt(std::max(0.0f, a * a - b * b))) / a;
        if (!std::isnormal(h)) {
            h = 0;
        }
        h_filter.push_back(h);
    }

    auto pad_x = PadTo2(r);
    qwqdsp_spectral::IppRealFFT fft;
    fft.Init(pad_x.size());
    std::vector<float> pad_x_fft; pad_x_fft.resize(pad_x.size() + 2);
    fft.FFT(pad_x.data(), pad_x_fft.data());

    for (int i = 0; i < pad_x_fft.size(); i += 2) {
        float g = h_filter[i / 2];
        pad_x_fft[i] *= g;
        pad_x_fft[i + 1] *= g;
    }

    fft.IFFT(pad_x_fft.data(), pad_x.data());
    return pad_x;
}

static void Test2() {
    AudioFile<float> file;
    file.load(R"(C:\Users\Kawai\Music\wormhole.wav)");
    auto& x = file.samples.front();

    std::vector<float> noise; noise.reserve(x.size());
    std::vector<float> r; r.reserve(x.size());
    qwqdsp_oscillator::WhiteNoise noise_gen;
    for (size_t i = 0; i < x.size(); ++i) {
        float n = noise_gen.Next() * noise_gain;
        noise.push_back(n);
        r.push_back(n + x[i]);
    }

    auto pad_r = PadTo2(r);
    auto pad_noise = PadTo2(noise);
    auto psd_r = PowerSpectrum(pad_r);
    auto psd_noise = PowerSpectrum(pad_noise);
    auto out = WienerFilter2(r, psd_r, psd_noise);
    
    file.setNumSamplesPerChannel(x.size());
    file.setNumChannels(1);
    file.samples.front() = x;
    file.save("2/x.wav");

    file.setNumSamplesPerChannel(noise.size());
    file.setNumChannels(1);
    file.samples.front() = noise;
    file.save("2/noise.wav");

    file.setNumSamplesPerChannel(out.size());
    file.setNumChannels(1);
    file.samples.front() = out;
    file.save("2/wiener_out.wav");

    file.setNumSamplesPerChannel(r.size());
    file.setNumChannels(1);
    file.samples.front() = r;
    file.save("2/wiener_in.wav");
}

int main() {
    Test1();
    Test2();
}
