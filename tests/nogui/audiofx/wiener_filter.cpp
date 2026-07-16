#include "AudioFile.h"
#include "work_dir.hpp"
#include <qwqdsp/fx/uniform_convolution.hpp>
#include <qwqdsp/oscillator/blep_coeff.hpp>
#include <qwqdsp/oscillator/noise.hpp>
#include <qwqdsp/oscillator/polyblep.hpp>
#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/window/blackman.hpp>
#include <qwqdsp/window/hann.hpp>

constexpr float noise_gain = 1e-1f;

static std::vector<float> PowerSpectrum(const std::vector<float>& x) {
    const size_t n = x.size();
    qwqdsp_spectral::RealFFT fft;
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

// ----------------------------------------
// test1
// standard wiener filter use Noise and Voice
// ----------------------------------------

static std::vector<float> WienerFilter(const std::vector<float>& x, const std::vector<float>& psd_x,
                                       const std::vector<float>& psd_noise) {
    std::vector<float> h_filter;
    h_filter.reserve(psd_x.size());
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
    std::vector<float> pad_x_fft;
    pad_x_fft.resize(pad_x.size() + 2);
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
    file.load(qwqdsp_support::WormholeWav());
    auto& x = file.samples.front();

    std::vector<float> noise;
    noise.reserve(x.size());
    std::vector<float> mix;
    mix.reserve(x.size());
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

    file.setNumSamplesPerChannel(noise.size());
    file.setNumChannels(1);
    file.samples.front() = noise;
    file.save(qwqdsp_support::OutputFile("wiener_1_noise.wav"));

    file.setNumSamplesPerChannel(out.size());
    file.setNumChannels(1);
    file.samples.front() = out;
    file.save(qwqdsp_support::OutputFile("wiener_1_out.wav"));

    file.setNumSamplesPerChannel(mix.size());
    file.setNumChannels(1);
    file.samples.front() = mix;
    file.save(qwqdsp_support::OutputFile("wiener_1_in.wav"));
}

// ----------------------------------------
// test2
// Noise And Dirty signal
// ----------------------------------------

static std::vector<float> WienerFilter2(const std::vector<float>& r, const std::vector<float>& psd_r,
                                        const std::vector<float>& psd_noise) {
    std::vector<float> h_filter;
    h_filter.reserve(psd_r.size());
    for (size_t i = 0; i < psd_r.size(); ++i) {
        float a = psd_r[i];
        float b = psd_noise[i];
        // float h = (std::sqrt(std::max(0.0f, a * a - b * b))) / a;
        float h = (a - b) / (a + b);
        if (!std::isnormal(h)) {
            h = 0;
        }
        h_filter.push_back(h);
    }

    auto pad_x = PadTo2(r);
    qwqdsp_spectral::IppRealFFT fft;
    fft.Init(pad_x.size());
    std::vector<float> pad_x_fft;
    pad_x_fft.resize(pad_x.size() + 2);
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
    file.load(qwqdsp_support::WormholeWav());
    auto& x = file.samples.front();

    std::vector<float> noise;
    noise.reserve(x.size());
    std::vector<float> r;
    r.reserve(x.size());
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

    file.setNumSamplesPerChannel(noise.size());
    file.setNumChannels(1);
    file.samples.front() = noise;
    file.save(qwqdsp_support::OutputFile("wiener_2_noise.wav"));

    file.setNumSamplesPerChannel(out.size());
    file.setNumChannels(1);
    file.samples.front() = out;
    file.save(qwqdsp_support::OutputFile("wiener_2_out.wav"));

    file.setNumSamplesPerChannel(r.size());
    file.setNumChannels(1);
    file.samples.front() = r;
    file.save(qwqdsp_support::OutputFile("wiener_2_in.wav"));
}

// ----------------------------------------
// test3
// Wiener Filter For Vocoding
// result: Strange high frequency carry and Chirps
// ----------------------------------------

static void Test3() {
    AudioFile<float> file;
    file.load(qwqdsp_support::WormholeWav());
    auto& x = file.samples.front();

    constexpr int block_size = 2048;
    constexpr int hop_size = 1024;
    qwqdsp_oscillator::PolyBlep<qwqdsp_oscillator::blep_coeff::Triangle> osc;
    osc.SetFreq(88.0f, file.getSampleRate());

    qwqdsp_oscillator::WhiteNoise noise;

    int size = x.size();

    int wpos = 0;
    std::vector<float> synthsis;
    synthsis.reserve(x.size());
    std::vector<float> ola_buffer;
    ola_buffer.resize(block_size);
    std::vector<float> modu_buffer;
    modu_buffer.resize(block_size);
    std::vector<float> carr_buffer;
    carr_buffer.resize(block_size);
    while (size != 0) {
        for (int i = 0; i < (block_size - hop_size); ++i) {
            modu_buffer[i] = modu_buffer[i + hop_size];
        }
        for (int i = 0; i < (block_size - hop_size); ++i) {
            carr_buffer[i] = carr_buffer[i + hop_size];
        }

        std::fill_n(modu_buffer.begin() + (block_size - hop_size), hop_size, 0);
        int can_read = std::min(hop_size, size);
        std::copy_n(x.begin() + wpos, can_read, modu_buffer.begin() + (block_size - hop_size));
        wpos += can_read;
        size -= can_read;

        for (int i = 0; i < hop_size; ++i) {
            carr_buffer[i + (block_size - hop_size)] = osc.Sawtooth();
            // carr_buffer[i + (block_size - hop_size)] = noise.Next();
        }

        auto modu_win = modu_buffer;
        auto carr_win = carr_buffer;
        qwqdsp_window::Hann::ApplyWindow(modu_win, true);
        qwqdsp_window::Hann::ApplyWindow(carr_win, true);

        auto psd_modu = PowerSpectrum(modu_win);
        auto psd_carry = PowerSpectrum(carr_win);
        auto out = WienerFilter(carr_buffer, psd_modu, psd_carry);
        qwqdsp_window::Hann::ApplyWindow(out, true);
        for (int i = 0; i < block_size; ++i) {
            ola_buffer[i] += out[i];
        }

        std::copy_n(ola_buffer.begin(), hop_size, std::back_inserter(synthsis));
        for (int i = 0; i < (block_size - hop_size); ++i) {
            ola_buffer[i] = ola_buffer[i + hop_size];
        }
        for (int i = 0; i < hop_size; ++i) {
            ola_buffer[i + (block_size - hop_size)] = 0;
        }
    }

    file.setNumSamplesPerChannel(synthsis.size());
    file.setNumChannels(1);
    file.samples.front() = synthsis;
    file.save(qwqdsp_support::OutputFile("wiener_3_out.wav"));
}

static void Test4() {
    AudioFile<float> file;
    file.load(qwqdsp_support::WormholeWav());
    auto x = file.samples.front();

    file.load(qwqdsp_support::InputFile("carry.wav"));
    auto n = file.samples.front();

    constexpr int block_size = 512;
    constexpr int hop_size = 128;
    int size = std::min(x.size(), n.size());
    int wpos = 0;
    std::vector<float> synthsis;
    synthsis.reserve(size);
    std::vector<float> ola_buffer;
    ola_buffer.resize(block_size);
    std::vector<float> modu_buffer;
    modu_buffer.resize(block_size);
    std::vector<float> carr_buffer;
    carr_buffer.resize(block_size);
    while (size != 0) {
        for (int i = 0; i < (block_size - hop_size); ++i) {
            modu_buffer[i] = modu_buffer[i + hop_size];
        }
        for (int i = 0; i < (block_size - hop_size); ++i) {
            carr_buffer[i] = carr_buffer[i + hop_size];
        }

        std::fill_n(modu_buffer.begin() + (block_size - hop_size), hop_size, 0);
        int can_read = std::min(hop_size, size);
        std::copy_n(x.begin() + wpos, can_read, modu_buffer.begin() + (block_size - hop_size));

        std::fill_n(carr_buffer.begin() + (block_size - hop_size), hop_size, 0);
        std::copy_n(n.begin() + wpos, can_read, carr_buffer.begin() + (block_size - hop_size));

        wpos += can_read;
        size -= can_read;

        auto modu_win = modu_buffer;
        auto carr_win = carr_buffer;
        qwqdsp_window::Hann::ApplyWindow(modu_win, true);
        qwqdsp_window::Hann::ApplyWindow(carr_win, true);

        auto psd_modu = PowerSpectrum(modu_win);
        auto psd_carry = PowerSpectrum(carr_win);
        auto out = WienerFilter(carr_buffer, psd_modu, psd_carry);
        qwqdsp_window::Hann::ApplyWindow(out, true);
        for (int i = 0; i < block_size; ++i) {
            ola_buffer[i] += out[i];
        }

        std::copy_n(ola_buffer.begin(), hop_size, std::back_inserter(synthsis));
        for (int i = 0; i < (block_size - hop_size); ++i) {
            ola_buffer[i] = ola_buffer[i + hop_size];
        }
        for (int i = 0; i < hop_size; ++i) {
            ola_buffer[i + (block_size - hop_size)] = 0;
        }
    }

    file.setNumSamplesPerChannel(synthsis.size());
    file.setNumChannels(1);
    file.samples.front() = synthsis;
    file.save(qwqdsp_support::OutputFile("wiener_3_out2.wav"));
}

int main() {
    Test1();
    Test2();
    Test3();
    Test4();
}
