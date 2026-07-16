#include "AudioFile.h"
#include "work_dir.hpp"
#include <iostream>
#include <qwqdsp/filter/fir.hpp>
#include <qwqdsp/fx/uniform_convolution.hpp>
#include <qwqdsp/oscillator/noise.hpp>
#include <qwqdsp/spectral/complex_fft_adv.hpp>
#include <qwqdsp/spectral/real_fft_adv.hpp>
#include <qwqdsp/window/hann.hpp>
#include <random>
#include <vector>

constexpr int NextPow2(int n) {
    int v = n - 1;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return v + 1;
}

template <int ir_size>
void MinPhase(std::span<float, ir_size> x) {
    constexpr int pad_ratio = 16;
    constexpr int pad_size = ir_size * pad_ratio;

    static float fir_pad[pad_size]{};
    std::copy_n(x.begin(), ir_size, fir_pad);

    qwqdsp_spectral::ComplexFftAdv fft;
    fft.Init(pad_size);
    constexpr size_t num_bins = fft.NumBins(pad_size);
    static float gains[num_bins];
    fft.FFTGainPhase(fir_pad, gains);

    static float log_gains[num_bins];
    for (size_t i = 0; i < num_bins; ++i) {
        log_gains[i] = std::log(gains[i] + 1e-18f);
    }

    static float phases[num_bins]{};
    static float ir[pad_size];
    fft.IFFT(ir, log_gains, phases);
    ir[0] = 0;
    ir[num_bins / 2] = 0;
    for (size_t i = num_bins / 2 + 1; i < num_bins; ++i) {
        ir[i] = -ir[i];
    }

    fft.FFT(ir, log_gains, phases);
    fft.IFFTGainPhase(ir, gains, phases);

    for (size_t i = 0; i < ir_size; ++i) {
        x[i] = ir[i];
    }
}

template <int ir_size>
void RandomPhaseIr(std::span<float, ir_size> x, float ripple_db = 0.1f) {
    qwqdsp_spectral::RealFftAdv fft;
    fft.Init(ir_size);

    constexpr int bins = fft.NumBins(ir_size);

    static float ir[ir_size];
    static float gains[bins];
    static float phases[bins];

    float max_phase_change = std::acos(2 * std::pow(10.0f, -ripple_db / 20.0f) - 1.0f);

    float bin_phase = 0.0f;
    static std::uniform_real_distribution<float> dist{-max_phase_change, max_phase_change};
    static std::default_random_engine rng{std::random_device{}()};
    for (int i = 0; i < bins; ++i) {
        gains[i] = 1.0f;
        phases[i] = bin_phase;
        bin_phase += dist(rng);
    }

    fft.IFFTGainPhase(ir, gains, phases);
    fft.TimeDomainShift(ir);
    qwqdsp_window::Hann::ApplyWindow(ir, false);

    // MinPhase<ir_size>(ir);

    fft.Init(ir_size * 2);
    static float pad[ir_size * 2]{};
    std::copy_n(ir, ir_size, pad);

    constexpr int pad_bins = fft.NumBins(ir_size * 2);
    static float pad_g[pad_bins];
    fft.FFTGainPhase(pad, pad_g);

    float max_g = *std::max_element(pad_g, pad_g + pad_bins);
    for (int i = 0; i < ir_size; ++i) {
        x[i] = max_g * ir[i];
    }
}

template <int ir_size>
void RandomPhaseIr2(std::span<float, ir_size> x, float ripple_db = 0.1f) {
    std::vector<float> ir;
    ir.resize(ir_size);
    std::vector<float> gains;
    gains.resize(ir_size);
    std::vector<float> phases;
    phases.resize(ir_size);

    float max_phase_change = std::acos(2 * std::pow(10.0f, -ripple_db / 20.0f) - 1.0f);

    float bin_phase = 0.0f;
    static std::uniform_real_distribution<float> dist{-max_phase_change, max_phase_change};
    static std::default_random_engine rng{std::random_device{}()};
    for (int i = 0; i < ir_size; ++i) {
        gains[i] = 1.0f;
        phases[i] = bin_phase;
        bin_phase += dist(rng);
    }

    for (int n = 0; n < ir_size; ++n) {
        ir[n] = 0.0f;
        for (int k = 0; k < ir_size; ++k) {
            ir[n] += gains[k] * std::cos(2.0f * std::numbers::pi_v<float> * k * n / ir_size + phases[k]);
        }
    }
    for (int i = 0; i < ir_size; ++i) {
        ir[i] /= ir_size;
    }
    {
        std::vector<float> tmp;
        tmp.resize(ir_size);
        std::copy_n(ir.begin(), ir_size / 2, tmp.begin());
        for (size_t i = 0; i < ir_size / 2; ++i) {
            ir[i] = ir[i + ir_size / 2];
        }
        std::copy_n(tmp.begin(), ir_size / 2, ir.begin() + ir_size / 2);
    }
    qwqdsp_window::Hann::ApplyWindow(ir, false);

    // norm
    int fft_size = NextPow2(ir_size) * 2;

    qwqdsp_spectral::RealFftAdv fft;
    fft.Init(fft_size);
    std::vector<float> pad;
    pad.resize(fft_size);
    std::copy_n(ir.begin(), ir_size, pad.begin());

    int pad_bins = fft.NumBins(fft_size);
    std::vector<float> pad_g;
    pad_g.resize(pad_bins);
    fft.FFTGainPhase(pad, pad_g);

    float max_g = *std::max_element(pad_g.begin(), pad_g.end());
    for (int i = 0; i < ir_size; ++i) {
        x[i] = ir[i] / max_g;
    }
}

static void RandomFeedback() {
    constexpr int block_size = 512;
    constexpr int ir_len = 16384;
    constexpr int audio_sr = 48000;
    constexpr int audio_sec = 20;
    constexpr int audio_samples = audio_sr * audio_sec;
    constexpr float ir_decay = 10;

    qwqdsp_fx::UniformConvolution conv;
    conv.Init(block_size);
    conv.Reset();

    static float ir[ir_len];
    RandomPhaseIr<ir_len>(ir, 10.0f);

    conv.SetIR(ir);

    static float input[audio_samples]{1.0f};
    float loop_size = block_size + ir_len / 2.0f;
    // float loop_size = block_size;
    float decay_size = ir_decay * audio_sr;
    float decay = std::pow(1e-3f, loop_size / decay_size);

    float lag = 0.0f;
    for (int i = 0; i < audio_samples; ++i) {
        float x = input[i];

        float temp[1]{x + lag * decay};
        conv.Process(temp);
        lag = temp[0];

        input[i] = temp[0];
    }

    AudioFile<float> file;
    file.setBitDepth(32);
    file.setNumChannels(1);
    file.setSampleRate(audio_sr);
    file.samples.resize(1);
    file.samples.front().resize(std::size(input));
    std::copy_n(input, std::size(input), file.samples.front().begin());
    file.save(qwqdsp_support::OutputFile("random_feedback.wav"));
}

static void RandomFeedback2() {
    constexpr int block_size = 512;
    constexpr int ir_len = 1145;
    constexpr int audio_sr = 48000;
    constexpr int audio_sec = 20;
    constexpr int audio_samples = audio_sr * audio_sec;
    constexpr float ir_decay = 10;

    qwqdsp_fx::UniformConvolution conv;
    conv.Init(block_size);
    conv.Reset();

    static float ir[ir_len];
    RandomPhaseIr2<ir_len>(ir, 1.0f);

    conv.SetIR(ir);

    static float input[audio_samples]{1.0f};
    float loop_size = block_size + ir_len / 2.0f;
    // float loop_size = block_size;
    float decay_size = ir_decay * audio_sr;
    float decay = std::pow(1e-3f, loop_size / decay_size);

    float lag = 0.0f;
    for (int i = 0; i < audio_samples; ++i) {
        float x = input[i];

        float temp[1]{x + lag * decay};
        conv.Process(temp);
        lag = temp[0];

        input[i] = temp[0];
    }

    AudioFile<float> file;
    file.setBitDepth(32);
    file.setNumChannels(1);
    file.setSampleRate(audio_sr);
    file.samples.resize(1);
    file.samples.front().resize(std::size(input));
    std::copy_n(input, std::size(input), file.samples.front().begin());
    file.save(qwqdsp_support::OutputFile("random_feedback2.wav"));
}

// -------------------- conv net test --------------------

template <int max_delay>
class Delay {
public:
    void Push(float x) {
        buffer_[wpos_] = x;
        wpos_ = Wrap(wpos_ + 1);
    }

    float GetBeforePush(int delay) {
        return buffer_[Wrap(wpos_ - delay)];
    }
private:
    static int Wrap(int pos) {
        pos %= max_delay;
        if (pos < 0) {
            pos += max_delay;
        }
        return pos;
    }

    float buffer_[max_delay]{};
    int wpos_{};
};

class ConvNet {
public:
    static constexpr int kOrder = 4;
    static constexpr int kBlockSize = 128;
    static constexpr float kRippleSize = 10.0f;

    void Init(float t60, int sr) {
        int base_delays[] = {1384, 1512, 1674, 1265};
        constexpr int ir_size[] = {1384 * 2, 1512 * 2, 1674 * 2, 1265 * 2};
        constexpr float riiple[] = {0.1f, 0.1f, 0.1f, 0.1f};

        for (int i = 0; i < kOrder; ++i) {
            int conv_delay = kBlockSize + ir_size[i] / 2;
            delay_val_[i] = std::max(base_delays[i] - conv_delay, 0);
        }

        for (int i = 0; i < kOrder; ++i) {
            int true_delay = delay_val_[i] + ir_size[i] / 2;
            decay_[i] = std::pow(0.001f, static_cast<float>(true_delay) / (t60 * sr));
        }

        for (int i = 0; i < kOrder; ++i) {
            conv_[i].Init(kBlockSize);
        }

        static float ir[ir_size[0]];
        RandomPhaseIr2<ir_size[0]>(ir, riiple[0]);
        conv_[0].SetIR(ir);

        static float ir2[ir_size[1]];
        RandomPhaseIr2<ir_size[1]>(ir2, riiple[1]);
        conv_[1].SetIR(ir2);

        static float ir3[ir_size[2]];
        RandomPhaseIr2<ir_size[2]>(ir3, riiple[2]);
        conv_[2].SetIR(ir3);

        static float ir4[ir_size[3]];
        RandomPhaseIr2<ir_size[3]>(ir4, riiple[3]);
        conv_[3].SetIR(ir4);
    }

    float Tick(float x) {
        float temp[kOrder]{x, x, x, x};
        for (int i = 0; i < kOrder; ++i) {
            temp[i] = x + delay_[i].GetBeforePush(delay_val_[i]) * decay_[i];
        }

        for (int i = 0; i < kOrder; ++i) {
            float tmp[1]{temp[i]};
            conv_[i].Process(tmp);
            temp[i] = tmp[0];
        }

        float y{};
        for (int i = 0; i < kOrder; ++i) {
            y += temp[i];
        }

        constexpr float kHadamardNorm = 0.5f;
        const float x0 = temp[0];
        const float x1 = temp[1];
        const float x2 = temp[2];
        const float x3 = temp[3];

        temp[0] = (x0 + x1 + x2 + x3) * kHadamardNorm;
        temp[1] = (x0 - x1 + x2 - x3) * kHadamardNorm;
        temp[2] = (x0 + x1 - x2 - x3) * kHadamardNorm;
        temp[3] = (x0 - x1 - x2 + x3) * kHadamardNorm;

        for (int i = 0; i < kOrder; ++i) {
            delay_[i].Push(temp[i]);
        }

        return y;
    }
private:
    Delay<4096> delay_[kOrder];
    qwqdsp_fx::UniformConvolution conv_[kOrder];

    std::array<int, kOrder> delay_val_;
    std::array<float, kOrder> decay_;
};

void RandomFDN() {
    constexpr int audio_sr = 48000;
    constexpr int audio_sec = 20;
    constexpr int audio_samples = audio_sr * audio_sec;
    constexpr float ir_decay = 10;

    ConvNet dsp;
    dsp.Init(ir_decay, audio_sr);

    static float input[audio_samples]{1.0f};
    for (int i = 0; i < audio_samples; ++i) {
        input[i] = dsp.Tick(input[i]);
    }

    AudioFile<float> file;
    file.setBitDepth(32);
    file.setNumChannels(1);
    file.setSampleRate(audio_sr);
    file.samples.resize(1);
    file.samples.front().resize(std::size(input));
    std::copy_n(input, std::size(input), file.samples.front().begin());
    if (file.save(qwqdsp_support::OutputFile("random_FDN.wav"))) {
        std::cout << "save file" << std::endl;
    }
    else {
        std::cout << "save failed" << std::endl;
    }
}

int main() {
    RandomFeedback();
    RandomFDN();
    RandomFeedback2();
}
