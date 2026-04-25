#include <qwqdsp/fx/uniform_convolution.hpp>
#include <qwqdsp/oscillator/noise.hpp>
#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/filter/fir.hpp>
#include "../../playing/AudioFile.h"

static void ShouldBeDelay32() {
    qwqdsp_fx::UniformConvolution conv;
    conv.Init(32);
    conv.Reset();

    float ir[512]{1.0f};
    conv.SetIR(ir);

    float input[2048]{1.0f};

    float lag = 0.0f;
    for (int i = 0; i < 2048; ++i) {
        float x = input[i];

        float temp[1]{x + lag * 0.99f};
        conv.Process(temp);
        lag = temp[0];

        input[i] = temp[0];
    }
}

static void ShouldBeDelay33() {
    qwqdsp_fx::UniformConvolution conv;
    conv.Init(32);
    conv.Reset();

    float ir[512]{0.0f, 1.0f};
    conv.SetIR(ir);

    float input[2048]{1.0f};

    float lag = 0.0f;
    for (int i = 0; i < 2048; ++i) {
        float x = input[i];

        float temp[1]{x + lag * 0.99f};
        conv.Process(temp);
        lag = temp[0];

        input[i] = temp[0];
    }
}

static void ShouldBeDelay32Plus34() {
    qwqdsp_fx::UniformConvolution conv;
    conv.Init(32);
    conv.Reset();

    float ir[512]{0.5f, 0.0f, 0.5f};
    conv.SetIR(ir);

    float input[2048]{1.0f};

    float lag = 0.0f;
    for (int i = 0; i < 2048; ++i) {
        float x = input[i];

        float temp[1]{x + lag * 0.99f};
        conv.Process(temp);
        lag = temp[0];

        input[i] = temp[0];
    }
}

static void RandomFeedback() {
    qwqdsp_fx::UniformConvolution conv;
    conv.Init(1024);
    conv.Reset();

    float ir[4096];
    qwqdsp_oscillator::WhiteNoise noise;
    for (int i = 0; i < std::size(ir); ++i) {
        float x = noise.Next();;
        ir[i] = x;
    }

    qwqdsp_spectral::RealFFT fft;
    fft.Init(std::size(ir));
    float gains[std::size(ir)/2 + 1];
    fft.FFTGainPhase(ir, gains);

    float max_g = *std::max_element(gains, gains + std::size(gains));
    for (int i = 0; i < std::size(ir); ++i) {
        ir[i] /= max_g;
    }

    conv.SetIR(ir);

    float input[8000*10]{1.0f};

    float lag = 0.0f;
    for (int i = 0; i < std::size(input); ++i) {
        float x = input[i];

        float temp[1]{x + lag * 0.99f};
        conv.Process(temp);
        lag = temp[0];

        input[i] = temp[0];
    }

    AudioFile<float> file;
    file.setBitDepth(32);
    file.setNumChannels(1);
    file.setSampleRate(8000);
    file.samples.resize(1);
    file.samples.front().resize(std::size(input));
    std::copy_n(input, std::size(input), file.samples.front().begin());
    [[maybe_unused]] bool succ = file.save("random_fb.wav");
}

class Delay {
public:
    static constexpr int buffer_size = 16384;

    void Push(float x) {
        buffer_[wpos_] = x;
        wpos_ = Wrap(wpos_ + 1);
    }

    float GetBeforePush(int delay) {
        return buffer_[Wrap(wpos_ - delay)];
    }
private:
    static int Wrap(int pos) {
        pos %= buffer_size;
        if (pos < 0) {
            pos += buffer_size;
        }
        return pos;
    }

    float buffer_[buffer_size]{};
    int wpos_{};
};

static void RandomPhase() {
    qwqdsp_fx::UniformConvolution conv;
    conv.Init(32);
    conv.Reset();

    constexpr int ir_size = 2048;
    constexpr int bins = ir_size / 2 + 1;
    float gains[bins];
    float phases[bins];
    float ir[ir_size];
    std::fill(gains, gains + bins, 1.0f);
    std::fill(phases, phases + bins, 1.0f);

    qwqdsp_spectral::RealFFT fft;
    fft.Init(ir_size);
    fft.IFFTGainPhase(ir, gains, phases);
    conv.SetIR(ir);

    constexpr int pad_mul = 2;
    float pad_ir[ir_size * pad_mul]{};
    std::copy(ir, ir + ir_size, pad_ir);
    fft.Init(ir_size * pad_mul);
    constexpr int pad_bins = ir_size * pad_mul / 2 + 1;
    float gains2[pad_bins];
    fft.FFTGainPhase(pad_ir, gains2);

    float max_g = *std::max_element(gains2, gains2 + pad_bins);
    for (int i = 0; i < ir_size; ++i) {
        ir[i] /= max_g;
    }

    qwqdsp_filter::FIRTranspose fir;
    fir.SetCoeff([&ir](std::vector<float>& coeff) {
        coeff.resize(ir_size);
        std::copy(ir, ir + ir_size, coeff.begin());
    });

    float input[2048]{1.0f};

    float lag = 0.0f;
    // Delay delay;
    for (int i = 0; i < 2048; ++i) {
        float x = input[i];

        float temp[1]{x + lag * 0.99f};
        // float temp[1]{x + delay.GetBeforePush(110) * 0.99f};
        conv.Process(temp);
        fir.Process(temp);
        lag = temp[0];
        // delay.Push(temp[0]);

        input[i] = temp[0];
    }
}

int main() {
    // ShouldBeDelay32();
    // ShouldBeDelay33();
    // ShouldBeDelay32Plus34();
    // RandomFeedback();
    RandomPhase();
}
