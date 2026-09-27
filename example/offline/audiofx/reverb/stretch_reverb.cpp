#include <AudioFile.h>
#include <audio_ops.hpp>
#include <complex>
#include <numbers>
#include <qwqdsp/convert.hpp>
#include <qwqdsp/segement/analyze_synthsis_offline2.hpp>
#include <qwqdsp/spectral/real_fft.hpp>
#include <random>
#include <vector>
#include <work_dir.hpp>

class Stretcher {
public:
    Stretcher(int in_size, int out_size) {
        fft_.Init(out_size);
        in_pad_.resize(out_size);
        reim_.resize(fft_.FFTSize() + 2);

        win_.resize(out_size, 1.0f);
        for (int i = 0; i < 256; ++i) {
            win_[i] = i / 256.0f;
            win_[out_size - i - 1] = i / 256.0f;
        }
    }

    void operator()(std::span<const float> input, std::span<float> output) {
        std::fill(in_pad_.begin(), in_pad_.end(), float{});
        std::copy(input.begin(), input.end(), in_pad_.begin());
        fft_.FFT(in_pad_.data(), reim_.data());

        int num_data = fft_.FFTSize() + 2;
        for (int i = 0; i < num_data; i += 2) {
            std::complex cpx{reim_[i], reim_[i + 1]};
            float g = std::abs(cpx);
            float random_phase = dist_(random_);
            float re = g * std::cos(random_phase);
            float im = g * std::sin(random_phase);
            reim_[i] = re;
            reim_[i + 1] = im;
        }

        fft_.IFFT(reim_.data(), output.data());
        ApplyDecay(output);
    }

    void ApplyDecay(std::span<float> x_vec) {
        float g = 1.0f;
        float mul = qwqdsp::convert::Ms2DecayDb(10000, 48000, -60);
        for (int i = 0; auto& x : x_vec) {
            x *= g * win_[i++];
            g *= mul;
        }
        // for (int i = 0; auto& x : x_vec) {
        //     x *= win_[i++];
        // }
    }
private:
    qwqdsp_spectral::RealFFT fft_;
    std::vector<float> in_pad_;
    std::vector<float> reim_;
    std::default_random_engine random_;
    std::uniform_real_distribution<float> dist_{-std::numbers::pi_v<float>, std::numbers::pi_v<float>};

    std::vector<float> win_;
};

int main() {
    constexpr int grain_size = 2048;
    constexpr int stretch_ratio = 64;
    constexpr int stretch_grain_size = stretch_ratio * grain_size;

    AudioFile<float> file{qwqdsp_support::WormholeWav()};
    auto const& x_vec = file.samples.front();

    qwqdsp_segement::AnalyzeSynthsisOffline2<false> as;
    as.SetInputHop(grain_size / 2);
    as.SetInputSize(grain_size);
    as.SetOutputSize(stretch_grain_size);
    as.SetOutputHop(grain_size / 2);

    std::vector<float> y_vec;
    Stretcher stretcher{grain_size, stretch_grain_size};
    as.Process(x_vec, y_vec, stretcher);
    qwqdsp_support::AudioOps::Normalize(y_vec);
    file.samples.front() = y_vec;

    y_vec.clear();
    as.Process(x_vec, y_vec, stretcher);
    qwqdsp_support::AudioOps::Normalize(y_vec);
    file.samples.push_back(y_vec);

    file.setNumSamplesPerChannel(y_vec.size());
    file.setNumChannels(2);
    file.save(qwqdsp_support::OutputFile("stretch_reverb.wav"));
}
