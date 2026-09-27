#include <AudioFile.h>
#include <audio_ops.hpp>
#include <cmath>
#include <numbers>
#include <qwqdsp/segement/analyze_synthsis_offline2.hpp>
#include <qwqdsp/spectral/real_fft.hpp>
#include <random>
#include <vector>
#include <work_dir.hpp>

class Stretcher {
public:
    Stretcher(int in_size, int out_size, float sample_rate, float rt60_ms) {
        fft_.Init(out_size);
        in_pad_.resize(out_size);
        reim_.resize(fft_.FFTSize() + 2);

        // 逐 bin 一阶平滑器的状态
        smoothed_mag_.resize(fft_.NumBins(), 0.0f);

        // 每 hop 的衰减因子：经过 rt60_ms 后幅度衰减 -60 dB
        float hop_size = in_size / 2.0f;
        float rt60_samples = rt60_ms * sample_rate / 1000.0f;
        decay_factor_ = std::pow(10.0f, -3.0f * hop_size / rt60_samples);

        // 淡入淡出窗，用于重叠相加
        win_.resize(out_size, 1.0f);
        for (int i = 0; i < 256; ++i) {
            win_[i] = i / 256.0f;
            win_[out_size - i - 1] = i / 256.0f;
        }
    }

    void operator()(std::span<const float> input, std::span<float> output) {
        // Zero-pad 输入
        std::fill(in_pad_.begin(), in_pad_.end(), float{});
        std::copy(input.begin(), input.end(), in_pad_.begin());

        // FFT
        fft_.FFT(in_pad_.data(), reim_.data());

        // 逐 bin：一阶平滑滤波器处理幅度，相位随机化
        int num_bins = static_cast<int>(fft_.NumBins());
        for (int i = 0; i < num_bins; ++i) {
            float re = reim_[2 * i];
            float im = reim_[2 * i + 1];
            float current_mag = std::sqrt(re * re + im * im);

            // Attack = 0（瞬时跟踪），Release 用一阶平滑衰减
            if (current_mag <= smoothed_mag_[i]) {
                smoothed_mag_[i] = smoothed_mag_[i] * decay_factor_ + current_mag * (1.0f - decay_factor_);
            }
            else {
                smoothed_mag_[i] = current_mag;
            }

            // 相位随机化
            float random_phase = dist_(random_) + std::atan2(im, re);
            reim_[2 * i] = smoothed_mag_[i] * std::cos(random_phase);
            reim_[2 * i + 1] = smoothed_mag_[i] * std::sin(random_phase);
        }

        // IFFT
        fft_.IFFT(reim_.data(), output.data());

        // 应用窗函数
        for (int i = 0; i < static_cast<int>(output.size()); ++i) {
            output[i] *= win_[i];
        }
    }
private:
    qwqdsp_spectral::RealFFT fft_;
    std::vector<float> in_pad_;
    std::vector<float> reim_;
    std::vector<float> smoothed_mag_;
    float decay_factor_;

    std::default_random_engine random_;
    std::uniform_real_distribution<float> dist_{-std::numbers::pi_v<float>, std::numbers::pi_v<float>};
    std::vector<float> win_;
};

int main() {
    constexpr int grain_size = 2048;
    constexpr int stretch_ratio = 32;
    constexpr int stretch_grain_size = stretch_ratio * grain_size;

    AudioFile<float> file{qwqdsp_support::WormholeWav()};
    auto const& x_vec = file.samples.front();

    qwqdsp_segement::AnalyzeSynthsisOffline2<false> as;
    as.SetInputHop(grain_size / 2);
    as.SetInputSize(grain_size);
    as.SetOutputSize(stretch_grain_size);
    as.SetOutputHop(grain_size / 2);

    std::vector<float> y_vec;
    Stretcher stretcher{grain_size, stretch_grain_size, 48000.0f, 1000.0f};
    as.Process(x_vec, y_vec, stretcher);
    qwqdsp_support::AudioOps::Normalize(y_vec);
    file.samples.front() = y_vec;

    y_vec.clear();
    as.Process(x_vec, y_vec, stretcher);
    qwqdsp_support::AudioOps::Normalize(y_vec);
    file.samples.push_back(y_vec);

    file.setNumSamplesPerChannel(y_vec.size());
    file.setNumChannels(2);
    file.save(qwqdsp_support::OutputFile("spectral_reverb.wav"));
}
