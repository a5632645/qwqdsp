#include <AudioFile.h>
#include <audio_ops.hpp>
#include <complex>
#include <qwqdsp/segement/analyze_synthsis_offline.hpp>
#include <qwqdsp/spectral/ipp_real_fft.hpp>
#include <qwqdsp/window/hann.hpp>
#include <qwqdsp/window/window.hpp>
#include <vector>
#include <work_dir.hpp>

// ------------------------------------------------------------
// FFT视作滤波器组，取共轭反转每个子带
// ------------------------------------------------------------

struct Piwarp1Class {
    Piwarp1Class(int fft_size) {
        window_.resize(fft_size);
        qwqdsp_window::BlackmanHarris::Window(window_, true);
        reim_.resize(fft_size + 2);
        fft_.Init(fft_size);
        buffer_.resize(fft_size);
    }

    void operator()(std::span<const float> in, std::span<float> out) noexcept {
        std::transform(window_.begin(), window_.end(), in.begin(), buffer_.begin(), std::multiplies<>{});
        // std::copy(in.begin(), in.end(), buffer_.begin());
        fft_.FFT(buffer_.data(), reim_.data());

        int reim_size = reim_.size();
        for (int i = 0; i < reim_size; i += 2) {
            reim_[i + 1] *= -2;
            reim_[i] *= 2;
        }

        fft_.IFFT(reim_.data(), out.data());
        // std::transform(window_.begin(), window_.end(), out.begin(), out.begin(), std::multiplies<>{});
    }

    std::vector<float> window_;
    std::vector<float> buffer_;
    qwqdsp_spectral::IppRealFFT fft_;
    std::vector<float> reim_;
};

static void Piwarp1() {
    AudioFile<float> file{qwqdsp_support::InputFile("carry.wav")};
    auto& x_vec = file.samples.front();

    constexpr int fft_size = 512;

    qwqdsp_segement::AnalyzeSynthsisOffline as;
    as.SetSize(fft_size);

    // 下采样倍数必须设置为弱COLA条件的倍数
    // 增加子带采样率会导致梳妆滤波
    as.SetInputHop(fft_size / 4);
    as.SetOutputHop(fft_size / 4);
    as.Reset();

    std::vector<float> y_vec;
    Piwarp1Class piwarp1{fft_size};
    as.Process(x_vec, y_vec, piwarp1);

    file.setNumSamplesPerChannel(y_vec.size());
    file.samples.front() = std::move(y_vec);
    file.setNumChannels(1);
    file.save(qwqdsp_support::OutputFile("piwarp1.wav"));
}

// ------------------------------------------------------------
// DFT视作滤波器组，取共轭反转每个子带
// ------------------------------------------------------------

struct Piwarp2Class {
    Piwarp2Class(int dft_size) {
        window_.resize(dft_size);
        qwqdsp_window::BlackmanNuttall::Window(window_, true);

        int fft_size = 1;
        while (fft_size < dft_size) {
            fft_size *= 2;
        }

        reim_.resize(fft_size + 2);
        fft_.Init(fft_size);
        buffer_.resize(fft_size);
    }

    void operator()(std::span<const float> in, std::span<float> out) noexcept {
        std::transform(window_.begin(), window_.end(), in.begin(), buffer_.begin(), std::multiplies<>{});
        // std::copy(in.begin(), in.end(), buffer_.begin());
        fft_.FFT(buffer_.data(), reim_.data());

        int reim_size = reim_.size();
        for (int i = 0; i < reim_size; i += 2) {
            reim_[i + 1] *= -2;
            reim_[i] *= 2;
        }

        fft_.IFFT(reim_.data(), out.data());
        // std::transform(window_.begin(), window_.end(), out.begin(), out.begin(), std::multiplies<>{});
    }

    std::vector<float> window_;
    std::vector<float> buffer_;
    qwqdsp_spectral::IppRealFFT fft_;
    std::vector<float> reim_;
};

static void Piwarp2() {
    AudioFile<float> file{qwqdsp_support::InputFile("carry.wav")};
    auto& x_vec = file.samples.front();

    constexpr int filter_banks = 256;
    constexpr int dft_size = filter_banks * 4 - 1;
    Piwarp2Class piwarp2{dft_size};

    qwqdsp_segement::AnalyzeSynthsisOffline as;
    as.SetSize(piwarp2.fft_.GetFFTSize());

    // 下采样倍数必须设置为弱COLA条件的倍数
    // 增加子带采样率会导致梳妆滤波
    as.SetInputHop(dft_size / 4);
    as.SetOutputHop(dft_size / 4);
    as.Reset();

    std::vector<float> y_vec;
    as.Process(x_vec, y_vec, piwarp2);

    file.setNumSamplesPerChannel(y_vec.size());
    file.samples.front() = std::move(y_vec);
    file.setNumChannels(1);
    file.save(qwqdsp_support::OutputFile("piwarp2.wav"));
}

// ------------------------------------------------------------
// DFT每个取共轭 等价于 实数信号反向时间
// ------------------------------------------------------------

struct Piwarp3Class {
    Piwarp3Class(int dft_size) {
        window_.resize(dft_size);
        qwqdsp_window::Hann::Window(window_, true);
    }

    void operator()(std::span<const float> in, std::span<float> out) noexcept {
        int len = in.size();
        for (int i = 0; i < in.size(); ++i) {
            out[i] = window_[i] * in[len - i - 1];
        }
    }

    std::vector<float> window_;
};

static void Piwarp3() {
    AudioFile<float> file{qwqdsp_support::InputFile("carry.wav")};
    // AudioFile<float> file{qwqdsp_support::SweepWav()};
    auto& x_vec = file.samples.front();

    constexpr int filter_banks = 50;

    qwqdsp_segement::AnalyzeSynthsisOffline as;
    as.SetSize(filter_banks);

    // 下采样倍数必须设置为弱COLA条件的倍数
    // 增加子带采样率会导致梳妆滤波
    as.SetInputHop(filter_banks / 2);
    as.SetOutputHop(filter_banks / 2);
    as.Reset();

    std::vector<float> y_vec;
    Piwarp3Class piwarp3{filter_banks};
    as.Process(x_vec, y_vec, piwarp3);

    file.setNumSamplesPerChannel(y_vec.size());
    file.samples.front() = std::move(y_vec);
    file.setNumChannels(1);
    file.save(qwqdsp_support::OutputFile("piwarp3.wav"));
}

int main() {
    Piwarp1();
    Piwarp2();
    Piwarp3();
}
