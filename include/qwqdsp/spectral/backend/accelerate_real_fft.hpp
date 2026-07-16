#ifdef QWQDSP_HAVE_ACCELERATE
#pragma once
#include <Accelerate/Accelerate.h>
#include <bit>
#include <cassert>
#include <vector>

namespace qwqdsp_spectral {

class AccelerateRealFFT {
public:
    ~AccelerateRealFFT() {
        if (setup_) {
            vDSP_destroy_fftsetup(setup_);
        }
    }

    void Init(size_t fft_size) {
        assert(std::has_single_bit(fft_size));
        fft_size_ = fft_size;
        log2n_ = static_cast<int>(std::countr_zero(fft_size));
        setup_ = vDSP_create_fftsetup(log2n_, FFT_RADIX2);
        realp_.resize(fft_size / 2);
        imagp_.resize(fft_size / 2);
    }

    /**
     * @brief 实数 FFT，输出 CCS 格式 [re, im] × num_bins
     * @param input  size = fft_size
     * @param output size = fft_size + 2
     */
    void FFT(const float* input, float* output) noexcept {
        DSPSplitComplex split = {realp_.data(), imagp_.data()};

        // Pack real input as even/odd into split complex
        vDSP_ctoz(reinterpret_cast<const DSPComplex*>(input), 1, &split, 1, fft_size_ / 2);

        // Forward real FFT
        vDSP_fft_zrip(setup_, &split, 1, log2n_, FFT_FORWARD);

        // Unpack split complex to CCS format (matching IPP convention)
        output[0] = realp_[0];
        output[1] = 0.0f;
        for (size_t i = 1; i < fft_size_ / 2; ++i) {
            output[2 * i] = realp_[i];
            output[2 * i + 1] = imagp_[i];
        }
        output[fft_size_] = imagp_[0];
        output[fft_size_ + 1] = 0.0f;
    }

    /**
     * @brief 逆实数 FFT
     * @param input  size = fft_size + 2, CCS 格式
     * @param output size = fft_size
     */
    void IFFT(const float* input, float* output) noexcept {
        DSPSplitComplex split = {realp_.data(), imagp_.data()};

        // CCS to split complex
        realp_[0] = input[0];
        imagp_[0] = input[fft_size_];
        for (size_t i = 1; i < fft_size_ / 2; ++i) {
            realp_[i] = input[2 * i];
            imagp_[i] = input[2 * i + 1];
        }

        // Inverse real FFT
        vDSP_fft_zrip(setup_, &split, 1, log2n_, FFT_INVERSE);

        // Unpack split complex to real signal (even/odd → interleaved)
        vDSP_ztoc(&split, 1, reinterpret_cast<DSPComplex*>(output), 1, fft_size_ / 2);

        // vDSP_fft_zrip round trip = 2N, scale to identity
        float scale = 1.0f / (2.0f * static_cast<float>(fft_size_));
        vDSP_vsmul(output, 1, &scale, output, 1, fft_size_);
    }

    size_t GetFFTSize() const noexcept {
        return fft_size_;
    }
private:
    size_t fft_size_{};
    int log2n_{};
    FFTSetup setup_{};
    std::vector<float> realp_;
    std::vector<float> imagp_;
};

} // namespace qwqdsp_spectral
#endif
