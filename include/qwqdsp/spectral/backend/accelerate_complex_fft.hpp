#ifdef QWQDSP_HAVE_ACCELERATE
#pragma once
#include <Accelerate/Accelerate.h>
#include <bit>
#include <cassert>
#include <vector>

namespace qwqdsp_spectral {

class AccelerateComplexFFT {
public:
    ~AccelerateComplexFFT() {
        if (setup_) {
            vDSP_destroy_fftsetup(setup_);
        }
    }

    void Init(size_t fft_size) {
        assert(std::has_single_bit(fft_size));
        fft_size_ = fft_size;
        log2n_ = static_cast<int>(std::countr_zero(fft_size));
        setup_ = vDSP_create_fftsetup(log2n_, FFT_RADIX2);
        realp_.resize(fft_size);
        imagp_.resize(fft_size);
    }

    /**
     * @brief 复数 FFT，频率输出 [0, 2π)
     * @param in_real  size = fft_size
     * @param in_imag  size = fft_size
     * @param out_real size = fft_size
     * @param out_imag size = fft_size
     */
    void FFT(const float* in_real, const float* in_imag, float* out_real, float* out_imag) noexcept {
        std::copy_n(in_real, fft_size_, realp_.data());
        std::copy_n(in_imag, fft_size_, imagp_.data());
        DSPSplitComplex split = {realp_.data(), imagp_.data()};
        vDSP_fft_zip(setup_, &split, 1, log2n_, FFT_FORWARD);
        std::copy_n(realp_.data(), fft_size_, out_real);
        std::copy_n(imagp_.data(), fft_size_, out_imag);
    }

    /**
     * @brief 逆复数 FFT
     * @param in_real  size = fft_size
     * @param in_imag  size = fft_size
     * @param out_real size = fft_size
     * @param out_imag size = fft_size
     */
    void IFFT(const float* in_real, const float* in_imag, float* out_real, float* out_imag) noexcept {
        std::copy_n(in_real, fft_size_, realp_.data());
        std::copy_n(in_imag, fft_size_, imagp_.data());
        DSPSplitComplex split = {realp_.data(), imagp_.data()};
        vDSP_fft_zip(setup_, &split, 1, log2n_, FFT_INVERSE);
        // vDSP_fft_zip 的 round trip = N, 缩放 1/N 以匹配 IPP
        float scale = 1.0f / static_cast<float>(fft_size_);
        vDSP_vsmul(realp_.data(), 1, &scale, out_real, 1, fft_size_);
        vDSP_vsmul(imagp_.data(), 1, &scale, out_imag, 1, fft_size_);
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
