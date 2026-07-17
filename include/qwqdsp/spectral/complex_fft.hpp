#pragma once

#if defined(QWQDSP_HAVE_IPP)
#include "backend/ipp_complex_fft.hpp"
#elif defined(QWQDSP_HAVE_ACCELERATE)
#include "backend/accelerate_complex_fft.hpp"
#else
#include "backend/oouras_complex_fft.hpp"
#endif

namespace qwqdsp_spectral {

/**
 * @brief 根据编译配置自动选择复数 FFT 后端（IPP > Accelerate > OOURA）
 *        提供裸 FFT/IFFT 接口，频率输出 [0, 2π)
 *        高级 span 封装请使用 ComplexFftAdv
 */
class ComplexFFT {
public:
    void Init(size_t fft_size) {
        backend_.Init(fft_size);
    }

    void FFT(const float* in_real, const float* in_imag, float* out_real, float* out_imag) noexcept {
        backend_.FFT(in_real, in_imag, out_real, out_imag);
    }

    void IFFT(const float* in_real, const float* in_imag, float* out_real, float* out_imag) noexcept {
        backend_.IFFT(in_real, in_imag, out_real, out_imag);
    }

    size_t NumBins() const noexcept {
        return fft_size_;
    }

    static constexpr size_t NumBins(size_t fft_size) noexcept {
        return fft_size;
    }

    size_t FFTSize() const noexcept {
        return fft_size_;
    }
private:
    size_t fft_size_{};
#if defined(QWQDSP_HAVE_IPP)
    IppComplexFFT backend_;
#elif defined(QWQDSP_HAVE_ACCELERATE)
    AccelerateComplexFFT backend_;
#else
    OourasComplexFFT backend_;
#endif
};

} // namespace qwqdsp_spectral
