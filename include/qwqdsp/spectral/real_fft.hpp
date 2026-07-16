#pragma once

#if defined(QWQDSP_HAVE_IPP)
#include "backend/ipp_real_fft.hpp"
#elif defined(QWQDSP_HAVE_ACCELERATE)
#include "backend/accelerate_real_fft.hpp"
#else
#include "backend/oouras_real_fft.hpp"
#endif

namespace qwqdsp_spectral {

/**
 * @brief 根据编译配置自动选择 FFT 后端（IPP > Accelerate > OOURA）
 *        提供裸 CCS 接口（Init / FFT / IFFT）
 *        高级 span 封装请使用 RealFftAdv
 */
class RealFFT {
public:
    void Init(size_t fft_size) {
        backend_.Init(fft_size);
    }

    void FFT(const float* input, float* output) noexcept {
        backend_.FFT(input, output);
    }

    void IFFT(const float* input, float* output) noexcept {
        backend_.IFFT(input, output);
    }

    size_t GetFFTSize() const noexcept {
        return backend_.GetFFTSize();
    }

    size_t NumBins() const noexcept {
        return GetFFTSize() / 2 + 1;
    }

    static constexpr size_t NumBins(size_t fft_size) noexcept {
        return fft_size / 2 + 1;
    }

    size_t FFTSize() const noexcept {
        return GetFFTSize();
    }

private:
#if defined(QWQDSP_HAVE_IPP)
    IppRealFFT backend_;
#elif defined(QWQDSP_HAVE_ACCELERATE)
    AccelerateRealFFT backend_;
#else
    OourasRealFFT backend_;
#endif
};

} // namespace qwqdsp_spectral
