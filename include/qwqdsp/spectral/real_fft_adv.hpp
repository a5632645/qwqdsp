#pragma once
#include "real_fft.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <span>
#include <vector>

namespace qwqdsp_spectral {

/**
 * @brief RealFFT 的高级封装，提供 span 接口、增益/相位、Hilbert 等便捷功能
 */
class RealFftAdv {
public:
    void Init(size_t fft_size) {
        fft_.Init(fft_size);
        buf_.resize(fft_size + 2);
        fft_size_ = fft_size;
    }

    // ------------------------------------------------------------
    // FFT
    // ------------------------------------------------------------

    void FFT(std::span<const float> time, std::span<std::complex<float>> spectral) noexcept {
        assert(time.size() == fft_size_);
        assert(spectral.size() == NumBins());
        fft_.FFT(time.data(), buf_.data());
        for (size_t i = 0; i < NumBins(); ++i) {
            spectral[i] = {buf_[2 * i], buf_[2 * i + 1]};
        }
    }

    void FFT(std::span<const float> time, std::span<float> real, std::span<float> imag) noexcept {
        assert(time.size() == fft_size_);
        assert(real.size() == NumBins());
        assert(imag.size() == NumBins());
        fft_.FFT(time.data(), buf_.data());
        for (size_t i = 0; i < NumBins(); ++i) {
            real[i] = buf_[2 * i];
            imag[i] = buf_[2 * i + 1];
        }
    }

    /**
     * @param phase 可选的，不需要请传入 {}
     */
    void FFTGainPhase(std::span<const float> time, std::span<float> gain, std::span<float> phase = {}) noexcept {
        assert(time.size() == fft_size_);
        assert(gain.size() == NumBins());
        fft_.FFT(time.data(), buf_.data());
        if (phase.empty()) {
            for (size_t i = 0; i < NumBins(); ++i) {
                float re = buf_[2 * i];
                float im = buf_[2 * i + 1];
                gain[i] = std::sqrt(re * re + im * im);
            }
        }
        else {
            assert(phase.size() == NumBins());
            for (size_t i = 0; i < NumBins(); ++i) {
                float re = buf_[2 * i];
                float im = buf_[2 * i + 1];
                gain[i] = std::sqrt(re * re + im * im);
                phase[i] = std::atan2(im, re);
            }
        }
    }

    // ------------------------------------------------------------
    // IFFT
    // ------------------------------------------------------------

    void IFFT(std::span<float> time, std::span<const std::complex<float>> spectral) noexcept {
        assert(time.size() == fft_size_);
        assert(spectral.size() == NumBins());
        for (size_t i = 0; i < NumBins(); ++i) {
            buf_[2 * i] = spectral[i].real();
            buf_[2 * i + 1] = spectral[i].imag();
        }
        fft_.IFFT(buf_.data(), time.data());
    }

    void IFFT(std::span<float> time, std::span<const float> real, std::span<const float> imag) noexcept {
        assert(time.size() == fft_size_);
        assert(real.size() == NumBins());
        assert(imag.size() == NumBins());
        for (size_t i = 0; i < NumBins(); ++i) {
            buf_[2 * i] = real[i];
            buf_[2 * i + 1] = imag[i];
        }
        fft_.IFFT(buf_.data(), time.data());
    }

    void IFFTGainPhase(std::span<float> time, std::span<const float> gain, std::span<const float> phase) noexcept {
        assert(time.size() == fft_size_);
        assert(gain.size() == NumBins());
        assert(phase.size() == NumBins());
        for (size_t i = 0; i < NumBins(); ++i) {
            auto const a = std::polar(gain[i], phase[i]);
            buf_[2 * i] = a.real();
            buf_[2 * i + 1] = a.imag();
        }
        fft_.IFFT(buf_.data(), time.data());
    }

    // ------------------------------------------------------------
    // Hilbert
    // ------------------------------------------------------------

    /**
     * @brief 通过 FFT 域 90° 相移计算 Hilbert 变换
     * @param clear_dc 是否清除直流分量
     */
    void Hilbert(std::span<const float> input, std::span<float> shift90, bool clear_dc) noexcept {
        assert(input.size() == fft_size_);
        assert(shift90.size() == fft_size_);
        fft_.FFT(input.data(), buf_.data());
        const size_t n = fft_size_ / 2;
        for (size_t i = 1; i < n; ++i) {
            // (Re, Im) → (-Im, Re) = 90° 相移
            float re = -buf_[2 * i + 1];
            float im = buf_[2 * i];
            buf_[2 * i] = re;
            buf_[2 * i + 1] = im;
        }
        if (clear_dc) {
            buf_[0] = 0.0f;
            buf_[1] = 0.0f;
            buf_[fft_size_] = 0.0f;
            buf_[fft_size_ + 1] = 0.0f;
        }
        fft_.IFFT(buf_.data(), shift90.data());
    }

    // ------------------------------------------------------------
    // 时域操作
    // ------------------------------------------------------------

    /**
     * @brief 0 ~ N → -N/2 ~ N/2（循环左移 N/2）
     */
    static void TimeDomainShift(std::span<float> block) noexcept {
        auto const n = block.size();
        assert(n > 0 && std::has_single_bit(n));
        std::rotate(block.begin(), block.begin() + static_cast<std::ptrdiff_t>(n / 2), block.end());
    }

    // ------------------------------------------------------------
    // 工具函数
    // ------------------------------------------------------------

    size_t NumBins() const noexcept {
        return fft_size_ / 2 + 1;
    }

    static constexpr size_t NumBins(size_t fft_size) noexcept {
        return fft_size / 2 + 1;
    }

    size_t FFTSize() const noexcept {
        return fft_size_;
    }
private:
    RealFFT fft_;
    std::vector<float> buf_;
    size_t fft_size_{};
};

} // namespace qwqdsp_spectral
