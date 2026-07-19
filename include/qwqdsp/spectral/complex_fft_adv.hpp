#pragma once
#include "complex_fft.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <span>
#include <vector>

namespace qwqdsp_spectral {

/**
 * @brief ComplexFFT 的高级封装，提供 span 接口、增益/相位、Hilbert 等便捷功能
 *
 * 默认输出 [0, 2π] 频率顺序。
 *           奈奎斯特     负频率        零            正频率
 *   -π ~ π     0      1 ~ n/2-1      n/2        n/2+1 ~ n-1
 *              零       正频率      奈奎斯特         负频率
 *   0 ~ 2π     0      1 ~ n/2-1      n/2        n/2+1 ~ n
 */
class ComplexFftAdv {
public:
    void Init(size_t fft_size) {
        fft_.Init(fft_size);
        real_buf_.resize(fft_size);
        imag_buf_.resize(fft_size);
        fft_size_ = fft_size;
    }

    // ------------------------------------------------------------
    // FFT
    // ------------------------------------------------------------

    void FFT(std::span<const float> time, std::span<std::complex<float>> spectral) noexcept {
        assert(time.size() == fft_size_);
        assert(spectral.size() == NumBins());
        std::copy(time.begin(), time.end(), real_buf_.begin());
        std::fill(imag_buf_.begin(), imag_buf_.end(), 0.0f);
        fft_.FFT(real_buf_.data(), imag_buf_.data(), real_buf_.data(), imag_buf_.data());
        for (size_t i = 0; i < fft_size_; ++i) {
            spectral[i] = {real_buf_[i], imag_buf_[i]};
        }
    }

    void FFT(std::span<const std::complex<float>> time, std::span<std::complex<float>> spectral) noexcept {
        assert(time.size() == fft_size_);
        assert(spectral.size() == NumBins());
        for (size_t i = 0; i < fft_size_; ++i) {
            real_buf_[i] = time[i].real();
            imag_buf_[i] = time[i].imag();
        }
        fft_.FFT(real_buf_.data(), imag_buf_.data(), real_buf_.data(), imag_buf_.data());
        for (size_t i = 0; i < fft_size_; ++i) {
            spectral[i] = {real_buf_[i], imag_buf_[i]};
        }
    }

    void FFT(std::span<const float> time, std::span<float> real, std::span<float> imag) noexcept {
        assert(time.size() == fft_size_);
        assert(real.size() == NumBins());
        assert(imag.size() == NumBins());
        std::copy(time.begin(), time.end(), real_buf_.begin());
        std::fill(imag_buf_.begin(), imag_buf_.end(), 0.0f);
        fft_.FFT(real_buf_.data(), imag_buf_.data(), real.data(), imag.data());
    }

    void FFT(std::span<const std::complex<float>> time, std::span<float> real, std::span<float> imag) noexcept {
        assert(time.size() == fft_size_);
        assert(real.size() == NumBins());
        assert(imag.size() == NumBins());
        for (size_t i = 0; i < fft_size_; ++i) {
            real_buf_[i] = time[i].real();
            imag_buf_[i] = time[i].imag();
        }
        fft_.FFT(real_buf_.data(), imag_buf_.data(), real.data(), imag.data());
    }

    /**
     * @param phase 可选的，不需要请传入 {}
     */
    void FFTGainPhase(std::span<const float> time, std::span<float> gain, std::span<float> phase = {}) noexcept {
        assert(time.size() == fft_size_);
        assert(gain.size() == NumBins());
        std::copy(time.begin(), time.end(), real_buf_.begin());
        std::fill(imag_buf_.begin(), imag_buf_.end(), 0.0f);
        fft_.FFT(real_buf_.data(), imag_buf_.data(), real_buf_.data(), imag_buf_.data());
        if (phase.empty()) {
            for (size_t i = 0; i < fft_size_; ++i) {
                gain[i] = std::sqrt(real_buf_[i] * real_buf_[i] + imag_buf_[i] * imag_buf_[i]);
            }
        }
        else {
            assert(phase.size() == NumBins());
            for (size_t i = 0; i < fft_size_; ++i) {
                gain[i] = std::sqrt(real_buf_[i] * real_buf_[i] + imag_buf_[i] * imag_buf_[i]);
                phase[i] = std::atan2(imag_buf_[i], real_buf_[i]);
            }
        }
    }

    /**
     * @param phase 可选的，不需要请传入 {}
     */
    void FFTGainPhase(std::span<const std::complex<float>> time, std::span<float> gain,
                      std::span<float> phase = {}) noexcept {
        assert(time.size() == fft_size_);
        assert(gain.size() == NumBins());
        for (size_t i = 0; i < fft_size_; ++i) {
            real_buf_[i] = time[i].real();
            imag_buf_[i] = time[i].imag();
        }
        fft_.FFT(real_buf_.data(), imag_buf_.data(), real_buf_.data(), imag_buf_.data());
        if (phase.empty()) {
            for (size_t i = 0; i < fft_size_; ++i) {
                gain[i] = std::sqrt(real_buf_[i] * real_buf_[i] + imag_buf_[i] * imag_buf_[i]);
            }
        }
        else {
            assert(phase.size() == NumBins());
            for (size_t i = 0; i < fft_size_; ++i) {
                gain[i] = std::sqrt(real_buf_[i] * real_buf_[i] + imag_buf_[i] * imag_buf_[i]);
                phase[i] = std::atan2(imag_buf_[i], real_buf_[i]);
            }
        }
    }

    // ------------------------------------------------------------
    // IFFT
    // ------------------------------------------------------------

    void IFFT(std::span<float> time, std::span<const std::complex<float>> spectral) noexcept {
        assert(time.size() == fft_size_);
        assert(spectral.size() == NumBins());
        for (size_t i = 0; i < fft_size_; ++i) {
            real_buf_[i] = spectral[i].real();
            imag_buf_[i] = spectral[i].imag();
        }
        fft_.IFFT(real_buf_.data(), imag_buf_.data(), real_buf_.data(), imag_buf_.data());
        std::copy_n(real_buf_.data(), fft_size_, time.data());
    }

    void IFFT(std::span<std::complex<float>> time, std::span<const std::complex<float>> spectral) noexcept {
        assert(time.size() == fft_size_);
        assert(spectral.size() == NumBins());
        for (size_t i = 0; i < fft_size_; ++i) {
            real_buf_[i] = spectral[i].real();
            imag_buf_[i] = spectral[i].imag();
        }
        fft_.IFFT(real_buf_.data(), imag_buf_.data(), real_buf_.data(), imag_buf_.data());
        for (size_t i = 0; i < fft_size_; ++i) {
            time[i] = {real_buf_[i], imag_buf_[i]};
        }
    }

    void IFFT(std::span<float> time, std::span<const float> real, std::span<const float> imag) noexcept {
        assert(time.size() == fft_size_);
        assert(real.size() == NumBins());
        assert(imag.size() == NumBins());
        fft_.IFFT(real.data(), imag.data(), real_buf_.data(), imag_buf_.data());
        std::copy_n(real_buf_.data(), fft_size_, time.data());
    }

    void IFFT(std::span<std::complex<float>> time, std::span<const float> real, std::span<const float> imag) noexcept {
        assert(time.size() == fft_size_);
        assert(real.size() == NumBins());
        assert(imag.size() == NumBins());
        fft_.IFFT(real.data(), imag.data(), real_buf_.data(), imag_buf_.data());
        for (size_t i = 0; i < fft_size_; ++i) {
            time[i] = {real_buf_[i], imag_buf_[i]};
        }
    }

    void IFFTGainPhase(std::span<float> time, std::span<const float> gain, std::span<const float> phase) noexcept {
        assert(time.size() == fft_size_);
        assert(gain.size() == NumBins());
        assert(phase.size() == NumBins());
        for (size_t i = 0; i < fft_size_; ++i) {
            auto const a = std::polar(gain[i], phase[i]);
            real_buf_[i] = a.real();
            imag_buf_[i] = a.imag();
        }
        fft_.IFFT(real_buf_.data(), imag_buf_.data(), real_buf_.data(), imag_buf_.data());
        std::copy_n(real_buf_.data(), fft_size_, time.data());
    }

    void IFFTGainPhase(std::span<std::complex<float>> time, std::span<const float> gain,
                       std::span<const float> phase) noexcept {
        assert(time.size() == fft_size_);
        assert(gain.size() == NumBins());
        assert(phase.size() == NumBins());
        for (size_t i = 0; i < fft_size_; ++i) {
            auto const a = std::polar(gain[i], phase[i]);
            real_buf_[i] = a.real();
            imag_buf_[i] = a.imag();
        }
        fft_.IFFT(real_buf_.data(), imag_buf_.data(), real_buf_.data(), imag_buf_.data());
        for (size_t i = 0; i < fft_size_; ++i) {
            time[i] = {real_buf_[i], imag_buf_[i]};
        }
    }

    // ------------------------------------------------------------
    // Hilbert
    // ------------------------------------------------------------

    /**
     * @brief 计算解析信号（Hilbert 变换）
     * @param clear_dc 是否清除直流和奈奎斯特分量
     */
    void Hilbert(std::span<const float> time, std::span<std::complex<float>> output, bool clear_dc) noexcept {
        assert(time.size() == fft_size_);
        assert(output.size() == fft_size_);
        std::copy(time.begin(), time.end(), real_buf_.begin());
        std::fill(imag_buf_.begin(), imag_buf_.end(), 0.0f);
        fft_.FFT(real_buf_.data(), imag_buf_.data(), real_buf_.data(), imag_buf_.data());

        size_t const nyquist = fft_size_ / 2;
        // 零化负频率 (nyquist+1 .. fft_size_-1)
        for (size_t i = nyquist + 1; i < fft_size_; ++i) {
            real_buf_[i] = 0.0f;
            imag_buf_[i] = 0.0f;
        }
        if (clear_dc) {
            real_buf_[0] = 0.0f;
            imag_buf_[0] = 0.0f;
            real_buf_[nyquist] = 0.0f;
            imag_buf_[nyquist] = 0.0f;
        }

        fft_.IFFT(real_buf_.data(), imag_buf_.data(), real_buf_.data(), imag_buf_.data());
        for (size_t i = 0; i < fft_size_; ++i) {
            output[i] = {real_buf_[i], imag_buf_[i]};
        }
    }

    /**
     * @brief 计算实数 Hilbert 变换结果（仅实部输出）
     *        解析信号为 x_a[n] = x[n] + j·H(x)[n]，此处返回 H(x)[n]（IFFT 的虚部）
     * @param clear_dc 是否清除直流分量
     */
    void Hilbert(std::span<const float> input, std::span<float> output90, bool clear_dc) noexcept {
        assert(input.size() == fft_size_);
        assert(output90.size() == fft_size_);
        std::copy(input.begin(), input.end(), real_buf_.begin());
        std::fill(imag_buf_.begin(), imag_buf_.end(), 0.0f);
        fft_.FFT(real_buf_.data(), imag_buf_.data(), real_buf_.data(), imag_buf_.data());

        size_t const nyquist = fft_size_ / 2;
        for (size_t i = nyquist + 1; i < fft_size_; ++i) {
            real_buf_[i] = 0.0f;
            imag_buf_[i] = 0.0f;
        }
        if (clear_dc) {
            real_buf_[0] = 0.0f;
            imag_buf_[0] = 0.0f;
            real_buf_[nyquist] = 0.0f;
            imag_buf_[nyquist] = 0.0f;
        }

        fft_.IFFT(real_buf_.data(), imag_buf_.data(), real_buf_.data(), imag_buf_.data());
        std::copy_n(imag_buf_.data(), fft_size_, output90.data());
    }

    // ------------------------------------------------------------
    // 工具函数
    // ------------------------------------------------------------

    size_t NumBins() const noexcept {
        return fft_size_;
    }

    static constexpr size_t NumBins(size_t fft_size) noexcept {
        return fft_size;
    }

    size_t FFTSize() const noexcept {
        return fft_size_;
    }

    /**
     * [0, 2π) → [-π, π) 或 [-π, π) → [0, 2π)
     */
    template <class TYPE>
    static void ShuffleFrequency(std::span<TYPE> buffer) {
        size_t const n = buffer.size();
        size_t neg_idx = n / 2;
        size_t pos_idx = 0;
        for (size_t i = 0; i < n / 2; ++i) {
            std::swap(buffer[pos_idx], buffer[neg_idx]);
            ++pos_idx;
            ++neg_idx;
        }
    }
private:
    ComplexFFT fft_;
    std::vector<float> real_buf_;
    std::vector<float> imag_buf_;
    size_t fft_size_{};
};

} // namespace qwqdsp_spectral
