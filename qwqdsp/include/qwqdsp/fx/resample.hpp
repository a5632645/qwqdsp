#pragma once
#include <vector>
#include <span>
#include "qwqdsp/filter/window_fir.hpp"
#include "qwqdsp/interpolation.hpp"
#include "qwqdsp/window/kaiser.hpp"

namespace qwqdsp::fx {
class Resample {
public:
    static constexpr size_t kOversample = 128;

    /**
     * @param atten >0
     * @param kernel_len 必须是奇数
     */
    void Init(float source_fs, float target_fs, float atten, size_t kernel_len) {
        assert(kernel_len % 2 == 1);

        phase_inc_ = source_fs / target_fs;
        kernel_len_ = kernel_len;

        kernel_.resize(kernel_len * kOversample + 2, 0.0f);
        float const beta = qwqdsp::window::Kaiser::Beta(atten);
        float const width = qwqdsp::window::Kaiser::MainLobeWidth(beta) * std::numbers::pi_v<float> * 2.0f / static_cast<float>(kernel_len);
        
        std::span kernel_block{kernel_.data(), kernel_len * kOversample};
        if (target_fs < source_fs) {
            float const cutoff = std::numbers::pi_v<float> * target_fs / source_fs - width;
            qwqdsp::filter::WindowFIR::Lowpass(kernel_block, cutoff / kOversample);
        }
        else {
            float const cutoff = std::numbers::pi_v<float> - width;
            qwqdsp::filter::WindowFIR::Lowpass(kernel_block, cutoff / kOversample);
        }

        for (auto& s : kernel_block) {
            s *= kOversample;
        }

        qwqdsp::window::Kaiser::ApplyWindow(kernel_block, beta, false);
        std::reverse(kernel_block.begin(), kernel_block.end());
    }

    std::vector<float> Process(std::span<float> x) noexcept {
        std::vector<float> r;

        float phase = 0.0f;
        int const ikernel_len = static_cast<int>(kernel_len_);
        int const xsize = static_cast<int>(x.size());
        int xrpos = -(static_cast<int>(kernel_len_) - 1) / 2;

        while (xrpos < xsize) {
            float const frac = 1.0f - phase;

            float sum{};
            int const begin = std::max(0, -xrpos);
            int const xibegin = std::max(0, frac == 0.0f ? xrpos : xrpos + 1);
            int const len = std::min(ikernel_len, xsize - xibegin);
            for (int i = begin; i < len; ++i) {
                float const krpos = (static_cast<float>(i) + frac) * kOversample;
                size_t const ikrpos = static_cast<size_t>(krpos);
                float const frac_krpos = krpos - std::floor(krpos);
                float const kernel_v = qwqdsp::Interpolation::Linear(kernel_[ikrpos], kernel_[ikrpos + 1], frac_krpos);
                sum += kernel_v * x[static_cast<size_t>(xibegin + i)];
            }
            r.push_back(sum);

            phase += phase_inc_;
            xrpos += static_cast<int>(std::floor(phase));
            phase = phase - std::floor(phase);
        }

        return r;
    }

private:
    float phase_inc_{};
    size_t kernel_len_{};
    std::vector<float> kernel_;
};
}