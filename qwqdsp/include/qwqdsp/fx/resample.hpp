#pragma once
#include <vector>
#include <span>
#include <cmath>
#include "qwqdsp/filter/window_fir.hpp"
#include "qwqdsp/interpolation.hpp"
#include "qwqdsp/window/kaiser.hpp"

// qwqfixme: 内核的位置不正确

namespace qwqdsp::fx {
class Resample {
public:
    /**
     * @param atten >0
     * @param kernel_len 必须是奇数
     */
    void Init(float source_fs, float target_fs, float atten, size_t kernel_len) noexcept {
        assert(kernel_len % 2 == 1);

        phase_inc_ = source_fs / target_fs;
        kernel_.resize(kernel_len + 3);
        float const beta = qwqdsp::window::Kaiser::Beta(atten);
        float const width = qwqdsp::window::Kaiser::MainLobeWidth(beta) * std::numbers::pi_v<float> * 2.0f / static_cast<float>(kernel_len);

        std::span kernel_org{kernel_.data() + 1, kernel_len};
        if (target_fs < source_fs) {
            float const cutoff = std::numbers::pi_v<float> * target_fs / source_fs - width;
            qwqdsp::filter::WindowFIR::Lowpass(kernel_org, cutoff);
        }
        else {
            float const cutoff = std::numbers::pi_v<float> - width;
            qwqdsp::filter::WindowFIR::Lowpass(kernel_org, cutoff);
        }

        qwqdsp::window::Kaiser::ApplyWindow(kernel_org, beta, false);
        std::reverse(kernel_org.begin(), kernel_org.end());
        kernel_[0] = kernel_org.front();
        kernel_[kernel_len + 1] = kernel_org[0];
        kernel_[kernel_len + 2] = kernel_org[1];
    }

    std::vector<float> Process(std::span<float> x) noexcept {
        std::vector<float> r;

        float phase = 0.0f;
        int const kernel_len = static_cast<int>(kernel_.size() - 3);
        int const xend = static_cast<int>(x.size());
        int const half_len = (kernel_len - 1) / 2;
        int rpos = -half_len;
        while (rpos < xend) {
            int const x_rbegin = std::max(0, rpos);
            int const x_rend = std::min(xend - 1, rpos + kernel_len);
            int const kernel_offset = x_rbegin - rpos;
            int const convo_len = x_rend - x_rbegin;
            float const frac = phase;
            float sum{};
            for (int i = 0; i < convo_len; ++i) {
                float const xx = x[static_cast<size_t>(x_rbegin + i)];
                size_t const kernel_rpos = static_cast<size_t>(kernel_offset + i + 1);
                float const yy = qwqdsp::Interpolation::PCHIP(kernel_[kernel_rpos - 1], kernel_[kernel_rpos], kernel_[kernel_rpos + 1], kernel_[kernel_rpos + 2], frac);
                sum += xx * yy;
            }
            r.push_back(sum);

            phase += phase_inc_;
            while (phase >= 1.0f) {
                phase -= 1.0f;
                ++rpos;
            }
        }

        return r;
    }

    std::span<const float> GetKernel() const {
        return {kernel_.data(), kernel_.size() - 3};
    }
private:
    float phase_inc_{};
    std::vector<float> kernel_;
};
}