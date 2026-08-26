#pragma once
#include "../interpolation.hpp"
#include "../window/kaiser.hpp"
#include <cassert>
#include <numbers>
#include <span>
#include <vector>

namespace qwqdsp_fx {
class Resample {
public:
    /**
     * @param atten (>0)dB 这决定了target_fs/2处的衰减值
     * @param kernel_len >=3 必须是奇数，越大过渡带越小，计算量越大
     * @param oversample >0 越大质量越好，内存需求越大
     */
    void Init(float source_fs, float target_fs, float atten, size_t kernel_len, size_t oversample) {
        assert(kernel_len % 2 == 1);
        assert(source_fs > 0.0f && target_fs > 0.0f && atten > 0.0f && oversample > 0);

        phase_inc_ = source_fs / target_fs;
        kernel_len_ = kernel_len;
        oversample_plus1_ = oversample + 1;

        size_t const lut_size = (kernel_len - 1) * oversample_plus1_ + 1;
        kernel_.resize(lut_size + 1, 0.0f);
        float const beta = qwqdsp_window::Kaiser::Beta(atten);
        float const width = qwqdsp_window::Kaiser::MainLobeWidth(beta) * std::numbers::pi_v<float>
                          * 2.0f / static_cast<float>(kernel_len);

        std::span kernel_block{kernel_.data(), lut_size};
        {
            float cutoff = 0.0f;
            cutoff = std::numbers::pi_v<float> * std::min(1.0f, target_fs / source_fs) * 0.97f;

            float const center = (static_cast<float>(kernel_block.size()) - 1.0f) / 2.0f;
            float const omega = cutoff / static_cast<float>(oversample_plus1_);
            for (size_t i = 0; i < kernel_block.size(); ++i) {
                float t = static_cast<float>(i) - center;
                [[unlikely]]
                if (t == 0.0f) {
                    kernel_block[i] = cutoff / std::numbers::pi_v<float>;
                }
                else {
                    kernel_block[i] =
                        std::sin(omega * t) * static_cast<float>(oversample_plus1_) / (std::numbers::pi_v<float> * t);
                }
            }
        }

        qwqdsp_window::Kaiser::ApplyWindow(kernel_block, beta, false);
    }

    std::vector<float> Process(std::span<float> x) {
        std::vector<float> r;

        double position = 0.0;
        int const half_len = static_cast<int>((kernel_len_ - 1) / 2);
        while (position < static_cast<double>(x.size())) {
            float sum{};
            float weight_sum{};
            int const base = static_cast<int>(std::floor(position));
            float const frac = static_cast<float>(position - static_cast<double>(base));
            for (int tap = -half_len; tap <= half_len; ++tap) {
                int const index = base + tap;
                if (index < 0 || index >= static_cast<int>(x.size())) continue;
                float const table_pos = static_cast<float>(half_len - tap) * static_cast<float>(oversample_plus1_)
                                      + frac * static_cast<float>(oversample_plus1_);
                size_t const table_index = std::min(static_cast<size_t>(table_pos), kernel_.size() - 2);
                float const table_frac = table_pos - static_cast<float>(table_index);
                float const weight =
                    qwqdsp::Interpolation::Linear(kernel_[table_index], kernel_[table_index + 1], table_frac);
                sum += weight * x[static_cast<size_t>(index)];
                weight_sum += weight;
            }
            r.push_back(std::abs(weight_sum) > 1.0e-12f ? sum / weight_sum : 0.0f);
            position += static_cast<double>(phase_inc_);
        }

        return r;
    }
private:
    float phase_inc_{};
    size_t oversample_plus1_{};
    size_t kernel_len_{};
    std::vector<float> kernel_;
};
} // namespace qwqdsp_fx
