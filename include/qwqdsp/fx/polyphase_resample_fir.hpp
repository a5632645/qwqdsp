#pragma once
#include <cassert>
#include <span>
#include <vector>

namespace qwqdsp_fx {
// ------------------------------------------------------------
// downsample
// ------------------------------------------------------------
class PolyphaseDownsamplerFir {
public:
    void Init(std::span<const float> coeff, int downsample) {
        downsample_ = downsample;

        int len_coeff = static_cast<int>(coeff.size());
        each_phase_len_ = len_coeff / downsample;
        if (len_coeff % downsample != 0) {
            ++each_phase_len_;
        }

        phase_coeffs_.resize(downsample * each_phase_len_);
        for (int i = 0; i < downsample; ++i) {
            // 反转子滤波器的求值顺序
            auto dst_it = phase_coeffs_.begin() + (downsample - i - 1) * each_phase_len_;
            for (int j = 0; j < each_phase_len_; ++j) {
                int idx = j * downsample + i;
                if (idx < coeff.size()) {
                    *dst_it = coeff[idx];
                }
                else {
                    *dst_it = 0;
                }
                ++dst_it;
            }
        }

        phase_state_.resize(phase_coeffs_.size());
    }

    void Reset() noexcept {
        std::fill(phase_state_.begin(), phase_state_.end(), 0);
    }

    float Tick(std::span<float> x) noexcept {
        assert(x.size() == downsample_);

        float y = 0;
        for (int i = 0; i < downsample_; ++i) {
            auto* lag_it = phase_state_.data() + i * each_phase_len_;
            auto* coeff_it = phase_coeffs_.data() + i * each_phase_len_;
            y += TransposeFir(lag_it, coeff_it, x[i], each_phase_len_);
        }
        return y;
    }

    std::vector<float> Process(std::span<float> x) {
        std::vector<float> out;

        int full_loop_count = x.size() / downsample_;
        int scalar_count = x.size() - full_loop_count * downsample_;
        out.reserve(scalar_count == 0 ? full_loop_count : full_loop_count + 1);

        for (int block = 0; block < full_loop_count; ++block) {
            int const block_start = block * downsample_;
            float y = 0;
            for (int j = 0; j < downsample_; ++j) {
                auto* lag_it = phase_state_.data() + j * each_phase_len_;
                auto* coeff_it = phase_coeffs_.data() + j * each_phase_len_;
                y += TransposeFir(lag_it, coeff_it, x[block_start + j], each_phase_len_);
            }
            out.push_back(y);
        }

        if (scalar_count != 0) {
            int const block_start = full_loop_count * downsample_;
            float y = 0;
            for (int j = 0; j < scalar_count; ++j) {
                auto* lag_it = phase_state_.data() + j * each_phase_len_;
                auto* coeff_it = phase_coeffs_.data() + j * each_phase_len_;
                y += TransposeFir(lag_it, coeff_it, x[block_start + j], each_phase_len_);
            }
            for (int j = scalar_count; j < downsample_; ++j) {
                auto* lag_it = phase_state_.data() + j * each_phase_len_;
                auto* coeff_it = phase_coeffs_.data() + j * each_phase_len_;
                y += TransposeFir(lag_it, coeff_it, 0, each_phase_len_);
            }
            out.push_back(y);
        }

        return out;
    }
private:
    static float TransposeFir(float* lag, float const* coeff, float x, int len) noexcept {
        float const y = coeff[0] * x + lag[0];
        for (size_t i = 0; i < len - 1; ++i) {
            lag[i] = lag[i + 1] + coeff[i + 1] * x;
        }
        lag[len - 1] = coeff[len - 1] * x;
        return y;
    }

    int downsample_{};
    int each_phase_len_{};
    std::vector<float> phase_coeffs_;
    std::vector<float> phase_state_;
};

// ------------------------------------------------------------
// upsample
// ------------------------------------------------------------
class PolyphaseUpsamplerFir {
public:
    void Init(std::span<const float> coeff, int upsample) {
        upsample_ = upsample;

        int len_coeff = static_cast<int>(coeff.size());
        each_phase_len_ = len_coeff / upsample;
        if (len_coeff % upsample != 0) {
            ++each_phase_len_;
        }

        phase_coeffs_.resize(upsample * each_phase_len_);
        for (int i = 0; i < upsample; ++i) {
            auto dst_it = phase_coeffs_.begin() + i * each_phase_len_;
            for (int j = 0; j < each_phase_len_; ++j) {
                int idx = j * upsample + i;
                if (idx < coeff.size()) {
                    *dst_it = coeff[idx] * upsample;
                }
                else {
                    *dst_it = 0;
                }
                ++dst_it;
            }
        }

        phase_state_.resize(each_phase_len_);
    }

    void Reset() noexcept {
        std::fill(phase_state_.begin(), phase_state_.end(), 0);
    }

    void Tick(float x, std::span<float> y) noexcept {
        assert(y.size() == upsample_);

        PushFirState(x);
        for (int i = 0; i < upsample_; ++i) {
            auto* coeff_it = phase_coeffs_.data() + i * each_phase_len_;
            y[i] = Fir(coeff_it);
        }
    }

    std::vector<float> Process(std::span<float> x) {
        std::vector<float> out;
        out.reserve(x.size() * upsample_);

        for (auto s : x) {
            PushFirState(s);
            for (int i = 0; i < upsample_; ++i) {
                auto* coeff_it = phase_coeffs_.data() + i * each_phase_len_;
                out.push_back(Fir(coeff_it));
            }
        }

        return out;
    }
private:
    float Fir(float const* coeff) noexcept {
        float sum = 0;
        for (int i = 0; i < each_phase_len_; ++i) {
            sum += phase_state_[i] * coeff[i];
        }
        return sum;
    }

    void PushFirState(float x) noexcept {
        for (int i = each_phase_len_ - 1; i > 0; --i) {
            phase_state_[i] = phase_state_[i - 1];
        }
        phase_state_[0] = x;
    }

    int upsample_{};
    int each_phase_len_{};
    std::vector<float> phase_coeffs_;
    std::vector<float> phase_state_;
};
} // namespace qwqdsp_fx
