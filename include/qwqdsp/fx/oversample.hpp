#pragma once
#include <algorithm>
#include <cassert>
#include <span>
#include <vector>

namespace qwqdsp_fx {
class Oversample {
public:
    void Init(std::span<const float> halfband_coeff, int states) {
        num_stages_ = states;
        coeff_.assign(halfband_coeff.begin(), halfband_coeff.end());

        int const len = static_cast<int>(coeff_.size());
        // 必须是奇数半带滤波器
        assert(len % 2 == 1);
        half_len_ = (len - 1) / 2;

        // 上采样每级状态大小 = even phase 长度 = half_len_ + 1
        // 下采样使用全速率转置 FIR，状态大小 = filter len
        up_state_.resize(num_stages_ * (half_len_ + 1));
        down_state_.resize(num_stages_ * len);

        Reset();
    }

    void Reset() noexcept {
        std::fill(up_state_.begin(), up_state_.end(), 0.0f);
        std::fill(down_state_.begin(), down_state_.end(), 0.0f);
    }

    void Upsample(std::span<float> x_vec, std::span<float> up_vec) noexcept {
        int const n_in = static_cast<int>(x_vec.size());
        if (num_stages_ == 0) {
            std::copy_n(x_vec.data(), n_in, up_vec.data());
            return;
        }

        // 中间级暂存 buffer（最大为第二大到最终级输出大小）
        int scratch_size = n_in;
        for (int s = 0; s < num_stages_ - 1; ++s) {
            scratch_size *= 2;
        }
        scratch_.resize(scratch_size);

        float const* src = x_vec.data();
        int src_n = n_in;
        float* dst = up_vec.data();

        for (int stage = 0; stage < num_stages_; ++stage) {
            int const dst_n = src_n * 2;

            float* state = up_state_.data() + stage * (half_len_ + 1);
            Upsample2x(src, src_n, dst, coeff_.data(), half_len_, state);

            src = dst;
            if (dst == up_vec.data()) {
                dst = scratch_.data();
            }
            else {
                dst = up_vec.data();
            }

            src_n = dst_n;
        }

        if (dst == up_vec.data()) {
            std::copy(scratch_.begin(), scratch_.end(), up_vec.begin());
        }
    }

    void Downsample(std::span<float> up_vec, std::span<float> y_vec) noexcept {
        int const n_in = static_cast<int>(up_vec.size());
        if (num_stages_ == 0) {
            std::copy_n(up_vec.data(), n_in, y_vec.data());
            return;
        }

        // 中间级暂存 buffer
        int scratch_size = n_in;
        scratch_.resize(scratch_size);
        scratch2_.resize(scratch_size);

        float const* src = up_vec.data();
        int src_n = n_in;
        float* dst = scratch_.data();

        for (int stage = 0; stage < num_stages_; ++stage) {
            int const dst_n = src_n / 2;

            float* state = down_state_.data() + stage * static_cast<int>(coeff_.size());
            Downsample2x(src, src_n, dst, coeff_.data(), half_len_, state);

            src = dst;
            src_n = dst_n;
            if (dst == scratch_.data()) {
                dst = scratch2_.data();
            }
            else {
                dst = scratch_.data();
            }
        }

        if (dst == scratch2_.data()) {
            std::copy_n(scratch_.begin(), src_n, y_vec.data());
        }
        else {
            std::copy_n(scratch2_.begin(), src_n, y_vec.data());
        }
    }

    /**
     * @return 整个超采样链路的总群延迟（输入采样单位）
     *
     * 每级 2x FIR 的群延迟为 half_len 个采样（以该级采样率计）。
     * 所有上采样级和降采样级各自贡献后折算到输入域：
     *   latency = 2 * half_len * (1 - 1/2^num_stages)
     *
     * 例如 num_stages=1 → half_len,  num_stages=2 → 1.5*half_len,
     * num_stages=3 → 1.75*half_len, 随级数增加趋近 2*half_len。
     */
    float Latency() const noexcept {
        if (num_stages_ == 0) return 0.0f;
        float const sum = 1.0f - 1.0f / static_cast<float>(1 << num_stages_);
        return 2.0f * static_cast<float>(half_len_) * sum;
    }
private:
    static void Upsample2x(float const* x, int n_in, float* y, float const* coeff, int half_len,
                           float* state) noexcept {
        int const even_len = half_len + 1;
        int const odd_len = half_len;

        for (int i = 0; i < n_in; ++i) {
            // 移位寄存器
            for (int j = even_len - 1; j > 0; --j) {
                state[j] = state[j - 1];
            }
            state[0] = x[i] * 2;

            // even phase
            float y0 = 0;
            for (int j = 0; j < even_len; ++j) {
                y0 += state[j] * coeff[2 * j];
            }
            y[2 * i] = y0;

            // odd phase
            float y1 = 0;
            for (int j = 0; j < odd_len; ++j) {
                y1 += state[j] * coeff[2 * j + 1];
            }
            y[2 * i + 1] = y1;
        }
    }

    static void Downsample2x(float const* x, int n_in, float* y, float const* coeff, int half_len,
                             float* state) noexcept {
        assert(n_in % 2 == 0);
        int const len = 2 * half_len + 1;
        int const n_out = n_in / 2;

        for (int i = 0; i < n_out; ++i) {
            // 推入 x[2*i]（直接型 FIR 延迟线右移）
            for (int j = len - 1; j > 0; --j)
                state[j] = state[j - 1];
            state[0] = x[2 * i];

            // 以 x[2*i] 为最新采样计算 FIR 输出
            // 对应 y[i] = sum_k coeff[k] * x[2*i - k]
            float y_filt = 0;
            for (int j = 0; j < len; ++j)
                y_filt += state[j] * coeff[j];
            y[i] = y_filt;

            // 推入 x[2*i + 1]（为下一轮迭代准备状态）
            for (int j = len - 1; j > 0; --j)
                state[j] = state[j - 1];
            state[0] = x[2 * i + 1];
        }
    }

    int num_stages_{};
    int half_len_{};
    std::vector<float> coeff_;
    std::vector<float> up_state_;
    std::vector<float> down_state_;
    std::vector<float> scratch_;
    std::vector<float> scratch2_;
};

} // namespace qwqdsp_fx
