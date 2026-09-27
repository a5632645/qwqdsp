#pragma once
#include <algorithm>
#include <cassert>
#include <span>
#include <vector>

namespace qwqdsp_adaptive {
class BurgLP {
public:
    void Init(size_t block_len) {
        eb_.resize(block_len);
    }

    /**
     * @brief 逐级 Burg 递推，求出各阶反射系数
     * @param x 输入信号。**会被就地修改**（递推的 upgoing 结果直接写在 x 上），
     *          需要保留原信号时调用方先自行拷贝
     * @param latticek 输出的反射系数，长度 = 想要的最大阶数
     * @note 需要 `x.size() >= latticek.size()`；内部缓冲由 `Init(block_len)` 分配，
     *       应保证 `block_len >= x.size()`
     *
     *                     +-----k----+
     *                     |          ↓
     *   x -------------------------> + -----> x
     *                     |    |
     *                     |    +--k-+
     *            +-----+  |         ↓
     *  eb -------|z^-1|-----------> + ------> eb
     *            +----+
     */
    void Process(std::span<float> x, std::span<float> latticek) noexcept {
        assert(eb_.size() >= x.size());
        assert(x.size() >= latticek.size());
        // std::copy(x.begin(), x.end(), eb_.begin());
        // for (auto& k : latticek) {
        //     float lag{};
        //     float up{};
        //     float down{};
        //     for (size_t i = 0; i < x.size(); ++i) {
        //         up += x[i] * lag;
        //         down += x[i] * x[i];
        //         down += lag * lag;
        //         lag = eb_[i];
        //     }
        //     k = -2.0f * up / down;

        //     lag = 0;
        //     for (size_t i = 0; i < x.size(); ++i) {
        //         float const upgo = x[i] + lag * k;
        //         float const downgo = lag + x[i] * k;
        //         lag = eb_[i];
        //         x[i] = upgo;
        //         eb_[i] = downgo;
        //     }
        // }

        std::copy(x.begin() + 1, x.end(), eb_.begin());
        for (size_t kidx = 0; auto& k : latticek) {
            ++kidx;

            float up{};
            float down{};
            for (size_t i = 0; i < x.size() - kidx; ++i) {
                up += x[i] * eb_[i];
                down += x[i] * x[i];
                down += eb_[i] * eb_[i];
            }
            k = -2.0f * up / down;

            for (size_t i = 0; i < x.size() - kidx - 1; ++i) {
                float const upgo = x[i] + eb_[i] * k;
                float const downgo = eb_[i + 1] + x[i + 1] * k;
                x[i] = upgo;
                eb_[i] = downgo;
            }
        }
    }

    /**
     * @note 执行之后k会被b复写, b.size() = k.size()
     * @param b sum b[i] * z^-i, i from 1 to k.size
     */
    static void Lattice2Tf(std::span<float> k, std::span<float> b) noexcept {
        for (size_t i = 0; i < k.size(); ++i) {
            for (size_t j = 0; j + 1 <= i; j++) {
                k[j] = b[j] + k[i] * b[i - j - 1];
            }

            for (size_t j = 0; j <= i; j++) {
                b[j] = k[j];
            }
        }
    }

    /**
     * @note 内部会先把两个多项式清零并从 0 阶 `A_0(z) = A_0(z^-1) = 1` 起步，
     *       调用方不需要（也不应该）预置初值
     * @note upgoing.size() = downgoing.size() = k.size() + 1
     * @param upgoing sum upgoing[i] * z^-i, i from 0 to k.size, 最小相位
     * @param downgoing sum downgoing[i] * z^-i, i from 0 to k.size, 最大相位
     */
    static void Lattice2Tf_KeepK(std::span<const float> k, std::span<float> upgoing,
                                 std::span<float> downgoing) noexcept {
        assert(upgoing.size() > k.size());
        assert(downgoing.size() > k.size());

        // 0 阶多项式：A_0(z) = A_0(z^-1) = 1
        std::fill(upgoing.begin(), upgoing.end(), 0.0f);
        std::fill(downgoing.begin(), downgoing.end(), 0.0f);
        upgoing[0] = 1.0f;
        downgoing[0] = 1.0f;

        for (size_t kidx = 0; kidx < k.size(); ++kidx) {
            for (size_t i = kidx + 1; i != 0; --i) {
                downgoing[i] = downgoing[i - 1];
            }
            downgoing[0] = 0;

            for (size_t i = 0; i < kidx + 2; ++i) {
                float up = upgoing[i] + k[kidx] * downgoing[i];
                float down = downgoing[i] + k[kidx] * upgoing[i];
                upgoing[i] = up;
                downgoing[i] = down;
            }
        }
    }
private:
    std::vector<float> eb_;
};
} // namespace qwqdsp_adaptive
