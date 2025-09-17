#pragma once
#include "lag_buffer.hpp"

namespace qwqdsp::adaptive {
template<size_t kSize, float kStep = 0.1f>
class NLMS {
public:
    void Reset() noexcept {
        w_.fill(float{});
    }

    /**
     * @return 预测
     */
    float Tick(float x, float target) noexcept {
        x_.Push(x);

        float pred{};
        float sum{};
        for (size_t i = 0; i < kSize; ++i) {
            pred += w_[i] * x_[i];
            sum += x_[i] * x_[i];
        }

        float const e = target - pred;
        for (size_t i = 0; i < kSize; ++i) {
            w_[i] += kStep * e * x_[i] / (sum + 1e-3f);
        }

        return pred;
    }
private:
    LagBufferStatic<kSize> x_;
    std::array<float, kSize> w_{};
};
}

// def nlms(x, d, N=4, mu=0.1):
//   nIters = min(len(x),len(d)) - N
//   u = np.zeros(N)
//   w = np.zeros(N)
//   e = np.zeros(nIters)
//   for n in range(nIters):
//     u[1:] = u[:-1]
//     u[0] = x[n]
//     e_n = d[n] - np.dot(u, w)
//     w = w + mu * e_n * u / (np.dot(u,u)+1e-3)
//     e[n] = e_n
//   return e
