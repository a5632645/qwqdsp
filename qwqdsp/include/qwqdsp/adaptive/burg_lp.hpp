#pragma once
#include <span>
#include <vector>

namespace qwqdsp::adaptive {
class BurgLP {
public:
    void Init(size_t block_len) {
        eb_.resize(block_len);
    }

    void Process(std::span<float> x, std::span<float> latticek) {
        std::copy(x.begin(), x.end(), eb_.begin());
        for (auto& k : latticek) {
            float lag{};
            float up{};
            float down{};
            for (size_t i = 0; i < x.size(); ++i) {
                up += x[i] * lag;
                down += x[i] * x[i];
                down += lag * lag;
                lag = eb_[i];
            }
            k = -2.0f * up / down;

            lag = 0;
            for (size_t i = 0; i < x.size(); ++i) {
                float const upgo = x[i] + lag * k;
                float const downgo = lag + x[i] * k;
                lag = eb_[i];
                x[i] = upgo;
                eb_[i] = downgo;
            }
        }
    }
private:
    std::vector<float> eb_;
};
}