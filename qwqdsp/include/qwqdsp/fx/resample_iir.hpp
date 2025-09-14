#pragma once
#include <vector>
#include <span>
#include "../source/elliptic_blep.hpp"

namespace qwqdsp::fx {
/**
 * @brief holters-parker IIR重采样器，使用Elliptic-blep库实现，移除了高通滤波器系数
 */
template<class TCoeff, size_t kPartialStep>
class ResampleIIR {
public:
    void Init(float source_fs, float target_fs) {
        blep_.Init(source_fs);
        blep_.SetCutoff(target_fs / 2 * TCoeff::fpass / TCoeff::fstop);
        phase_inc_ = source_fs / target_fs;
    }

    std::vector<float> Process(std::span<float> x) {
        std::vector<float> ret;

        blep_.Reset();

        float phase{};
        float max_one_{};
        size_t rpos{};
        blep_.Add(x[0]);
        while (rpos < x.size() - 1) {
            float const frac = phase;
            float const v = blep_.Get(frac);
            if (std::abs(v) > max_one_) {
                max_one_ = std::abs(v);
            }
            ret.push_back(v);

            phase += phase_inc_;
            size_t new_rpos = rpos + static_cast<size_t>(std::floor(phase));
            phase -= std::floor(phase);

            new_rpos = std::min(new_rpos, x.size() - 1);
            for (size_t i = rpos; i < new_rpos; ++i) {
                blep_.Step();
                blep_.Add(x[i + 1]);
            }
            rpos = new_rpos;
        }

        if (max_one_ == 0.0f) {
            max_one_ = 1.0f;
        }
        else {
            max_one_ = 0.9f / max_one_;
        }
        for (auto& s : ret) {
            s *= max_one_;
        }

        return ret;
    }
private:
    float phase_inc_{};
    signalsmith::blep::EllipticBlep<TCoeff, float, kPartialStep> blep_;
};
}