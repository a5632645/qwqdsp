#pragma once
#include <vector>
#include <span>
#include <memory>

namespace signalsmith::blep {
template<class Sample>
struct EllipticBlep;
}

namespace qwqdsp::fx {
/**
 * @brief holters-parker IIR重采样器，使用Elliptic-blep库实现，移除了高通滤波器系数
 */
class ResampleIIR {
public:
    ResampleIIR();
    ~ResampleIIR();

    /**
     * @param partial_step 越大质量越高
     */
    void Init(float source_fs, float target_fs, size_t partial_step = 255);
    std::vector<float> Process(std::span<float> x);
private:
    float phase_inc_{};
    std::unique_ptr<signalsmith::blep::EllipticBlep<float>> blep_;
};
}