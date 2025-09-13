#pragma once
#include <vector>
#include <span>
#include <memory>

namespace qwqdsp::fx {

namespace internal {
class ResampleIIRImpl;
}

class ResampleIIR {
public:
    ResampleIIR();
    ~ResampleIIR();

    void Init(float source_fs, float target_fs);
    std::vector<float> Process(std::span<float> x);
private:
    std::unique_ptr<internal::ResampleIIRImpl> impl_;
};
}