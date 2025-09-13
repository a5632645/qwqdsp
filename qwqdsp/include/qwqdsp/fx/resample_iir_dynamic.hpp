#pragma once
#include <memory>

namespace signalsmith::blep {
template<class Sample>
struct EllipticBlep;
}

namespace qwqdsp::fx {

class ResampleIIRDynamic {
public:
    ResampleIIRDynamic();
    ~ResampleIIRDynamic();

    void Init(float source_fs, float max_cutoff);
    void Reset() noexcept;
    void Set(float source_fs, float target_fs) noexcept;
    void Push(float x) noexcept;
    float Read() noexcept;
    bool IsReady() const noexcept;
private:
    std::unique_ptr<signalsmith::blep::EllipticBlep<float>> blep_;
    float buffer_[64]{};
    size_t wpos_{};
    size_t rpos_{};
    size_t buffer_size_{};

    float phase_inc_{};
    float phase_{};
    size_t need_{};
};
}