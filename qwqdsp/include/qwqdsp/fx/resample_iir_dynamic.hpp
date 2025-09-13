#pragma once
#include <memory>

namespace signalsmith::blep {
template<class Sample>
struct EllipticBlep;
}

namespace qwqdsp::fx {
/**
 * @brief holters-parker IIR重采样器，使用Elliptic-blep库实现，移除了高通滤波器系数
 */
class ResampleIIRDynamic {
public:
    ResampleIIRDynamic();
    ~ResampleIIRDynamic();

    void Init(float source_fs);
    void Reset() noexcept;

    void Set(float source_fs, float target_fs) noexcept;
    void SetRatio(float ratio) noexcept;
    void SetPitchShift(float shift) noexcept;

    void Push(float x) noexcept;
    float Read() noexcept;
    /**
     * @return false 请调用Push给予更多数据
     *          true 可以调用Read直到返回false
     */
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
    bool first_init_{};
    float source_fs_{};
};
}