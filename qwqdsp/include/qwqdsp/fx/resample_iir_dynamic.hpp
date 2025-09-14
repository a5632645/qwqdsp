#pragma once
#include <memory>
#include "../source/elliptic_blep.hpp"

namespace qwqdsp::fx {
/**
 * @brief holters-parker IIR重采样器，使用Elliptic-blep库实现，移除了高通滤波器系数
 */
template<class TCoeff, size_t kPartialSteps>
class ResampleIIRDynamic {
public:
    void Init(float source_fs) {
        source_fs_ = source_fs;
        blep_.Init(source_fs);
        Reset();
    }
    void Reset() noexcept {
        first_init_ = true;
        need_ = 1;
        wpos_ = 0;
        rpos_ = 0;
        phase_ = 0;
        buffer_size_ = 0;
        blep_.Reset();
    }

    void Set(float source_fs, float target_fs) noexcept {
        SetRatio(source_fs / target_fs);
    }
    void SetRatio(float ratio) noexcept {
        if (!first_init_) {
            float const last_phase = static_cast<float>(need_) + phase_ - phase_inc_;
            float new_phase = last_phase + ratio;
            need_ = static_cast<size_t>(std::floor(new_phase));
            phase_ = new_phase - std::floor(new_phase);
        }
        phase_inc_ = ratio;

        float cutoff = source_fs_ * 0.5f;
        if (ratio > 1.0f) {
            cutoff /= ratio;
        }
        blep_.SetCutoff(cutoff * TCoeff::fpass / TCoeff::fstop);
    }
    void SetPitchShift(float shift) noexcept {
        SetRatio(std::exp2(shift / 12.0f));
    }

    void Push(float x) noexcept {
        buffer_[wpos_++] = x;
        wpos_ &= 63;
        ++buffer_size_;
    }
    float Read() noexcept {
        first_init_ = false;

        for (size_t i = 0; i < need_; ++i) {
            blep_.Step();
            blep_.Add(buffer_[rpos_++], 0);
            rpos_ &= 63;
        }
        buffer_size_ -= need_;

        float const frac = phase_;
        phase_ += phase_inc_;
        need_ = static_cast<size_t>(std::floor(phase_));
        phase_ -= std::floor(phase_);

        return blep_.Get(frac);
    }
    /**
     * @return false 请调用Push给予更多数据
     *          true 可以调用Read直到返回false
     */
    bool IsReady() const noexcept {
        return buffer_size_ >= need_;
    }
private:
    signalsmith::blep::EllipticBlep<TCoeff, float, kPartialSteps> blep_;
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