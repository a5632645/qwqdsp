#include "qwqdsp/fx/resample_iir_dynamic.hpp"
#include "elliptic_blep.hpp"

namespace qwqdsp::fx {
ResampleIIRDynamic::ResampleIIRDynamic() = default;
ResampleIIRDynamic::~ResampleIIRDynamic() = default;

void ResampleIIRDynamic::Init(float source_fs, float max_cutoff) {
    blep_ = std::make_unique<signalsmith::blep::EllipticBlep<float>>(
        true, source_fs, max_cutoff, 255
    );
    Reset();
}

void ResampleIIRDynamic::Set(float source_fs, float target_fs) noexcept {
    phase_inc_ = source_fs / target_fs;
}

void ResampleIIRDynamic::Reset() noexcept {
    need_ = 1;
    wpos_ = 0;
    rpos_ = 0;
    phase_ = 0;
    buffer_size_ = 0;
    blep_->reset();
}

void ResampleIIRDynamic::Push(float x) noexcept {
    buffer_[wpos_++] = x;
    wpos_ &= 63;
    ++buffer_size_;
}

float ResampleIIRDynamic::Read() noexcept {
    for (size_t i = 0; i < need_; ++i) {
        blep_->step();
        blep_->add(buffer_[rpos_++], 0);
        rpos_ &= 63;
    }
    buffer_size_ -= need_;

    float const frac = phase_;
    phase_ += phase_inc_;
    need_ = static_cast<size_t>(std::floor(phase_));
    phase_ -= std::floor(phase_);

    return blep_->get(frac);
}

bool ResampleIIRDynamic::IsReady() const noexcept {
    return buffer_size_ >= need_;
}
}