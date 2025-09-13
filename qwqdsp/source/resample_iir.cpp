#include "qwqdsp/fx/resample_iir.hpp"
#include "elliptic_blep.hpp"

namespace qwqdsp::fx {

namespace internal {
class ResampleIIRImpl {
public:
    void Init(float source_fs, float target_fs) {
        blep_ = std::make_unique<signalsmith::blep::EllipticBlep<float>>(
            true, source_fs, target_fs * 20000.0f / 44100.0f, 255
        );
        phase_inc_ = source_fs / target_fs;
    }

    std::vector<float> Process(std::span<float> x) {
        std::vector<float> ret;

        blep_->reset();

        float phase{};
        size_t rpos{};
        blep_->add(x[0], 0);
        while (rpos < x.size() - 1) {
            float const frac = phase;
            ret.push_back(blep_->get(frac));

            phase += phase_inc_;
            size_t new_rpos = rpos + static_cast<size_t>(std::floor(phase));
            phase -= std::floor(phase);

            new_rpos = std::min(new_rpos, x.size() - 1);
            for (size_t i = rpos; i < new_rpos; ++i) {
                blep_->step();
                blep_->add(x[i + 1], 0);
            }
            rpos = new_rpos;
        }

        return ret;
    }
private:
    float phase_inc_{};
    std::unique_ptr<signalsmith::blep::EllipticBlep<float>> blep_;
};
} // qwqdsp::fx::internal

ResampleIIR::ResampleIIR() {
    impl_ = std::make_unique<internal::ResampleIIRImpl>();
}

ResampleIIR::~ResampleIIR() = default;

void ResampleIIR::Init(float source_fs, float target_fs) {
    impl_->Init(source_fs, target_fs);
}

std::vector<float> ResampleIIR::Process(std::span<float> x) {
    return impl_->Process(x);
}

}