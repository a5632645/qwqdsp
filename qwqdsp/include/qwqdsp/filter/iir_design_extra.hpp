#pragma once
#include "qwqdsp/filter/iir_design.hpp"

namespace qwqdsp {
class IIRDesignExtra {
public:
    using ZPK = IIRDesign::ZPK;
    static constexpr auto pi = IIRDesign::pi;

    /**
     * @param atten (0, 1)
     */
    static double ButterworthAttenGain(std::span<ZPK> ret, size_t num_filter, double atten) {
        return ButterworthAtten(ret, num_filter, (1.0f - atten * atten) / atten);
    }

    /**
     * @param atten >0
     */
    static double ButterworthAttenDb(std::span<ZPK> ret, size_t num_filter, double atten) {
        return ButterworthAtten(ret, num_filter, std::pow(10.0, atten / 10.0) - 1.0);
    }
private:
    static double ButterworthAtten(std::span<ZPK> ret, size_t num_filter, double square_epsi) {
        assert(ret.size() >= num_filter);

        double const g = 1.0 / std::pow(square_epsi, 0.25 / num_filter);
        size_t const n = 2 * num_filter;
        size_t i = 0;
        for (size_t k = 1; k <= num_filter; ++k) {
            double phi = (2.0 * k - 1.0) * pi / (2.0 * n);
            ret[i].p = g * std::complex{-std::sin(phi), std::cos(phi)};
            ++i;
        }
        return g * g;
    }
};
}