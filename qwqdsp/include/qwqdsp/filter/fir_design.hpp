#pragma once
#include <span>

namespace qwqdsp::filter {
struct FirDesign {
    /**
     * @param cutoff [0, pi]
     * @param transist_width [0, pi]
     * @param stop_band_weight [1, 420]
     */
    static void LowpassLeastSquare(
        std::span<float> x,
        float cutoff,
        float transist_width, float stop_band_weight
    );
};
}