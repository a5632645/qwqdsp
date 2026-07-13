#pragma once
#include <span>

namespace qwqdsp_support {

struct AudioOps {
    static void Normalize(std::span<float> x_vec) noexcept {
        float max{};
        for (float x : x_vec) {
            max = std::max(max, std::abs(x));
        }

        if (max < 1e-10f)
            return;

        float g = 1.0f / max;
        for (float& x : x_vec) {
            x *= g;
        }
    }
};

} // namespace qwqdsp_support
