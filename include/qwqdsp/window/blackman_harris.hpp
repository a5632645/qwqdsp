#pragma once
#include <cmath>
#include <numbers>
#include <span>

namespace qwqdsp_window {
/// @ref https://www.recordingblogs.com/wiki/blackman-harris-window
struct BlackmanHarris {
    // 和分析有关的
    // f = width / N
    static constexpr float kMainlobeWidth = 4.0f;
    static constexpr float kSidelobe = -92.0f;
    static constexpr float kSidelobeRolloff = -14.4f;

    template <bool period>
    static float Get(int n, int L) noexcept {
        constexpr float a0 = 0.35875f;
        constexpr float a1 = 0.48829f;
        constexpr float a2 = 0.14128f;
        constexpr float a3 = 0.01168f;
        constexpr float twopi = std::numbers::pi_v<float> * 2;

        float t;
        if constexpr (period) {
            t = static_cast<float>(n) / static_cast<float>(L);
        }
        else {
            t = static_cast<float>(n) / static_cast<float>(L - 1);
        }

        return a0 - a1 * std::cos(twopi * t) + a2 * std::cos(twopi * 2 * t) - a3 * std::cos(twopi * 3 * t);
    }

    static void Window(std::span<float> x, bool for_analyze_not_fir) noexcept {
        const size_t L = x.size();
        if (for_analyze_not_fir) {
            for (size_t n = 0; n < x.size(); ++n) {
                x[n] = Get<true>(n, L);
            }
        }
        else {
            for (size_t n = 0; n < x.size(); ++n) {
                const float t = n / (L - 1.0f);
                x[n] = Get<false>(n, L);
            }
        }
    }

    static void ApplyWindow(std::span<float> x, bool for_analyze_not_fir) noexcept {
        const size_t L = x.size();
        if (for_analyze_not_fir) {
            for (size_t n = 0; n < x.size(); ++n) {
                x[n] *= Get<true>(n, L);
            }
        }
        else {
            for (size_t n = 0; n < x.size(); ++n) {
                const float t = n / (L - 1.0f);
                x[n] *= Get<false>(n, L);
            }
        }
    }

    static void DWindow(std::span<float> x) noexcept {
        constexpr float a0 = 0.35875f;
        constexpr float a1 = 0.48829f;
        constexpr float a2 = 0.14128f;
        constexpr float a3 = 0.01168f;
        constexpr float twopi = std::numbers::pi_v<float> * 2;

        size_t L = x.size();
        for (size_t n = 0; n < L; ++n) {
            float t = static_cast<float>(n) / static_cast<float>(L);
            x[n] = a0 - twopi * a1 * std::cos(twopi * t) + twopi * 2 * a2 * std::cos(twopi * 2 * t)
                 - a3 * twopi * 3 * std::cos(twopi * 3 * t);
        }
    }
};
} // namespace qwqdsp_window
