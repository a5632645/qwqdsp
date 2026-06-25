#pragma once
#include <cmath>
#include <numbers>
#include <span>

namespace qwqdsp_window {
/// @ref https://www.recordingblogs.com/wiki/blackman-window
struct Blackman {
    // 和分析有关的
    // f = width / N
    static constexpr float kMainlobeWidth = 3.0f;
    static constexpr float kSidelobe = -58.2336f;
    static constexpr float kSidelobeRolloff = -18.0f;

    static void Window(std::span<float> x, bool for_analyze_not_fir) noexcept {
        const size_t N = x.size();
        if (for_analyze_not_fir) {
            for (size_t n = 0; n < x.size(); ++n) {
                const float t = n / static_cast<float>(N);
                x[n] = Get<true>(n, N);
            }
        }
        else {
            for (size_t n = 0; n < x.size(); ++n) {
                const float t = n / (N - 1.0f);
                x[n] = Get<false>(n, N);
            }
        }
    }

    static void ApplyWindow(std::span<float> x, bool for_analyze_not_fir) noexcept {
        const size_t N = x.size();
        constexpr float twopi = std::numbers::pi_v<float> * 2;
        if (for_analyze_not_fir) {
            for (size_t n = 0; n < x.size(); ++n) {
                const float t = n / static_cast<float>(N);
                x[n] *= Get<true>(n, N);
            }
        }
        else {
            for (size_t n = 0; n < x.size(); ++n) {
                const float t = n / (N - 1.0f);
                x[n] *= Get<false>(n, N);
            }
        }
    }

    static void DWindow(std::span<float> x) noexcept {
        const size_t N = x.size();
        constexpr float twopi = std::numbers::pi_v<float> * 2;
        constexpr float a1 = 0.496562f;
        constexpr float a2 = 0.076849f;
        for (size_t n = 0; n < N; ++n) {
            const float t = static_cast<float>(n) / N;
            x[n] = a1 * twopi * std::sin(twopi * t) - a2 * twopi * 2 * std::sin(twopi * 2 * t);
        }
    }
private:
    template <bool period>
    static float Get(int n, int L) noexcept {
        constexpr float a0 = 0.42659f;
        constexpr float a1 = 0.496562f;
        constexpr float a2 = 0.076849f;
        constexpr float twopi = std::numbers::pi_v<float> * 2;

        float t;
        if constexpr (period) {
            t = static_cast<float>(n) / static_cast<float>(L);
        }
        else {
            t = static_cast<float>(n) / static_cast<float>(L - 1);
        }

        return a0 - a1 * std::cos(twopi * t) + a2 * std::cos(twopi * 2 * t);
    }
};
} // namespace qwqdsp_window
