#pragma once
#include <cmath>
#include <numbers>
#include <span>

namespace qwqdsp_window {
/// @ref https://www.recordingblogs.com/wiki/hamming-window
struct Hamming {
    // 和分析有关的
    // f = width / N
    static constexpr float kMainlobeWidth = 2.0f;
    static constexpr float kSidelobe = -43.7547f;
    static constexpr float kSidelobeRolloff = -6.0f;

    static void Window(std::span<float> x, bool for_analyze_not_fir) noexcept {
        const size_t N = x.size();
        if (for_analyze_not_fir) {
            for (size_t n = 0; n < N; ++n) {
                const float t = n / static_cast<float>(N) - 0.5f;
                x[n] = 0.53836f + 0.46164f * std::cos(std::numbers::pi_v<float> * 2 * t);
            }
        }
        else {
            for (size_t n = 0; n < N; ++n) {
                const float t = n / (N - 1.0f) - 0.5f;
                x[n] = 0.53836f + 0.46164f * std::cos(std::numbers::pi_v<float> * 2 * t);
            }
        }
    }

    static void ApplyWindow(std::span<float> x, bool for_analyze_not_fir) noexcept {
        const size_t N = x.size();
        if (for_analyze_not_fir) {
            for (size_t n = 0; n < N; ++n) {
                const float t = n / static_cast<float>(N) - 0.5f;
                x[n] *= 0.53836f + 0.46164f * std::cos(std::numbers::pi_v<float> * 2 * t);
            }
        }
        else {
            for (size_t n = 0; n < N; ++n) {
                const float t = n / (N - 1.0f) - 0.5f;
                x[n] *= 0.53836f + 0.46164f * std::cos(std::numbers::pi_v<float> * 2 * t);
            }
        }
    }

    static void DWindow(std::span<float> x) noexcept {
        const size_t N = x.size();
        for (size_t n = 0; n < N; ++n) {
            const float t = n / static_cast<float>(N) - 0.5f;
            x[n] = -0.46164f * std::numbers::pi_v<float> * 2 * std::sin(std::numbers::pi_v<float> * 2 * t);
        }
    }
};
} // namespace qwqdsp_window
