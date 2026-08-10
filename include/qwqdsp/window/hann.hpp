#pragma once
#include <cmath>
#include <numbers>
#include <span>

namespace qwqdsp_window {
/// @ref https://www.recordingblogs.com/wiki/hann-window
struct Hann {
    // 和分析有关的
    // f = width / N
    static constexpr float kMainlobeWidth = 3.0f;
    static constexpr float kSidelobe = -31.5565f;
    static constexpr float kSidelobeRolloff = -18.0f;

    static void Window(std::span<float> x, bool for_analyze_not_fir) noexcept {
        const size_t N = x.size();
        if (for_analyze_not_fir) {
            for (size_t n = 0; n < N; ++n) {
                const float t = static_cast<float>(n) / static_cast<float>(N);
                x[n] = 0.5f * (1.0f - std::cos(2.0f * std::numbers::pi_v<float> * t));
            }
        }
        else {
            for (size_t n = 0; n < N; ++n) {
                const float t = static_cast<float>(n) / static_cast<float>(N - 1);
                x[n] = 0.5f * (1.0f - std::cos(2.0f * std::numbers::pi_v<float> * t));
            }
        }
    }

    static void ApplyWindow(std::span<float> x, bool for_analyze_not_fir) noexcept {
        const size_t N = x.size();
        if (for_analyze_not_fir) {
            for (size_t n = 0; n < N; ++n) {
                const float t = static_cast<float>(n) / static_cast<float>(N);
                x[n] *= 0.5f * (1.0f - std::cos(2.0f * std::numbers::pi_v<float> * t));
            }
        }
        else {
            for (size_t n = 0; n < N; ++n) {
                const float t = static_cast<float>(n) / static_cast<float>(N - 1);
                x[n] *= 0.5f * (1.0f - std::cos(2.0f * std::numbers::pi_v<float> * t));
            }
        }
    }

    static void DWindow(std::span<float> x) noexcept {
        const size_t N = x.size();
        for (size_t n = 0; n < N; ++n) {
            const float t = static_cast<float>(n) / static_cast<float>(N);
            x[n] = std::numbers::pi_v<float> * std::sin(std::numbers::pi_v<float> * 2 * t);
        }
    }
};
} // namespace qwqdsp_window
