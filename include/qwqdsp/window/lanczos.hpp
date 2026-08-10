#pragma once
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>

namespace qwqdsp_window {
/// @ref https://www.recordingblogs.com/wiki/lanczos-window
struct Lanczos {
    // 和分析有关的
    // f = width / N
    static constexpr float kMainlobeWidth = 1.625f;
    static constexpr float kSidelobe = -26.5935f;
    static constexpr float kSidelobeRolloff = -12.0f;

    static void Window(std::span<float> x, bool for_analyze_not_fir) noexcept {
        if (for_analyze_not_fir) {
            const size_t N = x.size();
            for (size_t i = 0; i < N; ++i) {
                x[i] = Sinc(2.0f * static_cast<float>(i) / static_cast<float>(N) - 1.0f);
            }
        }
        else {
            const size_t N = x.size();
            for (size_t i = 0; i < N; ++i) {
                x[i] = Sinc(2.0f * static_cast<float>(i) / static_cast<float>(N - 1) - 1.0f);
            }
        }
    }

    static void ApplyWindow(std::span<float> x, bool for_analyze_not_fir) noexcept {
        if (for_analyze_not_fir) {
            const size_t N = x.size();
            for (size_t i = 0; i < N; ++i) {
                x[i] *= Sinc(2.0f * static_cast<float>(i) / static_cast<float>(N) - 1.0f);
            }
        }
        else {
            const size_t N = x.size();
            for (size_t i = 0; i < N; ++i) {
                x[i] *= Sinc(2.0f * static_cast<float>(i) / static_cast<float>(N - 1) - 1.0f);
            }
        }
    }

    static void DWindow(std::span<float> x) noexcept {
        const size_t N = x.size();
        constexpr float pi = std::numbers::pi_v<float>;
        for (size_t i = 0; i < N; ++i) {
            float t = 2.0f * static_cast<float>(i) / static_cast<float>(N) - 1.0f;
            if (t == 0.0f) {
                x[i] = 0.0f;
            }
            else {
                x[i] = (pi * t * std::cos(pi * t) - std::sin(pi * t)) / (pi * t * t);
            }
        }
    }
private:
    static float Sinc(float x) noexcept {
        x *= std::numbers::pi_v<float>;
        if (x == 0.0f) {
            return 1.0f;
        }
        else {
            return std::sin(x) / x;
        }
    }
};
} // namespace qwqdsp_window
