#pragma once
#include <cstddef>
#include <numbers>
#include <span>
#include <cmath>
#include <cassert>

namespace qwqdsp::window {
struct Kaiser {
    // 和分析有关的
    // f = width / N
    static float MainLobeWidth(float beta) noexcept {
        float a = beta / std::numbers::pi_v<float>;
        return std::sqrt(1.0f + a * a);
    }
    static constexpr float kSidelobeRolloff = -6.0f;
    // 和滤波器设计有关的
    // 卷积之后第一个旁瓣的大小
    // qwqfixme: 可能需要补充信息
    // static constexpr float kStopband = -53.0f;
    // static constexpr float kTransmit = 3.3f;

    static void Window(std::span<float> window, float beta, bool for_analyze_not_fir) noexcept {
        const size_t N = window.size();
        auto inc = 2.0f / (static_cast<float>(N) - 1.0f);
        if (for_analyze_not_fir) {
            inc = 2.0f / static_cast<float>(N);
        }
        auto down = 1.0f / std::cyl_bessel_i(0.0f, beta);
        for (size_t i = 0; i < N; ++i) {
            auto t = -1.0f + static_cast<float>(i) * inc;
            auto arg = std::sqrt(1.0f - t * t);
            window[i] = std::cyl_bessel_i(0.0f, beta * arg) * down;
        }
    }

    static void ApplyWindow(std::span<float> x, float beta, bool for_analyze_not_fir) noexcept {
        const size_t N = x.size();
        auto inc = 2.0f / (static_cast<float>(N) - 1.0f);
        if (for_analyze_not_fir) {
            inc = 2.0f / static_cast<float>(N);
        }
        auto down = 1.0f / std::cyl_bessel_i(0.0f, beta);
        for (size_t i = 0; i < N; ++i) {
            auto t = -1.0f + static_cast<float>(i) * inc;
            auto arg = std::sqrt(1.0f - t * t);
            x[i] *= std::cyl_bessel_i(0.0f, beta * arg) * down;
        }
    }

    static void Window(std::span<float> window, std::span<float> dwindow, float beta) noexcept {
        constexpr auto kTimeDelta = 0.001f;
        const size_t N = window.size();
        auto inc = 2.0f / static_cast<float>(N);
        auto down = 1.0f / std::cyl_bessel_i(0.0f, beta);
        for (size_t i = 0; i < N; ++i) {
            auto t = -1.0f + static_cast<float>(i) * inc;
            auto arg = std::sqrt(1.0f - t * t);
            window[i] = std::cyl_bessel_i(0.0f, beta * arg) * down;
            if (i == 0) {
                dwindow.front() = (std::cyl_bessel_i(0.0f, beta * std::sqrt(1.0f - (t + kTimeDelta) * (t + kTimeDelta))) * down - window.front()) / kTimeDelta;
            }
            else if (i == N - 1) {
                dwindow.back() = (std::cyl_bessel_i(0.0f, beta * std::sqrt(1.0f - (t - kTimeDelta) * (t - kTimeDelta))) * down - window.back()) / -kTimeDelta;
            }
            else {
                dwindow[i] = std::cyl_bessel_i(1.0f, beta * arg) * beta * (-t / arg) * down;
            }
        }
    }

    /**
     * @param side_lobe >0
     * @ref https://ww2.mathworks.cn/help/signal/ref/kaiser.html
     */
    static float Beta(float side_lobe) noexcept {
        assert(side_lobe > 0);
        if (side_lobe < 21.0f) {
            return 0.0f;
        }
        else if (side_lobe <= 50.0f) {
            return 0.5842f * std::pow(side_lobe - 21.0f, 0.4f) 
                           + 0.07886f * (side_lobe - 21.0f);
        }
        else {
            return 0.1102f * (side_lobe - 8.7f);
        }
    }
};
}