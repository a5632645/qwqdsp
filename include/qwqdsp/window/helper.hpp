#pragma once
#include <cassert>
#include <numeric>
#include <span>

namespace qwqdsp_window {
struct Helper {
    static void Normalize(std::span<float> x) noexcept {
        float gain = NormalizeGain(x);
        for (auto& v : x) {
            v *= gain;
        }
    }

    static float NormalizeGain(std::span<const float> x) noexcept {
        float gain = 2.0f / std::accumulate(x.begin(), x.end(), 0.0f);
        return gain;
    }

    static void TWindow(std::span<float> buffer, std::span<const float> window) noexcept {
        assert(buffer.size() == window.size());
        const size_t N = window.size();
        float offset = 0.5f * (static_cast<float>(N) - 1.0f);
        for (size_t k = 0; k < N; ++k) {
            buffer[k] = window[k] * (static_cast<float>(k) - offset);
        }
    }

    static void ZeroPhasePad(std::span<float> output, std::span<const float> input) noexcept {
        assert(output.size() >= input.size());

        std::fill(output.begin(), output.end(), 0.0f);
        size_t La = input.size();
        for (size_t j = 0; j < La / 2; ++j) {
            output[j] = input[j + La / 2];
        }
        const size_t pad = output.size() - La / 2;
        for (size_t j = 0; j < La / 2; ++j) {
            output[pad + j] = input[j];
        }
    }

    static void ZeroPhaseApply(std::span<const float> input, std::span<const float> window,
                               std::span<float> output) noexcept {
        assert(input.size() == window.size());
        assert(output.size() >= window.size());

        std::fill(output.begin(), output.end(), 0.0f);
        size_t La = window.size();
        for (size_t j = 0; j < La / 2; ++j) {
            output[j] = input[j + La / 2] * window[j + La / 2];
        }
        const size_t pad = output.size() - La / 2;
        for (size_t j = 0; j < La / 2; ++j) {
            output[pad + j] = input[j] * window[j];
        }
    }

    static void ZeroPad(std::span<float> output, std::span<const float> input) noexcept {
        assert(output.size() >= input.size());
        auto it = std::copy(input.begin(), input.end(), output.begin());
        std::fill(it, output.end(), 0.0f);
    }
};
} // namespace qwqdsp_window
