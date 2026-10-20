/// @brief 自适应分析窗的 NSDF 音高跟踪 + TD-PSOLA autotune

#pragma once

#include "autotune_common.hpp"

namespace autotune_lab {

constexpr int kAdaptiveHop = 128;
constexpr int kAdaptiveMinWindow = 256;
constexpr int kAdaptiveMaxWindow = 2048;
constexpr int kAdaptiveInitialWindow = 1024;

struct AdaptiveInfo {
    int hop = kAdaptiveHop;
    int min_window = 0;
    int max_window = 0;
    double mean_window = 0.0;
    int frames = 0;
};

/// @brief 归一化平方差函数，返回周期（样本）和清晰度
inline bool EstimateNsdf(std::span<const float> block, float fs, float& period, float& clarity) {
    int const N = static_cast<int>(block.size());
    int const min_lag = std::max(2, static_cast<int>(std::lround(fs / kMaxF0)));
    int const max_lag = std::min(N / 2, static_cast<int>(std::lround(fs / kMinF0)));
    if (max_lag <= min_lag + 1) {
        return false;
    }

    std::vector<float> nsdf(static_cast<size_t>(max_lag) + 1, 0.0f);
    for (int lag = min_lag; lag <= max_lag; ++lag) {
        int const length = N - lag;
        float acf = 0.0f;
        float power = 0.0f;
        for (int i = 0; i < length; ++i) {
            float const a = block[static_cast<size_t>(i)];
            float const b = block[static_cast<size_t>(i + lag)];
            acf += a * b;
            power += a * a + b * b;
        }
        nsdf[static_cast<size_t>(lag)] = (power > 1.0e-12f) ? 2.0f * acf / power : 0.0f;
    }

    float best_peak = 0.0f;
    for (int lag = min_lag + 1; lag < max_lag; ++lag) {
        float const value = nsdf[static_cast<size_t>(lag)];
        if (value > nsdf[static_cast<size_t>(lag - 1)] && value >= nsdf[static_cast<size_t>(lag + 1)]) {
            best_peak = std::max(best_peak, value);
        }
    }
    if (best_peak <= 0.0f) {
        return false;
    }

    int best_lag = -1;
    for (int lag = min_lag + 1; lag < max_lag; ++lag) {
        float const value = nsdf[static_cast<size_t>(lag)];
        if (value > nsdf[static_cast<size_t>(lag - 1)] && value >= nsdf[static_cast<size_t>(lag + 1)]
            && value >= 0.9f * best_peak && value > 0.5f) {
            best_lag = lag;
            break;
        }
    }
    if (best_lag < 0) {
        return false;
    }

    float const y0 = nsdf[static_cast<size_t>(best_lag - 1)];
    float const y1 = nsdf[static_cast<size_t>(best_lag)];
    float const y2 = nsdf[static_cast<size_t>(best_lag + 1)];
    float const denom = y0 - 2.0f * y1 + y2;
    float delta = (std::abs(denom) > 1.0e-8f) ? 0.5f * (y0 - y2) / denom : 0.0f;
    delta = std::clamp(delta, -0.5f, 0.5f);

    period = static_cast<float>(best_lag) + delta;
    clarity = y1;
    return clarity > 0.5f;
}

inline PitchTrack TrackPitchAdaptive(std::span<const float> x, float fs, AdaptiveInfo& info) {
    size_t const N = x.size();
    std::vector<PitchSample> samples;
    float prev_period = 0.0f;
    int min_window = std::numeric_limits<int>::max();
    int max_window = 0;
    long long sum_window = 0;
    int frames = 0;

    for (size_t pos = 0; pos < N; pos += kAdaptiveHop) {
        int window = (prev_period > 0.0f)
            ? static_cast<int>(std::lround(2.5f * prev_period))
            : kAdaptiveInitialWindow;
        window = std::clamp(window, kAdaptiveMinWindow, kAdaptiveMaxWindow);

        float period = 0.0f;
        float clarity = 0.0f;
        bool found = false;
        while (window <= kAdaptiveMaxWindow && pos + static_cast<size_t>(window) <= N) {
            auto const block = x.subspan(pos, static_cast<size_t>(window));
            if (BlockRms(block) < kSilenceRms) {
                break;
            }
            if (EstimateNsdf(block, fs, period, clarity)) {
                found = true;
                break;
            }
            if (window == kAdaptiveMaxWindow) {
                break;
            }
            window = std::min(window * 2, kAdaptiveMaxWindow);
        }

        if (!found) {
            continue;
        }

        samples.push_back({pos + static_cast<size_t>(window) / 2, fs / period});
        prev_period = period;
        min_window = std::min(min_window, window);
        max_window = std::max(max_window, window);
        sum_window += window;
        ++frames;
    }

    info = {};
    info.hop = kAdaptiveHop;
    info.min_window = (frames > 0) ? min_window : 0;
    info.max_window = max_window;
    info.mean_window = (frames > 0) ? static_cast<double>(sum_window) / static_cast<double>(frames) : 0.0;
    info.frames = frames;

    PitchTrack track;
    track.f0 = BuildF0Track(N, std::move(samples), kAdaptiveHop);
    track.ratio = BuildRatio(track.f0, fs);
    return track;
}

inline bool RunAdaptive(std::span<const float> input, int fs, std::string const& name) {
    float const fs_f = static_cast<float>(fs);
    AdaptiveInfo info;
    auto const track = TrackPitchAdaptive(input, fs_f, info);
    auto const marks = BuildMarks(input, track.f0, fs_f);
    std::vector<float> output = Synthesize(input, marks, track.ratio, fs_f);

    if (info.frames > 0) {
        std::cout << std::format(
            "  adaptive analysis: hop {} ({:.1f} ms), window {}..{} (mean {:.0f}, {:.1f} ms half)\n",
            info.hop, 1000.0 * static_cast<double>(info.hop) / static_cast<double>(fs),
            info.min_window, info.max_window, info.mean_window,
            500.0 * info.mean_window / static_cast<double>(fs));
    }
    else {
        std::cout << "  adaptive analysis: no voiced frames\n";
    }
    return SaveAndReport(input, output, fs, name);
}

} // namespace autotune_lab
