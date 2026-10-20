/// @brief autotune lab 的共享参数、音高轨迹后处理与 TD-PSOLA 合成

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <iostream>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

#include "AudioFile.h"
#include "work_dir.hpp"

#include <qwqdsp/pitch/yin.hpp>

namespace autotune_lab {

constexpr float kMinF0 = 90.0f;
constexpr float kMaxF0 = 1000.0f;
constexpr float kSilenceRms = 1.0e-4f;
constexpr float kCorrection = 1.0f;
constexpr float kRetuneSeconds = 0.02f;
constexpr float kMaxShiftSemitones = 12.0f;

// C 大调: C D E F G A B
constexpr uint16_t kScaleMask = 0x0AB5;

struct PitchTrack {
    std::vector<float> f0;
    std::vector<float> ratio;
};

struct PitchSample {
    size_t center;
    float f0;
};

struct Mark {
    size_t pos;
    float period;
};

struct PitchStats {
    float median_f0 = 0.0f;
    float median_abs_cents = 0.0f;
    int voiced = 0;
};

inline float BlockRms(std::span<const float> x) {
    double sum = 0.0;
    for (float v : x) {
        sum += static_cast<double>(v) * static_cast<double>(v);
    }
    return static_cast<float>(std::sqrt(sum / static_cast<double>(x.size())));
}

inline int FindNearestNote(float midi, uint16_t mask) {
    int const base = static_cast<int>(std::lround(midi));
    for (int d = 0; d <= 6; ++d) {
        int const up = base + d;
        int const down = base - d;
        int const up_pc = ((up % 12) + 12) % 12;
        int const down_pc = ((down % 12) + 12) % 12;
        if (mask & (1u << up_pc)) {
            return up;
        }
        if (d != 0 && (mask & (1u << down_pc))) {
            return down;
        }
    }
    return base;
}

inline float ClampPeriod(float period, float fs) {
    return std::clamp(period, fs / kMaxF0, fs / kMinF0);
}

inline float PeriodAt(std::span<const float> f0, size_t pos, float fs) {
    float const period = (pos < f0.size() && f0[pos] > 0.0f) ? fs / f0[pos] : fs / 200.0f;
    return ClampPeriod(period, fs);
}

inline std::vector<float> BuildF0Track(size_t num_samples, std::vector<PitchSample> samples, int hop) {
    std::vector<float> f0(num_samples, 0.0f);
    if (samples.empty()) {
        return f0;
    }

    std::sort(samples.begin(), samples.end(), [](PitchSample const& a, PitchSample const& b) {
        return a.center < b.center;
    });
    std::vector<PitchSample> unique;
    unique.reserve(samples.size());
    for (auto const& s : samples) {
        if (!unique.empty() && s.center == unique.back().center) {
            unique.back() = s;
        }
        else {
            unique.push_back(s);
        }
    }

    for (auto const& s : unique) {
        f0[s.center] = s.f0;
    }
    for (size_t k = 0; k + 1 < unique.size(); ++k) {
        if (unique[k + 1].center - unique[k].center > static_cast<size_t>(8 * hop)) {
            continue;
        }
        for (size_t i = unique[k].center + 1; i < unique[k + 1].center; ++i) {
            float const t = static_cast<float>(i - unique[k].center)
                          / static_cast<float>(unique[k + 1].center - unique[k].center);
            f0[i] = unique[k].f0 + t * (unique[k + 1].f0 - unique[k].f0);
        }
    }
    return f0;
}

inline std::vector<float> BuildRatio(std::span<const float> f0, float fs) {
    std::vector<float> ratio(f0.size(), 1.0f);
    float const alpha = 1.0f - std::exp(-1.0f / (fs * kRetuneSeconds));
    float smoothed = 1.0f;
    float const max_ratio = std::exp2(kMaxShiftSemitones / 12.0f);
    for (size_t i = 0; i < f0.size(); ++i) {
        float target = 1.0f;
        if (f0[i] > 0.0f) {
            float const midi = 69.0f + 12.0f * std::log2(f0[i] / 440.0f);
            int const nearest = FindNearestNote(midi, kScaleMask);
            float const target_midi = midi + kCorrection * (static_cast<float>(nearest) - midi);
            target = std::exp2((target_midi - midi) / 12.0f);
            target = std::clamp(target, 1.0f / max_ratio, max_ratio);
        }
        smoothed += alpha * (target - smoothed);
        ratio[i] = smoothed;
    }
    return ratio;
}

inline size_t SnapPeak(std::span<const float> x, size_t center, size_t radius, float sign) {
    if (x.size() < 3) {
        return center;
    }
    size_t const lo = (center > radius) ? center - radius : 1;
    size_t const hi = std::min(x.size() - 2, center + radius);
    size_t best = std::clamp(center, lo, hi);
    float best_value = -std::numeric_limits<float>::infinity();
    for (size_t i = lo; i <= hi; ++i) {
        float const value = sign * x[i];
        if (value > best_value) {
            best_value = value;
            best = i;
        }
    }
    return best;
}

inline std::vector<Mark> BuildMarks(std::span<const float> x, std::span<const float> f0, float fs) {
    std::vector<Mark> marks;
    size_t const N = x.size();
    if (N < 3) {
        return marks;
    }

    size_t start = 0;
    while (start < N && f0[start] <= 0.0f) {
        ++start;
    }
    if (start >= N) {
        return marks;
    }

    float const sign = x[start] >= 0.0f ? 1.0f : -1.0f;
    float const default_period = ClampPeriod(fs / 200.0f, fs);
    size_t const default_radius = std::max<size_t>(2, static_cast<size_t>(std::lround(default_period * 0.3f)));
    size_t anchor = SnapPeak(x, start, default_radius, sign);

    std::vector<Mark> head;
    size_t pos = anchor;
    while (true) {
        float const period = PeriodAt(f0, pos, fs);
        float const predicted = static_cast<float>(pos) - period;
        if (predicted < 0.0f) {
            break;
        }
        size_t const radius = std::max<size_t>(2, static_cast<size_t>(std::lround(period * 0.3f)));
        size_t const mark = SnapPeak(x, static_cast<size_t>(std::lround(predicted)), radius, sign);
        if (mark >= pos) {
            break;
        }
        head.push_back({mark, period});
        pos = mark;
    }
    std::reverse(head.begin(), head.end());
    marks = std::move(head);

    marks.push_back({anchor, PeriodAt(f0, anchor, fs)});
    pos = anchor;
    float const max_period = fs / kMinF0;
    while (true) {
        float const period = PeriodAt(f0, pos, fs);
        float const predicted = static_cast<float>(pos) + period;
        if (predicted + max_period >= static_cast<float>(N)) {
            break;
        }
        size_t const radius = std::max<size_t>(2, static_cast<size_t>(std::lround(period * 0.3f)));
        size_t mark = SnapPeak(x, static_cast<size_t>(std::lround(predicted)), radius, sign);
        if (mark <= pos) {
            mark = pos + 1;
        }
        if (mark >= N) {
            break;
        }
        marks.push_back({mark, period});
        pos = mark;
    }
    return marks;
}

inline void AddGrain(std::span<const float> x, size_t src, float left, float right,
                     std::span<float> out, std::span<float> norm, size_t dst) {
    int const l = std::max(1, static_cast<int>(std::lround(left)));
    int const r = std::max(1, static_cast<int>(std::lround(right)));
    for (int i = 0; i < l; ++i) {
        int const si = static_cast<int>(src) - l + i;
        int const di = static_cast<int>(dst) - l + i;
        if (si < 0 || si >= static_cast<int>(x.size()) || di < 0 || di >= static_cast<int>(out.size())) {
            continue;
        }
        float const w = 0.5f * (1.0f - std::cos(std::numbers::pi_v<float> * static_cast<float>(i) / static_cast<float>(l)));
        out[di] += x[si] * w;
        norm[di] += w;
    }
    for (int i = 0; i < r; ++i) {
        int const si = static_cast<int>(src) + i;
        int const di = static_cast<int>(dst) + i;
        if (si < 0 || si >= static_cast<int>(x.size()) || di < 0 || di >= static_cast<int>(out.size())) {
            continue;
        }
        float const w = 0.5f * (1.0f + std::cos(std::numbers::pi_v<float> * static_cast<float>(i) / static_cast<float>(r)));
        out[di] += x[si] * w;
        norm[di] += w;
    }
}

inline std::vector<float> Synthesize(std::span<const float> x, std::span<const Mark> marks,
                                     std::span<const float> ratio, float fs) {
    size_t const N = x.size();
    std::vector<float> out(N, 0.0f);
    std::vector<float> norm(N, 0.0f);
    if (marks.size() < 2) {
        return std::vector<float>(x.begin(), x.end());
    }

    float const min_period = fs / kMaxF0;
    float const max_period = fs / kMinF0;
    size_t const last = marks.size() - 1;

    float syn = 0.0f;
    size_t cursor = 0;
    while (syn < static_cast<float>(N)) {
        float const src_time = syn;
        while (cursor + 1 < marks.size() && static_cast<float>(marks[cursor + 1].pos) <= src_time) {
            ++cursor;
        }
        size_t best = cursor;
        if (cursor + 1 < marks.size()) {
            float const d0 = std::abs(static_cast<float>(marks[cursor].pos) - src_time);
            float const d1 = std::abs(static_cast<float>(marks[cursor + 1].pos) - src_time);
            if (d1 < d0) {
                best = cursor + 1;
            }
        }

        float left = (best > 0) ? static_cast<float>(marks[best].pos - marks[best - 1].pos) : marks[best].period;
        float right = (best < last) ? static_cast<float>(marks[best + 1].pos - marks[best].pos) : marks[best].period;
        left = std::clamp(left, min_period, max_period * 2.0f);
        right = std::clamp(right, min_period, max_period * 2.0f);

        AddGrain(x, marks[best].pos, left, right, out, norm,
                 static_cast<size_t>(std::lround(syn)));

        size_t const ri = std::min(static_cast<size_t>(std::lround(src_time)), N - 1);
        float const r = std::max(ratio[ri], 0.01f);
        float step = marks[best].period / r;
        step = std::clamp(step, min_period * 0.75f, max_period * 1.25f);
        syn += step;
    }

    for (size_t i = 0; i < N; ++i) {
        if (norm[i] > 1.0e-8f) {
            out[i] /= norm[i];
        }
    }
    return out;
}

inline PitchStats AnalyzePitch(std::span<const float> x, float fs) {
    constexpr int kReportBlock = 1024;
    constexpr int kReportHop = 256;
    PitchStats stats;
    if (x.size() < static_cast<size_t>(kReportBlock)) {
        return stats;
    }

    qwqdsp_pitch::Yin yin;
    yin.Init(fs, kReportBlock);
    yin.SetMinPitch(kMinF0);
    yin.SetMaxPitch(kMaxF0);
    yin.SetThreshold(0.2f);

    std::vector<float> block(kReportBlock);
    std::vector<float> pitches;
    std::vector<float> cents;
    for (size_t pos = 0; pos + static_cast<size_t>(kReportBlock) <= x.size(); pos += kReportHop) {
        std::copy_n(x.begin() + static_cast<std::ptrdiff_t>(pos), kReportBlock, block.begin());
        if (BlockRms(block) < kSilenceRms) {
            continue;
        }

        yin.Process(block);
        auto const pitch = yin.GetPitch();
        if (pitch.pitch_hz < kMinF0 || pitch.pitch_hz > kMaxF0
            || pitch.non_period_ratio >= 0.25f) {
            continue;
        }

        float const midi = 69.0f + 12.0f * std::log2(pitch.pitch_hz / 440.0f);
        int const nearest = FindNearestNote(midi, kScaleMask);
        pitches.push_back(pitch.pitch_hz);
        cents.push_back(100.0f * std::abs(midi - static_cast<float>(nearest)));
    }

    if (pitches.empty()) {
        return stats;
    }
    std::sort(pitches.begin(), pitches.end());
    std::sort(cents.begin(), cents.end());
    stats.median_f0 = pitches[pitches.size() / 2];
    stats.median_abs_cents = cents[cents.size() / 2];
    stats.voiced = static_cast<int>(pitches.size());
    return stats;
}

inline bool LoadMono(std::string const& path, std::vector<float>& out, int& fs) {
    AudioFile<float> file;
    if (!file.load(path)) {
        return false;
    }
    fs = static_cast<int>(file.getSampleRate());
    int const num_channels = file.getNumChannels();
    int const num_samples = file.getNumSamplesPerChannel();
    out.assign(static_cast<size_t>(num_samples), 0.0f);
    if (num_channels == 1) {
        out = file.samples[0];
    }
    else {
        float const inv = 1.0f / static_cast<float>(num_channels);
        for (int i = 0; i < num_samples; ++i) {
            float sum = 0.0f;
            for (int ch = 0; ch < num_channels; ++ch) {
                sum += file.samples[ch][i];
            }
            out[static_cast<size_t>(i)] = sum * inv;
        }
    }
    return true;
}

inline bool SaveMono(std::string const& path, std::span<const float> x, int fs) {
    AudioFile<float> file;
    file.setBitDepth(32);
    file.setNumChannels(1);
    file.setNumSamplesPerChannel(static_cast<int>(x.size()));
    file.setSampleRate(static_cast<uint32_t>(fs));
    file.samples[0].assign(x.begin(), x.end());
    return file.save(path);
}

inline bool CheckOutput(std::span<const float> input, std::span<const float> output,
                        int fs, std::string const& name) {
    bool ok = output.size() == input.size();
    float peak = 0.0f;
    double sum = 0.0;
    bool finite = true;
    for (float v : output) {
        if (!std::isfinite(v)) {
            finite = false;
        }
        peak = std::max(peak, std::abs(v));
        sum += static_cast<double>(v) * static_cast<double>(v);
    }
    float const rms = output.empty() ? 0.0f : static_cast<float>(std::sqrt(sum / static_cast<double>(output.size())));
    if (!finite || peak < 1.0e-6f) {
        ok = false;
    }

    auto const in_stats = AnalyzePitch(input, static_cast<float>(fs));
    auto const out_stats = AnalyzePitch(output, static_cast<float>(fs));
    std::cout << std::format("{}: len {} -> {}, peak {:.3f}, rms {:.4f}, voiced {} -> {}\n",
                             name, input.size(), output.size(), peak, rms,
                             in_stats.voiced, out_stats.voiced);
    std::cout << std::format("  f0 median {:.2f} -> {:.2f} Hz, |cents| median {:.2f} -> {:.2f}\n",
                             in_stats.median_f0, out_stats.median_f0,
                             in_stats.median_abs_cents, out_stats.median_abs_cents);
    return ok;
}

inline bool SaveAndReport(std::span<const float> input, std::span<const float> output,
                          int fs, std::string const& name) {
    bool const ok = CheckOutput(input, output, fs, name);
    std::filesystem::create_directories(qwqdsp_support::GetOutputDir());
    std::string const path = qwqdsp_support::OutputFile(name);
    if (!SaveMono(path, output, fs)) {
        std::cout << std::format("  failed to save {}\n", path);
        return false;
    }
    std::cout << std::format("  saved {}\n", path);
    return ok;
}

} // namespace autotune_lab
