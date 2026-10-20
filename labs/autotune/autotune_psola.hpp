/// @brief 固定分析窗的 YIN + TD-PSOLA autotune（基线）

#pragma once

#include "autotune_common.hpp"

namespace autotune_lab {

constexpr int kPsolaBlock = 1024;
constexpr int kPsolaHop = 256;

inline PitchTrack TrackPitchPsola(std::span<const float> x, float fs) {
    size_t const N = x.size();
    qwqdsp_pitch::Yin yin;
    yin.Init(fs, kPsolaBlock);
    yin.SetMinPitch(kMinF0);
    yin.SetMaxPitch(kMaxF0);
    yin.SetThreshold(0.2f);

    std::vector<float> block(kPsolaBlock);
    std::vector<PitchSample> samples;
    for (size_t pos = 0; pos + static_cast<size_t>(kPsolaBlock) <= N; pos += kPsolaHop) {
        std::copy_n(x.begin() + static_cast<std::ptrdiff_t>(pos), kPsolaBlock, block.begin());
        if (BlockRms(block) < kSilenceRms) {
            continue;
        }

        yin.Process(block);
        auto const pitch = yin.GetPitch();
        if (pitch.pitch_hz < kMinF0 || pitch.pitch_hz > kMaxF0
            || pitch.non_period_ratio >= 0.25f) {
            continue;
        }
        samples.push_back({pos + static_cast<size_t>(kPsolaBlock) / 2, pitch.pitch_hz});
    }

    PitchTrack track;
    track.f0 = BuildF0Track(N, std::move(samples), kPsolaHop);
    track.ratio = BuildRatio(track.f0, fs);
    return track;
}

inline bool RunPsola(std::span<const float> input, int fs, std::string const& name) {
    float const fs_f = static_cast<float>(fs);
    auto const track = TrackPitchPsola(input, fs_f);
    auto const marks = BuildMarks(input, track.f0, fs_f);
    std::vector<float> output = Synthesize(input, marks, track.ratio, fs_f);

    std::cout << std::format("  psola analysis: block {} ({:.1f} ms), hop {} ({:.1f} ms)\n",
                             kPsolaBlock, 1000.0f * static_cast<float>(kPsolaBlock) / fs_f,
                             kPsolaHop, 1000.0f * static_cast<float>(kPsolaHop) / fs_f);
    return SaveAndReport(input, output, fs, name);
}

} // namespace autotune_lab
