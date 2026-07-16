#include "AudioFile.h"
#include <qwqdsp/pitch/pmpm.hpp>
#include <qwqdsp/pitch/pyin.hpp>
#include <work_dir.hpp>

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

// ------------------------------------------------------------
// Test info
// ------------------------------------------------------------
struct DegradedTest {
    const char* filename;
    const char* label;
    float expected_hz; // E3 ~ 164.81
};

static constexpr float kE3 = 164.81f;

static const DegradedTest kTests[] = {
    {"Viola-deg0.wav", "deg0(clean)",                 kE3},
    {"Viola-deg1.wav", "deg1(mild_noise)",            kE3},
    {"Viola-deg2.wav", "deg2(noise+distortion)",      kE3},
    {"Viola-deg3.wav", "deg3(noise+dist+background)", kE3},
    {"Viola-deg4.wav", "deg4(mp3_roundtrip)",         kE3},
    {"Viola-deg5.wav", "deg5(horrible)",              kE3},
};

// ------------------------------------------------------------
// Process a single file with a pitch estimator
// ------------------------------------------------------------
template <typename Estimator>
static float processFile(Estimator& est, const std::vector<float>& audio, int sample_rate, int block_size,
                         int hop_size) {
    est.Init(static_cast<float>(sample_rate), block_size, hop_size);
    est.SetMinPitch(kE3 * 0.5f);
    est.SetMaxPitch(kE3 * 2.0f);

    std::vector<float> block(block_size, 0.0f);
    for (size_t pos = 0; pos + static_cast<size_t>(block_size) <= audio.size(); pos += static_cast<size_t>(hop_size)) {
        std::copy_n(audio.begin() + static_cast<long>(pos), block_size, block.begin());
        est.Process(block);
    }

    auto track = est.GetPitchTrack();
    if (track.empty())
        return 0.0f;

    float sum = 0.0f;
    int count = 0;
    for (auto p : track) {
        if (p > 0.0f) {
            sum += p;
            ++count;
        }
    }
    return (count > 0) ? (sum / static_cast<float>(count)) : 0.0f;
}

// ------------------------------------------------------------
// main
// ------------------------------------------------------------
int main() {
    std::string dir = (qwqdsp_support::GetInputDir() / "pitch").string();

    std::cout << std::format("pMPM vs pYIN on degraded viola (E3={:.1f} Hz)\n", kE3);
    std::cout << "================================================================\n";
    std::cout << std::format("  {:26s}  {:>10s}  {:>10s}  {:>10s}  {:>10s}\n", "file", "pmpm(Hz)", "pyin(Hz)",
                             "pmpm_err%", "pyin_err%");
    std::cout << "  -----------------------------------------------------------------\n";

    int failed_pmpm = 0;
    int failed_pyin = 0;

    for (auto const& t : kTests) {
        std::string path = dir + "/" + t.filename;

        AudioFile<float> af;
        if (!af.load(path)) {
            std::cout << std::format("  {:26s}  {:>10s}  {:>10s}  {:>10s}  {:>10s}\n", t.label, "ERR", "ERR", "-", "-");
            ++failed_pmpm;
            ++failed_pyin;
            continue;
        }

        int sample_rate = static_cast<int>(af.getSampleRate());
        int num_channels = af.getNumChannels();
        int num_samples = af.getNumSamplesPerChannel();

        // Mix down to mono if needed
        std::vector<float> audio(num_samples, 0.0f);
        if (num_channels == 1) {
            for (int i = 0; i < num_samples; ++i)
                audio[i] = af.samples[0][i];
        }
        else {
            float inv_ch = 1.0f / static_cast<float>(num_channels);
            for (int i = 0; i < num_samples; ++i) {
                float sum = 0.0f;
                for (int ch = 0; ch < num_channels; ++ch)
                    sum += af.samples[ch][i];
                audio[i] = sum * inv_ch;
            }
        }

        int block_size = 2048;
        int hop_size = 512;

        qwqdsp_pitch::Pmpm pmpm;
        float pmpm_pitch = processFile(pmpm, audio, sample_rate, block_size, hop_size);

        qwqdsp_pitch::Pyin pyin;
        pyin.SetThresholdDistribution(2);
        float pyin_pitch = processFile(pyin, audio, sample_rate, block_size, hop_size);

        float pmpm_err = (pmpm_pitch > 0.0f) ? std::abs(pmpm_pitch - kE3) / kE3 * 100.0f : -1.0f;
        float pyin_err = (pyin_pitch > 0.0f) ? std::abs(pyin_pitch - kE3) / kE3 * 100.0f : -1.0f;

        constexpr float kErrTol = 10.0f;
        if (pmpm_err > kErrTol || pmpm_pitch <= 0.0f)
            ++failed_pmpm;
        if (pyin_err > kErrTol || pyin_pitch <= 0.0f)
            ++failed_pyin;

        std::cout << std::format("  {:26s}  {:8.1f}   {:8.1f}   {:8.2f}   {:8.2f}", t.label, pmpm_pitch, pyin_pitch,
                                 pmpm_err, pyin_err);

        if (pmpm_err > kErrTol || pmpm_pitch <= 0.0f)
            std::cout << "  !pmpm";
        if (pyin_err > kErrTol || pyin_pitch <= 0.0f)
            std::cout << "  !pyin";
        std::cout << "\n";
    }

    std::cout << "  ================================================================\n";
    std::cout << std::format("  pMPM fails: {}, pYIN fails: {}\n", failed_pmpm, failed_pyin);
    return failed_pmpm + failed_pyin;
}
