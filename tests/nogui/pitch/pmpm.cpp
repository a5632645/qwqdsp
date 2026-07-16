#include <qwqdsp/pitch/pmpm.hpp>
#include <qwqdsp/pitch/pmpm/pmpm_core.hpp>

#include <cmath>
#include <cstdio>
#include <numbers>
#include <vector>

// ------------------------------------------------------------
// 辅助: 生成正弦波
// ------------------------------------------------------------
static void generateSine(std::vector<float>& buffer, float freq, float sample_rate, int num_samples,
                         float amplitude = 0.5f) noexcept {
    for (int i = 0; i < num_samples; ++i) {
        buffer[i] = amplitude * std::sin(2.0f * std::numbers::pi_v<float> * freq * i / sample_rate);
    }
}

// ------------------------------------------------------------
// 辅助: 检查候选列表中是否包含目标频率
// ------------------------------------------------------------
static bool hasFrequency(std::span<const qwqdsp_pitch::PyinCandidate> candidates, float target_hz,
                         float tolerance = 0.05f) noexcept {
    for (auto const& c : candidates) {
        if (std::abs(c.pitch_hz - target_hz) / target_hz < tolerance) {
            return true;
        }
    }
    return false;
}

// ------------------------------------------------------------
// 测试 PmpmCore: 正确频率应出现在候选列表中
// ------------------------------------------------------------
static int testPmpmCore() noexcept {
    std::printf("  [PmpmCore] ... ");

    float const fs = 44100.0f;
    int const block_size = 2048;
    float const test_freq = 440.0f;

    std::vector<float> buffer(block_size);
    generateSine(buffer, test_freq, fs, block_size);

    qwqdsp_pitch::PmpmCore core;
    core.Init(fs, block_size);
    core.SetMinPitch(test_freq * 0.5f);
    core.SetMaxPitch(test_freq * 2.0f);

    auto candidates = core.Process(buffer, 10, 0.001f);

    if (candidates.empty()) {
        std::printf("FAIL: no candidates\n");
        return 1;
    }

    bool found = hasFrequency(candidates, test_freq);
    std::printf("%zu candidates, 440 Hz %s", candidates.size(), found ? "found" : "MISSING");
    for (auto const& c : candidates) {
        std::printf("  [%.1f Hz, %.4f]", c.pitch_hz, c.probability);
    }

    if (!found) {
        std::printf("  FAIL\n");
        return 1;
    }

    // 概率排序检查
    for (size_t i = 1; i < candidates.size(); ++i) {
        if (candidates[i].probability > candidates[i - 1].probability) {
            std::printf("  FAIL: not sorted by probability\n");
            return 1;
        }
    }

    std::printf("  OK\n");
    return 0;
}

// ------------------------------------------------------------
// 测试多种频率
// ------------------------------------------------------------
static int testMultipleFrequencies() noexcept {
    std::printf("  [PmpmCore Multi-Freq] ... ");

    float const fs = 44100.0f;
    int const block_size = 2048;
    struct {
        float freq;
        float min_tolerance;
    } tests[] = {
        {110.0f, 0.10f},
        {220.0f, 0.08f},
        {440.0f, 0.05f},
        {880.0f, 0.05f}
    };

    qwqdsp_pitch::PmpmCore core;
    core.Init(fs, block_size);

    std::vector<float> buffer(block_size);

    for (auto const& t : tests) {
        core.SetMinPitch(t.freq * 0.5f);
        core.SetMaxPitch(t.freq * 2.0f);

        generateSine(buffer, t.freq, fs, block_size);
        auto candidates = core.Process(buffer, 5, 0.001f);

        if (candidates.empty()) {
            std::printf("FAIL (%.0f Hz): no candidates\n", t.freq);
            return 1;
        }

        float best = candidates[0].pitch_hz;
        float err = std::abs(best - t.freq) / t.freq;

        if (err > t.min_tolerance) {
            std::printf("FAIL (%.0f Hz): got %.1f Hz (err=%.2f%%)\n", t.freq, best, err * 100.0f);
            return 1;
        }
    }

    std::printf("all OK\n");
    return 0;
}

// ------------------------------------------------------------
// 测试 PmpmCore + MonoPitchHmm 完整流程
// ------------------------------------------------------------
static int testPmpmFull() noexcept {
    std::printf("  [Pmpm Full Pipeline] ... ");

    float const fs = 44100.0f;
    int const block_size = 2048;
    int const step_size = 512;
    float const test_freq = 261.63f; // 避开 440 的整数倍子谐波
    int const num_frames = 30;

    qwqdsp_pitch::Pmpm pmpm;
    pmpm.Init(fs, block_size, step_size);
    pmpm.SetMinPitch(200.0f);
    pmpm.SetMaxPitch(500.0f);

    std::vector<float> buffer(block_size);
    for (int f = 0; f < num_frames; ++f) {
        int offset = f * step_size;
        for (int i = 0; i < block_size; ++i) {
            buffer[i] = 0.5f * std::sin(2.0f * std::numbers::pi_v<float> * test_freq * (offset + i) / fs);
        }
        pmpm.Process(buffer);
    }

    auto pitch_track = pmpm.GetPitchTrack();

    if (pitch_track.empty()) {
        std::printf("FAIL: empty pitch track\n");
        return 1;
    }

    int voiced_count = 0;
    float pitch_sum = 0.0f;
    for (auto p : pitch_track) {
        if (p > 0.0f) {
            voiced_count++;
            pitch_sum += p;
        }
    }

    if (voiced_count == 0) {
        std::printf("FAIL: no voiced frames\n");
        return 1;
    }

    float avg_pitch = pitch_sum / voiced_count;
    float pitch_error = std::abs(avg_pitch - test_freq) / test_freq;

    std::printf("voiced=%d/%d, avg=%.1f Hz (err=%.2f%%)", voiced_count, num_frames, avg_pitch, pitch_error * 100.0f);

    if (pitch_error > 0.15f) {
        std::printf("  FAIL\n");
        return 1;
    }

    std::printf("  OK\n");
    return 0;
}

// ------------------------------------------------------------
int main() {
    std::printf("pMPM 测试\n");
    std::printf("============================\n\n");

    int failed = 0;

    std::printf("--- PmpmCore ---\n");
    failed += testPmpmCore();
    failed += testMultipleFrequencies();
    std::printf("\n");

    std::printf("--- 完整流程 ---\n");
    failed += testPmpmFull();
    std::printf("\n");

    std::printf("============================\n");
    if (failed == 0) {
        std::printf("全部通过!\n");
    }
    else {
        std::printf("失败 %d 个测试\n", failed);
    }

    return failed;
}
