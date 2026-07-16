#include <qwqdsp/pitch/pyin.hpp>
#include <qwqdsp/pitch/pyin/pyin_core.hpp>
#include <qwqdsp/pitch/pyin/pyin_hmm.hpp>

#include <cmath>
#include <format>
#include <iostream>
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
// 测试 PyinCore: 正确频率应出现在候选列表中
// ------------------------------------------------------------
static int testPyinCore() noexcept {
    std::cout << "  [PyinCore] ... ";

    float const fs = 44100.0f;
    int const block_size = 2048;
    float const test_freq = 440.0f;

    std::vector<float> buffer(block_size);
    generateSine(buffer, test_freq, fs, block_size);

    qwqdsp_pitch::PyinCore core;
    core.Init(fs, block_size);
    core.SetThresholdDistribution(2);

    auto candidates = core.Process(buffer, 10, 0.001f);

    if (candidates.empty()) {
        std::cout << "FAIL: no candidates\n";
        return 1;
    }

    bool found = hasFrequency(candidates, test_freq);
    std::cout << std::format("{} candidates, 440 Hz {}", candidates.size(), found ? "found" : "MISSING");
    for (auto const& c : candidates) {
        std::cout << std::format("  [{:.1f} Hz, {:.4f}]", c.pitch_hz, c.probability);
    }

    if (!found) {
        std::cout << "  FAIL\n";
        return 1;
    }

    // 概率排序检查
    for (size_t i = 1; i < candidates.size(); ++i) {
        if (candidates[i].probability > candidates[i - 1].probability) {
            std::cout << "  FAIL: not sorted by probability\n";
            return 1;
        }
    }

    std::cout << "  OK\n";
    return 0;
}

// ------------------------------------------------------------
// 测试多种频率
// ------------------------------------------------------------
static int testMultipleFrequencies() noexcept {
    std::cout << "  [PyinCore Multi-Freq] ... ";

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

    qwqdsp_pitch::PyinCore core;
    core.Init(fs, block_size);
    core.SetThresholdDistribution(2);

    std::vector<float> buffer(block_size);

    for (auto const& t : tests) {
        // 调整 min/max pitch 范围使候选聚焦
        core.SetMinPitch(t.freq * 0.5f);
        core.SetMaxPitch(t.freq * 2.0f);

        generateSine(buffer, t.freq, fs, block_size);
        auto candidates = core.Process(buffer, 5, 0.001f);

        if (candidates.empty()) {
            std::cout << std::format("FAIL ({:.0f} Hz): no candidates\n", t.freq);
            return 1;
        }

        float best = candidates[0].pitch_hz;
        float err = std::abs(best - t.freq) / t.freq;

        if (err > t.min_tolerance) {
            std::cout << std::format("FAIL ({:.0f} Hz): got {:.1f} Hz (err={:.2f}%)\n", t.freq, best, err * 100.0f);
            return 1;
        }
    }

    std::cout << "all OK\n";
    return 0;
}

// ------------------------------------------------------------
// 测试 PyinCore + MonoPitchHmm 完整流程
// ------------------------------------------------------------
static int testPyinFull() noexcept {
    std::cout << "  [Pyin Full Pipeline] ... ";

    float const fs = 44100.0f;
    int const block_size = 2048;
    int const step_size = 512;
    float const test_freq = 261.63f; // 避开 440 的整数倍子谐波
    int const num_frames = 30;

    qwqdsp_pitch::Pyin pyin;
    pyin.Init(fs, block_size, step_size);
    pyin.SetThresholdDistribution(2);
    pyin.SetMinPitch(200.0f);
    pyin.SetMaxPitch(500.0f);

    std::vector<float> buffer(block_size);
    for (int f = 0; f < num_frames; ++f) {
        int offset = f * step_size;
        for (int i = 0; i < block_size; ++i) {
            buffer[i] = 0.5f * std::sin(2.0f * std::numbers::pi_v<float> * test_freq * (offset + i) / fs);
        }
        pyin.Process(buffer);
    }

    auto pitch_track = pyin.GetPitchTrack();

    if (pitch_track.empty()) {
        std::cout << "FAIL: empty pitch track\n";
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
        std::cout << "FAIL: no voiced frames\n";
        return 1;
    }

    float avg_pitch = pitch_sum / voiced_count;
    float pitch_error = std::abs(avg_pitch - test_freq) / test_freq;

    std::cout << std::format("voiced={}/{}, avg={:.1f} Hz (err={:.2f}%)", voiced_count, num_frames, avg_pitch,
                             pitch_error * 100.0f);

    if (pitch_error > 0.15f) {
        std::cout << "  FAIL\n";
        return 1;
    }

    std::cout << "  OK\n";
    return 0;
}

// ------------------------------------------------------------
// 测试 SparseHmm 维特比解码
// ------------------------------------------------------------
static int testViterbi() noexcept {
    std::cout << "  [MonoPitchHmm Viterbi] ... ";

    qwqdsp_pitch::MonoPitchHmm hmm;
    hmm.Init();

    // 构造模拟观测序列: 440 Hz 持续 5 帧
    size_t const n_frame = 5;
    std::vector<std::vector<double>> obs_prob(n_frame);

    for (size_t f = 0; f < n_frame; ++f) {
        std::vector<std::pair<double, double>> pitch_prob;
        pitch_prob.emplace_back(69.0, 0.8); // 69 MIDI = 440 Hz
        obs_prob[f] = hmm.CalculateObsProb(pitch_prob);
    }

    auto path = hmm.DecodeViterbi(obs_prob);

    if (path.size() != n_frame) {
        std::cout << "FAIL: wrong path length\n";
        return 1;
    }

    // 检查路径是否包含有声状态 (正频率)
    bool has_voiced = false;
    for (auto s : path) {
        if (hmm.GetFrequency(s) > 0) {
            has_voiced = true;
            break;
        }
    }

    if (!has_voiced) {
        std::cout << "FAIL: no voiced states in path\n";
        return 1;
    }

    std::cout << "OK\n";
    return 0;
}

// ------------------------------------------------------------
int main() {
    std::cout << "pYIN test\n";
    std::cout << "============================\n\n";

    int failed = 0;

    std::cout << "--- PyinCore ---\n";
    failed += testPyinCore();
    failed += testMultipleFrequencies();
    std::cout << "\n";

    std::cout << "--- HMM ---\n";
    failed += testViterbi();
    std::cout << "\n";

    std::cout << "--- Full Pipeline ---\n";
    failed += testPyinFull();
    std::cout << "\n";

    std::cout << "============================\n";
    if (failed == 0) {
        std::cout << "All passed!\n";
    }
    else {
        std::cout << std::format("Failed {} test(s)\n", failed);
    }

    return failed;
}
