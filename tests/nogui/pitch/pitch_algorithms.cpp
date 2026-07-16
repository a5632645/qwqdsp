#include <qwqdsp/pitch/fast_yin.hpp>
#include <qwqdsp/pitch/mpm.hpp>

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
// 测试 FastYin
// ------------------------------------------------------------
static int testFastYin() noexcept {
    std::printf("  [FastYin] ... ");
    int failed = 0;

    float const fs = 44100.0f;
    int const block_size = 2048;
    float test_freqs[] = {110.0f, 220.0f, 440.0f, 880.0f};

    qwqdsp_pitch::FastYin yin;
    yin.Init(fs, block_size * 2);

    std::vector<float> buffer(block_size * 2);

    for (auto test_freq : test_freqs) {
        // FastYin 的 block 需要是 Init 时 size 的一半
        // 因为 delta_corr_ 是 size/2, 而 Process 里 max_tal = num_samples/2
        // 所以 block 大小应为 Init 时传入的 size
        generateSine(buffer, test_freq, fs, block_size * 2);

        // 重置并配置频率范围
        yin.Init(fs, block_size * 2);
        yin.SetMinPitch(test_freq * 0.5f);
        yin.SetMaxPitch(test_freq * 2.0f);
        yin.SetThreshold(0.15f);

        yin.Process(buffer);
        auto pitch = yin.GetPitch();

        float err = (pitch.pitch_hz > 0) ? std::abs(pitch.pitch_hz - test_freq) / test_freq : 1.0f;

        if (err > 0.05f) {
            std::printf("\n    FAIL (%.0f Hz): got %.1f Hz (err=%.2f%%)", test_freq, pitch.pitch_hz, err * 100.0f);
            failed++;
        }
        else {
            std::printf(" %.0f=%.1f(%.1f%%)", test_freq, pitch.pitch_hz, err * 100.0f);
        }
    }

    if (failed == 0)
        std::printf("  OK\n");
    return failed;
}

// ------------------------------------------------------------
// 测试 MPM
// ------------------------------------------------------------
static int testMPM() noexcept {
    std::printf("  [MPM] ... ");
    int failed = 0;

    float const fs = 44100.0f;
    int const block_size = 2048;
    float test_freqs[] = {110.0f, 220.0f, 440.0f, 880.0f};
    float tolerances[] = {0.08f, 0.06f, 0.05f, 0.05f}; // 低频允许更大误差

    std::vector<float> buffer(block_size);
    qwqdsp_pitch::MPM mpm;

    for (int i = 0; i < 4; ++i) {
        float test_freq = test_freqs[i];
        float tolerance = tolerances[i];

        mpm.Init(fs, block_size);
        mpm.SetMinPitch(test_freq * 0.5f);
        mpm.SetMaxPitch(test_freq * 2.0f);

        generateSine(buffer, test_freq, fs, block_size);

        mpm.Process(buffer);
        auto pitch = mpm.GetPitch();

        float err = (pitch.pitch_hz > 0) ? std::abs(pitch.pitch_hz - test_freq) / test_freq : 1.0f;

        if (err > tolerance) {
            std::printf("\n    FAIL (%.0f Hz): got %.1f Hz (err=%.2f%%)", test_freq, pitch.pitch_hz, err * 100.0f);
            failed++;
        }
        else {
            std::printf(" %.0f=%.1f(%.1f%%)", test_freq, pitch.pitch_hz, err * 100.0f);
        }
    }

    if (failed == 0)
        std::printf("  OK\n");
    return failed;
}

// ------------------------------------------------------------
// 测试非周期性检测 (噪声输入应返回高 non_period_ratio)
// ------------------------------------------------------------
static int testNonPeriodic() noexcept {
    std::printf("  [NonPeriodic] ... ");
    int failed = 0;

    float const fs = 44100.0f;
    int const block_size = 2048;

    std::vector<float> noise(block_size);
    for (auto& s : noise) {
        s = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 2.0f - 1.0f;
    }

    // FastYin
    {
        qwqdsp_pitch::FastYin yin;
        yin.Init(fs, block_size);
        yin.Process(noise);
        auto p = yin.GetPitch();
        if (p.non_period_ratio < 0.5f) {
            std::printf("\n    FastYin non_period=%.3f (expected >0.5)", p.non_period_ratio);
            failed++;
        }
    }

    // MPM
    {
        qwqdsp_pitch::MPM mpm;
        mpm.Init(fs, block_size);
        mpm.Process(noise);
        auto p = mpm.GetPitch();
        if (p.non_period_ratio < 0.5f) {
            std::printf("\n    MPM non_period=%.3f (expected >0.5)", p.non_period_ratio);
            failed++;
        }
    }

    if (failed == 0)
        std::printf("OK\n");
    return failed;
}

// ------------------------------------------------------------
int main() {
    std::printf("音高检测算法测试\n");
    std::printf("============================\n\n");

    int failed = 0;

    std::printf("--- FastYin ---\n");
    failed += testFastYin();
    std::printf("\n");

    std::printf("--- MPM ---\n");
    failed += testMPM();
    std::printf("\n");

    std::printf("--- 非周期性检测 ---\n");
    failed += testNonPeriodic();
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
