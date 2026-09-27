#include <AudioFile.h>
#include <audio_ops.hpp>
#include <cmath>
#include <random>
#include <vector>
#include <work_dir.hpp>

// ------------------------------------------------------------
// velvet1
// ------------------------------------------------------------

/**
 * @brief 生成 velvet noise 脉冲响应，用于混响
 * @param ir           输出 IR 数组
 * @param density      每秒脉冲数（典型值 500～2000）
 * @param rt60_samples RT60 对应的采样点数
 * @param sample_rate  采样率
 * @param rng          随机数引擎
 */
static void GenerateVelvetIR(std::vector<float>& ir, float density, int rt60_samples, float sample_rate,
                             std::default_random_engine& rng) {
    ir.assign(rt60_samples, 0.0f);

    int num_pulses = static_cast<int>(density * rt60_samples / sample_rate);
    std::uniform_int_distribution<int> pos_dist(0, rt60_samples - 1);
    std::uniform_int_distribution<int> sign_dist(0, 1);

    // 指数衰减包络：每采样点的衰减因子，RT60 后 -60 dB
    float decay_per_sample = std::pow(10.0f, -3.0f / rt60_samples);

    for (int i = 0; i < num_pulses; ++i) {
        int pos = pos_dist(rng);
        float sign = sign_dist(rng) == 0 ? 1.0f : -1.0f;
        float envelope = std::pow(decay_per_sample, static_cast<float>(pos));
        ir[pos] += sign * envelope;
    }

    // 归一化
    float peak = 0.0f;
    for (auto x : ir) {
        peak = std::max(peak, std::abs(x));
    }
    if (peak > 0.0f) {
        for (auto& x : ir) {
            x /= peak;
        }
    }
}

/**
 * @brief 稀疏卷积：用 velvet noise IR 对输入做直接卷积
 */
static void SparseConvolve(std::span<const float> input, std::span<const float> ir, std::vector<float>& output) {
    int out_len = static_cast<int>(input.size() + ir.size() - 1);
    output.assign(out_len, 0.0f);

    // 收集非零脉冲位置
    struct Impulse {
        int pos;
        float amp;
    };
    std::vector<Impulse> impulses;
    impulses.reserve(2048);
    for (int i = 0; i < static_cast<int>(ir.size()); ++i) {
        if (ir[i] != 0.0f) {
            impulses.push_back({i, ir[i]});
        }
    }

    // 稀疏直接卷积
    for (auto [pos, amp] : impulses) {
        for (int i = 0; i < static_cast<int>(input.size()); ++i) {
            output[i + pos] += input[i] * amp;
        }
    }
}

int main() {
    constexpr float sample_rate = 48000.0f;
    constexpr float rt60_ms = 2000.0f; // 混响时长 2 秒
    constexpr float density = 1000.0f; // 脉冲密度 1000/sec

    AudioFile<float> file{qwqdsp_support::InputFile("drumloop.wav")};
    auto x_vec = file.samples.front();

    // 生成 velvet noise IR
    int rt60_samples = static_cast<int>(rt60_ms / 1000.0f * sample_rate);
    std::vector<float> ir;
    std::default_random_engine rng{42};
    GenerateVelvetIR(ir, density, rt60_samples, sample_rate, rng);

    // 稀疏卷积
    std::vector<float> y_vec;
    SparseConvolve(x_vec, ir, y_vec);

    qwqdsp_support::AudioOps::Normalize(y_vec);
    file.samples.front() = y_vec;
    file.setNumSamplesPerChannel(static_cast<int>(y_vec.size()));
    file.setNumChannels(1);
    file.save(qwqdsp_support::OutputFile("velvet_reverb.wav"));
}
