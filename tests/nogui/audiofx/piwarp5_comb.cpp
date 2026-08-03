// ------------------------------------------------------------
//  piwarp5_comb.cpp
//  Piwarp5 (Kaiser 窗 + 帧内时间反转) 参数扫描 — 观察 comb filter
//
//  参数矩阵 (len = M × len_mult, M=32):
//    len*   beta
//     4     32/16/8
//     2     16/8/4
//     1      8/4/2
//
//  主流程复制自 piwarp.cpp 的 Piwarp5 (hop = M-1 = 31),
//  输入可选: sweep.wav (kUseSweep=true) 或 WhiteNoise 纯白噪声
// ------------------------------------------------------------

#include <algorithm>
#include <AudioFile.h>
#include <format>
#include <iostream>
#include <qwqdsp/oscillator/noise.hpp>
#include <qwqdsp/segement/analyze_synthsis_offline.hpp>
#include <qwqdsp/window/kaiser.hpp>
#include <string>
#include <utility>
#include <vector>
#include <work_dir.hpp>

// ------------------------------------------------------------
//  Piwarp5Class (复制自 piwarp.cpp)
// ------------------------------------------------------------
struct Piwarp5Class {
    Piwarp5Class(int dft_size, float kaiser_beta) {
        window_.resize(dft_size);
        // 仅 Kaiser 窗 (窗函数本身, 非低通原型滤波器系数)
        qwqdsp_window::Kaiser::Window(window_, kaiser_beta, false);
    }

    void operator()(std::span<const float> in, std::span<float> out) noexcept {
        int len = in.size();
        for (int i = 0; i < in.size(); ++i) {
            out[i] = 2 * window_[i] * in[len - i - 1];
        }
    }

    std::vector<float> window_;
};

// ------------------------------------------------------------
//  参数矩阵
// ------------------------------------------------------------
struct CombCase {
    int len_mult;   // len = M × len_mult
    float beta;     // Kaiser 窗 β
};

static constexpr CombCase kCases[] = {
    // len*4
    {4, 32.0f}, {4, 16.0f}, {4, 8.0f},
    // len*2
    {2, 16.0f}, {2, 8.0f}, {2, 4.0f},
    // len*1
    {1, 8.0f}, {1, 4.0f}, {1, 2.0f},
};

// 输入选择: true = sweep.wav, false = WhiteNoise 纯白噪声
static constexpr bool kUseSweep = true;

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------
int main() {
    constexpr int kSampleRate = 48000;

    // 加载/生成输入
    std::vector<float> x_vec;
    if constexpr (kUseSweep) {
        AudioFile<float> in_file{qwqdsp_support::SweepWav()};
        x_vec.assign(in_file.samples.front().begin(), in_file.samples.front().end());
        std::cout << std::format("input: sweep.wav ({}) samples\n", x_vec.size());
    }
    else {
        // 纯白噪声 (noise.wav 非纯噪声, 会污染自相关/comb 分析)
        constexpr int kNumSamples = kSampleRate * 2;  // 2 秒
        x_vec.resize(static_cast<size_t>(kNumSamples));
        qwqdsp_oscillator::WhiteNoise noise;
        noise.SetSeed(0x12345678u);
        for (float& x : x_vec) {
            x = noise.Next();
        }
        std::cout << std::format("input: WhiteNoise ({}) samples\n", x_vec.size());
    }

    constexpr int filter_banks = 32;        // M — 子带数
    constexpr int hop = filter_banks - 1;   // R = M-1 (复制自 Piwarp5)

    for (const CombCase& c : kCases) {
        const int len = filter_banks * c.len_mult;

        qwqdsp_segement::AnalyzeSynthsisOffline as;
        as.SetSize(len);
        as.SetInputHop(hop);
        as.SetOutputHop(hop);
        as.Reset();

        std::vector<float> y_vec;
        Piwarp5Class piwarp5{len, c.beta};
        as.Process(x_vec, y_vec, piwarp5);

        // 归一化输出到 0.9 峰值 — 帧反转 + 重叠相加可能放大信号, 避免 clip
        {
            float peak = 0.0f;
            for (const float v : y_vec) {
                peak = std::max(peak, std::abs(v));
            }
            if (peak > 1e-9f) {
                const float g = 0.9f / peak;
                for (float& v : y_vec) {
                    v *= g;
                }
            }
        }

        AudioFile<float> out_file;
        out_file.setNumChannels(1);
        out_file.setBitDepth(32);
        out_file.setSampleRate(kSampleRate);
        out_file.setNumSamplesPerChannel(y_vec.size());
        std::copy(y_vec.begin(), y_vec.end(), out_file.samples[0].begin());
        const std::string prefix = kUseSweep ? "piwarp5_comb_sweep_" : "piwarp5_comb_";
        const std::string name =
            std::format("{}len{}x_b{}.wav", prefix, c.len_mult, static_cast<int>(c.beta));
        out_file.save(qwqdsp_support::OutputFile(name));
        std::cout << std::format("len{}x beta{:2} -> {}\n", c.len_mult, static_cast<int>(c.beta), name);
    }
    return 0;
}
