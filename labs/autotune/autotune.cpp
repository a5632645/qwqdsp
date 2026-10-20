/// @brief autotune lab 集中运行入口，对比固定窗与自适应窗两种算法

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <format>
#include <iostream>
#include <numbers>
#include <span>
#include <string>
#include <vector>

#include "autotune_adaptive.hpp"
#include "autotune_common.hpp"
#include "autotune_pghi.hpp"
#include "autotune_psola.hpp"

namespace {

using namespace autotune_lab;

constexpr float kWormholeSeconds = 3.8f;

std::vector<float> MakeChirp(int fs) {
    // C3 -> C5 的指数扫频正弦, 用来听跨音阶的量化台阶。
    constexpr float kStartHz = 130.8128f;
    constexpr float kEndHz = 523.2511f;
    constexpr float kDuration = 5.0f;
    size_t const num_samples = static_cast<size_t>(static_cast<float>(fs) * kDuration);
    std::vector<float> x(num_samples, 0.0f);
    float phase = 0.0f;
    float const octaves = std::log2(kEndHz / kStartHz);
    for (size_t i = 0; i < num_samples; ++i) {
        float const t = static_cast<float>(i) / static_cast<float>(fs);
        float const u = std::min(t / kDuration, 1.0f);
        float const f = kStartHz * std::exp2(octaves * u);
        phase += 2.0f * std::numbers::pi_v<float> * f / static_cast<float>(fs);
        x[i] = 0.5f * std::sin(phase);
    }

    size_t const fade = static_cast<size_t>(fs) * 20 / 1000;
    for (size_t i = 0; i < fade && i < num_samples; ++i) {
        float const w = static_cast<float>(i) / static_cast<float>(fade);
        x[i] *= w;
        x[num_samples - 1 - i] *= w;
    }
    return x;
}

bool RunCase(std::string const& label, std::span<const float> input, int fs, std::string const& base) {
    std::cout << std::format("== {} ==\n", label);
    bool ok = true;
    ok &= RunPsola(input, fs, base + "_psola.wav");
    ok &= RunAdaptive(input, fs, base + "_adaptive.wav");
    ok &= RunPghi(input, fs, base + "_pghi.wav");
    std::cout << "\n";
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    {
        std::vector<float> input;
        int fs = 0;
        std::string const path = qwqdsp_support::InputFile("wormhole.wav");
        if (LoadMono(path, input, fs)) {
            size_t const max_samples = static_cast<size_t>(static_cast<float>(fs) * kWormholeSeconds);
            if (input.size() > max_samples) {
                input.resize(max_samples);
            }
            std::string const label = std::format("wormhole first {:.1f} s", kWormholeSeconds);
            ok &= RunCase(label, input, fs, "autotune_wormhole");
        }
        else {
            std::cout << std::format("failed to load {}\n", path);
            ok = false;
        }
    }
    {
        auto const input = MakeChirp(44100);
        std::filesystem::create_directories(qwqdsp_support::GetOutputDir());
        if (!SaveMono(qwqdsp_support::OutputFile("autotune_chirp_dry.wav"), input, 44100)) {
            std::cout << "failed to save autotune_chirp_dry.wav\n";
            ok = false;
        }
        ok &= RunCase("chirp C3 -> C5", input, 44100, "autotune_chirp");
    }
    return ok ? 0 : 1;
}
