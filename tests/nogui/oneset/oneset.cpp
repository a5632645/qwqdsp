#include <AudioFile.h>
#include <audio_ops.hpp>
#include <cmath>
#include <format>
#include <iostream>
#include <span>
#include <vector>
#include <work_dir.hpp>

#include "baseline.hpp"
#include "complex_domain.hpp"
#include "fdk_attack.hpp"
#include "hpss.hpp"
#include "oneset_common.hpp"
#include "opus_transient.hpp"
#include "spectral_flux.hpp"
#include "superflux.hpp"

using qwqdsp_test::BaselineTransientDetector;
using qwqdsp_test::ComplexDomainDetector;
using qwqdsp_test::FdkAttackDetector;
using qwqdsp_test::HpssPercussiveOnsetDetector;
using qwqdsp_test::OnsetResult;
using qwqdsp_test::OpusTransientDetector;
using qwqdsp_test::SpectralFluxDetector;
using qwqdsp_test::SuperFluxDetector;

// ------------------------------------------------------------
// 标记叠加
// ------------------------------------------------------------
/**
 * @brief 在标记声道（初始静音）的 onset 采样位置生成短促衰减点击。
 *
 * 点击为 5ms 指数衰减正弦（约 2kHz），幅度 0.9，避免爆音。
 * 输出为立体声：第一声道原始音频，第二声道标记点击。
 */
static void AddOnsetMarkers(std::vector<float>& x, float sample_rate, std::span<const size_t> onset_samples) {
    constexpr double kFreq = 2000.0;
    constexpr double kDecay = 3000.0; // 指数衰减率 (1/s)
    const size_t click_len = static_cast<size_t>(std::lround(0.005 * sample_rate));
    const double omega = 2.0 * 3.14159265358979323846 * kFreq / static_cast<double>(sample_rate);

    for (size_t pos : onset_samples) {
        if (pos >= x.size())
            continue;
        for (size_t i = 0; i < click_len && pos + i < x.size(); ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(sample_rate);
            const float click = static_cast<float>(0.9 * std::exp(-kDecay * t) * std::sin(omega * static_cast<double>(i)));
            x[pos + i] += click;
        }
    }
}

// ------------------------------------------------------------
// 算法运行
// ------------------------------------------------------------
/**
 * @brief 运行一个检测器，叠加标记并保存立体声 wav。
 *
 * 第一声道为原始音频，第二声道为 onset 位置处的标记点击（纯脉冲），
 * 便于试听对比各算法的检测位置。
 *
 * @param algo_name  算法名（用于输出文件名与控制台）
 * @param input_name 输入名（drumloop / wormhole）
 * @param result     检测结果
 * @param fs         采样率
 * @param input      原始输入（第一声道，不污染）
 */
static void RunOnsetDetector(const char* algo_name, const char* input_name, const OnsetResult& result, float fs,
                             const std::vector<float>& input) {
    auto ch0 = input;
    std::vector<float> ch1(ch0.size(), 0.0f); // 标记声道（初始静音）
    AddOnsetMarkers(ch1, fs, result.onset_samples);

    AudioFile<float> out;
    out.setNumChannels(2);
    out.setNumSamplesPerChannel(ch0.size());
    out.samples = {std::move(ch0), std::move(ch1)};
    out.setSampleRate(static_cast<int>(fs));
    out.save(qwqdsp_support::OutputFile(std::format("oneset_{}_{}.wav", algo_name, input_name)));
    std::cout << std::format("  {} {}: {} onsets\n", algo_name, input_name, result.onset_samples.size()) << std::flush;
}

int main() {
    // --------------------------------------------------------
    // 加载输入
    // --------------------------------------------------------
    const std::vector<std::pair<const char*, std::string>> inputs = {
        {"drumloop", qwqdsp_support::InputFile("drumloop.wav")},
        {"wormhole", qwqdsp_support::InputFile("wormhole.wav")},
    };

    for (const auto& [input_name, path] : inputs) {
        AudioFile<float> file{path};
        const auto x_vec = file.samples.front(); // 每算法共享的原始输入副本
        const float fs = static_cast<float>(file.getSampleRate());
        std::cout << std::format("=== {}: {} samples, {} Hz ===\n", input_name, x_vec.size(), file.getSampleRate())
                  << std::flush;

        // 1. 基线（当前 transient_vocoder 逻辑）
        {
            BaselineTransientDetector det;
            det.SetFrameSize(2048);
            det.SetHopSize(512);
            auto r = det.Detect(x_vec, fs);
            RunOnsetDetector("baseline", input_name, r, fs, x_vec);
        }

        // 2. 谱通量（librosa）
        {
            SpectralFluxDetector det;
            det.SetFrameSize(2048);
            det.SetHopSize(512);
            auto r = det.Detect(x_vec, fs);
            RunOnsetDetector("flux", input_name, r, fs, x_vec);
        }

        // 3. 复合域（Dixon 2006）
        {
            ComplexDomainDetector det;
            det.SetFrameSize(2048);
            det.SetHopSize(512);
            auto r = det.Detect(x_vec, fs);
            RunOnsetDetector("complex", input_name, r, fs, x_vec);
        }

        // 4. SuperFlux（Boeck 2013）
        {
            SuperFluxDetector det;
            det.SetFrameSize(2048);
            det.SetHopSize(512);
            auto r = det.Detect(x_vec, fs);
            RunOnsetDetector("superflux", input_name, r, fs, x_vec);
        }

        // 5. Opus CELT 瞬态（工业级，时域）
        {
            OpusTransientDetector det;
            det.SetFrameSize(960); // 48kHz 20ms 帧
            auto r = det.Detect(x_vec, fs);
            RunOnsetDetector("opus", input_name, r, fs, x_vec);
        }

        // 6. FDK-AAC 块切换攻击（工业级，时域）
        {
            FdkAttackDetector det;
            det.SetGranuleSize(1024);
            auto r = det.Detect(x_vec, fs);
            RunOnsetDetector("fdk", input_name, r, fs, x_vec);
        }

        // 7. HPSS 打击乐分离（librosa）
        {
            HpssPercussiveOnsetDetector det;
            det.SetFrameSize(2048);
            det.SetHopSize(512);
            auto r = det.Detect(x_vec, fs);
            RunOnsetDetector("hpss", input_name, r, fs, x_vec);
        }

        std::cout << "\n" << std::flush;
    }

    std::cout << std::format("all done\n") << std::flush;
    return 0;
}
