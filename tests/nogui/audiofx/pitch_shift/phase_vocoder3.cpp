#include <AudioFile.h>
#include <audio_ops.hpp>
#include <format>
#include <iostream>
#include <qwqdsp/fx/phase_vocoder2.hpp>
#include <qwqdsp/segement/analyze_synthsis_offline2.hpp>
#include <qwqdsp/window/helper.hpp>
#include <work_dir.hpp>

static std::vector<float> LinearResample(std::span<const float> in, float ratio) {
    const size_t out_len = static_cast<size_t>(std::round(static_cast<float>(in.size()) / ratio));
    if (out_len < 2)
        return {};
    std::vector<float> out(out_len);
    const float step = static_cast<float>(in.size() - 1) / static_cast<float>(out_len - 1);
    for (size_t i = 0; i < out_len; ++i) {
        const float pos = step * static_cast<float>(i);
        const size_t idx = static_cast<size_t>(pos);
        const float frac = pos - static_cast<float>(idx);
        const size_t nxt = std::min(idx + 1, in.size() - 1);
        out[i] = in[idx] + frac * (in[nxt] - in[idx]);
    }
    return out;
}

static void Process(const char* name, float kt, float kp) {
    auto wav_path = qwqdsp_support::SweepWav();
    AudioFile<float> file{wav_path};
    auto& x_vec = file.samples.front();

    std::cout << std::format("{} (kt={:.2f}, kp={:.2f}): {} -> ", name, kt, kp, x_vec.size()) << std::flush;

    qwqdsp_fx::PhaseGradientVocoder dsp;
    dsp.Init(4096, 8192, 1024);
    dsp.Reset();
    int hop_analyze = dsp.SetTimeStretch(kt * kp);

    qwqdsp_segement::AnalyzeSynthsisOffline2 as;
    as.SetInputHop(hop_analyze);
    as.SetOutputHop(1024);
    as.SetInputSize(4096);
    as.SetOutputSize(4096);

    std::vector<float> out;
    as.Process(x_vec, out,
               [&dsp](std::span<const float> input, std::span<float> output) { dsp.Process(input, output); });

    out = LinearResample(out, kp);

    std::cout << std::format("{}\n", out.size()) << std::flush;

    qwqdsp_support::AudioOps::Normalize(out);
    file.setNumSamplesPerChannel(out.size());
    file.samples[0] = out;
    file.setNumChannels(1);
    file.save(qwqdsp_support::OutputFile(name));
    std::cout << std::format("saved {}\n\n", name) << std::flush;
}

// ------------------------------------------------------------
// Delta 测试
// ------------------------------------------------------------
/**
 * @brief 使用 delta 8192 处为 1.0）作为输入测试声码器。
 */
static void TestDelta(const char* name, float kt, float kp) {
    constexpr size_t kLen = 65536;
    constexpr size_t kDeltaPos = 8192;
    std::vector<float> x(kLen, 0.0f);
    x[kDeltaPos] = 1.0f;

    std::cout << std::format("{} (kt={:.2f}, kp={:.2f}): delta@{} -> ", name, kt, kp, kDeltaPos) << std::flush;

    qwqdsp_fx::PhaseGradientVocoder dsp;
    dsp.Init(4096, 8192, 1024);
    dsp.Reset();
    int hop_analyze = dsp.SetTimeStretch(kt * kp);

    qwqdsp_segement::AnalyzeSynthsisOffline2 as;
    as.SetInputHop(hop_analyze);
    as.SetOutputHop(1024);
    as.SetInputSize(4096);
    as.SetOutputSize(4096);

    std::vector<float> out;
    as.Process(x, out, [&dsp](std::span<const float> input, std::span<float> output) { dsp.Process(input, output); });
    std::cout << std::format("{}\n", out.size()) << std::flush;

    out = LinearResample(out, kp);

    qwqdsp_support::AudioOps::Normalize(out);
    AudioFile<float> file;
    file.setNumSamplesPerChannel(out.size());
    file.samples.resize(1);
    file.samples[0] = out;
    file.setNumChannels(1);
    file.save(qwqdsp_support::OutputFile(name));
    std::cout << std::format("saved {}\n\n", name) << std::flush;
}

int main() {
    Process("PV3_ts_1.5x.wav", 1.5f, 1.0f);
    Process("PV3_ps_1.5x.wav", 1.0f, 1.5f);
    Process("PV3_ts1.5_ps1.5.wav", 1.5f, 1.5f);

    TestDelta("PV3_delta_identity.wav", 1.0f, 1.0f);
    TestDelta("PV3_delta_ts1.5.wav", 1.5f, 1.0f);
    TestDelta("PV3_delta_ps1.5.wav", 1.0f, 1.5f);

    std::cout << std::format("all done\n") << std::flush;
}
