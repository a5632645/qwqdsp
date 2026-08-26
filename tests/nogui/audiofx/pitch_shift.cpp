#include <AudioFile.h>
#include <audio_ops.hpp>
#include <format>
#include <iostream>
#include <vector>
#include <work_dir.hpp>

#include "pitch_shift/phase_vocoder2.hpp"
#include "pitch_shift/phase_vocoder3.hpp"
#include "pitch_shift/psola.hpp"
#include "pitch_shift/wsola.hpp"

using qwqdsp_test::PhaseGradientVocoder;
using qwqdsp_test::Psola;
using qwqdsp_test::RunPhaseVocoder3;
using qwqdsp_test::RunWsola;

// ------------------------------------------------------------
// WSOLA
// ------------------------------------------------------------
/**
 * @brief WSOLA 时间拉伸，输入 wormhole.wav，输出 wsola.wav。
 */
static void RunWsolaTest() {
    AudioFile<float> infile;
    infile.load(qwqdsp_support::WormholeWav());
    auto& data = infile.samples.front();

    std::cout << std::format("wsola (stretch=2.0): {} -> ", data.size()) << std::flush;

    std::vector<float> out = RunWsola(data, 2.0f);

    std::cout << std::format("{}\n", out.size()) << std::flush;

    AudioFile<float>::AudioBuffer buf;
    buf.push_back(std::move(out));
    infile.setAudioBuffer(buf);
    infile.save(qwqdsp_support::OutputFile("wsola.wav"));
    std::cout << std::format("saved wsola.wav\n\n") << std::flush;
}

// ------------------------------------------------------------
// PSOLA
// ------------------------------------------------------------
/**
 * @brief PSOLA 音高/共振峰移动，输入 wormhole.wav，输出 4 组 wav。
 */
static void RunPsolaTest() {
    const auto wav_path = qwqdsp_support::WormholeWav();
    std::cout << std::format("loading {}\n", wav_path) << std::flush;

    AudioFile<float> file{wav_path};
    auto& x_vec = file.samples.front();
    const float fs = file.getSampleRate();
    std::cout << std::format("sample_rate={}, len={}\n", fs, x_vec.size()) << std::flush;

    auto save_wav = [fs](const std::string& name, const std::vector<float>& out) {
        AudioFile<float> of;
        of.setBitDepth(32);
        of.setNumSamplesPerChannel(out.size());
        of.samples[0] = out;
        of.setNumChannels(1);
        of.setSampleRate(static_cast<int>(fs));
        of.save(qwqdsp_support::OutputFile(name));
        std::cout << std::format("saved {}\n", name) << std::flush;
    };

    // --------------------------------------------------------
    // 1. 仅音高偏移（共振峰不变）
    // --------------------------------------------------------
    {
        Psola psola;
        psola.sample_rate_ = fs;
        psola.pitch_shift_semitones_ = 7.0f; // 升纯五度
        psola.formant_shift_ = 1.0f;

        auto out = psola.Process(x_vec);
        std::cout << std::format("pitch_only: out_len={}\n", out.size()) << std::flush;
        save_wav("psola_pitch_only.wav", out);
    }

    // --------------------------------------------------------
    // 2. 仅共振峰偏移（音高不变）
    // --------------------------------------------------------
    {
        Psola psola;
        psola.sample_rate_ = fs;
        psola.pitch_shift_semitones_ = 0.0f; // 音高不变
        psola.formant_shift_ = 1.3f;         // 共振峰升高 1.3 倍

        auto out = psola.Process(x_vec);
        std::cout << std::format("formant_only: out_len={}\n", out.size()) << std::flush;
        save_wav("psola_formant_only.wav", out);
    }

    // --------------------------------------------------------
    // 3. 音高升高 + 共振峰降低（变矮胖声音）
    // --------------------------------------------------------
    {
        Psola psola;
        psola.sample_rate_ = fs;
        psola.pitch_shift_semitones_ = 5.0f; // 音高升
        psola.formant_shift_ = 0.7f;         // 共振峰降低

        auto out = psola.Process(x_vec);
        std::cout << std::format("pitch_up_formant_down: out_len={}\n", out.size()) << std::flush;
        save_wav("psola_pitch_up_formant_down.wav", out);
    }

    // --------------------------------------------------------
    // 4. 音高降低 + 共振峰升高（花栗鼠效果反转）
    // --------------------------------------------------------
    {
        Psola psola;
        psola.sample_rate_ = fs;
        psola.pitch_shift_semitones_ = -5.0f; // 音高降
        psola.formant_shift_ = 1.4f;          // 共振峰升高

        auto out = psola.Process(x_vec);
        std::cout << std::format("pitch_down_formant_up: out_len={}\n", out.size()) << std::flush;
        save_wav("psola_pitch_down_formant_up.wav", out);
    }

    std::cout << "\n" << std::flush;
}

// ------------------------------------------------------------
// Phase Vocoder 2（本地 PGHI 实现）
// ------------------------------------------------------------
/** 本地相位梯度声码器：正常音频处理并保存 */
static void ProcessPV2(const char* name, float kt, float kp) {
    auto wav_path = qwqdsp_support::SweepWav();
    AudioFile<float> file{wav_path};
    auto& x_vec = file.samples.front();

    std::cout << std::format("{} (kt={:.2f}, kp={:.2f}): {} -> ", name, kt, kp, x_vec.size()) << std::flush;

    PhaseGradientVocoder dsp;
    dsp.SetFrameSize(4096);
    dsp.SetOverSample(2);
    dsp.SetTimeStretch(kt);
    dsp.SetPitchShift(kp);

    auto out = dsp.Process(x_vec);
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
 * @brief 使用 delta（8192 处为 1.0）作为输入测试声码器。
 */
static void TestDeltaPV2(const char* name, float kt, float kp) {
    constexpr size_t kLen = 65536;
    constexpr size_t kDeltaPos = 8192;
    std::vector<float> x(kLen, 0.0f);
    x[kDeltaPos] = 1.0f;

    std::cout << std::format("{} (kt={:.2f}, kp={:.2f}): delta@{} -> ", name, kt, kp, kDeltaPos) << std::flush;

    PhaseGradientVocoder dsp;
    dsp.SetFrameSize(4096);
    dsp.SetOverSample(2);
    dsp.SetTimeStretch(kt);
    dsp.SetPitchShift(kp);

    auto out = dsp.Process(x);
    std::cout << std::format("{}\n", out.size()) << std::flush;

    qwqdsp_support::AudioOps::Normalize(out);
    AudioFile<float> file;
    file.setNumSamplesPerChannel(out.size());
    file.samples.resize(1);
    file.samples[0] = out;
    file.setNumChannels(1);
    file.save(qwqdsp_support::OutputFile(name));
    std::cout << std::format("saved {}\n\n", name) << std::flush;
}

// ------------------------------------------------------------
// Phase Vocoder 3（库版 PGHI + 离线分帧）
// ------------------------------------------------------------
/** 库版相位梯度声码器：正常音频处理并保存 */
static void ProcessPV3(const char* name, float kt, float kp) {
    auto wav_path = qwqdsp_support::SweepWav();
    AudioFile<float> file{wav_path};
    auto& x_vec = file.samples.front();

    std::cout << std::format("{} (kt={:.2f}, kp={:.2f}): {} -> ", name, kt, kp, x_vec.size()) << std::flush;

    auto out = RunPhaseVocoder3(x_vec, kt, kp);

    std::cout << std::format("{}\n", out.size()) << std::flush;

    qwqdsp_support::AudioOps::Normalize(out);
    file.setNumSamplesPerChannel(out.size());
    file.samples[0] = out;
    file.setNumChannels(1);
    file.save(qwqdsp_support::OutputFile(name));
    std::cout << std::format("saved {}\n\n", name) << std::flush;
}

/**
 * @brief 使用 delta（8192 处为 1.0）作为输入测试库版声码器。
 */
static void TestDeltaPV3(const char* name, float kt, float kp) {
    constexpr size_t kLen = 65536;
    constexpr size_t kDeltaPos = 8192;
    std::vector<float> x(kLen, 0.0f);
    x[kDeltaPos] = 1.0f;

    std::cout << std::format("{} (kt={:.2f}, kp={:.2f}): delta@{} -> ", name, kt, kp, kDeltaPos) << std::flush;

    auto out = RunPhaseVocoder3(x, kt, kp);
    std::cout << std::format("{}\n", out.size()) << std::flush;

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
    RunWsolaTest();

    RunPsolaTest();

    ProcessPV2("PV2_ts_1.5x.wav", 1.5f, 1.0f);
    ProcessPV2("PV2_ps_1.5x.wav", 1.0f, 1.5f);
    ProcessPV2("PV2_ts1.5_ps1.5.wav", 1.5f, 1.5f);

    TestDeltaPV2("PV2_delta_identity.wav", 1.0f, 1.0f);
    TestDeltaPV2("PV2_delta_ts1.5.wav", 1.5f, 1.0f);
    TestDeltaPV2("PV2_delta_ps1.5.wav", 1.0f, 1.5f);

    ProcessPV3("PV3_ts_1.5x.wav", 1.5f, 1.0f);
    ProcessPV3("PV3_ps_1.5x.wav", 1.0f, 1.5f);
    ProcessPV3("PV3_ts1.5_ps1.5.wav", 1.5f, 1.5f);

    TestDeltaPV3("PV3_delta_identity.wav", 1.0f, 1.0f);
    TestDeltaPV3("PV3_delta_ts1.5.wav", 1.5f, 1.0f);
    TestDeltaPV3("PV3_delta_ps1.5.wav", 1.0f, 1.5f);

    std::cout << std::format("all done\n") << std::flush;
}
