#include <AudioFile.h>
#include <audio_ops.hpp>
#include <format>
#include <iostream>
#include <vector>
#include <work_dir.hpp>

#include "pitch_shift/peak_map_pitch_shifter.hpp"
#include "pitch_shift/phase_gradient_transient.hpp"
#include "pitch_shift/phase_locked_vocoder.hpp"
#include "pitch_shift/phase_vocoder2.hpp"
#include "pitch_shift/phase_vocoder3.hpp"
#include "pitch_shift/pitcher_wsola.hpp"
#include "pitch_shift/psola.hpp"
#include "pitch_shift/transient_vocoder.hpp"
#include "pitch_shift/wsola.hpp"

using qwqdsp_test::OdfType;
using qwqdsp_test::PeakMapPitchShifter;
using qwqdsp_test::PhaseGradientTransientVocoder;
using qwqdsp_test::PhaseGradientVocoder;
using qwqdsp_test::PhaseLockedVocoder;
using qwqdsp_test::PitcherWsola;
using qwqdsp_test::Psola;
using qwqdsp_test::RunPhaseVocoder3;
using qwqdsp_test::RunWsola;
using qwqdsp_test::TransientVocoder;

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
    auto wav_path = qwqdsp_support::InputFile("sine.wav");
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
// Phase Vocoder 3（库版 PGHI + 离线分帧）
// ------------------------------------------------------------
/** 库版相位梯度声码器：正常音频处理并保存 */
static void ProcessPV3(const char* name, float kt, float kp) {
    auto wav_path = qwqdsp_support::InputFile("sine.wav");
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

// ------------------------------------------------------------
// Phase Locked Vocoder（相位锁定声码器）
// ------------------------------------------------------------
/**
 * @brief 相位锁定声码器（Laroche & Dolson identity phase locking），
 *        恢复谐波相位相干，消除普通声码器的 "phasiness"。
 */
static void ProcessPVL(const char* name, float kt, float kp) {
    auto wav_path = qwqdsp_support::InputFile("sine.wav");
    AudioFile<float> file{wav_path};
    auto& x_vec = file.samples.front();

    std::cout << std::format("{} (kt={:.2f}, kp={:.2f}): {} -> ", name, kt, kp, x_vec.size()) << std::flush;

    PhaseLockedVocoder dsp;
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
// Transient Vocoder（瞬态感知声码器）
// ------------------------------------------------------------
/**
 * @brief 瞬态感知相位锁定声码器（Röbel 2003），谱通量检测瞬态，
 *        起始处重置相位保留攻击锐度。
 */
static void ProcessTRV(const char* name, float kt, float kp) {
    auto wav_path = qwqdsp_support::InputFile("sine.wav");
    AudioFile<float> file{wav_path};
    auto& x_vec = file.samples.front();

    std::cout << std::format("{} (kt={:.2f}, kp={:.2f}): {} -> ", name, kt, kp, x_vec.size()) << std::flush;

    TransientVocoder dsp;
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
// Phase Gradient + Flux/SuperFlux 瞬态（PGHI 结合瞬态检测）
// ------------------------------------------------------------
/**
 * @brief PGHI 相位声码器 + flux/superflux 瞬态检测（瞬态帧相位重置）。
 *
 * @tparam Odf 瞬态 ODF 类型（Flux 或 SuperFlux）
 */
template <qwqdsp_test::OdfType Odf>
static void ProcessPGT(const char* name, float kt, float kp) {
    auto wav_path = qwqdsp_support::InputFile("sine.wav");
    AudioFile<float> file{wav_path};
    auto& x_vec = file.samples.front();

    std::cout << std::format("{} (kt={:.2f}, kp={:.2f}): {} -> ", name, kt, kp, x_vec.size()) << std::flush;

    PhaseGradientTransientVocoder<Odf> dsp;
    dsp.SetFrameSize(4096);
    dsp.SetOverSample(2);
    dsp.SetTimeStretch(kt);
    dsp.SetPitchShift(kp);
    dsp.SetSampleRate(static_cast<float>(file.getSampleRate()));

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
// PitcherWsola（WSOLA 音高移位器，移植自 qdelay Pitcher）
// ------------------------------------------------------------
/**
 * @brief WSOLA 环形缓冲音高移位器：半音输入，输出等长信号。
 */
static void ProcessPitcher(const char* name, float semitones, PitcherWsola::WindowMode mode) {
    auto wav_path = qwqdsp_support::InputFile("sine.wav");
    AudioFile<float> file{wav_path};
    auto& x_vec = file.samples.front();

    std::cout << std::format("{} (semitones={:.1f}): {} -> ", name, semitones, x_vec.size()) << std::flush;

    PitcherWsola dsp;
    dsp.SetWindowMode(mode);
    dsp.SetSemitones(semitones);

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
// PeakMapPitchShifter（峰域映射 + Phase Lock 音高移位器）
// ------------------------------------------------------------
/**
 * @brief 峰域映射音高移位器（参考 pitch_quantizer4）：按 kp 比例移动谱峰。
 */
static void ProcessPMS(const char* name, float kp) {
    auto wav_path = qwqdsp_support::InputFile("sine.wav");
    AudioFile<float> file{wav_path};
    auto& x_vec = file.samples.front();

    std::cout << std::format("{} (kp={:.2f}): {} -> ", name, kp, x_vec.size()) << std::flush;

    PeakMapPitchShifter dsp;
    dsp.Init(static_cast<float>(file.getSampleRate()), 512, 2048);
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

int main() {
    // RunWsolaTest();

    // RunPsolaTest();

    ProcessPV2("PV2_ts_1.5x.wav", 1.5f, 1.0f);
    ProcessPV2("PV2_ps_1.5x.wav", 1.0f, 2.0f);
    ProcessPV2("PV2_ts1.5_ps1.5.wav", 1.5f, 1.5f);

    ProcessPV3("PV3_ts_1.5x.wav", 1.5f, 1.0f);
    ProcessPV3("PV3_ps_1.5x.wav", 1.0f, 2.0f);
    ProcessPV3("PV3_ts1.5_ps1.5.wav", 1.5f, 1.5f);

    ProcessPVL("pvl_ts_1.5x.wav", 1.5f, 1.0f);
    ProcessPVL("pvl_ps_1.5x.wav", 1.0f, 2.0f);
    ProcessPVL("pvl_ts1.5_ps1.5.wav", 1.5f, 1.5f);

    ProcessTRV("trv_ts_1.5x.wav", 1.5f, 1.0f);
    ProcessTRV("trv_ps_1.5x.wav", 1.0f, 2.0f);
    ProcessTRV("trv_ts1.5_ps1.5.wav", 1.5f, 1.5f);

    ProcessPGT<qwqdsp_test::OdfType::Flux>("pgt_flux_ts_1.5x.wav", 1.5f, 1.0f);
    ProcessPGT<qwqdsp_test::OdfType::Flux>("pgt_flux_ps_1.5x.wav", 1.0f, 1.5f);
    ProcessPGT<qwqdsp_test::OdfType::Flux>("pgt_flux_ts1.5_ps1.5.wav", 1.5f, 1.5f);

    ProcessPGT<qwqdsp_test::OdfType::SuperFlux>("pgt_superflux_ts_1.5x.wav", 1.5f, 1.0f);
    ProcessPGT<qwqdsp_test::OdfType::SuperFlux>("pgt_superflux_ps_1.5x.wav", 1.0f, 1.5f);
    ProcessPGT<qwqdsp_test::OdfType::SuperFlux>("pgt_superflux_ts1.5_ps1.5.wav", 1.5f, 1.5f);

    ProcessPitcher("pitcher_up7_kmedium.wav", 7.0f, PitcherWsola::WindowMode::kMedium);
    ProcessPitcher("pitcher_down5_kmedium.wav", -5.0f, PitcherWsola::WindowMode::kMedium);
    ProcessPitcher("pitcher_up12_klarge.wav", 12.0f, PitcherWsola::WindowMode::kLarge);

    ProcessPMS("pms_ps_1.5x.wav", 1.5f);
    ProcessPMS("pms_ps_2.0x.wav", 2.0f);
    ProcessPMS("pms_ps_0.5x.wav", 0.5f);

    std::cout << std::format("all done\n") << std::flush;
}
