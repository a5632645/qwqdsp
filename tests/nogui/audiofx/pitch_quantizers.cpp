#include <AudioFile.h>
#include <audio_ops.hpp>
#include <work_dir.hpp>
#include <iostream>
#include <format>

#include "pitch_quantize/scale_helper.hpp"
#include "pitch_quantize/pitch_quantizer1.hpp"
#include "pitch_quantize/pitch_quantizer2.hpp"
#include "pitch_quantize/pitch_quantizer3.hpp"
#include "pitch_quantize/pitch_quantizer4.hpp"
#include "pitch_quantize/pitch_quantizer6.hpp"

// ------------------------------------------------------------
// 辅助：运行单个量化器并保存结果
// ------------------------------------------------------------
template <typename T>
static void RunAndSave(T& pq, std::span<const float> input, float fs,
                       const std::string& label, const std::string& out_name) {
    std::cout << std::format("\n===== {} =====\n", label) << std::flush;
    pq.SetKey(0, ScaleHelper::Type::kMajor);
    std::cout << std::format("key: {}\n", pq.KeyDescription()) << std::flush;

    auto out = pq.Process(input);
    std::cout << std::format("output len={}\n", out.size()) << std::flush;

    qwqdsp_support::AudioOps::Normalize(out);

    AudioFile<float> file;
    file.setBitDepth(32);
    file.setNumSamplesPerChannel(out.size());
    file.samples[0] = out;
    file.setNumChannels(1);
    file.setSampleRate(static_cast<int>(fs));
    file.save(qwqdsp_support::OutputFile(out_name));

    std::cout << std::format("saved {}\n", out_name) << std::flush;
}

// ------------------------------------------------------------
// main
// ------------------------------------------------------------
int main() {
    // ---- 输入 1: WormholeWav (人声) ----
    {
        const auto wav_path = qwqdsp_support::WormholeWav();
        std::cout << std::format("loading {}\n", wav_path) << std::flush;

        AudioFile<float> file{wav_path};
        auto& x_vec = file.samples.front();
        const float fs = file.getSampleRate();
        std::cout << std::format("sample_rate={}, len={}\n", fs, x_vec.size()) << std::flush;

        constexpr size_t kFftSize = 2048;

        // 版本 1: 朴素 PV
        {
            PitchQuantizer<false> pq;
            pq.Init2(fs, 3000, kFftSize);
            RunAndSave(pq, x_vec, fs, "PitchQuantizer (v1, naive PV)", "pitch_quantizer1_vocal.wav");
        }

        // 版本 2: 逐 bin 映射 + Puckette phase lock
        {
            PitchQuantizer2<true> pq;
            pq.Init(fs, kFftSize / 4, kFftSize);
            RunAndSave(pq, x_vec, fs, "PitchQuantizer2 (v2, per-bin mapping)", "pitch_quantizer2_vocal.wav");
        }

        // 版本 3: 峰域映射 + PGHI phase lock
        {
            PitchQuantizer3<true> pq;
            pq.Init(fs, kFftSize / 4, kFftSize);
            RunAndSave(pq, x_vec, fs, "PitchQuantizer3 (v3, peak-domain + PGHI phase lock)", "pitch_quantizer3_vocal.wav");
        }

        // 版本 4: 峰域映射 + 原始分析相位 phase lock
        {
            PitchQuantizer4<true> pq;
            pq.Init(fs, kFftSize / 4, kFftSize);
            RunAndSave(pq, x_vec, fs, "PitchQuantizer4 (v4, peak-domain + raw phase lock)", "pitch_quantizer4_vocal.wav");
        }

        // 版本 6: 峰域映射 + 纯 PGHI
        {
            PitchQuantizer6<true> pq;
            pq.Init(fs, kFftSize / 4, kFftSize);
            RunAndSave(pq, x_vec, fs, "PitchQuantizer6 (v6, peak-domain + pure PGHI)", "pitch_quantizer6_vocal.wav");
        }
    }

    // ---- 输入 2: drumloop.wav ----
    {
        const auto wav_path = qwqdsp_support::InputFile("drumloop.wav");
        std::cout << std::format("\nloading {}\n", wav_path) << std::flush;

        AudioFile<float> file{wav_path};
        auto& x_vec = file.samples.front();
        const float fs_drum = file.getSampleRate();
        std::cout << std::format("sample_rate={}, len={}\n", fs_drum, x_vec.size()) << std::flush;

        constexpr size_t kFftSize = 2048;

        // 版本 1: 朴素 PV (drumloop)
        {
            PitchQuantizer<false> pq;
            pq.Init2(fs_drum, 3000, kFftSize);
            RunAndSave(pq, x_vec, fs_drum, "PitchQuantizer (v1, drumloop)", "pitch_quantizer1_drum.wav");
        }

        // 版本 2: 逐 bin 映射 + Puckette phase lock (drumloop)
        {
            PitchQuantizer2<true> pq;
            pq.Init(fs_drum, kFftSize / 4, kFftSize);
            RunAndSave(pq, x_vec, fs_drum, "PitchQuantizer2 (v2, drumloop)", "pitch_quantizer2_drum.wav");
        }

        // 版本 3: 峰域映射 + PGHI phase lock (drumloop)
        {
            PitchQuantizer3<true> pq;
            pq.Init(fs_drum, kFftSize / 4, kFftSize);
            RunAndSave(pq, x_vec, fs_drum, "PitchQuantizer3 (v3, drumloop)", "pitch_quantizer3_drum.wav");
        }

        // 版本 4: 峰域映射 + 原始分析相位 phase lock (drumloop)
        {
            PitchQuantizer4<true> pq;
            pq.Init(fs_drum, kFftSize / 4, kFftSize);
            RunAndSave(pq, x_vec, fs_drum, "PitchQuantizer4 (v4, drumloop)", "pitch_quantizer4_drum.wav");
        }

        // 版本 6: 峰域映射 + 纯 PGHI (drumloop)
        {
            PitchQuantizer6<true> pq;
            pq.Init(fs_drum, kFftSize / 4, kFftSize);
            RunAndSave(pq, x_vec, fs_drum, "PitchQuantizer6 (v6, drumloop)", "pitch_quantizer6_drum.wav");
        }
    }

    std::cout << "\nAll pitch quantizers done.\n" << std::flush;
    return 0;
}
