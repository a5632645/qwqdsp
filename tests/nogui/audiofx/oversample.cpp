#include <cstdio>
#include <vector>

#include "AudioFile.h"
#include <qwqdsp/algebraic_waveshaper.hpp>
#include <qwqdsp/fx/oversample.hpp>
#include <work_dir.hpp>

// kHalfbandCoeffs (N=65)
static constexpr std::array<float, 65> kHalfbandCoeffs = {
    -0.0000003856f, -0.0003044707f, 0.0000010085f,  0.0005326581f,  -0.0000017348f, -0.0009739876f, 0.0000027020f,
    0.0016310629f,  -0.0000039780f, -0.0025679886f, 0.0000054613f,  0.0038614894f,  -0.0000072091f, -0.0056054843f,
    0.0000090709f,  0.0079210177f,  -0.0000111919f, -0.0109748173f, 0.0000132931f,  0.0150179768f,  -0.0000152761f,
    -0.0204685818f, 0.0000171944f,  0.0281074181f,  -0.0000190393f, -0.0396091882f, 0.0000203860f,  0.0593530469f,
    -0.0000213413f, -0.1034683564f, 0.0000220592f,  0.3174232112f,  0.4999779612f,  0.3174232112f,  0.0000220592f,
    -0.1034683564f, -0.0000213413f, 0.0593530469f,  0.0000203860f,  -0.0396091882f, -0.0000190393f, 0.0281074181f,
    0.0000171944f,  -0.0204685818f, -0.0000152761f, 0.0150179768f,  0.0000132931f,  -0.0109748173f, -0.0000111919f,
    0.0079210177f,  0.0000090709f,  -0.0056054843f, -0.0000072091f, 0.0038614894f,  0.0000054613f,  -0.0025679886f,
    -0.0000039780f, 0.0016310629f,  0.0000027020f,  -0.0009739876f, -0.0000017348f, 0.0005326581f,  0.0000010085f,
    -0.0003044707f, -0.0000003856f,
};

// ------------------------------------------------------------
// 测试 7：1kHz 正弦波失真对比（有无超采样）
// ------------------------------------------------------------
static void TestDistortionSweep() {
    std::printf("--- TestDistortionSweep ---\n");

    AudioFile<float> file;
    if (!file.load(qwqdsp_support::InputFile("sweep.wav"))) {
        std::printf("  SKIP (cannot load sweep.wav)\n");
        return;
    }

    auto const& input_ref = file.samples[0];
    int const n_samp = static_cast<int>(input_ref.size());
    float const kFs = file.getSampleRate();
    std::vector<float> input(input_ref.begin(), input_ref.end()); // 可变副本

    // 失真器：AlgebraicWaveshaper Naive
    qwqdsp::AlgebraicWaveshaper ws;

    // ---- 无超采样直接失真 ----
    std::vector<float> direct(n_samp);
    for (int i = 0; i < n_samp; ++i)
        direct[i] = ws.Naive(input[i]);

    // ---- 2x 超采样后失真 ----
    ws.Reset();
    int const stages = 1; // 2x
    qwqdsp_fx::Oversample ov;
    ov.Init(kHalfbandCoeffs, stages);

    int const up_n = n_samp * (1 << stages);
    std::vector<float> up(up_n);
    std::vector<float> processed(up_n);
    std::vector<float> down(n_samp);

    ov.Reset();
    ov.Upsample(input, up);
    for (int i = 0; i < up_n; ++i)
        processed[i] = ws.Naive(20 * up[i]);
    ov.Downsample(processed, down);

    // ---- 导出 ----
    auto save = [&](std::span<float> data, char const* name, int fs) {
        AudioFile<float> f;
        f.setNumChannels(1);
        f.setBitDepth(32);
        f.setSampleRate(fs);
        f.setNumSamplesPerChannel(data.size());
        std::copy(data.begin(), data.end(), f.samples[0].begin());
        f.save(qwqdsp_support::OutputFile(name));
    };

    save(down, "sweep_oversample_8x.wav", static_cast<int>(kFs));

    std::printf("  PASS (exported sweep_*.wav)\n");
}

// ------------------------------------------------------------
// AlgebraicWaveshaper 变体对比
// ------------------------------------------------------------
static void TestWaveshaperVariants() {
    std::printf("--- TestWaveshaperVariants ---\n");

    AudioFile<float> file;
    if (!file.load(qwqdsp_support::InputFile("sweep.wav"))) {
        std::printf("  SKIP (cannot load sweep.wav)\n");
        return;
    }

    auto const& input_ref = file.samples[0];
    int const n_samp = static_cast<int>(input_ref.size());
    float const kFs = file.getSampleRate();

    auto save = [&](std::span<float> data, char const* name) {
        AudioFile<float> f;
        f.setNumChannels(1);
        f.setBitDepth(32);
        f.setSampleRate(static_cast<int>(kFs));
        f.setNumSamplesPerChannel(data.size());
        std::copy(data.begin(), data.end(), f.samples[0].begin());
        f.save(qwqdsp_support::OutputFile(name));
    };

    // Naive — 无状态
    {
        std::vector<float> out(n_samp);
        for (int i = 0; i < n_samp; ++i)
            out[i] = qwqdsp::AlgebraicWaveshaper::Naive(20 * input_ref[i]);
        save(out, "ws_naive.wav");
    }

    // ADAA — 0.5 sample latency
    {
        qwqdsp::AlgebraicWaveshaper ws;
        ws.Reset();
        std::vector<float> out(n_samp);
        for (int i = 0; i < n_samp; ++i)
            out[i] = ws.ADAA(20 * input_ref[i]);
        save(out, "ws_adaa.wav");
    }

    // ADAA_MV — 1 sample latency
    {
        qwqdsp::AlgebraicWaveshaper ws;
        ws.Reset();
        std::vector<float> out(n_samp);
        for (int i = 0; i < n_samp; ++i)
            out[i] = ws.ADAA_MV(20 * input_ref[i]);
        save(out, "ws_adaa_mv.wav");
    }

    // ADAA_MV_Compensation — >=1 sample latency
    {
        qwqdsp::AlgebraicWaveshaper ws;
        ws.Reset();
        std::vector<float> out(n_samp);
        for (int i = 0; i < n_samp; ++i)
            out[i] = ws.ADAA_MV_Compensation(20 * input_ref[i]);
        save(out, "ws_adaa_mv_comp.wav");
    }

    std::printf("  PASS (exported ws_*.wav)\n");
}

// ------------------------------------------------------------
// 过采样 + ADAA_MV_Compensation 联合测试
// ------------------------------------------------------------
static void TestOversampleADAA() {
    std::printf("--- TestOversampleADAA ---\n");

    AudioFile<float> file;
    if (!file.load(qwqdsp_support::InputFile("sweep.wav"))) {
        std::printf("  SKIP (cannot load sweep.wav)\n");
        return;
    }

    auto const& input_ref = file.samples[0];
    int const n_samp = static_cast<int>(input_ref.size());
    float const kFs = file.getSampleRate();
    std::vector<float> input(input_ref.begin(), input_ref.end());

    auto save = [&](std::span<float> data, char const* name) {
        AudioFile<float> f;
        f.setNumChannels(1);
        f.setBitDepth(32);
        f.setSampleRate(static_cast<int>(kFs));
        f.setNumSamplesPerChannel(data.size());
        std::copy(data.begin(), data.end(), f.samples[0].begin());
        f.save(qwqdsp_support::OutputFile(name));
    };

    auto constexpr kGain = 20.0f;

    // 无超采样 + ADAA_MV_Compensation
    {
        qwqdsp::AlgebraicWaveshaper ws;
        ws.Reset();
        std::vector<float> out(n_samp);
        for (int i = 0; i < n_samp; ++i)
            out[i] = ws.ADAA_MV_Compensation(kGain * input[i]);
        save(out, "ov_adaa_direct.wav");
    }

    // 2x 超采样 -> ADAA_MV_Compensation -> 降采样
    int const stages = 1;
    qwqdsp_fx::Oversample ov;
    ov.Init(kHalfbandCoeffs, stages);

    int const up_n = n_samp * (1 << stages);
    std::vector<float> up(up_n);
    std::vector<float> processed(up_n);
    std::vector<float> down(n_samp);

    ov.Reset();
    ov.Upsample(input, up);
    {
        qwqdsp::AlgebraicWaveshaper ws;
        ws.Reset();
        for (int i = 0; i < up_n; ++i)
            processed[i] = ws.ADAA_MV_Compensation(kGain * up[i]);
    }
    ov.Downsample(processed, down);
    save(down, "ov_adaa_oversample_8x.wav");

    std::printf("  PASS (exported ov_adaa_*.wav)\n");
}

// ------------------------------------------------------------
// main
// ------------------------------------------------------------
int main() {
    TestDistortionSweep();
    TestWaveshaperVariants();
    TestOversampleADAA();
}
