#include <AudioFile.h>
#include <cmath>
#include <qwqdsp/pitch/simple_voicing_detector.hpp>
#include <vector>
#include <work_dir.hpp>

// ------------------------------------------------------------
//  vad — 使用 FixedVoicingDetector 将 wormhole.wav 分离为
//        浊音/清音两部分并分别导出
// ------------------------------------------------------------

int main() {
    // ── 加载音频 ──
    AudioFile<float> file{qwqdsp_support::WormholeWav()};
    auto const& x = file.samples.front();
    float const fs = file.getSampleRate();
    size_t const n = x.size();

    // ── 检测器 ──
    qwqdsp_pitch::SimpleVoicingDetector detector;
    detector.Init(fs);

    // ── 逐采样处理 ──
    std::vector<float> voiced(n, 0.0f);
    std::vector<float> unvoiced(n, 0.0f);

    for (size_t i = 0; i < n; ++i) {
        detector.ProcessSample(x[i]);

        float prob = detector.FrameResult();

        if (prob >= qwqdsp_pitch::SimpleVoicingDetector::kThreshold) {
            voiced[i] = x[i];   // 浊音保留
            unvoiced[i] = 0.0f; // 清音静音
        }
        else {
            voiced[i] = 0.0f;   // 浊音静音
            unvoiced[i] = x[i]; // 清音保留
        }
    }

    // ── 导出 ──
    auto save = [&](std::string const& name, std::vector<float> const& data) {
        AudioFile<float> out;
        out.setNumChannels(1);
        out.setNumSamplesPerChannel(data.size());
        out.samples[0] = data;
        out.setSampleRate(static_cast<int>(fs));
        out.save(qwqdsp_support::OutputFile(name));
    };

    save("vad_voiced.wav", voiced);
    save("vad_unvoiced.wav", unvoiced);
}
