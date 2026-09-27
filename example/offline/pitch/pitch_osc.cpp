#include <AudioFile.h>
#include <qwqdsp/convert.hpp>
#include <qwqdsp/filter/svf.hpp>
#include <qwqdsp/oscillator/noise.hpp>
#include <qwqdsp/oscillator/vic_sine_osc.hpp>
#include <qwqdsp/pitch/yin.hpp>
#include <qwqdsp/segement/analyze_auto.hpp>
#include <work_dir.hpp>

class PitchOsc {
public:
    void Init(float fs, int hop, int size) {
        fs_ = fs;
        hop_ = hop;

        yin_.SetThreshold(0.2f);
        yin_.SetMinPitch(50.0f);
        yin_.SetMaxPitch(500.0f);
        yin_.Init(fs, size);
    }

    void operator()(std::span<const float> block) noexcept {
        yin_.Process(block);

        auto pitch = yin_.GetPitch();
        float impluse_g = pitch.non_period_ratio < 0.3f ? 1.0f - pitch.non_period_ratio : 0.0f;
        float noise_g = (1.0f - impluse_g) * 0.1f;
        sine_.SetFreq(pitch.pitch_hz, fs_);
        for (int i = 0; i < hop_; ++i) {
            y_.push_back(sine_.Tick() * impluse_g + noise_g * noise_.Next());
        }
    }

    std::vector<float> GetOutput() {
        return y_;
    }
private:
    float fs_;
    int hop_;
    std::vector<float> y_;
    qwqdsp_oscillator::VicSineOsc sine_;
    qwqdsp_oscillator::WhiteNoise noise_;
    qwqdsp_pitch::Yin yin_;
};

int main() {
    AudioFile<float> file{qwqdsp_support::WormholeWav()};

    // qwqdsp_filter::SVF lpf;
    // lpf.MakeLowpass(qwqdsp::convert::Freq2W(500.0f, file.getSampleRate()), std::numbers::sqrt2_v<float> / 2);
    // for (auto& x : file.samples.front()) {
    //     x = lpf.Tick(x);
    // }

    qwqdsp_segement::AnalyzeAuto<false> aa;
    aa.SetHop(512);
    aa.SetSize(1024);

    PitchOsc pitch_osc;
    pitch_osc.Init(file.getSampleRate(), 512, 1024);
    aa.Process(file.samples.front(), pitch_osc);

    file.samples.front() = pitch_osc.GetOutput();
    file.setNumSamplesPerChannel(file.samples.front().size());
    file.save(qwqdsp_support::OutputFile("pitch_osc.wav"));
}
