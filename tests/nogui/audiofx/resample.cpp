#include "AudioFile.h"
#include "work_dir.hpp"
#include <qwqdsp/filter/window_fir.hpp>
#include <qwqdsp/fx/polyphase_resample_fir.hpp>
#include <qwqdsp/fx/resample.hpp>
#include <qwqdsp/fx/resample_coeffs.h>
#include <qwqdsp/fx/resample_iir.hpp>
#include <qwqdsp/window/blackman_nuttall.hpp>


constexpr float kTargetFs = 88200.0f;

static void ResampleIir() {
    AudioFile<float> infile;
    infile.load(qwqdsp_support::InputFile("sweep.wav"));
    auto& sweep = infile.samples.front();

    qwqdsp_fx::ResampleIIR<qwqdsp_fx::coeff::BestCoeffs<float>, 128> resample;
    resample.Init(infile.getSampleRate(), kTargetFs);
    auto sweep_resample = resample.Process<float>(sweep);

    AudioFile<float> outfile;
    outfile.setNumChannels(1);
    outfile.setBitDepth(32);
    outfile.setSampleRate(kTargetFs);
    outfile.setNumSamplesPerChannel(sweep_resample.size());
    outfile.samples.front() = std::move(sweep_resample);
    outfile.save(qwqdsp_support::OutputFile("sweep_resample_iir.wav"));
}

static void ResampleFir() {
    AudioFile<float> infile;
    infile.load(qwqdsp_support::InputFile("sweep.wav"));
    auto& sweep = infile.samples.front();

    qwqdsp_fx::Resample resample;
    resample.Init(infile.getSampleRate(), kTargetFs, 180, 127, 128);
    auto sweep_resample = resample.Process(sweep);

    AudioFile<float> outfile;
    outfile.setNumChannels(1);
    outfile.setBitDepth(32);
    outfile.setSampleRate(kTargetFs);
    outfile.setNumSamplesPerChannel(sweep_resample.size());
    outfile.samples.front() = std::move(sweep_resample);
    outfile.save(qwqdsp_support::OutputFile("sweep_resample_fir.wav"));
}

static void ResamplePolyphaseFir() {
    constexpr int coeff_len = 1024;
    constexpr int ratio = 2;
    float coeff[coeff_len];
    qwqdsp_filter::WindowFIR::Lowpass(coeff, std::numbers::pi_v<float> / ratio);
    qwqdsp_window::BlackmanNuttall::ApplyWindow(coeff, false);

    AudioFile<float> infile;
    infile.load(qwqdsp_support::InputFile("sweep.wav"));
    auto& sweep = infile.samples.front();

    {
        qwqdsp_fx::PolyphaseDownsamplerFir resample;
        resample.Init(coeff, ratio);
        resample.Reset();
        auto sweep_resample = resample.Process(sweep);

        AudioFile<float> outfile;
        outfile.setNumChannels(1);
        outfile.setBitDepth(32);
        outfile.setSampleRate(infile.getSampleRate() / ratio);
        outfile.setNumSamplesPerChannel(sweep_resample.size());
        outfile.samples.front() = std::move(sweep_resample);
        outfile.save(qwqdsp_support::OutputFile("sweep_polyphase_resample_fir_down.wav"));
    }
    {
        qwqdsp_fx::PolyphaseUpsamplerFir resample;
        resample.Init(coeff, ratio);
        resample.Reset();
        auto sweep_resample = resample.Process(sweep);

        AudioFile<float> outfile;
        outfile.setNumChannels(1);
        outfile.setBitDepth(32);
        outfile.setSampleRate(infile.getSampleRate() * ratio);
        outfile.setNumSamplesPerChannel(sweep_resample.size());
        outfile.samples.front() = std::move(sweep_resample);
        outfile.save(qwqdsp_support::OutputFile("sweep_polyphase_resample_fir_up.wav"));
    }
}

int main() {
    ResampleIir();
    ResampleFir();
    ResamplePolyphaseFir();
}
