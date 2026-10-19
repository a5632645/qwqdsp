#include <AudioFile.h>
#include <cmath>
#include <cstdlib>
#include <format>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <work_dir.hpp>

#include "pitch_shifter_probe.hpp"

using pitch_shift_rt::Algorithm;
using pitch_shift_rt::RealtimePitchShifter;
using pitch_shift_rt::WindowType;

namespace {

using Shifter = RealtimePitchShifter<WindowType::BlackmanHarrisThreeTerm, 512>;

struct ProbeCase {
    std::string name;
    Shifter::ProbeConfig config;
};

void PrintInterestingTrace(const Shifter& shifter) {
    for (const auto& frame : shifter.getDetectorTrace()) {
        if (!frame.rise && frame.guard && !frame.transient) {
            continue;
        }
        std::cout << std::format(
            "  trace frame={:4} t={:8.4f}s hop={:4} flux={:.6f} peak={:.6f} median={:.6f} thr={:.6f} rise={} "
            "guard={} cooldown={} reset={}\n",
            frame.frame_index, frame.time_seconds, frame.analysis_hop, frame.flux, frame.flux_peak, frame.median,
            frame.threshold, frame.rise ? 1 : 0, frame.guard ? 1 : 0, frame.cooldown_before,
            frame.transient ? 1 : 0);
    }
}

double Rms(const std::vector<float>& x, double sample_rate, double t0, double t1) {
    const std::size_t begin = static_cast<std::size_t>(t0 * sample_rate);
    const std::size_t end = std::min(x.size(), static_cast<std::size_t>(t1 * sample_rate));
    double sum = 0.0;
    for (std::size_t i = begin; i < end; ++i) {
        sum += static_cast<double>(x[i]) * x[i];
    }
    return end > begin ? std::sqrt(sum / static_cast<double>(end - begin)) : 0.0;
}

double AttackRatio(const std::vector<float>& x, double sample_rate, double t0, double t1) {
    const std::size_t window = std::max<std::size_t>(1, static_cast<std::size_t>(0.005 * sample_rate));
    const std::size_t begin = static_cast<std::size_t>(t0 * sample_rate);
    const std::size_t end = std::min(x.size(), static_cast<std::size_t>(t1 * sample_rate));
    double max_rms = 0.0;
    for (std::size_t start = begin; start + window <= end; ++start) {
        double sum = 0.0;
        for (std::size_t i = 0; i < window; ++i) {
            sum += static_cast<double>(x[start + i]) * x[start + i];
        }
        max_rms = std::max(max_rms, std::sqrt(sum / static_cast<double>(window)));
    }
    const double steady = Rms(x, sample_rate, 0.45, 0.9);
    return steady > 0.0 ? max_rms / steady : 0.0;
}

void RunCase(const std::vector<float>& input, float sample_rate, float semitones, const ProbeCase& probe_case,
             bool print_trace) {
    auto shifter = std::make_unique<Shifter>();
    shifter->init(sample_rate);
    shifter->setAlgorithm(Algorithm::PghiDSPark);
    shifter->setPitchShift(semitones);
    shifter->setProbeConfig(probe_case.config);

    std::vector<float> output(input.size(), 0.0f);
    shifter->process(input.data(), output.data(), input.size());

    const auto& events = shifter->getTransientEvents();
    const auto& applied = shifter->getAppliedJumps();
    const auto& applied_bins = shifter->getAppliedBinCounts();
    std::cout << std::format("\n=== {} | semitones={:+.1f}: resets={} ===\n", probe_case.name, semitones,
                             events.size());
    for (std::size_t i = 0; i < events.size(); ++i) {
        const auto& event = events[i];
        const double gap = i == 0 ? 0.0 : event.time_seconds - events[i - 1].time_seconds;
        const float applied_jump = i < applied.size() ? applied[i] : 0.0f;
        const std::size_t applied_bin_count = i < applied_bins.size() ? applied_bins[i] : 0;
        std::cout << std::format(
            "reset #{:<2} frame={:4} t={:8.4f}s gap={:7.4f}s flux={:.6f} median={:.6f} rising={:3} "
            "max_rise={:.4f} rise=[{:7.1f},{:7.1f}]Hz potential_jump={:.3f}rad applied_jump={:.3f}rad "
            "reset_bins={}\n",
            i + 1, event.frame_index, event.time_seconds, gap, event.flux, event.median, event.rising_band_count,
            event.max_band_rise, event.rise_f0, event.rise_f1, event.phase_jump, applied_jump, applied_bin_count);
    }

    if (print_trace) {
        PrintInterestingTrace(*shifter);
    }

    double peak = 0.0;
    std::size_t peak_index = 0;
    for (std::size_t i = 0; i < output.size(); ++i) {
        const double value = std::fabs(static_cast<double>(output[i]));
        if (value > peak) {
            peak = value;
            peak_index = i;
        }
    }
    std::cout << std::format("output peak={:.4f} at t={:.4f}s attack_ratio={:.3f}\n", peak,
                             static_cast<double>(peak_index) / sample_rate,
                             AttackRatio(output, sample_rate, 0.2, 0.4));

    AudioFile<float> file;
    file.setBitDepth(32);
    file.setSampleRate(static_cast<int>(sample_rate));
    file.setNumChannels(1);
    file.setNumSamplesPerChannel(static_cast<int>(output.size()));
    file.samples[0] = std::move(output);
    const auto name = std::format("pitch_shift_transient_{}_{:+.0f}st.wav", probe_case.name, semitones);
    if (!file.save(qwqdsp_support::OutputFile(name))) {
        std::cout << "save failed\n";
    }
}

} // namespace

int main(int argc, char** argv) {
    AudioFile<float> file;
    if (!file.load(qwqdsp_support::InputFile("sine_with_drum.wav"))) {
        std::cout << "cannot load sine_with_drum.wav\n";
        return 1;
    }
    const float sample_rate = static_cast<float>(file.getSampleRate());
    const auto& input = file.samples.front();
    std::cout << std::format("file: {} samples, {} Hz, {} ch\n", input.size(), file.getSampleRate(),
                             file.getNumChannels());

    const float semitones = argc > 1 ? std::strtof(argv[1], nullptr) : 0.0f;

    ProbeCase baseline;
    baseline.name = "baseline_full";
    baseline.config.reset_mode = Shifter::ResetMode::kFull;

    ProbeCase weighted;
    weighted.name = "weighted";
    weighted.config.reset_mode = Shifter::ResetMode::kWeighted;

    ProbeCase no_reset;
    no_reset.name = "no_reset";
    no_reset.config.reset_mode = Shifter::ResetMode::kNone;

    ProbeCase peak_seeded;
    peak_seeded.name = "peak_seeded";
    peak_seeded.config.reset_mode = Shifter::ResetMode::kPeakSeeded;

    ProbeCase peak_05;
    peak_05.name = "peak_05";
    peak_05.config.reset_mode = Shifter::ResetMode::kPeakSeeded;
    peak_05.config.peak_rise_floor_ratio = 0.05f;

    ProbeCase peak_02;
    peak_02.name = "peak_02";
    peak_02.config.reset_mode = Shifter::ResetMode::kPeakSeeded;
    peak_02.config.peak_rise_floor_ratio = 0.02f;

    ProbeCase peak_00;
    peak_00.name = "peak_00";
    peak_00.config.reset_mode = Shifter::ResetMode::kPeakSeeded;
    peak_00.config.peak_rise_floor_ratio = 0.0f;

    ProbeCase guarded_peak_seeded;
    guarded_peak_seeded.name = "guarded_peak_seeded";
    guarded_peak_seeded.config.reset_mode = Shifter::ResetMode::kPeakSeeded;
    guarded_peak_seeded.config.flux_peak_ratio = 0.20f;
    guarded_peak_seeded.config.flux_peak_decay = 0.98f;

    ProbeCase opt_final;
    opt_final.name = "opt_final";
    opt_final.config.reset_mode = Shifter::ResetMode::kPeakSeeded;
    opt_final.config.flux_peak_ratio = 0.20f;
    opt_final.config.flux_peak_decay = 0.98f;
    opt_final.config.peak_rise_floor_ratio = 0.05f;

    ProbeCase guarded_full;
    guarded_full.name = "guarded_full";
    guarded_full.config.reset_mode = Shifter::ResetMode::kFull;
    guarded_full.config.flux_peak_ratio = 0.20f;
    guarded_full.config.flux_peak_decay = 0.98f;

    ProbeCase guarded_weighted;
    guarded_weighted.name = "guarded_weighted";
    guarded_weighted.config.reset_mode = Shifter::ResetMode::kWeighted;
    guarded_weighted.config.flux_peak_ratio = 0.20f;
    guarded_weighted.config.flux_peak_decay = 0.98f;

    ProbeCase low_band_guard;
    low_band_guard.name = "low_band_guard";
    low_band_guard.config.reset_mode = Shifter::ResetMode::kFull;
    low_band_guard.config.require_low_band_rise = true;
    low_band_guard.config.low_band_limit_hz = 2000.0f;

    RunCase(input, sample_rate, semitones, baseline, true);
    RunCase(input, sample_rate, semitones, no_reset, false);
    RunCase(input, sample_rate, semitones, weighted, false);
    RunCase(input, sample_rate, semitones, peak_seeded, false);
    RunCase(input, sample_rate, semitones, peak_05, false);
    RunCase(input, sample_rate, semitones, peak_02, false);
    RunCase(input, sample_rate, semitones, peak_00, false);
    RunCase(input, sample_rate, semitones, guarded_peak_seeded, false);
    RunCase(input, sample_rate, semitones, opt_final, false);
    RunCase(input, sample_rate, semitones, guarded_full, true);
    RunCase(input, sample_rate, semitones, guarded_weighted, false);
    RunCase(input, sample_rate, semitones, low_band_guard, false);
    return 0;
}
