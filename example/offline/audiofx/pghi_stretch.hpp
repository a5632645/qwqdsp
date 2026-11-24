#pragma once

#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <numbers>
#include <queue>
#include <qwqdsp/fx/resample.hpp>
#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/window/hann.hpp>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace qwqdsp_example {

struct PghiSettings {
    size_t window_size = 2048;
    bool transients = true;
    bool chirp = false;
};

struct PghiResult {
    std::vector<float> samples;
    size_t onset_count = 0;
    size_t anchor_count = 0;
    size_t chirp_frames = 0;
};

namespace detail {

static constexpr double kTwoPi = 2 * std::numbers::pi;
static constexpr float kUnassigned = std::numeric_limits<float>::quiet_NaN();

static inline float WrapPhase(double phase) {
    return static_cast<float>(std::remainder(phase, kTwoPi));
}

class CenteredFft {
public:
    explicit CenteredFft(size_t size) : frame_(size), spectrum_(size + 2) {
        fft_.Init(size);
    }

    void Analyze(std::span<const float> samples, size_t channels, size_t channel,
                 int64_t center, std::span<const float> window) {
        std::fill(frame_.begin(), frame_.end(), 0.0f);
        const auto half = static_cast<int64_t>(window.size() / 2);
        const auto frames = static_cast<int64_t>(samples.size() / channels);
        const auto start = center - half;
        const auto first = std::max<int64_t>(0, -start);
        const auto end = std::min<int64_t>(static_cast<int64_t>(window.size()), frames - start);
        for (auto i = first; i < end; ++i) {
            const auto offset = static_cast<size_t>(i - half) & (frame_.size() - 1);
            frame_[offset] = window[static_cast<size_t>(i)] * samples[static_cast<size_t>(start + i) * channels + channel];
        }
        fft_.FFT(frame_.data(), spectrum_.data());
    }

    std::complex<double> Bin(size_t bin) const {
        return {spectrum_[2 * bin], spectrum_[2 * bin + 1]};
    }

    std::span<float> Spectrum() { return spectrum_; }

    std::span<const float> Inverse() {
        fft_.IFFT(spectrum_.data(), frame_.data());
        return frame_;
    }

private:
    qwqdsp_spectral::RealFFT fft_;
    std::vector<float> frame_;
    std::vector<float> spectrum_;
};

static inline void MeasureMagnitude(CenteredFft& fft, std::span<const float> input, size_t channels,
                                    int64_t center, std::span<const float> window, std::span<float> magnitude) {
    std::fill(magnitude.begin(), magnitude.end(), 0.0f);
    for (size_t channel = 0; channel < channels; ++channel) {
        fft.Analyze(input, channels, channel, center, window);
        const auto spectrum = fft.Spectrum();
        for (size_t k = 0; k < magnitude.size(); ++k) {
            magnitude[k] = std::max(magnitude[k], std::hypot(spectrum[2 * k], spectrum[2 * k + 1]));
        }
    }
}

static inline std::vector<int64_t> DetectOnsets(std::span<const float> input, size_t channels, unsigned sample_rate,
                                                std::span<const float> window) {
    const size_t fft_size = std::bit_ceil(window.size());
    const size_t bins = fft_size / 2 + 1;
    const auto half = static_cast<int64_t>(window.size() / 2);
    const auto hop = half / 2;
    const auto length = static_cast<int64_t>(input.size() / channels);
    CenteredFft fft{fft_size};
    std::vector<float> current(bins), previous(bins);
    float peak_magnitude = 0;
    for (int64_t center = 0; center < length + half; center += hop) {
        MeasureMagnitude(fft, input, channels, center, window, current);
        peak_magnitude = std::max(peak_magnitude, *std::max_element(current.begin(), current.end()));
    }
    if (peak_magnitude == 0) {
        return {};
    }
    const float threshold = 0.001f * peak_magnitude;
    std::vector<double> flux;
    for (int64_t center = 0; center < length + half; center += hop) {
        MeasureMagnitude(fft, input, channels, center, window, current);
        size_t occupied = 0;
        double sum = 0;
        for (size_t k = 1; k + 1 < bins; ++k) {
            const float rise = current[k] - std::max({previous[k - 1], previous[k], previous[k + 1]});
            if (rise > threshold) {
                ++occupied;
                sum += rise;
            }
        }
        flux.push_back(8 * occupied >= bins ? sum : 0);
        previous.swap(current);
    }

    const double minimum = 0.08 * *std::max_element(flux.begin(), flux.end());
    const auto cooldown = static_cast<int64_t>(std::lround(0.012 * sample_rate));
    std::vector<size_t> peaks;
    for (size_t i = 0; i < flux.size(); ++i) {
        if (flux[i] <= minimum || (i > 0 && flux[i] <= flux[i - 1])
            || (i + 1 < flux.size() && flux[i] < flux[i + 1])) {
            continue;
        }
        if (peaks.empty() || static_cast<int64_t>(i - peaks.back()) * hop >= cooldown) {
            peaks.push_back(i);
        } else if (flux[i] > flux[peaks.back()]) {
            peaks.back() = i;
        }
    }

    const auto attack = [&](int64_t frame) {
        float strongest = 0;
        const auto base = static_cast<size_t>(frame) * channels;
        for (size_t channel = 0; channel < channels; ++channel) {
            const float prior = frame > 0 ? input[base + channel - channels] : 0;
            strongest = std::max(strongest, std::abs(input[base + channel] - prior));
        }
        return strongest;
    };
    std::vector<int64_t> onsets;
    for (const size_t peak : peaks) {
        const auto center = static_cast<int64_t>(peak) * hop;
        const auto first = std::max<int64_t>(0, center - half);
        const auto end = std::min(length, center + half + 1);
        float strongest = 0;
        for (auto frame = first; frame < end; ++frame) {
            strongest = std::max(strongest, attack(frame));
        }
        for (auto frame = first; strongest > 0 && frame < end; ++frame) {
            if (attack(frame) >= 0.7f * strongest) {
                if (onsets.empty() || frame > onsets.back()) {
                    onsets.push_back(frame);
                }
                break;
            }
        }
    }
    return onsets;
}

struct ChirpCorrection {
    double center_bin = 0;
    double phase_offset = 0;
    double phase_curve = 0;
    double log_gain = 0;
    double gain_curve = 0;

    double Phase(size_t bin) const {
        const double offset = static_cast<double>(bin) - center_bin;
        return phase_offset + offset * offset * phase_curve;
    }

    float Gain(size_t bin) const {
        const double offset = static_cast<double>(bin) - center_bin;
        return static_cast<float>(std::exp(std::clamp(log_gain + offset * offset * gain_curve, -12.0, std::log(4.0))));
    }
};

struct ChirpFit {
    double center_bin = 0;
    double beta = 0;
    double confidence = 0;

    ChirpCorrection ForStretch(double scale, double stretch) const {
        const double beta_target = beta / stretch;
        const double a = 1 + beta * beta;
        const double b = 1 + beta_target * beta_target;
        return {
            center_bin,
            0.5 * confidence * (std::atan(beta_target) - stretch * std::atan(beta)),
            0.5 * confidence * scale * scale * (stretch * beta / a - beta_target / b),
            0.25 * confidence * std::log(a / b),
            0.5 * confidence * scale * scale * (1 / a - 1 / b),
        };
    }
};

static inline ChirpFit FitChirp(CenteredFft& fft, std::span<const float> input,
                                int64_t center, std::span<const float> magnitude,
                                std::span<const float> moment1, std::span<const float> moment2,
                                double scale) {
    const size_t peak = static_cast<size_t>(std::max_element(magnitude.begin(), magnitude.end()) - magnitude.begin());
    if (peak < 3 || peak + 3 >= magnitude.size() || magnitude[peak] <= 0) {
        return {};
    }
    // Fit only one isolated, dominant ridge in this frame.
    for (size_t k = 1; k + 1 < magnitude.size(); ++k) {
        if (k != peak && (k + 1 < peak || k > peak + 1) && magnitude[k] >= 0.12f * magnitude[peak]
            && magnitude[k] > magnitude[k - 1] && magnitude[k] >= magnitude[k + 1]) {
            return {};
        }
    }
    const auto c0 = fft.Bin(peak);
    fft.Analyze(input, 1, 0, center, moment1);
    const auto mu = fft.Bin(peak) / c0;
    fft.Analyze(input, 1, 0, center, moment2);
    const auto h = mu * mu - fft.Bin(peak) / c0;
    if (std::abs(h) < 1e-6) {
        return {};
    }
    const auto inverse = 1.0 / h;
    const auto detuning = std::complex<double>{0, -1} * mu * inverse;
    const double beta = inverse.imag();
    const double error = std::max(std::abs(inverse.real() + 1), std::abs(detuning.imag()));
    if (error >= 0.15 || std::abs(beta) > 12) {
        return {};
    }
    const double confidence = std::clamp((0.15 - error) / 0.10, 0.0, 1.0);
    return {static_cast<double>(peak) - detuning.real() / scale, beta, confidence};
}

struct PghiBin {
    float magnitude = 0;
    float time = 0;       // Phase advance minus the bin carrier, modulo 2π.
    float frequency = 0;  // Phase advance to the next bin, modulo 2π.
};

using NodeIndex = uint32_t;

struct Anchor {
    NodeIndex index;
    int64_t onset;
};

struct PhaseLink {
    NodeIndex parent;
    NodeIndex child;
    float phase_step;
};

struct PghiPaths {
    std::vector<NodeIndex> roots;
    std::vector<PhaseLink> links;
};

class PghiGrid {
public:
    PghiGrid(size_t bins, size_t rows) : bins_(bins) {
        if (rows > std::min<size_t>(std::numeric_limits<NodeIndex>::max(), nodes_.max_size()) / bins) {
            throw std::length_error("Audio is too long for the offline STFT grid.");
        }
        nodes_.resize(rows * bins);
    }

    std::span<PghiBin> Row(size_t row) { return {nodes_.data() + row * bins_, bins_}; }

    PghiPaths BuildPaths(std::span<const Anchor> anchors) const {
        float peak = 0;
        for (const auto& node : nodes_) {
            peak = std::max(peak, node.magnitude);
        }
        const float threshold = peak * 1e-6f;
        enum class State : uint8_t { Unvisited, Searching, Linked };
        std::vector<State> state(nodes_.size(), State::Unvisited);
        using Entry = std::pair<float, NodeIndex>;
        std::priority_queue<Entry> heap;
        PghiPaths paths;
        for (const auto& anchor : anchors) {
            state[anchor.index] = State::Linked;
            heap.emplace(nodes_[anchor.index].magnitude, anchor.index);
        }
        const auto propagate = [&] {
            while (!heap.empty()) {
                const auto index = heap.top().second;
                heap.pop();
                VisitNeighbors(index, [&](NodeIndex next, float edge) {
                    if (nodes_[next].magnitude > threshold && state[next] != State::Linked) {
                        state[next] = State::Linked;
                        paths.links.push_back({index, next, edge});
                        heap.emplace(nodes_[next].magnitude, next);
                    }
                });
            }
        };
        propagate();

        // Find each detached component's strongest bin without sorting the whole grid.
        std::deque<NodeIndex> pending;
        for (NodeIndex i = 0; i < nodes_.size(); ++i) {
            if (nodes_[i].magnitude <= threshold || state[i] != State::Unvisited) {
                continue;
            }
            NodeIndex seed = i;
            pending.push_back(i);
            state[i] = State::Searching;
            while (!pending.empty()) {
                const NodeIndex index = pending.front();
                pending.pop_front();
                if (nodes_[index].magnitude > nodes_[seed].magnitude) {
                    seed = index;
                }
                VisitNeighbors(index, [&](NodeIndex next, float) {
                    if (nodes_[next].magnitude > threshold && state[next] == State::Unvisited) {
                        state[next] = State::Searching;
                        pending.push_back(next);
                    }
                });
            }
            state[seed] = State::Linked;
            paths.roots.push_back(seed);
            heap.emplace(nodes_[seed].magnitude, seed);
            propagate();
        }
        return paths;
    }

private:
    template <class Visitor>
    void VisitNeighbors(size_t index, Visitor visit) const {
        const size_t bin = index % bins_;
        if (bin > 0) {
            visit(static_cast<NodeIndex>(index - 1), -nodes_[index - 1].frequency);
        }
        if (bin + 1 < bins_) {
            visit(static_cast<NodeIndex>(index + 1), nodes_[index].frequency);
        }
        if (index >= bins_) {
            visit(static_cast<NodeIndex>(index - bins_), -nodes_[index - bins_].time);
        }
        if (index + bins_ < nodes_.size()) {
            visit(static_cast<NodeIndex>(index + bins_), nodes_[index].time);
        }
    }

    size_t bins_;
    std::vector<PghiBin> nodes_;
};

static inline std::vector<Anchor> PickAnchors(PghiGrid& grid, CenteredFft& fft, std::span<const float> input,
                                               size_t channels, std::span<const float> window,
                                               std::span<const int64_t> onsets, int64_t analysis_hop, size_t bins) {
    std::vector<float> before(bins), event(bins), best(bins);
    std::vector<size_t> selected(bins);
    size_t last_row = std::numeric_limits<size_t>::max();
    std::vector<Anchor> anchors;
    for (const auto onset : onsets) {
        const auto row = static_cast<size_t>((onset + analysis_hop / 2) / analysis_hop);
        if (row != last_row) {
            std::fill(best.begin(), best.end(), 0.0f);
            last_row = row;
        }
        MeasureMagnitude(fft, input, channels, onset - static_cast<int64_t>(window.size()), window, before);
        MeasureMagnitude(fft, input, channels, onset, window, event);
        const auto nodes = grid.Row(row);
        float peak = 0;
        for (size_t k = 0; k < bins; ++k) {
            if (event[k] <= 2 * before[k] || event[k] < 0.5f * nodes[k].magnitude || nodes[k].magnitude <= 0) {
                event[k] = 0;
            }
            peak = std::max(peak, event[k]);
        }
        for (size_t k = 0; k < bins; ++k) {
            if (event[k] > 0.15f * peak && event[k] > best[k]) {
                if (best[k] == 0) {
                    selected[k] = anchors.size();
                    anchors.push_back({static_cast<NodeIndex>(row * bins + k), onset});
                } else {
                    anchors[selected[k]].onset = onset;
                }
                best[k] = event[k];
            }
        }
    }
    return anchors;
}

} // namespace detail

class PghiAnalysis {
public:
    PghiAnalysis(std::vector<float> input, size_t channels, unsigned sample_rate, const PghiSettings& settings)
        : input_(std::move(input)), channels_(channels) {
        if (channels == 0 || input_.empty() || input_.size() % channels != 0 || sample_rate == 0
            || !std::has_single_bit(settings.window_size) || settings.window_size < 256 || settings.window_size > 8192
            || (settings.chirp && channels != 1)) {
            throw std::invalid_argument("Invalid PGHI settings or input format.");
        }
        const size_t length = settings.window_size;
        fft_size_ = 2 * length;
        bins_ = fft_size_ / 2 + 1;
        half_window_ = static_cast<int64_t>(length / 2);
        bin_width_ = detail::kTwoPi / fft_size_;
        sigma_ = length / 14.0;
        // The maximum internal stretch is 8x, including the pitch compensation.
        analysis_hop_ = static_cast<int64_t>(length / (settings.chirp ? 32 : 16));
        const auto last_input = static_cast<int64_t>(input_.size() / channels - 1);
        rows_ = static_cast<size_t>(1 + last_input / analysis_hop_ + (last_input % analysis_hop_ != 0));
        window_.resize(length);
        Analyze(sample_rate, settings);
    }

    size_t FrameCount() const { return rows_; }
    int64_t AnalysisHop() const { return analysis_hop_; }
    bool HasCachedSpectra() const { return !spectra_.empty(); }

    PghiResult Render(double duration_ratio, double semitones) const {
        if (!std::isfinite(duration_ratio) || duration_ratio < 0.25 || duration_ratio > 4) {
            throw std::invalid_argument("Time stretch must be between 0.25 and 4.");
        }
        const float pitch_step = PitchStep(semitones);
        const double internal_stretch = duration_ratio * pitch_step;
        const size_t wanted = FramesFor(duration_ratio);
        PghiResult result;
        if (internal_stretch == 1) {
            result.samples = pitch_step == 1 ? input_ : ResamplePitch(input_, pitch_step, wanted);
        } else {
            result.samples = RenderWola(internal_stretch, FramesFor(internal_stretch));
            result.onset_count = onset_count_;
            result.anchor_count = anchors_.size();
            result.chirp_frames = chirp_frames_;
            if (pitch_step != 1) {
                result.samples = ResamplePitch(result.samples, pitch_step, wanted);
            }
        }
        return result;
    }

    PghiResult RenderResampled(double semitones) const {
        const float step = PitchStep(semitones);
        return {step == 1 ? input_ : ResamplePitch(input_, step, FramesFor(1.0 / step))};
    }

private:
    static float PitchStep(double semitones) {
        if (!std::isfinite(semitones) || semitones < -12 || semitones > 12) {
            throw std::invalid_argument("Pitch shift must be between -12 and +12 semitones.");
        }
        return static_cast<float>(std::exp2(semitones / 12));
    }

    size_t FramesFor(double stretch) const {
        return static_cast<size_t>(std::max<int64_t>(1, std::llround((input_.size() / channels_) * stretch)));
    }

    size_t SampleCount(size_t frames) const {
        if (frames > std::vector<float>().max_size() / channels_) {
            throw std::length_error("Output audio is too long.");
        }
        return frames * channels_;
    }

    std::vector<float> ResamplePitch(std::span<const float> interleaved, float step, size_t frames) const {
        const size_t input_frames = interleaved.size() / channels_;
        const size_t needed = static_cast<size_t>(std::floor((frames - 1) * static_cast<double>(step))) + 1;
        const size_t padded = std::max(input_frames, needed);
        if (padded > static_cast<size_t>(std::numeric_limits<int>::max()) - 127) {
            throw std::length_error("Audio is too long for the pitch resampler.");
        }
        qwqdsp_fx::Resample resampler;
        resampler.Init(step, 1, 100, 127, 128);
        std::vector<float> input(padded), output(SampleCount(frames));
        for (size_t channel = 0; channel < channels_; ++channel) {
            for (size_t i = 0; i < input_frames; ++i) {
                input[i] = interleaved[i * channels_ + channel];
            }
            const auto resampled = resampler.Process(input);
            if (resampled.size() < frames) {
                throw std::runtime_error("The resampler produced fewer audio frames than expected.");
            }
            for (size_t i = 0; i < frames; ++i) {
                output[i * channels_ + channel] = resampled[i];
            }
        }
        return output;
    }

    size_t SpectrumOffset(size_t row, size_t channel) const {
        return (row * channels_ + channel) * (fft_size_ + 2);
    }

    void Analyze(unsigned sample_rate, const PghiSettings& settings) {
        using namespace detail;
        const size_t length = window_.size();
        const auto input_frames = static_cast<int64_t>(input_.size() / channels_);
        std::vector<float> moment1, moment2;
        if (settings.chirp) {
            moment1.resize(length);
            moment2.resize(length);
            chirp_.resize(rows_);
            for (size_t i = 0; i < length; ++i) {
                const double u = (static_cast<double>(i) - half_window_) / sigma_;
                window_[i] = static_cast<float>(std::exp(-0.5 * u * u));
                moment1[i] = static_cast<float>(u * window_[i]);
                moment2[i] = static_cast<float>(u * u * window_[i]);
            }
        } else {
            qwqdsp_window::Hann::Window(window_, true);
        }

        const size_t short_size = std::min(length, 2 * std::max<size_t>(16, std::lround(0.0027 * sample_rate)));
        std::vector<float> short_window(short_size);
        qwqdsp_window::Hann::Window(short_window, true);
        const auto onsets = settings.transients ? DetectOnsets(input_, channels_, sample_rate, short_window)
                                                 : std::vector<int64_t>{};
        onset_count_ = onsets.size();

        PghiGrid grid{bins_, rows_};
        CenteredFft fft{fft_size_};
        constexpr size_t kSpectrumBudget = 128 * 1024 * 1024;
        const size_t row_values = (fft_size_ + 2) * channels_;
        if (row_values <= kSpectrumBudget / sizeof(float) / rows_) {
            spectra_.resize(rows_ * row_values);
        }
        std::vector<float> previous_phase(channels_ * bins_), previous_magnitude(channels_ * bins_);
        std::vector<float> time_support(bins_), frequency_support(bins_), frame_magnitude(bins_);
        for (size_t row = 0; row < rows_; ++row) {
            const auto center = static_cast<int64_t>(row) * analysis_hop_;
            const auto nodes = grid.Row(row);
            const auto previous_nodes = row ? grid.Row(row - 1) : std::span<PghiBin>{};
            std::fill(time_support.begin(), time_support.end(), -1.0f);
            std::fill(frequency_support.begin(), frequency_support.end(), -1.0f);
            for (size_t channel = 0; channel < channels_; ++channel) {
                fft.Analyze(input_, channels_, channel, center, window_);
                const auto spectrum = fft.Spectrum();
                if (!spectra_.empty()) {
                    std::copy(spectrum.begin(), spectrum.end(), spectra_.data() + SpectrumOffset(row, channel));
                }
                float preceding_phase = 0;
                float preceding_magnitude = 0;
                for (size_t k = 0; k < bins_; ++k) {
                    const float re = spectrum[2 * k];
                    const float im = spectrum[2 * k + 1];
                    const float magnitude = std::hypot(re, im);
                    const float phase = std::atan2(im, re);
                    const size_t index = channel * bins_ + k;
                    nodes[k].magnitude = std::max(nodes[k].magnitude, magnitude);
                    const float time_magnitude = std::min(magnitude, previous_magnitude[index]);
                    if (row && time_magnitude > time_support[k]) {
                        const double carrier = bin_width_ * analysis_hop_ * k;
                        previous_nodes[k].time = WrapPhase(phase - previous_phase[index] - carrier);
                        time_support[k] = time_magnitude;
                    }
                    const float frequency_magnitude = std::min(magnitude, preceding_magnitude);
                    if (k && frequency_magnitude > frequency_support[k]) {
                        nodes[k - 1].frequency = WrapPhase(phase - preceding_phase);
                        frequency_support[k] = frequency_magnitude;
                    }
                    previous_phase[index] = phase;
                    previous_magnitude[index] = magnitude;
                    preceding_phase = phase;
                    preceding_magnitude = magnitude;
                    frame_magnitude[k] = magnitude;
                }
            }
            if (settings.chirp && center >= half_window_ && center + half_window_ <= input_frames) {
                const auto onset = std::lower_bound(onsets.begin(), onsets.end(), center - half_window_);
                if (onset == onsets.end() || *onset > center + half_window_) {
                    chirp_[row] = FitChirp(fft, input_, center, frame_magnitude, moment1, moment2, sigma_ * bin_width_);
                    chirp_frames_ += chirp_[row].confidence > 0;
                }
            }
        }
        anchors_ = PickAnchors(grid, fft, input_, channels_, short_window, onsets, analysis_hop_, bins_);
        paths_ = grid.BuildPaths(anchors_);
    }

    std::vector<float> Integrate(double stretch, std::span<const int64_t> centers,
                                  std::span<const detail::ChirpCorrection> chirp) const {
        using namespace detail;
        std::vector<float> phase(rows_ * bins_, kUnassigned);
        std::vector<double> time_scale(rows_ - 1);
        for (size_t row = 0; row + 1 < rows_; ++row) {
            time_scale[row] = static_cast<double>(centers[row + 1] - centers[row]) / analysis_hop_ - 1;
        }
        for (const auto& anchor : anchors_) {
            const size_t row = anchor.index / bins_;
            const size_t bin = anchor.index % bins_;
            const auto displacement = std::llround(stretch * anchor.onset) - centers[row]
                                      - (anchor.onset - static_cast<int64_t>(row) * analysis_hop_);
            const double potential = chirp.empty() ? 0 : chirp[row].Phase(bin);
            phase[anchor.index] = WrapPhase(-bin_width_ * bin * displacement - potential);
        }
        for (const auto seed : paths_.roots) {
            phase[seed] = chirp.empty() ? 0 : WrapPhase(-chirp[seed / bins_].Phase(seed % bins_));
        }
        for (const auto& link : paths_.links) {
            const auto offset = static_cast<int64_t>(link.child) - link.parent;
            double step = (stretch - 1) * link.phase_step;
            if (std::abs(offset) == static_cast<int64_t>(bins_)) {
                const size_t row = std::min(link.parent, link.child) / bins_;
                const double carrier = bin_width_ * analysis_hop_ * (link.child % bins_);
                step = time_scale[row] * (link.phase_step + (offset > 0 ? carrier : -carrier));
            }
            phase[link.child] = WrapPhase(phase[link.parent] + step);
        }
        return phase;
    }

    std::vector<float> RenderWola(double stretch, size_t output_frames) const {
        using namespace detail;
        std::vector<int64_t> centers(rows_);
        for (size_t row = 0; row < rows_; ++row) {
            centers[row] = std::llround(stretch * static_cast<double>(row) * analysis_hop_);
        }
        std::vector<ChirpCorrection> chirp;
        chirp.reserve(chirp_.size());
        for (const auto& fit : chirp_) {
            chirp.push_back(fit.ForStretch(sigma_ * bin_width_, stretch));
        }
        const auto phase = Integrate(stretch, centers, chirp);
        std::vector<float> samples(SampleCount(output_frames)), normalization(output_frames);
        std::vector<std::complex<float>> rotation(bins_);
        CenteredFft fft{fft_size_};
        for (size_t row = 0; row < rows_; ++row) {
            const auto start = centers[row] - half_window_;
            const auto first = std::max<int64_t>(0, -start);
            const auto end = std::min<int64_t>(window_.size(), static_cast<int64_t>(output_frames) - start);
            for (size_t k = 0; k < bins_; ++k) {
                const float correction = phase[row * bins_ + k];
                const bool applied = std::isfinite(correction);
                const double potential = !chirp.empty() && applied ? chirp[row].Phase(k) : 0;
                const float gain = !chirp.empty() && applied ? chirp[row].Gain(k) : 1;
                rotation[k] = std::polar(gain, applied ? WrapPhase(correction + potential) : 0.0f);
            }
            for (auto i = first; i < end; ++i) {
                normalization[static_cast<size_t>(start + i)] += window_[i] * window_[i];
            }
            for (size_t channel = 0; channel < channels_; ++channel) {
                const auto spectrum = fft.Spectrum();
                if (spectra_.empty()) {
                    fft.Analyze(input_, channels_, channel, static_cast<int64_t>(row) * analysis_hop_, window_);
                } else {
                    std::copy_n(spectra_.data() + SpectrumOffset(row, channel), spectrum.size(), spectrum.data());
                }
                for (size_t k = 0; k < bins_; ++k) {
                    const auto coefficient = std::complex<float>{spectrum[2 * k], spectrum[2 * k + 1]} * rotation[k];
                    spectrum[2 * k] = coefficient.real();
                    spectrum[2 * k + 1] = coefficient.imag();
                }
                spectrum[1] = 0;
                spectrum.back() = 0;
                const auto frame = fft.Inverse();
                for (auto i = first; i < end; ++i) {
                    const auto offset = static_cast<size_t>(i - half_window_) & (fft_size_ - 1);
                    samples[static_cast<size_t>(start + i) * channels_ + channel] += window_[i] * frame[offset];
                }
            }
        }
        for (size_t frame = 0; frame < output_frames; ++frame) {
            if (normalization[frame] > 0) {
                for (size_t channel = 0; channel < channels_; ++channel) {
                    samples[frame * channels_ + channel] /= normalization[frame];
                }
            }
        }
        return samples;
    }

    std::vector<float> input_;
    size_t channels_;
    size_t fft_size_ = 0;
    size_t bins_ = 0;
    size_t rows_ = 0;
    int64_t half_window_ = 0;
    int64_t analysis_hop_ = 0;
    double bin_width_ = 0;
    double sigma_ = 0;
    size_t onset_count_ = 0;
    size_t chirp_frames_ = 0;
    std::vector<float> window_;
    std::vector<float> spectra_;
    std::vector<detail::ChirpFit> chirp_;
    std::vector<detail::Anchor> anchors_;
    detail::PghiPaths paths_;
};

} // namespace qwqdsp_example
