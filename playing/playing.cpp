#include "AudioFile.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <numbers>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

static std::vector<float> MakeSine(size_t count, float cycles_per_sample) {
	std::vector<float> x(count);
	for (size_t i = 0; i < x.size(); ++i) {
		x[i] = std::sin(2.0f * std::numbers::pi_v<float> * cycles_per_sample * static_cast<float>(i));
	}
	return x;
}

class PolyphaseFirAnalysisFilterBank {
public:
	void Init(size_t num_bands, size_t taps_per_phase) {
		if (num_bands == 0 || taps_per_phase == 0) {
			throw std::invalid_argument("num_bands and taps_per_phase must be non-zero");
		}

		num_bands_ = num_bands;
		taps_per_phase_ = taps_per_phase;

		BuildPrototype();
		BuildPolyphaseBranches();
		BuildModulationTable();
	}

	std::vector<std::vector<std::complex<float>>> Process(std::span<const float> x) const {
		assert(num_bands_ != 0);
		assert(taps_per_phase_ != 0);

		size_t const frame_count = x.size() / num_bands_;
		std::vector<std::vector<std::complex<float>>> bands(num_bands_);
		for (auto& band : bands) {
			band.reserve(frame_count);
		}

		std::vector<float> phase_sums(num_bands_);
		for (size_t frame = 0; frame < frame_count; ++frame) {
			size_t const output_time = frame * num_bands_;
			std::fill(phase_sums.begin(), phase_sums.end(), 0.0f);

			for (size_t phase = 0; phase < num_bands_; ++phase) {
				float sum = 0.0f;
				auto const coeffs = PhaseCoeffs(phase);
				for (size_t tap = 0; tap < taps_per_phase_; ++tap) {
					size_t const delay = phase + tap * num_bands_;
					if (delay <= output_time) {
						sum += coeffs[tap] * x[output_time - delay];
					}
				}
				phase_sums[phase] = sum;
			}

			for (size_t band = 0; band < num_bands_; ++band) {
				std::complex<float> y{};
				auto const modulation = ModulationCoeffs(band);
				for (size_t phase = 0; phase < num_bands_; ++phase) {
					y += phase_sums[phase] * modulation[phase];
				}
				bands[band].push_back(y);
			}
		}

		return bands;
	}

	size_t NumBands() const noexcept {
		return num_bands_;
	}

private:
	void BuildPrototype() {
		size_t const filter_len = num_bands_ * taps_per_phase_;
		prototype_.resize(filter_len);

		float const cutoff = std::numbers::pi_v<float> / static_cast<float>(num_bands_);
		float const center = (static_cast<float>(filter_len) - 1.0f) * 0.5f;
		for (size_t i = 0; i < filter_len; ++i) {
			float const t = static_cast<float>(i) - center;
			float sinc = cutoff / std::numbers::pi_v<float>;
			if (t != 0.0f) {
				sinc = std::sin(cutoff * t) / (std::numbers::pi_v<float> * t);
			}

			float const window = 0.5f - 0.5f * std::cos(
				2.0f * std::numbers::pi_v<float> * static_cast<float>(i) / static_cast<float>(filter_len - 1));
			prototype_[i] = sinc * window;
		}

		float const gain = std::accumulate(prototype_.begin(), prototype_.end(), 0.0f);
		if (gain != 0.0f) {
			for (auto& coeff : prototype_) {
				coeff /= gain;
			}
		}
	}

	void BuildPolyphaseBranches() {
		phases_.assign(num_bands_ * taps_per_phase_, 0.0f);
		for (size_t phase = 0; phase < num_bands_; ++phase) {
			for (size_t tap = 0; tap < taps_per_phase_; ++tap) {
				phases_[phase * taps_per_phase_ + tap] = prototype_[phase + tap * num_bands_];
			}
		}
	}

	void BuildModulationTable() {
		modulation_.assign(num_bands_ * num_bands_, {});
		for (size_t band = 0; band < num_bands_; ++band) {
			for (size_t phase = 0; phase < num_bands_; ++phase) {
				float const angle = -2.0f * std::numbers::pi_v<float> *
									static_cast<float>(band * phase) / static_cast<float>(num_bands_);
				modulation_[band * num_bands_ + phase] = std::complex<float>(std::cos(angle), std::sin(angle));
			}
		}
	}

	std::span<const float> PhaseCoeffs(size_t phase) const noexcept {
		return {phases_.data() + phase * taps_per_phase_, taps_per_phase_};
	}

	std::span<const std::complex<float>> ModulationCoeffs(size_t band) const noexcept {
		return {modulation_.data() + band * num_bands_, num_bands_};
	}

	size_t num_bands_{};
	size_t taps_per_phase_{};
	std::vector<float> prototype_;
	std::vector<float> phases_;
	std::vector<std::complex<float>> modulation_;
};

static float RmsMagnitude(std::span<const std::complex<float>> x) {
	if (x.empty()) {
		return 0.0f;
	}

	float sum = 0.0f;
	for (auto const& sample : x) {
		sum += std::norm(sample);
	}
	return std::sqrt(sum / static_cast<float>(x.size()));
}

static void RunFilterBankSelfTest() {
	PolyphaseFirAnalysisFilterBank bank;
	bank.Init(8, 8);
	auto bands = bank.Process(MakeSine(1024, 0.02f));
	assert(bands.size() == 8);
	assert(!bands[0].empty());
	for (auto const& band : bands) {
		assert(band.size() == 128);
	}
}

static void SaveBandAsAudio(
	std::span<const std::complex<float>> band,
	std::string const& base_path,
	size_t band_idx,
	int sample_rate)
{
	std::string out_path = base_path;
	size_t dot_pos = out_path.rfind('.');
	if (dot_pos != std::string::npos) {
		out_path = out_path.substr(0, dot_pos);
	}
	out_path += "_band" + std::to_string(band_idx) + ".wav";

	std::vector<float> real_samples(band.size());
	std::ranges::transform(band, real_samples.begin(),
		[](std::complex<float> const& s) { return s.real(); });

	AudioFile<float> out_file;
	out_file.setNumChannels(1);
	out_file.setNumSamplesPerChannel(static_cast<int>(real_samples.size()));
	out_file.setSampleRate(sample_rate);
	out_file.setBitDepth(32);
	out_file.samples[0] = std::move(real_samples);
	out_file.save(out_path);

	std::cout << "  saved: " << out_path << '\n';
}

int main() {
	RunFilterBankSelfTest();

	AudioFile<float> file;
	std::string const path = R"(C:\Users\Kawai\Music\modu.wav)";
	std::vector<float> input;
	int sample_rate = 48000;

	if (file.load(path) && !file.samples.empty() && !file.samples.front().empty()) {
		input = file.samples.front();
		sample_rate = file.getSampleRate();
		std::cout << "loaded: " << path << '\n';
	}
	else {
		input = MakeSine(sample_rate, 0.02f);
		std::cout << "audio file unavailable, using generated sine fallback\n";
	}

	constexpr size_t kNumBands = 8;
	constexpr size_t kTapsPerPhase = 12;

	PolyphaseFirAnalysisFilterBank bank;
	bank.Init(kNumBands, kTapsPerPhase);
	auto bands = bank.Process(input);

	int const decimated_fs = sample_rate / static_cast<int>(kNumBands);

	std::cout << "source sample rate: " << sample_rate << " Hz\n";
	std::cout << "source samples: " << input.size() << '\n';
	std::cout << "bands: " << bank.NumBands() << '\n';
	std::cout << "decimated sample rate: " << decimated_fs << " Hz\n";
	std::cout << "saving bands...\n";

	for (size_t band = 0; band < bands.size(); ++band) {
		std::cout << "band " << band
				  << ": frames=" << bands[band].size()
				  << ", rms=" << RmsMagnitude(bands[band]) << '\n';
		SaveBandAsAudio(bands[band], path, band, decimated_fs);
	}

	return 0;
}
