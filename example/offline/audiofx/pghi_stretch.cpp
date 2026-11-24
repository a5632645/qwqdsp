#include <algorithm>
#include <array>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <format>
#include <iostream>
#include <limits>
#include <miniaudio.h>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif

#include "pghi_stretch.hpp"

namespace qwqdsp_example::detail {

namespace fs = std::filesystem;

struct Audio {
    unsigned sample_rate = 0;
    unsigned channels = 0;
    std::vector<float> samples;
};

template <class String>
void TrimLine(String& line) {
    using Char = typename String::value_type;
    const String spaces{Char(' '), Char('\t'), Char('\r'), Char('\n')};
    const auto first = line.find_first_not_of(spaces);
    if (first == String::npos) {
        line.clear();
        return;
    }
    line = line.substr(first, line.find_last_not_of(spaces) - first + 1);
}

std::string ReadLine(std::string_view prompt) {
    std::cout << prompt << std::flush;
    std::string line;
    if (!std::getline(std::cin, line)) {
        throw std::runtime_error("End of input.");
    }
    TrimLine(line);
    return line;
}

template <class String>
void Unquote(String& line) {
    if (line.size() >= 2 && line.front() == '"' && line.back() == '"') {
        line = line.substr(1, line.size() - 2);
    }
}

fs::path ReadPath() {
    constexpr std::string_view kPrompt = "Input audio (WAV, FLAC or MP3): ";
#ifdef _WIN32
    const HANDLE input = GetStdHandle(STD_INPUT_HANDLE);
    DWORD mode = 0;
    if (GetConsoleMode(input, &mode)) {
        std::cout << kPrompt << std::flush;
        std::wstring line;
        std::array<wchar_t, 512> buffer{};
        do {
            DWORD count = 0;
            if (!ReadConsoleW(input, buffer.data(), static_cast<DWORD>(buffer.size()), &count, nullptr) || !count) {
                throw std::runtime_error("Could not read the input path.");
            }
            line.append(buffer.data(), count);
        } while (line.back() != L'\n');
        TrimLine(line);
        Unquote(line);
        return fs::path(line);
    }
#endif
    auto line = ReadLine(kPrompt);
    Unquote(line);
    return fs::path(std::u8string(line.begin(), line.end()));
}

template <class T>
T ParseNumber(std::string_view line, T default_value) {
    if (line.empty()) {
        return default_value;
    }
    if (line.starts_with('+') && line.size() > 1 && line[1] != '-') {
        line.remove_prefix(1);
    }
    T result{};
    const auto [end, error] = std::from_chars(line.data(), line.data() + line.size(), result);
    if (error != std::errc{} || end != line.data() + line.size()) {
        throw std::invalid_argument("Expected a number.");
    }
    return result;
}

template <class T>
T ReadNumber(std::string_view prompt, T default_value) {
    return ParseNumber(ReadLine(prompt), default_value);
}

bool ReadSwitch(std::string_view prompt, bool default_value) {
    const auto line = ReadLine(prompt);
    if (line.empty()) {
        return default_value;
    }
    if (line == "y" || line == "Y") {
        return true;
    }
    if (line == "n" || line == "N") {
        return false;
    }
    throw std::invalid_argument("Expected y or n.");
}

class Decoder {
public:
    explicit Decoder(const fs::path& path) {
        const auto config = ma_decoder_config_init(ma_format_f32, 0, 0);
#ifdef _WIN32
        const auto result = ma_decoder_init_file_w(path.c_str(), &config, &decoder_);
#else
        const auto result = ma_decoder_init_file(path.c_str(), &config, &decoder_);
#endif
        if (result != MA_SUCCESS) {
            throw std::runtime_error(std::format("Could not open input audio ({})", ma_result_description(result)));
        }
    }

    Decoder(const Decoder&) = delete;
    Decoder& operator=(const Decoder&) = delete;
    ~Decoder() { ma_decoder_uninit(&decoder_); }

    Audio Read() {
        Audio audio;
        const auto format = ma_decoder_get_data_format(&decoder_, nullptr, &audio.channels, &audio.sample_rate, nullptr, 0);
        if (format != MA_SUCCESS || audio.channels == 0 || audio.sample_rate == 0) {
            throw std::runtime_error("Could not read the audio format.");
        }
        constexpr size_t kBlock = 4096;
        std::vector<float> buffer(kBlock * audio.channels);
        while (true) {
            ma_uint64 read = 0;
            const auto result = ma_decoder_read_pcm_frames(&decoder_, buffer.data(), kBlock, &read);
            if (result != MA_SUCCESS && result != MA_AT_END) {
                throw std::runtime_error(std::format("Audio decoding failed ({})", ma_result_description(result)));
            }
            audio.samples.insert(audio.samples.end(), buffer.begin(), buffer.begin() + read * audio.channels);
            if (result == MA_AT_END || read == 0) {
                break;
            }
        }
        if (audio.samples.empty() || !std::all_of(audio.samples.begin(), audio.samples.end(), [](float sample) {
                return std::isfinite(sample);
            })) {
            throw std::runtime_error("The input needs a nonempty stream of finite samples.");
        }
        return audio;
    }

private:
    ma_decoder decoder_{};
};

std::string PathText(const fs::path& path) {
    const auto text = path.u8string();
    return {reinterpret_cast<const char*>(text.data()), text.size()};
}

fs::path ReserveOutputPath(const fs::path& input, const std::string& suffix) {
    auto base = input.stem();
    base += suffix;
    const auto cwd = fs::current_path();
    for (size_t number = 1;; ++number) {
        auto name = base;
        if (number > 1) {
            name += std::format("_{}", number);
        }
        name += ".wav";
        const auto output = cwd / name;
        if (fs::exists(fs::symlink_status(output))) {
            continue;
        }
#ifdef _WIN32
        const HANDLE handle = CreateFileW(output.c_str(), GENERIC_WRITE, 0, nullptr, CREATE_NEW,
                                         FILE_ATTRIBUTE_NORMAL, nullptr);
        if (handle != INVALID_HANDLE_VALUE) {
            CloseHandle(handle);
            return output;
        }
        const auto error = GetLastError();
        if (error != ERROR_FILE_EXISTS && error != ERROR_ALREADY_EXISTS) {
            throw std::system_error(static_cast<int>(error), std::system_category(), "Could not create the output WAV");
        }
#else
        if (auto* file = std::fopen(output.c_str(), "wbx")) {
            std::fclose(file);
            return output;
        }
        const int error = errno;
        if (error != EEXIST) {
            throw std::system_error(error, std::generic_category(), "Could not create the output WAV");
        }
#endif
    }
}

fs::path WriteWav(const fs::path& input_path, const std::string& suffix, const Audio& audio, const std::vector<float>& output) {
    if (output.size() > (std::numeric_limits<uint32_t>::max() - 256) / sizeof(float)) {
        throw std::length_error("The float WAV would exceed the RIFF size limit.");
    }
    if (!std::all_of(output.begin(), output.end(), [](float sample) { return std::isfinite(sample); })) {
        throw std::runtime_error("PGHI produced non-finite samples.");
    }
    const auto path = ReserveOutputPath(input_path, suffix);
    const auto config = ma_encoder_config_init(ma_encoding_format_wav, ma_format_f32, audio.channels, audio.sample_rate);
    ma_encoder encoder{};
#ifdef _WIN32
    const auto opened = ma_encoder_init_file_w(path.c_str(), &config, &encoder);
#else
    const auto opened = ma_encoder_init_file(path.c_str(), &config, &encoder);
#endif
    if (opened != MA_SUCCESS) {
        std::error_code ignored;
        fs::remove(path, ignored);
        throw std::runtime_error(std::format("Could not open the output WAV ({})", ma_result_description(opened)));
    }
    const auto frames = output.size() / audio.channels;
    ma_uint64 written = 0;
    const auto result = ma_encoder_write_pcm_frames(&encoder, output.data(), frames, &written);
    ma_encoder_uninit(&encoder);
    if (result != MA_SUCCESS || written != frames) {
        std::error_code ignored;
        fs::remove(path, ignored);
        throw std::runtime_error("Could not write the complete output WAV.");
    }
    return path;
}

} // namespace qwqdsp_example::detail

int main() {
    using namespace qwqdsp_example::detail;
    try {
        const auto input_path = ReadPath();
        if (input_path.empty()) {
            throw std::invalid_argument("Enter a path to an audio file.");
        }
        auto audio = Decoder(input_path).Read();
        const size_t input_frames = audio.samples.size() / audio.channels;
        std::cout << std::format("{} Hz, {} channel(s), {:.3f} s\n", audio.sample_rate, audio.channels,
                                 static_cast<double>(input_frames) / audio.sample_rate);

        qwqdsp_example::PghiSettings settings;
        settings.window_size = ReadNumber("Window size (power of two, 256-8192) [2048]: ", settings.window_size);
        settings.transients = ReadSwitch("Transient anchors [Y/n]: ", settings.transients);
        settings.chirp = ReadSwitch("Experimental Gaussian chirp correction (mono) [y/N]: ", settings.chirp);
        if (settings.chirp && audio.channels != 1) {
            throw std::invalid_argument("Gaussian chirp correction currently requires mono input. Use n for stereo or multichannel.");
        }
        if (!std::has_single_bit(settings.window_size) || settings.window_size < 256 || settings.window_size > 8192) {
            throw std::invalid_argument("Window size must be a power of two between 256 and 8192.");
        }

        std::cout << "Analyzing PGHI once...\n" << std::flush;
        const qwqdsp_example::PghiAnalysis analysis(std::move(audio.samples), audio.channels, audio.sample_rate, settings);
        std::cout << std::format("Ready: {} STFT frames, {} sample analysis hop; spectra {}.\n",
                                 analysis.FrameCount(), analysis.AnalysisHop(),
                                 analysis.HasCachedSpectra() ? "cached" : "recomputed as needed for synthesis");

        double last_duration = 2;
        double last_pitch = 0;
        while (true) {
            std::cout << std::format("\nTime stretch 0.25-4 [{:.6g}], r = resample only, q = quit: ", last_duration) << std::flush;
            std::string command;
            if (!std::getline(std::cin, command)) {
                break;
            }
            TrimLine(command);
            if (command == "q" || command == "Q") {
                break;
            }
            try {
                const bool resample_only = command == "r" || command == "R";
                const double duration = resample_only ? last_duration : ParseNumber(command, last_duration);
                if (!std::isfinite(duration) || duration < 0.25 || duration > 4) {
                    throw std::invalid_argument("Time stretch must be between 0.25 and 4.");
                }
                const double pitch = ReadNumber(std::format("Pitch shift (-12 to +12 semitones) [{:+.6g}]: ", last_pitch), last_pitch);
                std::cout << (resample_only ? "Resampling...\n" : "Rendering from cached PGHI analysis...\n") << std::flush;
                const auto result = resample_only ? analysis.RenderResampled(pitch) : analysis.Render(duration, pitch);
                auto suffix = resample_only ? std::string("_resampled") : std::format("_pghi_{:.6g}x", duration);
                if (pitch != 0) {
                    suffix += std::format("_p{:+.6g}st", pitch);
                }
                const auto output_path = WriteWav(input_path, suffix, audio, result.samples);
                std::cout << std::format("{} frames -> {} frames, {} onset(s), {} phase anchor(s), {} chirp frame(s)\nSaved {}\n",
                                         input_frames, result.samples.size() / audio.channels, result.onset_count,
                                         result.anchor_count, result.chirp_frames, PathText(output_path));
                last_duration = duration;
                last_pitch = pitch;
            } catch (const std::exception& error) {
                if (std::cin.eof()) {
                    break;
                }
                std::cerr << std::format("PGHI: {}\n", error.what());
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << std::format("PGHI: {}\n", error.what());
        return 1;
    }
}
