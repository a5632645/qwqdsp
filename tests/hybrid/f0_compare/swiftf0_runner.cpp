#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>
#include "../../gui/pitch/swift_f0_rt/decimator.hpp"
#include "../../gui/pitch/swift_f0_rt/inference_fast.hpp"
#include "../../gui/pitch/swift_f0_rt/pitch_frontend.hpp"

namespace {

struct WavData {
    int sample_rate{};
    int channels{};
    int bits_per_sample{};
    std::uint16_t format{};
    std::vector<float> samples;
};

/**
 * @brief 读取 little-endian 无符号 16 位整数。
 * @param data WAV 字节缓冲区。
 * @param offset 整数起始偏移。
 * @return 解码后的整数。
 */
std::uint16_t ReadU16(const std::vector<std::uint8_t>& data, size_t offset) {
    if (offset + 2 > data.size()) {
        throw std::runtime_error("WAV 文件截断");
    }
    return static_cast<std::uint16_t>(data[offset]) |
           (static_cast<std::uint16_t>(data[offset + 1]) << 8U);
}

/**
 * @brief 读取 little-endian 无符号 32 位整数。
 * @param data WAV 字节缓冲区。
 * @param offset 整数起始偏移。
 * @return 解码后的整数。
 */
std::uint32_t ReadU32(const std::vector<std::uint8_t>& data, size_t offset) {
    if (offset + 4 > data.size()) {
        throw std::runtime_error("WAV 文件截断");
    }
    return static_cast<std::uint32_t>(data[offset]) |
           (static_cast<std::uint32_t>(data[offset + 1]) << 8U) |
           (static_cast<std::uint32_t>(data[offset + 2]) << 16U) |
           (static_cast<std::uint32_t>(data[offset + 3]) << 24U);
}

/**
 * @brief 检查 WAV 字节缓冲区中的四字节标签。
 * @param data WAV 字节缓冲区。
 * @param offset 标签起始偏移。
 * @param tag 期望的四字节标签。
 * @return 标签相同返回 true。
 */
bool HasTag(const std::vector<std::uint8_t>& data, size_t offset, const char (&tag)[5]) {
    return offset + 4 <= data.size() && data[offset] == static_cast<std::uint8_t>(tag[0]) &&
           data[offset + 1] == static_cast<std::uint8_t>(tag[1]) &&
           data[offset + 2] == static_cast<std::uint8_t>(tag[2]) &&
           data[offset + 3] == static_cast<std::uint8_t>(tag[3]);
}

/**
 * @brief 从 WAV 文件读取单声道浮点样本。
 * @param path WAV 文件路径。
 * @return 音频格式与单声道样本。
 */
WavData ReadWav(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("无法打开 WAV: " + path);
    }
    std::vector<std::uint8_t> data((std::istreambuf_iterator<char>(file)), {});
    if (data.size() < 12 || !HasTag(data, 0, "RIFF") || !HasTag(data, 8, "WAVE")) {
        throw std::runtime_error("不是 RIFF/WAVE 文件: " + path);
    }

    WavData wav;
    size_t fmt_offset = 0;
    size_t fmt_size = 0;
    size_t audio_offset = 0;
    size_t audio_size = 0;
    size_t offset = 12;
    while (offset + 8 <= data.size()) {
        std::uint32_t const chunk_size = ReadU32(data, offset + 4);
        size_t const chunk_begin = offset + 8;
        size_t const chunk_end = chunk_begin + static_cast<size_t>(chunk_size);
        if (chunk_end > data.size()) {
            throw std::runtime_error("WAV chunk 越界");
        }
        if (HasTag(data, offset, "fmt ")) {
            fmt_offset = chunk_begin;
            fmt_size = chunk_size;
        }
        else if (HasTag(data, offset, "data")) {
            audio_offset = chunk_begin;
            audio_size = chunk_size;
        }
        offset = chunk_end + (chunk_size & 1U);
    }

    if (fmt_offset == 0 || fmt_size < 16 || audio_offset == 0) {
        throw std::runtime_error("WAV 缺少 fmt 或 data chunk");
    }
    wav.format = ReadU16(data, fmt_offset);
    wav.channels = static_cast<int>(ReadU16(data, fmt_offset + 2));
    wav.sample_rate = static_cast<int>(ReadU32(data, fmt_offset + 4));
    wav.bits_per_sample = static_cast<int>(ReadU16(data, fmt_offset + 14));
    if ((wav.format != 1 && wav.format != 3) || wav.channels <= 0 || wav.sample_rate <= 0 ||
        (wav.format == 1 && wav.bits_per_sample != 16 && wav.bits_per_sample != 24 &&
         wav.bits_per_sample != 32) ||
        (wav.format == 3 && wav.bits_per_sample != 32)) {
        throw std::runtime_error("swiftf0 helper 只支持 PCM 16/24/32 位或 IEEE float WAV");
    }

    size_t const bytes_per_sample = static_cast<size_t>(wav.bits_per_sample / 8);
    size_t const bytes_per_frame = bytes_per_sample * static_cast<size_t>(wav.channels);
    if (bytes_per_frame == 0 || audio_size % bytes_per_frame != 0) {
        throw std::runtime_error("WAV data chunk 大小无效");
    }
    size_t const frame_count = audio_size / bytes_per_frame;
    wav.samples.resize(frame_count, 0.0f);
    for (size_t frame = 0; frame < frame_count; ++frame) {
        float sum = 0.0f;
        for (int channel = 0; channel < wav.channels; ++channel) {
            size_t const sample_offset = audio_offset + frame * bytes_per_frame +
                                          static_cast<size_t>(channel) * bytes_per_sample;
            float sample = 0.0f;
            if (wav.format == 3) {
                std::uint32_t bits = ReadU32(data, sample_offset);
                static_assert(sizeof(float) == sizeof(bits));
                std::memcpy(&sample, &bits, sizeof(sample));
            }
            else if (wav.bits_per_sample == 16) {
                std::int16_t const value = static_cast<std::int16_t>(ReadU16(data, sample_offset));
                sample = static_cast<float>(value) / 32768.0f;
            }
            else if (wav.bits_per_sample == 24) {
                std::uint32_t raw = static_cast<std::uint32_t>(data[sample_offset]) |
                                    (static_cast<std::uint32_t>(data[sample_offset + 1]) << 8U) |
                                    (static_cast<std::uint32_t>(data[sample_offset + 2]) << 16U);
                if ((raw & 0x00800000U) != 0) {
                    raw |= 0xFF000000U;
                }
                sample = static_cast<float>(static_cast<std::int32_t>(raw)) / 8388608.0f;
            }
            else {
                std::int32_t const value = static_cast<std::int32_t>(ReadU32(data, sample_offset));
                sample = static_cast<float>(value) / 2147483648.0f;
            }
            sum += sample;
        }
        wav.samples[frame] = sum / static_cast<float>(wav.channels);
    }
    return wav;
}

/**
 * @brief 对整段音频执行仓库内 swiftf0 前端和推理。
 * @param wav 输入 WAV 数据。
 * @return 每行依次为分析中心时间（秒）、F0（Hz）和置信度。
 */
std::vector<std::array<float, 3>> AnalyzeSwiftF0(const WavData& wav) {
    if (wav.sample_rate != 48000) {
        throw std::runtime_error("swiftf0 helper 要求输入采样率为 48000 Hz");
    }

    swift_f0_rt::Decimator decimator;
    decimator.Init(3, 80.0f, 32);
    std::vector<float> decimated;
    decimator.Process(wav.samples, decimated);

    swift_f0_rt::LogMagStft stft;
    stft.Init();
    std::vector<float> log_mag;
    std::vector<float> center_samples;
    stft.Process(decimated, [&](std::span<const float> frame, std::int64_t center_16k) {
        log_mag.insert(log_mag.end(), frame.begin(), frame.end());
        center_samples.push_back(static_cast<float>(center_16k * 3 - decimator.DelaySamples()));
    });
    if (center_samples.empty()) {
        throw std::runtime_error("swiftf0 没有生成分析帧");
    }

    swift_f0_rt::FastInference inference;
    inference.Init();
    std::vector<float> pitch(center_samples.size());
    std::vector<float> confidence(center_samples.size());
    inference.Process(log_mag.data(), static_cast<int>(center_samples.size()), pitch.data(),
                      confidence.data());

    std::vector<std::array<float, 3>> result;
    result.reserve(center_samples.size());
    for (size_t i = 0; i < center_samples.size(); ++i) {
        result.push_back({center_samples[i] / 48000.0f, pitch[i], confidence[i]});
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            std::cerr << "用法: swiftf0_runner.exe <wav>\n";
            return 2;
        }
        WavData const wav = ReadWav(argv[1]);
        std::vector<std::array<float, 3>> const result = AnalyzeSwiftF0(wav);
        std::cout << "# time_s,f0_hz,confidence\n";
        for (auto const& row : result) {
            std::cout << row[0] << ',' << row[1] << ',' << row[2] << '\n';
        }
        return 0;
    }
    catch (std::exception const& error) {
        std::cerr << "swiftf0 helper 错误: " << error.what() << '\n';
        return 1;
    }
}
