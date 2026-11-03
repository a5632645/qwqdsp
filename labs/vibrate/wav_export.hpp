#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace qwqdsp_vibrate {

inline void WriteU16(std::ofstream& file, uint16_t value) {
    const char bytes[2] = {static_cast<char>(value & 0xFF),
                           static_cast<char>((value >> 8) & 0xFF)};
    file.write(bytes, 2);
}

inline void WriteU32(std::ofstream& file, uint32_t value) {
    const char bytes[4] = {
        static_cast<char>(value & 0xFF), static_cast<char>((value >> 8) & 0xFF),
        static_cast<char>((value >> 16) & 0xFF),
        static_cast<char>((value >> 24) & 0xFF)};
    file.write(bytes, 4);
}

// 单声道 16-bit PCM；自动把峰值归一化到 0.9。
inline bool WriteWavMono16(const std::string& path,
                           const std::vector<float>& samples,
                           int sample_rate) {
    if (samples.empty() || sample_rate <= 0) {
        return false;
    }
    float peak = 0.0f;
    for (float value : samples) {
        peak = std::max(peak, std::abs(value));
    }
    const float gain = peak > 1e-9f ? 0.9f / peak : 0.0f;

    std::ofstream file(path, std::ios::binary);
    if (!file) {
        return false;
    }
    const uint32_t data_bytes =
        static_cast<uint32_t>(samples.size() * sizeof(int16_t));
    file.write("RIFF", 4);
    WriteU32(file, 36 + data_bytes);
    file.write("WAVE", 4);
    file.write("fmt ", 4);
    WriteU32(file, 16);
    WriteU16(file, 1);  // PCM
    WriteU16(file, 1);  // mono
    WriteU32(file, static_cast<uint32_t>(sample_rate));
    WriteU32(file, static_cast<uint32_t>(sample_rate * 2));
    WriteU16(file, 2);  // block align
    WriteU16(file, 16); // bits per sample
    file.write("data", 4);
    WriteU32(file, data_bytes);
    for (float value : samples) {
        const long clamped = std::lround(
            static_cast<double>(std::clamp(value * gain, -1.0f, 1.0f)) * 32767.0);
        WriteU16(file, static_cast<uint16_t>(static_cast<int16_t>(clamped)));
    }
    return file.good();
}

inline std::vector<float> ResampleLinear(const std::vector<float>& input,
                                         double input_rate, int output_rate,
                                         std::size_t output_count) {
    std::vector<float> output(output_count, 0.0f);
    if (input.empty() || input_rate <= 0.0 || output_rate <= 0) {
        return output;
    }
    const double step = input_rate / static_cast<double>(output_rate);
    for (std::size_t index = 0; index < output_count; ++index) {
        const double position = static_cast<double>(index) * step;
        const std::size_t i = std::min(
            static_cast<std::size_t>(position), input.size() - 1);
        const std::size_t j = std::min(i + 1, input.size() - 1);
        const float fraction =
            static_cast<float>(position - static_cast<double>(i));
        output[index] = input[i] + (input[j] - input[i]) * fraction;
    }
    return output;
}

// 数据按 8 kHz 生成，但 WAV 头直接写成 96 kHz：
// 播放时相当于 12 倍升调；96 s 物理时间 = 768k 采样点 = 96 kHz 下 8 s。
constexpr int kExportSampleRate = 8000;
constexpr int kExportPlaybackRate = 96000;
// 96 s 物理时间。
constexpr float kExportSeconds = 96.0f;
// 导出用低阻尼，避免升调后衰减太快。
constexpr float kExportDamping = 0.03f;
constexpr int kExportStepsPerSecond1D = 16000;
constexpr int kExportStepsPerSecond2D = 1000;

// 从一个初始状态开始运行模型，记录探针，再重采样到 8 kHz 并写出 WAV。
template <typename Model, typename Step, typename Sample>
bool ExportProbeWav(const std::string& path, Model model, int steps_per_second,
                    float seconds, Step step, Sample sample) {
    if (steps_per_second <= 0 || seconds <= 0.0f) {
        return false;
    }
    const int step_count =
        static_cast<int>(std::lround(static_cast<double>(steps_per_second) * seconds));
    std::vector<float> values;
    values.reserve(static_cast<std::size_t>(step_count));
    const float dt = 1.0f / static_cast<float>(steps_per_second);
    for (int i = 0; i < step_count; ++i) {
        values.push_back(sample(model));
        step(model, dt);
    }

    const std::size_t output_count = static_cast<std::size_t>(
        std::lround(static_cast<double>(kExportSampleRate) * seconds));
    const std::vector<float> audio = ResampleLinear(
        values, steps_per_second, kExportSampleRate, output_count);
    return WriteWavMono16(path, audio, kExportPlaybackRate);
}

}  // namespace qwqdsp_vibrate
