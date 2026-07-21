#include "AudioFile.h"
#include "stb_image_write.h"
#include "work_dir.hpp"
#include <qwqdsp/pitch/swift_f0.hpp>
#include <qwqdsp/pitch/hide/swift_f0_model.hpp>
#include <qwqdsp/spectral/real_fft_adv.hpp>
#include <qwqdsp/window/hann.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <format>
#include <iostream>
#include <numbers>
#include <vector>

#include <qwqdsp/colormap/magma.hpp>

// ------------------------------------------------------------
// 模型参数常量
// ------------------------------------------------------------
constexpr int kSampleRate = 16000;
constexpr int kNFFT = 1024;
constexpr int kHop = 256;
constexpr int kPad = (kNFFT - kHop) / 2;
constexpr float kEps = 1.0e-8f;

// ------------------------------------------------------------
// swift_f0 基频检测 + STFT 谱渲染为 PNG
// ------------------------------------------------------------
int main() {
    // ---- 读取 WAV ----
    auto input_path = qwqdsp_support::GetInputDir() / "wormhole.wav";
    std::cout << std::format("Reading: {}\n", input_path.string());

    AudioFile<float> audio;
    if (!audio.load(input_path.string())) {
        std::cout << "FAIL: cannot load WAV file\n";
        return 1;
    }

    float fs = audio.getSampleRate();
    int const num_channels = audio.getNumChannels();
    int const num_samples = audio.getNumSamplesPerChannel();
    auto const& samples = audio.samples[0];

    std::cout << std::format("  fs={:.0f} Hz, ch={}, samples={}\n", fs, num_channels, num_samples);

    // ---- 重采样到 16kHz (简单线性插值) ----
    std::vector<float> audio_16k;
    int num_samples_16k;
    if (fs != kSampleRate) {
        float ratio = kSampleRate / fs;
        num_samples_16k = static_cast<int>(num_samples * ratio) + 1;
        audio_16k.resize(num_samples_16k);
        for (int i = 0; i < num_samples_16k; ++i) {
            float src_pos = i / ratio;
            int idx = static_cast<int>(src_pos);
            float frac = src_pos - idx;
            idx = std::min(idx, num_samples - 2);
            audio_16k[i] = samples[idx] * (1.0f - frac) + samples[idx + 1] * frac;
        }
        std::cout << std::format("  Resampled: {} -> {} samples\n", num_samples, num_samples_16k);
    }
    else {
        audio_16k.assign(samples.begin(), samples.end());
        num_samples_16k = num_samples;
    }

    // ---- swift_f0 模型参数对应的帧数 ----
    int const padded_len = num_samples_16k + 2 * kPad;
    int const num_frames = (padded_len - kNFFT) / kHop + 1;

    // ---- 图像参数 ----
    int const img_w = num_frames;
    int const img_h = 400;
    int const bpp = 4;
    std::vector<uint8_t> image(static_cast<size_t>(img_w * img_h * bpp), 0);

    // ---- 频谱频率范围 (对数, 适配 16kHz 采样率) ----
    float const freq_min = 55.0f;
    float const freq_max = 7000.0f;

    auto freqToY = [&](float hz) -> int {
        float norm = (std::log2(hz) - std::log2(freq_min)) / (std::log2(freq_max) - std::log2(freq_min));
        norm = std::clamp(norm, 0.0f, 1.0f);
        return static_cast<int>((1.0f - norm) * (img_h - 1));
    };

    // ---- STFT 预备 ----
    qwqdsp_spectral::RealFftAdv stft_fft;
    stft_fft.Init(kNFFT);
    std::vector<float> frame(kNFFT);
    std::vector<float> gain(stft_fft.NumBins()); // 513 bins

    // ---- STFT 能量统计 ----
    float mag_max = -1e9f;
    std::vector<std::vector<float>> frame_mags(num_frames);

    std::cout << std::format("STFT analysis of {} frames ...\n", num_frames);

    for (int f = 0; f < num_frames; ++f) {
        int offset = f * kHop;

        // 零填充 + 取帧
        for (int i = 0; i < kNFFT; ++i) {
            int src_idx = offset + i - kPad;
            if (src_idx >= 0 && src_idx < num_samples_16k)
                frame[i] = audio_16k[src_idx];
            else
                frame[i] = 0.0f;
        }

        // Hann 窗
        qwqdsp_window::Hann::ApplyWindow(frame, true);

        // FFT → 幅度谱
        stft_fft.FFTGainPhase(frame, gain);

        // dB 转换
        auto& mags = frame_mags[f];
        mags.resize(gain.size());
        for (size_t b = 0; b < mags.size(); ++b) {
            mags[b] = (gain[b] > 0.0f) ? 20.0f * std::log10(gain[b]) : -120.0f;
            if (mags[b] > mag_max)
                mag_max = mags[b];
        }
    }

    float mag_min = mag_max - 80.0f;

    std::cout << std::format("  Spectral dynamic range: {:.0f} - {:.0f} dB\n", mag_min, mag_max);

    // ---- 绘制 STFT 谱 ----
    for (int x = 0; x < img_w; ++x) {
        auto const& mags = frame_mags[x];
        for (int y = 0; y < img_h; ++y) {
            float norm_y = 1.0f - static_cast<float>(y) / (img_h - 1);
            float freq = freq_min * std::pow(2.0f, norm_y * std::log2(freq_max / freq_min));

            int bin = static_cast<int>(freq * kNFFT / kSampleRate);
            bin = std::clamp(bin, 0, static_cast<int>(mags.size()) - 1);

            float dB = mags[bin];
            float norm = (dB - mag_min) / (mag_max - mag_min);
            norm = std::clamp(norm, 0.0f, 1.0f);
            int ci = static_cast<int>(norm * 255.0f);
            ci = std::clamp(ci, 0, 255);

            size_t idx = static_cast<size_t>(y * img_w * bpp + x * bpp);
            image[idx + 0] = qwqdsp_colormap::Magma::kTable[ci][0];
            image[idx + 1] = qwqdsp_colormap::Magma::kTable[ci][1];
            image[idx + 2] = qwqdsp_colormap::Magma::kTable[ci][2];
            image[idx + 3] = 255;
        }
    }

    // ---- 绘制参考线 (半音) ----
    float note_freqs[] = {65.41f,  73.42f,  82.41f,  87.31f,  98.00f,  110.0f,  123.47f, 130.81f, 146.83f, 164.81f,
                          174.61f, 196.0f,  220.0f,  246.94f, 261.63f, 293.66f, 329.63f, 349.23f, 392.0f,  440.0f,
                          493.88f, 523.25f, 587.33f, 659.25f, 698.46f, 783.99f, 880.0f,  987.77f, 1046.5f, 1174.7f,
                          1318.5f, 1396.9f, 1568.0f, 1760.0f, 1975.5f, 2093.0f, 2349.3f, 2637.0f, 2793.8f, 3136.0f,
                          3520.0f, 3951.1f, 4186.0f, 4698.6f, 5274.0f, 5587.7f, 6271.9f, 7040.0f};
    for (auto f : note_freqs) {
        if (f < freq_min || f > freq_max)
            continue;
        int y = freqToY(f);
        for (int x = 0; x < img_w; ++x) {
            size_t idx = static_cast<size_t>(y * img_w * bpp + x * bpp);
            image[idx + 0] = static_cast<uint8_t>(image[idx + 0] * 0.6f + 180 * 0.4f);
            image[idx + 1] = static_cast<uint8_t>(image[idx + 1] * 0.6f + 180 * 0.4f);
            image[idx + 2] = static_cast<uint8_t>(image[idx + 2] * 0.6f + 180 * 0.4f);
        }
    }

    // ---- swift_f0 推理 ----
    std::cout << "swift_f0 inference ...\n";

    // 构建 log-mag 频谱 [T, 132]
    Eigen::Tensor<float, 2> log_mag(num_frames, qwqdsp_swift_f0::kNumMelBins);
    for (int t = 0; t < num_frames; ++t) {
        auto const& mags = frame_mags[t];
        for (int f = 0; f < qwqdsp_swift_f0::kNumMelBins; ++f) {
            int bin = f + qwqdsp_swift_f0::kSliceStart; // 3:135
            float m = (bin < static_cast<int>(mags.size())) ? std::pow(10.0f, mags[bin] / 20.0f) : 0.0f;
            log_mag(t, f) = std::log(m + kEps);
        }
    }

    // 推理
    Eigen::VectorXf pitch_hz(num_frames);
    Eigen::VectorXf confidence(num_frames);

    qwqdsp_swift_f0::SwiftF0Inference inference;
    inference.Process(log_mag, pitch_hz, confidence);

    std::cout << std::format("  Done, {} frames\n", num_frames);

    // ---- 绘制基频轨迹 (亮青色，叠加在 STFT 上) ----
    int prev_y = -1;
    for (int x = 0; x < img_w; ++x) {
        float hz = pitch_hz(x);
        float conf = confidence(x);

        // 仅显示置信度 > 0.5 的有声帧
        if (conf < 0.5f || hz < freq_min || hz > freq_max) {
            prev_y = -1;
            continue;
        }

        int y = freqToY(hz);
        if (y < 0 || y >= img_h) {
            prev_y = -1;
            continue;
        }

        // 连线
        if (prev_y >= 0) {
            int y1 = std::min(prev_y, y);
            int y2 = std::max(prev_y, y);
            for (int ly = y1; ly <= y2; ++ly) {
                size_t idx = static_cast<size_t>(ly * img_w * bpp + x * bpp);
                image[idx + 0] = 0;
                image[idx + 1] = 255;
                image[idx + 2] = 200;
                image[idx + 3] = 255;
            }
        }

        // 3×3 点，颜色根据置信度亮度
        int brightness = static_cast<int>(180 + 75 * conf);
        for (int dx = -1; dx <= 1; ++dx) {
            int px = x + dx;
            if (px < 0 || px >= img_w)
                continue;
            for (int dy = -1; dy <= 1; ++dy) {
                int py = y + dy;
                if (py < 0 || py >= img_h)
                    continue;
                size_t idx = static_cast<size_t>(py * img_w * bpp + px * bpp);
                image[idx + 0] = 0;
                image[idx + 1] = static_cast<uint8_t>(brightness);
                image[idx + 2] = static_cast<uint8_t>(brightness * 200 / 255);
                image[idx + 3] = 255;
            }
        }
        prev_y = y;
    }

    // ---- 保存 PNG ----
    auto out_dir = qwqdsp_support::GetOutputDir();
    std::filesystem::create_directories(out_dir);
    auto out_path = out_dir / "wormhole_swift_f0.png";

    int result = stbi_write_png(out_path.string().c_str(), img_w, img_h, bpp, image.data(), img_w * bpp);

    if (result) {
        std::cout << std::format("Output: {} ({}x{})\n", out_path.string(), img_w, img_h);
    }
    else {
        std::cout << "FAIL: cannot write PNG\n";
        return 1;
    }

    // ---- 统计 ----
    int voiced = 0;
    for (int i = 0; i < num_frames; ++i) {
        if (confidence(i) > 0.5f)
            voiced++;
    }
    std::cout << std::format("Voiced frames: {}/{}\n", voiced, num_frames);
    std::cout << "Done!\n";
    return 0;
}
