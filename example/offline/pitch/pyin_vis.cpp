#include "AudioFile.h"
#include "stb_image_write.h"
#include "work_dir.hpp"
#include <qwqdsp/pitch/pyin.hpp>
#include <qwqdsp/pitch/hide/pyin_core.hpp>
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
// 将 pYIN 基频检测 + STFT 谱渲染为 PNG 图像
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

    std::cout << std::format("  fs={} Hz, ch={}, samples={}\n", static_cast<int>(fs), num_channels, num_samples);

    // ---- 参数 ----
    int const block_size = 2048;
    int const step_size = 512;
    int const num_frames = (num_samples - block_size) / step_size + 1;

    // ---- pYIN ----
    qwqdsp_pitch::Pyin pyin;
    pyin.Init(fs, block_size, step_size);
    pyin.SetMinPitch(60.0f);
    pyin.SetMaxPitch(900.0f);

    // ---- STFT 预备 ----
    qwqdsp_spectral::RealFftAdv stft_fft;
    stft_fft.Init(block_size);
    std::vector<float> window(block_size);
    std::vector<float> fft_frame(block_size);
    std::vector<std::complex<float>> stft_spec(stft_fft.NumBins());

    std::vector<float> block(block_size);

    // ---- 图像参数 ----
    int const img_w = num_frames;
    int const img_h = 400;
    int const bpp = 4;
    std::vector<uint8_t> image(static_cast<size_t>(img_w * img_h * bpp), 0);

    // ---- 频谱频率范围 (对数) ----
    float const freq_min = 55.0f; // A1
    float const freq_max = 2000.0f;

    auto freqToY = [&](float hz) -> int {
        float norm = (std::log2(hz) - std::log2(freq_min)) / (std::log2(freq_max) - std::log2(freq_min));
        norm = std::clamp(norm, 0.0f, 1.0f);
        return static_cast<int>((1.0f - norm) * (img_h - 1));
    };

    // ---- STFT 能量统计 (用于动态范围) ----
    float mag_max = -1e9f;
    std::vector<std::vector<float>> frame_mags(num_frames);

    std::cout << std::format("STFT analysis of {} frames ...\n", num_frames);

    for (int f = 0; f < num_frames; ++f) {
        int offset = f * step_size;
        for (int i = 0; i < block_size; ++i)
            block[i] = samples[offset + i];

        // ---- STFT ----
        fft_frame = block;
        qwqdsp_window::Hann::ApplyWindow(fft_frame, true);
        std::vector<float> gain(stft_fft.NumBins());
        stft_fft.FFTGainPhase(fft_frame, gain);

        auto& mags = frame_mags[f];
        mags.resize(gain.size());
        for (size_t b = 0; b < mags.size(); ++b) {
            mags[b] = gain[b];
            if (mags[b] > 0.0f) {
                mags[b] = 20.0f * std::log10(mags[b]);
            }
            else {
                mags[b] = -120.0f;
            }
            if (mags[b] > mag_max)
                mag_max = mags[b];
        }

        // ---- pYIN ----
        pyin.Process(block);
    }

    float mag_min = mag_max - 80.0f; // 80 dB 动态范围

    std::cout << std::format("  Spectral dynamic range: {:.0f} - {:.0f} dB\n", mag_min, mag_max);

    // ---- 绘制 STFT 谱 ----
    for (int x = 0; x < img_w; ++x) {
        auto const& mags = frame_mags[x];
        for (int y = 0; y < img_h; ++y) {
            // y 坐标映射到频率
            float norm_y = 1.0f - static_cast<float>(y) / (img_h - 1);
            float freq = freq_min * std::pow(2.0f, norm_y * std::log2(freq_max / freq_min));

            // 找最近的 FFT bin
            int bin = static_cast<int>(freq * block_size / fs);
            if (bin < 0)
                bin = 0;
            if (bin >= static_cast<int>(mags.size()))
                bin = static_cast<int>(mags.size()) - 1;

            // dB → colormap index
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
    float note_freqs[] = {65.41f,  73.42f,  82.41f,  87.31f,  98.00f,  110.0f,  123.47f, 130.81f, 146.83f,
                          164.81f, 174.61f, 196.0f,  220.0f,  246.94f, 261.63f, 293.66f, 329.63f, 349.23f,
                          392.0f,  440.0f,  493.88f, 523.25f, 587.33f, 659.25f, 698.46f, 783.99f, 880.0f,
                          987.77f, 1046.5f, 1174.7f, 1318.5f, 1396.9f, 1568.0f, 1760.0f, 1975.5f};
    for (auto f : note_freqs) {
        if (f < freq_min || f > freq_max)
            continue;
        int y = freqToY(f);
        for (int x = 0; x < img_w; ++x) {
            size_t idx = static_cast<size_t>(y * img_w * bpp + x * bpp);
            // 半透明灰线
            image[idx + 0] = static_cast<uint8_t>(image[idx + 0] * 0.6f + 180 * 0.4f);
            image[idx + 1] = static_cast<uint8_t>(image[idx + 1] * 0.6f + 180 * 0.4f);
            image[idx + 2] = static_cast<uint8_t>(image[idx + 2] * 0.6f + 180 * 0.4f);
            // alpha 不变
        }
    }

    // ---- HMM 解码 ----
    auto pitch_track = pyin.GetPitchTrack();
    std::cout << std::format("Decoding complete, {} frames\n", pitch_track.size());

    // ---- 绘制基频轨迹 (亮青色，叠加在 STFT 上) ----
    int prev_y = -1;
    for (int x = 0; x < img_w && static_cast<size_t>(x) < pitch_track.size(); ++x) {
        float hz = pitch_track[x];
        if (hz <= 0.0f) {
            prev_y = -1;
            continue;
        }

        int y = freqToY(hz);
        if (y < 0 || y >= img_h) {
            prev_y = -1;
            continue;
        }

        if (prev_y >= 0) {
            int y1 = std::min(prev_y, y);
            int y2 = std::max(prev_y, y);
            for (int ly = y1; ly <= y2; ++ly) {
                size_t idx = static_cast<size_t>(ly * img_w * bpp + x * bpp);
                image[idx + 0] = 0;   // R
                image[idx + 1] = 255; // G
                image[idx + 2] = 200; // B
                image[idx + 3] = 255;
            }
        }

        // 3px 点
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
                image[idx + 1] = 255;
                image[idx + 2] = 200;
                image[idx + 3] = 255;
            }
        }
        prev_y = y;
    }

    // ---- 保存 PNG ----
    auto out_dir = qwqdsp_support::GetOutputDir();
    std::filesystem::create_directories(out_dir);
    auto out_path = out_dir / "wormhole_pyin.png";

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
    for (auto p : pitch_track)
        if (p > 0.0f)
            voiced++;
    std::cout << std::format("Voiced frames: {}/{}\n", voiced, pitch_track.size());
    std::cout << "Done!\n";
    return 0;
}
