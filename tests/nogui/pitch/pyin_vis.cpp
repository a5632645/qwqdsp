#include "AudioFile.h"
#include "stb_image_write.h"
#include "work_dir.hpp"
#include <qwqdsp/pitch/pyin.hpp>
#include <qwqdsp/pitch/pyin/pyin_core.hpp>
#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/window/hann.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <numbers>
#include <vector>

// ------------------------------------------------------------
// Magma colormap (matplotlib), 256 entries, RGBA 8-bit
// ------------------------------------------------------------
static constexpr uint8_t kMagma[256][4] = {
    {0,   0,   4,   255},
    {1,   0,   6,   255},
    {2,   0,   10,  255},
    {3,   0,   14,  255},
    {5,   0,   18,  255},
    {7,   0,   22,  255},
    {9,   0,   26,  255},
    {11,  1,   30,  255},
    {13,  1,   34,  255},
    {15,  2,   38,  255},
    {17,  2,   42,  255},
    {19,  3,   46,  255},
    {21,  4,   49,  255},
    {24,  5,   53,  255},
    {26,  6,   56,  255},
    {28,  7,   60,  255},
    {30,  9,   63,  255},
    {32,  10,  66,  255},
    {34,  12,  69,  255},
    {36,  13,  72,  255},
    {38,  15,  75,  255},
    {40,  17,  78,  255},
    {42,  18,  81,  255},
    {44,  20,  83,  255},
    {46,  22,  86,  255},
    {48,  24,  88,  255},
    {50,  26,  91,  255},
    {52,  28,  93,  255},
    {53,  30,  95,  255},
    {55,  32,  97,  255},
    {57,  34,  99,  255},
    {59,  36,  101, 255},
    {60,  38,  103, 255},
    {62,  40,  105, 255},
    {64,  42,  107, 255},
    {65,  44,  108, 255},
    {67,  46,  110, 255},
    {68,  48,  111, 255},
    {70,  50,  113, 255},
    {71,  52,  114, 255},
    {73,  54,  115, 255},
    {74,  56,  116, 255},
    {76,  58,  117, 255},
    {77,  60,  118, 255},
    {78,  62,  119, 255},
    {80,  64,  120, 255},
    {81,  66,  121, 255},
    {82,  68,  122, 255},
    {84,  70,  122, 255},
    {85,  72,  123, 255},
    {86,  74,  123, 255},
    {87,  76,  124, 255},
    {89,  78,  124, 255},
    {90,  80,  125, 255},
    {91,  82,  125, 255},
    {92,  84,  125, 255},
    {93,  86,  126, 255},
    {94,  87,  126, 255},
    {95,  89,  126, 255},
    {96,  91,  126, 255},
    {97,  93,  126, 255},
    {98,  95,  126, 255},
    {99,  97,  126, 255},
    {100, 99,  126, 255},
    {101, 100, 126, 255},
    {102, 102, 126, 255},
    {103, 104, 126, 255},
    {104, 106, 126, 255},
    {105, 108, 125, 255},
    {106, 109, 125, 255},
    {107, 111, 125, 255},
    {108, 113, 124, 255},
    {109, 115, 124, 255},
    {110, 116, 123, 255},
    {111, 118, 123, 255},
    {112, 120, 122, 255},
    {113, 121, 122, 255},
    {114, 123, 121, 255},
    {115, 125, 120, 255},
    {116, 126, 120, 255},
    {117, 128, 119, 255},
    {118, 130, 118, 255},
    {119, 131, 117, 255},
    {120, 133, 116, 255},
    {121, 135, 115, 255},
    {122, 136, 114, 255},
    {123, 138, 113, 255},
    {124, 139, 112, 255},
    {125, 141, 111, 255},
    {126, 142, 110, 255},
    {128, 144, 109, 255},
    {129, 146, 108, 255},
    {130, 147, 107, 255},
    {131, 149, 105, 255},
    {132, 150, 104, 255},
    {133, 152, 103, 255},
    {134, 153, 101, 255},
    {135, 155, 100, 255},
    {136, 156, 99,  255},
    {137, 158, 97,  255},
    {138, 159, 96,  255},
    {140, 161, 94,  255},
    {141, 162, 93,  255},
    {142, 164, 91,  255},
    {143, 165, 90,  255},
    {144, 167, 88,  255},
    {145, 168, 86,  255},
    {146, 170, 85,  255},
    {147, 171, 83,  255},
    {148, 173, 81,  255},
    {149, 174, 80,  255},
    {151, 176, 78,  255},
    {152, 177, 76,  255},
    {153, 179, 74,  255},
    {154, 180, 72,  255},
    {155, 181, 71,  255},
    {156, 183, 69,  255},
    {157, 184, 67,  255},
    {158, 186, 65,  255},
    {160, 187, 63,  255},
    {161, 189, 61,  255},
    {162, 190, 59,  255},
    {163, 191, 57,  255},
    {164, 193, 55,  255},
    {165, 194, 53,  255},
    {166, 196, 51,  255},
    {167, 197, 49,  255},
    {168, 198, 47,  255},
    {170, 200, 45,  255},
    {171, 201, 43,  255},
    {172, 203, 41,  255},
    {173, 204, 39,  255},
    {174, 205, 37,  255},
    {175, 207, 34,  255},
    {176, 208, 32,  255},
    {177, 209, 30,  255},
    {179, 211, 28,  255},
    {180, 212, 26,  255},
    {181, 213, 23,  255},
    {182, 215, 21,  255},
    {183, 216, 19,  255},
    {184, 217, 16,  255},
    {185, 219, 14,  255},
    {186, 220, 12,  255},
    {188, 221, 9,   255},
    {189, 223, 7,   255},
    {190, 224, 5,   255},
    {191, 225, 3,   255},
    {192, 227, 1,   255},
    {193, 228, 0,   255},
    {194, 229, 0,   255},
    {196, 231, 1,   255},
    {197, 232, 3,   255},
    {198, 233, 6,   255},
    {199, 235, 10,  255},
    {200, 236, 14,  255},
    {202, 237, 19,  255},
    {203, 239, 24,  255},
    {204, 240, 29,  255},
    {205, 241, 35,  255},
    {207, 242, 41,  255},
    {208, 244, 47,  255},
    {209, 245, 53,  255},
    {210, 246, 59,  255},
    {212, 247, 65,  255},
    {213, 248, 71,  255},
    {214, 249, 77,  255},
    {216, 250, 83,  255},
    {217, 251, 89,  255},
    {218, 252, 95,  255},
    {219, 253, 102, 255},
    {221, 254, 108, 255},
    {222, 254, 114, 255},
    {223, 255, 121, 255},
    {224, 255, 127, 255},
    {226, 255, 134, 255},
    {227, 255, 140, 255},
    {228, 255, 147, 255},
    {229, 255, 153, 255},
    {231, 255, 160, 255},
    {232, 255, 167, 255},
    {233, 255, 173, 255},
    {234, 255, 180, 255},
    {236, 255, 187, 255},
    {237, 255, 193, 255},
    {238, 255, 200, 255},
    {239, 255, 207, 255},
    {241, 254, 214, 255},
    {242, 254, 221, 255},
    {243, 253, 228, 255},
    {244, 253, 235, 255},
    {246, 252, 242, 255},
    {247, 252, 249, 255},
    {248, 251, 255, 255},
    {249, 250, 255, 255},
    {251, 250, 255, 255},
    {252, 249, 255, 255},
    {253, 248, 255, 255},
    {254, 247, 255, 255},
    {255, 247, 255, 255}
};

// ------------------------------------------------------------
// 将 pYIN 基频检测 + STFT 谱渲染为 PNG 图像
// ------------------------------------------------------------
int main() {
    // ---- 读取 WAV ----
    auto input_path = qwqdsp_support::GetInputDir() / "wormhole.wav";
    std::printf("读取: %s\n", input_path.string().c_str());

    AudioFile<float> audio;
    if (!audio.load(input_path.string())) {
        std::printf("FAIL: 无法加载 WAV 文件\n");
        return 1;
    }

    float fs = audio.getSampleRate();
    int const num_channels = audio.getNumChannels();
    int const num_samples = audio.getNumSamplesPerChannel();
    auto const& samples = audio.samples[0];

    std::printf("  fs=%d Hz, ch=%d, samples=%d\n", static_cast<int>(fs), num_channels, num_samples);

    // ---- 参数 ----
    int const block_size = 2048;
    int const step_size = 512;
    int const num_frames = (num_samples - block_size) / step_size + 1;

    // ---- pYIN ----
    qwqdsp_pitch::Pyin pyin;
    pyin.Init(fs, block_size, step_size);
    pyin.SetThresholdDistribution(2);
    pyin.SetMinPitch(60.0f);
    pyin.SetMaxPitch(900.0f);

    // ---- STFT 预备 ----
    qwqdsp_spectral::RealFFT stft_fft;
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

    std::printf("STFT 分析 %d 帧 ...\n", num_frames);

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

    std::printf("  频谱动态: %.0f–%.0f dB\n", mag_min, mag_max);

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
            image[idx + 0] = kMagma[ci][0];
            image[idx + 1] = kMagma[ci][1];
            image[idx + 2] = kMagma[ci][2];
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
    std::printf("解码完成, %zu 帧\n", pitch_track.size());

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
        std::printf("输出: %s (%dx%d)\n", out_path.string().c_str(), img_w, img_h);
    }
    else {
        std::printf("FAIL: 写入 PNG 失败\n");
        return 1;
    }

    // ---- 统计 ----
    int voiced = 0;
    for (auto p : pitch_track)
        if (p > 0.0f)
            voiced++;
    std::printf("有声帧: %d/%zu\n", voiced, pitch_track.size());
    std::printf("完成!\n");
    return 0;
}
