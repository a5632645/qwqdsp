/// @brief 实时 LPC 声码器演示界面
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <format>
#include <iostream>
#include <string>

#include "miniaudio.h"
#include "raylib.h"
#include "slider.hpp"

#include "lpc_vocoder_dsp.hpp"

namespace lpc_vocoder {

static constexpr float kSampleRate = 48000.0f;
static constexpr int kWindowWidth = 640;
static constexpr int kWindowHeight = 400;
static constexpr const char* kWindowTitle = "LPC Vocoder — Analysis / Synthesis";
static constexpr int kDefaultOrder = 40;
static constexpr std::size_t kDisplaySize = 8192;
static constexpr std::size_t kDisplaySamples = 2048;

static detail::LpcVocoder g_vocoder;
static std::atomic<int> g_order{kDefaultOrder};
static std::atomic<float> g_mix{1.0f};
static std::atomic<float> g_noise_mix{0.1f};
static std::atomic<float> g_voicing_threshold{0.6f};
static std::atomic<float> g_f0_hz{0.0f};
static std::atomic<float> g_voicing{0.0f};
static std::atomic<float> g_input_peak{0.0f};
static std::atomic<float> g_output_peak{0.0f};
static std::array<float, kDisplaySize> g_display_input{};
static std::array<float, kDisplaySize> g_display_output{};
static std::atomic<std::size_t> g_display_pos{0};

// ------------------------------------------------------------
// 音频回调
// ------------------------------------------------------------

/**
 * @brief 响应 miniaudio 的全双工单声道数据请求
 * @param[in] device 发起回调的音频设备
 * @param[out] output_buffer 单声道浮点输出缓冲
 * @param[in] input_buffer 单声道浮点输入缓冲
 * @param[in] frame_count 本次请求的音频帧数
 */
extern "C" void MaCallback(ma_device* device, void* output_buffer, const void* input_buffer,
                           ma_uint32 frame_count) {
    (void)device;

    auto* output = static_cast<float*>(output_buffer);
    const auto* input = static_cast<const float*>(input_buffer);
    if (output == nullptr) {
        return;
    }
    if (input == nullptr) {
        std::fill_n(output, frame_count, 0.0f);
        return;
    }

    const auto stats = g_vocoder.Process(
        input, output, frame_count, g_order.load(std::memory_order_relaxed),
        g_mix.load(std::memory_order_relaxed), g_noise_mix.load(std::memory_order_relaxed),
        g_voicing_threshold.load(std::memory_order_relaxed));

    g_f0_hz.store(stats.f0_hz, std::memory_order_relaxed);
    g_voicing.store(stats.voicing, std::memory_order_relaxed);
    g_input_peak.store(std::max(stats.input_peak, g_input_peak.load(std::memory_order_relaxed) * 0.82f),
                       std::memory_order_relaxed);
    g_output_peak.store(std::max(stats.output_peak, g_output_peak.load(std::memory_order_relaxed) * 0.82f),
                        std::memory_order_relaxed);

    std::size_t position = g_display_pos.load(std::memory_order_relaxed);
    for (ma_uint32 i = 0; i < frame_count; ++i) {
        g_display_input[position] = input[i];
        g_display_output[position] = output[i];
        position = (position + 1) % kDisplaySize;
    }
    g_display_pos.store(position, std::memory_order_release);
}

// ------------------------------------------------------------
// 绘制
// ------------------------------------------------------------

/**
 * @brief 绘制一段按列取最小最大值压缩后的波形
 * @param[in] samples 等间隔采样的波形数据
 * @param[in] x 绘图区左侧坐标
 * @param[in] y 绘图区顶部坐标
 * @param[in] width 绘图区宽度
 * @param[in] height 绘图区高度
 * @param[in] color 波形颜色
 */
static void DrawSignal(const std::array<float, kDisplaySamples>& samples, int x, int y, int width, int height,
                       Color color) {
    const int center_y = y + height / 2;
    const float half_height = static_cast<float>(height) * 0.45f;
    for (int column = 0; column < width; ++column) {
        const std::size_t begin = static_cast<std::size_t>(column) * kDisplaySamples / static_cast<std::size_t>(width);
        const std::size_t end =
            std::max(begin + 1, static_cast<std::size_t>(column + 1) * kDisplaySamples / static_cast<std::size_t>(width));
        float low = samples[begin];
        float high = samples[begin];
        for (std::size_t i = begin + 1; i < end; ++i) {
            low = std::min(low, samples[i]);
            high = std::max(high, samples[i]);
        }
        const float top = static_cast<float>(center_y) - std::clamp(high, -1.0f, 1.0f) * half_height;
        const float bottom = static_cast<float>(center_y) - std::clamp(low, -1.0f, 1.0f) * half_height;
        DrawLineV({static_cast<float>(x + column), top}, {static_cast<float>(x + column), bottom}, color);
    }
}

/**
 * @brief 绘制一个紧凑的水平峰值电平表
 * @param[in] x 电平表左侧坐标
 * @param[in] y 电平表顶部坐标
 * @param[in] width 电平表宽度
 * @param[in] value 线性峰值幅度
 * @param[in] label 电平表标签
 * @param[in] color 填充颜色
 */
static void DrawMeter(int x, int y, int width, float value, const char* label, Color color) {
    const float normalized = std::clamp(value, 0.0f, 1.0f);
    const int bar_width = width - 26;
    const int fill_width = static_cast<int>(static_cast<float>(bar_width) * normalized);
    DrawText(label, x, y - 2, 10, Color{126, 132, 139, 255});
    DrawRectangle(x + 26, y, bar_width, 7, Color{45, 48, 52, 255});
    DrawRectangle(x + 26, y, fill_width, 7, color);
}

/**
 * @brief 绘制波形面板与清浊音、基频状态
 * @param[in] input_view 最近的输入波形
 * @param[in] output_view 最近的输出波形
 * @param[in] audio_running 音频设备是否成功启动
 */
static void DrawPanel(const std::array<float, kDisplaySamples>& input_view,
                      const std::array<float, kDisplaySamples>& output_view, bool audio_running) {
    static constexpr Color kBackground{18, 20, 23, 255};
    static constexpr Color kPanel{24, 26, 29, 255};
    static constexpr Color kText{232, 233, 235, 255};
    static constexpr Color kMuted{126, 132, 139, 255};
    static constexpr Color kAccent{76, 201, 151, 255};
    static constexpr Color kWarning{225, 172, 74, 255};
    static constexpr Color kInputWave{92, 132, 178, 255};

    ClearBackground(kBackground);
    DrawRectangle(0, 0, kWindowWidth, 72, kPanel);
    DrawRectangle(0, 71, kWindowWidth, 1, Color{48, 52, 57, 255});
    DrawText("LPC VOCODER", 24, 18, 25, kText);
    DrawText("BURG ANALYSIS / ROSENBERG + MIXED EXCITATION", 25, 48, 10, kMuted);

    const Color status_color = audio_running ? kAccent : kWarning;
    DrawCircle(500, 30, 5.0f, status_color);
    DrawText(audio_running ? "AUDIO ONLINE" : "AUDIO OFFLINE", 513, 24, 12, status_color);

    const float f0 = g_f0_hz.load(std::memory_order_relaxed);
    const float voicing = g_voicing.load(std::memory_order_relaxed);
    const std::string pitch_text = f0 > 0.0f ? std::format("F0  {:6.1f} Hz", f0) : std::string{"F0       -- Hz"};
    DrawText(pitch_text.c_str(), 24, 236, 13, f0 > 0.0f ? kText : kMuted);

    const bool voiced = voicing > 0.5f;
    DrawText(voiced ? "VOICED" : "UNVOICED", 176, 236, 13, voiced ? kAccent : kMuted);
    DrawMeter(286, 242, 120, g_input_peak.load(std::memory_order_relaxed), "IN", kInputWave);
    DrawMeter(430, 242, 120, g_output_peak.load(std::memory_order_relaxed), "OUT", kAccent);

    g_input_peak.store(g_input_peak.load(std::memory_order_relaxed) * 0.95f, std::memory_order_relaxed);
    g_output_peak.store(g_output_peak.load(std::memory_order_relaxed) * 0.95f, std::memory_order_relaxed);

    static constexpr int kPanelX = 24;
    static constexpr int kPanelY = 260;
    static constexpr int kPanelW = 592;
    static constexpr int kPanelH = 128;
    static constexpr int kWaveY = 280;
    static constexpr int kWaveH = 100;
    DrawRectangle(kPanelX, kPanelY, kPanelW, kPanelH, Color{20, 22, 26, 255});
    DrawRectangleLinesEx(Rectangle{static_cast<float>(kPanelX), static_cast<float>(kPanelY),
                                   static_cast<float>(kPanelW), static_cast<float>(kPanelH)},
                         1.0f, Color{45, 48, 52, 255});

    DrawRectangle(kPanelX + 8, kPanelY + 6, 8, 8, kInputWave);
    DrawText("INPUT", kPanelX + 21, kPanelY + 5, 10, kMuted);
    DrawRectangle(kPanelX + 72, kPanelY + 6, 8, 8, kAccent);
    DrawText("SYNTHESIS", kPanelX + 85, kPanelY + 5, 10, kMuted);

    DrawRectangle(kPanelX + 1, kWaveY + kWaveH / 2, kPanelW - 2, 1, Color{42, 45, 50, 255});
    DrawSignal(input_view, kPanelX + 2, kWaveY, kPanelW - 4, kWaveH, kInputWave);
    DrawSignal(output_view, kPanelX + 2, kWaveY, kPanelW - 4, kWaveH, kAccent);
}

} // namespace lpc_vocoder

/**
 * @brief 启动实时 LPC 声码器演示程序
 * @return 正常退出返回零，音频设备初始化失败时仍保留界面并返回零
 */
int main() {
    using namespace lpc_vocoder;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, kWindowTitle);
    SetTargetFPS(60);

    g_vocoder.Init(kSampleRate);

    Knob order_knob;
    order_knob.set_bound(20, 90, 136, 132)
        .set_title("LPC Order")
        .set_range(8.0f, 80.0f, 1.0f, static_cast<float>(kDefaultOrder))
        .set_value(static_cast<float>(kDefaultOrder))
        .set_sensitivity(2)
        .set_name_font_size(17)
        .set_number_font_size(13)
        .set_fore_color(Color{235, 235, 238, 255})
        .set_bg_color(Color{24, 26, 29, 255});
    order_knob.value_to_text_function = [](float value) { return std::format("{:.0f}", value); };
    order_knob.on_value_change = [](float value) {
        g_order.store(static_cast<int>(std::lround(value)), std::memory_order_relaxed);
    };
    g_order.store(static_cast<int>(std::lround(order_knob.get_value())), std::memory_order_relaxed);

    Knob balance_knob;
    balance_knob.set_bound(168, 90, 136, 132)
        .set_title("Pulse / Noise")
        .set_range(0.0f, 1.0f, 0.01f, 0.1f)
        .set_value(0.1f)
        .set_sensitivity(2)
        .set_name_font_size(17)
        .set_number_font_size(13)
        .set_fore_color(Color{235, 235, 238, 255})
        .set_bg_color(Color{24, 26, 29, 255});
    balance_knob.value_to_text_function = [](float value) { return std::format("{:.0f}%", value * 100.0f); };
    balance_knob.on_value_change = [](float value) { g_noise_mix.store(value, std::memory_order_relaxed); };
    g_noise_mix.store(balance_knob.get_value(), std::memory_order_relaxed);

    Knob threshold_knob;
    threshold_knob.set_bound(316, 90, 136, 132)
        .set_title("Voicing Thresh")
        .set_range(0.1f, 0.9f, 0.01f, 0.6f)
        .set_value(0.6f)
        .set_sensitivity(2)
        .set_name_font_size(17)
        .set_number_font_size(13)
        .set_fore_color(Color{235, 235, 238, 255})
        .set_bg_color(Color{24, 26, 29, 255});
    threshold_knob.value_to_text_function = [](float value) { return std::format("{:.2f}", value); };
    threshold_knob.on_value_change = [](float value) { g_voicing_threshold.store(value, std::memory_order_relaxed); };
    g_voicing_threshold.store(threshold_knob.get_value(), std::memory_order_relaxed);

    Knob mix_knob;
    mix_knob.set_bound(464, 90, 136, 132)
        .set_title("Mix")
        .set_range(0.0f, 1.0f, 0.01f, 1.0f)
        .set_value(1.0f)
        .set_sensitivity(2)
        .set_name_font_size(17)
        .set_number_font_size(13)
        .set_fore_color(Color{235, 235, 238, 255})
        .set_bg_color(Color{24, 26, 29, 255});
    mix_knob.value_to_text_function = [](float value) { return std::format("{:.0f}%", value * 100.0f); };
    mix_knob.on_value_change = [](float value) { g_mix.store(value, std::memory_order_relaxed); };
    g_mix.store(mix_knob.get_value(), std::memory_order_relaxed);

    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.playback.format = ma_format_f32;
    config.playback.channels = 1;
    config.sampleRate = static_cast<ma_uint32>(kSampleRate);
    config.dataCallback = MaCallback;
    config.pUserData = nullptr;
    config.periodSizeInMilliseconds = 10;

    ma_device device{};
    const ma_result init_result = ma_device_init(nullptr, &config, &device);
    bool audio_running = init_result == MA_SUCCESS;
    if (audio_running && ma_device_start(&device) != MA_SUCCESS) {
        audio_running = false;
    }
    if (!audio_running) {
        std::cout << std::format("miniaudio: 设备启动失败 ({})\n", static_cast<int>(init_result));
    }

    while (!WindowShouldClose()) {
        std::array<float, kDisplaySamples> input_view{};
        std::array<float, kDisplaySamples> output_view{};
        const std::size_t position = g_display_pos.load(std::memory_order_acquire);
        const std::size_t start = (position + kDisplaySize - kDisplaySamples) % kDisplaySize;
        for (std::size_t i = 0; i < kDisplaySamples; ++i) {
            const std::size_t index = (start + i) % kDisplaySize;
            input_view[i] = g_display_input[index];
            output_view[i] = g_display_output[index];
        }

        BeginDrawing();
        DrawPanel(input_view, output_view, audio_running);
        order_knob.display();
        balance_knob.display();
        threshold_knob.display();
        mix_knob.display();
        EndDrawing();
    }

    if (init_result == MA_SUCCESS) {
        ma_device_stop(&device);
        ma_device_uninit(&device);
    }
    CloseWindow();
    return 0;
}
