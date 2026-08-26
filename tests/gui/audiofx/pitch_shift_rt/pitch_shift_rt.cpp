#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <format>
#include <iostream>
#include <string>

#include "miniaudio.h"
#include "raylib.h"
#include "slider.hpp"

#include "pitch_shifter_rt.hpp"

namespace pitch_shift_rt {

static constexpr float kSampleRate = 48000.0f;
static constexpr int kWindowWidth = 760;
static constexpr int kWindowHeight = 260;
static constexpr const char* kWindowTitle = "Spectral Pitch Shifter - Real-time";
static constexpr std::array<const char*, 5> kAlgorithmNames{"PGHI", "PGHI + Flux", "PGHI + SuperFlux", "Phase Lock",
                                                            "Transient Vocoder"};

class PitchShiftDemo {
public:
    /**
     * @brief 初始化实时处理器和音高旋钮
     */
    void init() {
        processor_.init(kSampleRate);
        pitch_knob_.set_bound((kWindowWidth - 92) / 2, 126, 92, 98)
            .set_title("Pitch Shift")
            .set_range(-12.0f, 12.0f, 1.0f, 0.0f)
            .set_value(0.0f)
            .set_sensitivity(2)
            .set_name_font_size(15)
            .set_number_font_size(13)
            .set_fore_color(Color{235, 235, 238, 255})
            .set_bg_color(Color{24, 26, 29, 255});
        pitch_knob_.value_to_text_function = [](float value) { return std::format("{:+.0f} st", value); };
        pitch_knob_.on_value_change = [this](float value) { pitch_shift_.store(value, std::memory_order_relaxed); };
    }

    /**
     * @brief 处理音频设备提供的一段单声道采样
     * @param[in] input 单声道输入缓冲
     * @param[out] output 单声道输出缓冲
     * @param[in] frame_count 缓冲中的采样数
     */
    void processAudio(const float* input, float* output, std::size_t frame_count) {
        processor_.setAlgorithm(static_cast<Algorithm>(algorithm_.load(std::memory_order_relaxed)));
        processor_.setPitchShift(pitch_shift_.load(std::memory_order_relaxed));
        processor_.process(input, output, frame_count);
    }

    /**
     * @brief 绘制算法选择器、音高旋钮和设备状态
     * @param[in] audio_running 音频设备是否成功启动
     */
    void draw(bool audio_running) {
        static constexpr Color kBackground{18, 20, 23, 255};
        static constexpr Color kPanel{24, 26, 29, 255};
        static constexpr Color kText{232, 233, 235, 255};
        static constexpr Color kMuted{126, 132, 139, 255};
        static constexpr Color kAccent{76, 201, 151, 255};
        static constexpr Color kWarning{225, 172, 74, 255};

        ClearBackground(kBackground);
        DrawRectangle(0, 0, kWindowWidth, 64, kPanel);
        DrawRectangle(0, 63, kWindowWidth, 1, Color{48, 52, 57, 255});
        DrawText("SPECTRAL PITCH SHIFTER", 22, 15, 23, kText);
        DrawText("REAL-TIME / MONO", 23, 42, 10, kMuted);

        const Color status_color = audio_running ? kAccent : kWarning;
        DrawCircle(kWindowWidth - 126, 30, 5.0f, status_color);
        DrawText(audio_running ? "AUDIO ONLINE" : "AUDIO OFFLINE", kWindowWidth - 114, 24, 12, status_color);

        drawAlgorithmSelector();
        pitch_knob_.display();
        DrawText("right click: reset", kWindowWidth - 118, kWindowHeight - 17, 10, kMuted);
    }

    /**
     * @brief 响应 miniaudio 的全双工单声道数据请求
     * @param[in] device 发起回调的音频设备
     * @param[out] output_buffer 单声道浮点输出缓冲
     * @param[in] input_buffer 单声道浮点输入缓冲
     * @param[in] frame_count 本次请求的音频帧数
     */
    static void audioCallback(ma_device* device, void* output_buffer, const void* input_buffer, ma_uint32 frame_count) {
        auto* output = static_cast<float*>(output_buffer);
        const auto* input = static_cast<const float*>(input_buffer);
        auto* demo = static_cast<PitchShiftDemo*>(device->pUserData);
        if (output == nullptr) {
            return;
        }
        if (input == nullptr || demo == nullptr) {
            std::fill_n(output, frame_count, 0.0f);
            return;
        }
        demo->processAudio(input, output, static_cast<std::size_t>(frame_count));
    }
private:
    void drawAlgorithmSelector() {
        constexpr float kLeft = 20.0f;
        constexpr float kTop = 82.0f;
        constexpr float kGap = 5.0f;
        constexpr float kHeight = 31.0f;
        constexpr float kWidth = (static_cast<float>(kWindowWidth) - 2.0f * kLeft - 4.0f * kGap) / 5.0f;
        const Vector2 mouse = GetMousePosition();
        const int selected = algorithm_.load(std::memory_order_relaxed);

        for (std::size_t i = 0; i < kAlgorithmNames.size(); ++i) {
            const Rectangle bounds{kLeft + static_cast<float>(i) * (kWidth + kGap), kTop, kWidth, kHeight};
            const bool hovered = CheckCollisionPointRec(mouse, bounds);
            const bool active = selected == static_cast<int>(i);
            const Color border = hovered ? Color{190, 194, 198, 255} : Color{75, 80, 85, 255};
            const Color text = active ? Color{18, 20, 23, 255} : Color{178, 182, 187, 255};
            if (active) {
                DrawRectangleRec(bounds, Color{232, 233, 235, 255});
            }
            else {
                DrawRectangleLinesEx(bounds, 1.0f, border);
            }
            const int font_size = 12;
            const int text_width = MeasureText(kAlgorithmNames[i], font_size);
            DrawText(kAlgorithmNames[i], static_cast<int>(bounds.x + (bounds.width - text_width) * 0.5f),
                     static_cast<int>(bounds.y + 9.0f), font_size, text);
            if (hovered && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
                algorithm_.store(static_cast<int>(i), std::memory_order_relaxed);
            }
        }
    }

    RealtimePitchShifter<WindowType::BlackmanHarrisThreeTerm, 512> processor_;
    Knob pitch_knob_;
    std::atomic<float> pitch_shift_{0.0f};
    std::atomic<int> algorithm_{0};
};

} // namespace pitch_shift_rt

/**
 * @brief 启动实时频域移调演示程序
 * @return 正常退出返回零，音频设备不可用时仍保留界面
 */
int main() {
    using namespace pitch_shift_rt;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, kWindowTitle);
    SetTargetFPS(60);

    static PitchShiftDemo demo;
    demo.init();

    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.playback.format = ma_format_f32;
    config.playback.channels = 1;
    config.sampleRate = static_cast<ma_uint32>(kSampleRate);
    config.dataCallback = PitchShiftDemo::audioCallback;
    config.pUserData = &demo;
    config.periodSizeInMilliseconds = 5;

    ma_device device{};
    const ma_result init_result = ma_device_init(nullptr, &config, &device);
    ma_result audio_result = init_result;
    bool audio_running = false;
    if (init_result == MA_SUCCESS) {
        audio_result = ma_device_start(&device);
        audio_running = audio_result == MA_SUCCESS;
    }
    if (!audio_running) {
        std::cout << std::format("miniaudio: device start failed ({})\n", static_cast<int>(audio_result));
    }

    while (!WindowShouldClose()) {
        BeginDrawing();
        demo.draw(audio_running);
        EndDrawing();
    }

    if (init_result == MA_SUCCESS) {
        if (audio_running) {
            ma_device_stop(&device);
        }
        ma_device_uninit(&device);
    }
    CloseWindow();
    return 0;
}
