#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <span>

#include "miniaudio.h"
#include "raylib.h"
#include "slider.hpp"

#include "qwqdsp/window/hann.hpp"

namespace formant_shifter {

static constexpr float kSampleRate = 48000.0f;
static constexpr int kWindowWidth = 480;
static constexpr int kWindowHeight = 240;
static constexpr const char* kWindowTitle = "Formant Shifter — Granular";
static constexpr float kMinFormantSemitones = -12.0f;
static constexpr float kMaxFormantSemitones = 12.0f;
static constexpr float kRatioMin = 0.5f;
static constexpr float kRatioMax = 2.0f;

namespace detail {

// ------------------------------------------------------------
//  粒子式共振峰移动器
// ------------------------------------------------------------
struct FormantShifter {
    // ----- 常量 -----
    static constexpr std::size_t kHop = 256;             // 粒子间隔
    static constexpr std::size_t kGrainSize = 2 * kHop;  // 粒子长度
    static constexpr std::size_t kLatency = kGrainSize;  // 输出延迟
    static constexpr std::size_t kRingSize = 8192;       // 环形缓冲长度（2 的幂）
    static constexpr std::size_t kRingMask = kRingSize - 1;
    static_assert((kRingSize & (kRingSize - 1)) == 0, "环形缓冲长度必须为二的幂");

    /**
     * @brief 将半音偏移量换算为频率倍率
     * @param[in] semitones 半音偏移量
     * @return 对应的频率倍率
     */
    static float semitonesToRatio(float semitones) noexcept {
        return std::exp2(semitones / 12.0f);
    }

    /**
     * @brief 构造函数，使用 qwqdsp 的 Hann 窗预计算窗权重
     * @note 周期窗（for_analyze_not_fir=true）保证相邻粒子窗平方和稳定
     */
    FormantShifter() noexcept {
        qwqdsp_window::Hann::Window(std::span<float>{window_}, true);
    }

    /**
     * @brief 初始化环形缓冲与采样计数
     */
    void init() noexcept {
        input_ring_.fill(0.0f);
        ola_ring_.fill(0.0f);
        weight_ring_.fill(0.0f);
        input_position_ = 0;
    }

    /**
     * @brief 处理一段单声道实时音频
     * @param[in] input 输入采样缓冲
     * @param[out] output 输出采样缓冲
     * @param[in] frame_count 缓冲中的采样数
     * @param[in] formant_ratio 共振峰频率倍率
     */
    void process(const float* input, float* output, std::size_t frame_count, float formant_ratio,
                 bool reverse_grain) noexcept {
        const float safe_ratio = std::clamp(formant_ratio, kRatioMin, kRatioMax);
        for (std::size_t i = 0; i < frame_count; ++i) {
            output[i] = processSample(input[i], safe_ratio, reverse_grain);
        }
    }

private:
    /**
     * @brief 将绝对采样位置映射到环形缓冲索引
     * @param[in] position 绝对采样位置
     * @return 环形缓冲中的数组索引
     */
    static std::size_t ringIndex(std::uint64_t position) noexcept {
        return static_cast<std::size_t>(position) & kRingMask;
    }

    /**
     * @brief 处理单个采样并输出延迟结果
     * @param[in] input_sample 当前输入采样
     * @param[in] ratio 共振峰频率倍率
     * @return 当前输出采样
     */
    float processSample(float input_sample, float ratio, bool reverse_grain) noexcept {
        input_ring_[ringIndex(input_position_)] = input_sample;

        // 先读取延迟位置的 OLA 结果，避免被本次新粒子覆盖
        float out = 0.0f;
        if (input_position_ >= kLatency) {
            const std::size_t index = ringIndex(input_position_ - kLatency);
            const float weight = weight_ring_[index];
            if (weight > 0.0001f) {
                out = ola_ring_[index] / weight;
            }
            ola_ring_[index] = 0.0f;
            weight_ring_[index] = 0.0f;
        }

        // 每 hop 个采样合成一个新粒子
        if (input_position_ >= kGrainSize && (input_position_ % kHop) == 0) {
            addGrain(input_position_ - kGrainSize, ratio, reverse_grain);
        }

        ++input_position_;
        return std::clamp(out, -1.0f, 1.0f);
    }

    /**
     * @brief 读取一段粒子、重采样并 OLA 叠加到输出环形缓冲
     * @param[in] grain_start 粒子起始绝对采样位置
     * @param[in] ratio 共振峰频率倍率
     * @note 流程：读 size 个采样 → Hann 窗 →（可选）时间反转 → 按倍率重采样 → 截断/补零到 size → 再 Hann 窗 → OLA
     */
    void addGrain(std::uint64_t grain_start, float ratio, bool reverse_grain) noexcept {
        // 1. 读取 size 个采样并施加 Hann 窗
        std::array<float, kGrainSize> grain;
        for (std::size_t i = 0; i < kGrainSize; ++i) {
            grain[i] = input_ring_[ringIndex(grain_start + i)] * window_[i];
        }

        // 1.5 分析窗后时间反转（可选）
        if (reverse_grain) {
            std::reverse(grain.begin(), grain.end());
        }

        // 2. 按倍率重采样到 new_len 长度
        const std::size_t new_len = std::clamp(
            static_cast<std::size_t>(std::lround(static_cast<double>(kGrainSize) / static_cast<double>(ratio))),
            std::size_t{2}, 2 * kGrainSize);

        // 3. 截断或补零到 size 长度（线性插值）
        std::array<float, kGrainSize> resampled;
        resampled.fill(0.0f);
        const std::size_t keep = std::min(new_len, kGrainSize);
        const double scale = static_cast<double>(kGrainSize - 1) / static_cast<double>(new_len - 1);
        for (std::size_t i = 0; i < keep; ++i) {
            const double pos = static_cast<double>(i) * scale;
            const std::size_t i0 = static_cast<std::size_t>(pos);
            const std::size_t i1 = std::min(i0 + 1, kGrainSize - 1);
            const float fraction = static_cast<float>(pos - static_cast<double>(i0));
            resampled[i] = std::lerp(grain[i0], grain[i1], fraction);
        }

        // 4. 施加 Hann 窗并按 hop 间隔 OLA 叠加（以窗平方和归一化）
        for (std::size_t i = 0; i < kGrainSize; ++i) {
            const float window = window_[i];
            const std::size_t index = ringIndex(grain_start + i);
            ola_ring_[index] += resampled[i] * window;
            weight_ring_[index] += window * window;
        }
    }

    std::array<float, kGrainSize> window_{};
    std::array<float, kRingSize> input_ring_{};
    std::array<float, kRingSize> ola_ring_{};
    std::array<float, kRingSize> weight_ring_{};
    std::uint64_t input_position_ = 0;
};

} // namespace detail

static detail::FormantShifter s_formant_shifter;
static std::atomic<float> s_formant_semitones{0.0f};
static std::atomic<bool> s_reverse_grain{false};

// ------------------------------------------------------------
//  miniaudio 回调
// ------------------------------------------------------------
extern "C" void MaCallback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;

    auto* input = static_cast<const float*>(pInput);
    auto* output = static_cast<float*>(pOutput);

    if (input == nullptr || output == nullptr)
        return;

    const float ratio = detail::FormantShifter::semitonesToRatio(s_formant_semitones.load(std::memory_order_relaxed));
    const bool reverse = s_reverse_grain.load(std::memory_order_relaxed);
    s_formant_shifter.process(input, output, static_cast<std::size_t>(frameCount), ratio, reverse);
}

} // namespace formant_shifter

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------
int main(void) {
    using namespace formant_shifter;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, kWindowTitle);
    SetTargetFPS(60);

    // ── 初始化 DSP ──
    s_formant_shifter.init();

    // ── miniaudio full-duplex（单声道入 → 单声道出）──
    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.playback.format = ma_format_f32;
    config.playback.channels = 1;
    config.sampleRate = (ma_uint32)kSampleRate;
    config.dataCallback = MaCallback;
    config.pUserData = nullptr;
    config.periodSizeInMilliseconds = 10;

    ma_device device;
    ma_result result = ma_device_init(nullptr, &config, &device);
    bool audio_ok = (result == MA_SUCCESS);
    if (audio_ok) {
        ma_device_start(&device);
    }

    // ── Formant 旋钮 ──
    Knob formant_knob;
    formant_knob.set_bound(175, 78, 130, 132)
        .set_title("Formant")
        .set_range(kMinFormantSemitones, kMaxFormantSemitones, 0.1f, 0.0f)
        .set_value(0.0f)
        .set_sensitivity(2)
        .set_name_font_size(17)
        .set_number_font_size(13)
        .set_fore_color(Color{235, 235, 238, 255})
        .set_bg_color(Color{24, 26, 29, 255});
    formant_knob.value_to_text_function = [](float value) { return std::format("{:+.1f} st", value); };
    formant_knob.on_value_change = [](float value) {
        s_formant_semitones.store(value, std::memory_order_relaxed);
    };
    s_formant_semitones.store(formant_knob.get_value(), std::memory_order_relaxed);

    // ── main loop ──
    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(Color{18, 20, 23, 255});

        // 顶栏
        DrawRectangle(0, 0, kWindowWidth, 72, Color{24, 26, 29, 255});
        DrawRectangle(0, 71, kWindowWidth, 1, Color{48, 52, 57, 255});
        DrawText("FORMANT SHIFTER", 24, 18, 25, Color{232, 233, 235, 255});
        DrawText("REAL-TIME GRANULAR PROCESSOR", 25, 48, 10, Color{126, 132, 139, 255});

        // 音频状态
        constexpr int status_begin_x = 300;
        const Color status_color = audio_ok ? Color{76, 201, 151, 255} : Color{225, 172, 74, 255};
        DrawCircle(status_begin_x, 30, 5.0f, status_color);
        DrawText(audio_ok ? "AUDIO ONLINE" : "AUDIO OFFLINE", status_begin_x + 13, 24, 12, status_color);

        // 旋钮
        formant_knob.display();

        // ── REVERSE toggle ──
        {
            static constexpr Rectangle kToggleRect{170, 213, 140, 22};
            const Vector2 mouse_pos = GetMousePosition();
            if (CheckCollisionPointRec(mouse_pos, kToggleRect) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
                s_reverse_grain.store(!s_reverse_grain.load(std::memory_order_relaxed), std::memory_order_relaxed);
            }
            const bool reverse = s_reverse_grain.load(std::memory_order_relaxed);
            DrawRectangleRec(kToggleRect, reverse ? Color{76, 201, 151, 255} : Color{40, 42, 46, 255});
            DrawRectangleLinesEx(kToggleRect, 1.0f, Color{60, 64, 70, 255});
            DrawText(reverse ? "REVERSE: ON" : "REVERSE: OFF", static_cast<int>(kToggleRect.x) + 32,
                     static_cast<int>(kToggleRect.y) + 6, 12,
                     reverse ? Color{10, 10, 12, 255} : Color{150, 155, 162, 255});
        }

        EndDrawing();
    }

    // ── cleanup ──
    if (audio_ok) {
        ma_device_stop(&device);
        ma_device_uninit(&device);
    }
    CloseWindow();
    return 0;
}
