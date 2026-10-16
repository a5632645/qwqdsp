#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <format>

#include "miniaudio.h"
#include "raylib.h"

#include "echo/ima_int16_delay_line.hpp"
#include "echo/naive_echo_gui.hpp"
#include "echo/quantized_delay_line.hpp"

// ------------------------------------------------------------
//  EchoGuiConcept — echo GUI 包装器必须满足的接口
//
//  任何用作存储格式槽位的类型都需要提供以下成员：
//   - static constexpr int kWindowWidth
//   - static constexpr int kWindowHeight
//   - static constexpr const char* kStorageName
//   - static constexpr const char* kStorageDescription
//   - void Init(float sample_rate)
//   - void Process(const float* input, float* output, size_t frame_count) noexcept
//   - void CreateKnobs()
//   - void DisplayKnobs()
// ------------------------------------------------------------
template <typename T>
concept EchoGuiConcept = requires(T& obj, const float* in, float* out, size_t n, float sr) {
    { T::kWindowWidth } -> std::convertible_to<int>;
    { T::kWindowHeight } -> std::convertible_to<int>;
    { T::kStorageName } -> std::convertible_to<const char*>;
    { T::kStorageDescription } -> std::convertible_to<const char*>;

    { obj.Init(sr) } -> std::same_as<void>;
    { obj.Process(in, out, n) } noexcept -> std::same_as<void>;
    { obj.CreateKnobs() } -> std::same_as<void>;
    { obj.DisplayKnobs() } -> std::same_as<void>;
};

// ------------------------------------------------------------
//  存储格式槽位 — 新增一种格式：加一个实例，再往 kAlgorithms 里加一项
// ------------------------------------------------------------
static NaiveEchoGui<echo_demo::FloatStorage> s_f32;
static NaiveEchoGui<echo_demo::Q8_23Storage> s_q8_23;
static NaiveEchoGui<echo_demo::Q3_20Storage> s_q3_20;
static NaiveEchoGui<echo_demo::Q3_12Storage> s_q3_12;
static NaiveEchoGui<echo_demo::Adpcm8Storage> s_adpcm;
static NaiveEchoGui<echo_demo::ImaInt16Storage> s_ima16;

static_assert(EchoGuiConcept<NaiveEchoGui<echo_demo::FloatStorage>>,
              "NaiveEchoGui 必须满足 EchoGuiConcept：需要 kWindowWidth/kWindowHeight/kStorageName/"
              "kStorageDescription 静态常量，以及 Init/Process/CreateKnobs/DisplayKnobs 方法");
static_assert(EchoGuiConcept<NaiveEchoGui<echo_demo::Q8_23Storage>>,
              "NaiveEchoGui 必须满足 EchoGuiConcept：需要 kWindowWidth/kWindowHeight/kStorageName/"
              "kStorageDescription 静态常量，以及 Init/Process/CreateKnobs/DisplayKnobs 方法");
static_assert(EchoGuiConcept<NaiveEchoGui<echo_demo::Q3_20Storage>>,
              "NaiveEchoGui 必须满足 EchoGuiConcept：需要 kWindowWidth/kWindowHeight/kStorageName/"
              "kStorageDescription 静态常量，以及 Init/Process/CreateKnobs/DisplayKnobs 方法");
static_assert(EchoGuiConcept<NaiveEchoGui<echo_demo::Q3_12Storage>>,
              "NaiveEchoGui 必须满足 EchoGuiConcept：需要 kWindowWidth/kWindowHeight/kStorageName/"
              "kStorageDescription 静态常量，以及 Init/Process/CreateKnobs/DisplayKnobs 方法");
static_assert(EchoGuiConcept<NaiveEchoGui<echo_demo::Adpcm8Storage>>,
              "NaiveEchoGui 必须满足 EchoGuiConcept：需要 kWindowWidth/kWindowHeight/kStorageName/"
              "kStorageDescription 静态常量，以及 Init/Process/CreateKnobs/DisplayKnobs 方法");
static_assert(EchoGuiConcept<NaiveEchoGui<echo_demo::ImaInt16Storage>>,
              "NaiveEchoGui 必须满足 EchoGuiConcept：需要 kWindowWidth/kWindowHeight/kStorageName/"
              "kStorageDescription 静态常量，以及 Init/Process/CreateKnobs/DisplayKnobs 方法");

struct EchoAlgorithm {
    const char* name;
    const char* description;
    void (*Init)(float sample_rate);
    void (*Process)(const float* input, float* output, size_t frame_count) noexcept;
    void (*CreateKnobs)();
    void (*DisplayKnobs)();
};

static const std::array<EchoAlgorithm, 6> kAlgorithms{{
    {
        NaiveEchoGui<echo_demo::FloatStorage>::kStorageName,
        NaiveEchoGui<echo_demo::FloatStorage>::kStorageDescription,
        [](float sample_rate) { s_f32.Init(sample_rate); },
        [](const float* in, float* out, size_t n) noexcept { s_f32.Process(in, out, n); },
        [] { s_f32.CreateKnobs(); },
        [] { s_f32.DisplayKnobs(); },
    },
    {
        NaiveEchoGui<echo_demo::Q8_23Storage>::kStorageName,
        NaiveEchoGui<echo_demo::Q8_23Storage>::kStorageDescription,
        [](float sample_rate) { s_q8_23.Init(sample_rate); },
        [](const float* in, float* out, size_t n) noexcept { s_q8_23.Process(in, out, n); },
        [] { s_q8_23.CreateKnobs(); },
        [] { s_q8_23.DisplayKnobs(); },
    },
    {
        NaiveEchoGui<echo_demo::Q3_20Storage>::kStorageName,
        NaiveEchoGui<echo_demo::Q3_20Storage>::kStorageDescription,
        [](float sample_rate) { s_q3_20.Init(sample_rate); },
        [](const float* in, float* out, size_t n) noexcept { s_q3_20.Process(in, out, n); },
        [] { s_q3_20.CreateKnobs(); },
        [] { s_q3_20.DisplayKnobs(); },
    },
    {
        NaiveEchoGui<echo_demo::Q3_12Storage>::kStorageName,
        NaiveEchoGui<echo_demo::Q3_12Storage>::kStorageDescription,
        [](float sample_rate) { s_q3_12.Init(sample_rate); },
        [](const float* in, float* out, size_t n) noexcept { s_q3_12.Process(in, out, n); },
        [] { s_q3_12.CreateKnobs(); },
        [] { s_q3_12.DisplayKnobs(); },
    },
    {
        NaiveEchoGui<echo_demo::Adpcm8Storage>::kStorageName,
        NaiveEchoGui<echo_demo::Adpcm8Storage>::kStorageDescription,
        [](float sample_rate) { s_adpcm.Init(sample_rate); },
        [](const float* in, float* out, size_t n) noexcept { s_adpcm.Process(in, out, n); },
        [] { s_adpcm.CreateKnobs(); },
        [] { s_adpcm.DisplayKnobs(); },
    },
    {
        NaiveEchoGui<echo_demo::ImaInt16Storage>::kStorageName,
        NaiveEchoGui<echo_demo::ImaInt16Storage>::kStorageDescription,
        [](float sample_rate) { s_ima16.Init(sample_rate); },
        [](const float* in, float* out, size_t n) noexcept { s_ima16.Process(in, out, n); },
        [] { s_ima16.CreateKnobs(); },
        [] { s_ima16.DisplayKnobs(); },
    },
}};

// 当前选中的存储格式下标，UI 线程写入、音频线程读取
static std::atomic<int> s_active{0};

static constexpr int kWindowWidth = NaiveEchoGui<echo_demo::FloatStorage>::kWindowWidth;
static constexpr int kWindowHeight = NaiveEchoGui<echo_demo::FloatStorage>::kWindowHeight;

// ------------------------------------------------------------
//  PeakMeter — 单声道输入/输出的峰值电平
// ------------------------------------------------------------
struct PeakMeter {
    std::atomic<float> raw_in{0.0f};
    std::atomic<float> raw_out{0.0f};

    void MeasureInput(const float* input, size_t frame_count) noexcept {
        raw_in.store(Peak(input, frame_count), std::memory_order_relaxed);
    }

    void MeasureOutput(const float* output, size_t frame_count) noexcept {
        raw_out.store(Peak(output, frame_count), std::memory_order_relaxed);
    }

    // --------------------------------------------------------
    //  Update — UI 线程每帧调用，instant attack + 指数 hold 衰减
    // --------------------------------------------------------
    void Update(float decay) noexcept {
        UpdateOne(disp_in_, raw_in, decay);
        UpdateOne(disp_out_, raw_out, decay);
    }

    // --------------------------------------------------------
    //  Draw — 绘制输入/输出两条电平条
    // --------------------------------------------------------
    void Draw(float x, float y, float total_h, float bar_w, float gap) const {
        auto bar_color = [](float norm) -> Color {
            if (norm < 0.667f) {
                float const t = norm / 0.667f;
                return Color{static_cast<unsigned char>(t * 255.0f), 255, 0, 255};
            }
            float const t = (norm - 0.667f) / 0.333f;
            return Color{255, static_cast<unsigned char>((1.0f - t) * 255.0f), 0, 255};
        };

        auto draw_one = [&](float bx, float norm, const char* label) {
            DrawRectangle(static_cast<int>(bx), static_cast<int>(y), static_cast<int>(bar_w),
                          static_cast<int>(total_h), Color{25, 25, 30, 255});
            DrawRectangleLines(static_cast<int>(bx), static_cast<int>(y), static_cast<int>(bar_w),
                               static_cast<int>(total_h), Color{60, 60, 70, 255});
            if (norm > 0.001f) {
                int const fill_h = static_cast<int>(norm * total_h);
                DrawRectangle(static_cast<int>(bx) + 1, static_cast<int>(y + total_h - fill_h),
                              static_cast<int>(bar_w) - 2, fill_h, bar_color(norm));
            }
            DrawText(label, static_cast<int>(bx), static_cast<int>(y + total_h + 4), 9, Color{140, 140, 150, 255});
        };

        draw_one(x, ToNorm(disp_in_), "in");
        draw_one(x + bar_w + gap, ToNorm(disp_out_), "out");
    }
private:
    float disp_in_ = 0.0f;
    float disp_out_ = 0.0f;

    static float Peak(const float* buffer, size_t frame_count) noexcept {
        float peak = 0.0f;
        for (size_t i = 0; i < frame_count; ++i) {
            peak = std::max(peak, std::abs(buffer[i]));
        }
        return peak;
    }

    static void UpdateOne(float& disp, std::atomic<float>& raw, float decay) noexcept {
        float const r = raw.load(std::memory_order_relaxed);
        disp = (r > disp) ? r : disp * decay;
    }

    /// @brief 线性幅度转 [-60dB, 0dB] 归一化显示值
    static float ToNorm(float amp) noexcept {
        float const db = 20.0f * std::log10(std::max(amp, 1e-10f));
        return std::clamp((db + 60.0f) / 60.0f, 0.0f, 1.0f);
    }
};

static PeakMeter s_peak_meter;

// ------------------------------------------------------------
//  顶部格式按钮
// ------------------------------------------------------------
namespace algorithm_bar {
// 占用窗口顶部空白区：旋钮自 y=42 起、电平表自 x=298 起，按钮行均不与两者重叠。
// 格式名长短不一，按钮宽度按文字实测宽度自适应，整行从 x 起依次排开。
static constexpr float kX = 14.0f;
static constexpr float kY = 6.0f;
static constexpr float kHeight = 26.0f;
static constexpr float kGap = 6.0f;
static constexpr float kPadding = 6.0f;
static constexpr int kFontSize = 10;

/// @brief 第 index 个按钮的宽度（文字宽度 + 两侧留白）
static int WidthOf(size_t index) {
    return MeasureText(kAlgorithms[index].name, kFontSize) + static_cast<int>(2 * kPadding);
}

/// @brief 第 index 个按钮的矩形
static Rectangle Bounds(size_t index) {
    float x = kX;
    for (size_t i = 0; i < index; ++i) {
        x += static_cast<float>(WidthOf(i)) + kGap;
    }
    return Rectangle{x, kY, static_cast<float>(WidthOf(index)), kHeight};
}

/// @brief 返回鼠标位置命中的格式下标，未命中返回 -1
static int HitTest(Vector2 mouse) {
    for (size_t i = 0; i < kAlgorithms.size(); ++i) {
        if (CheckCollisionPointRec(mouse, Bounds(i))) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

/// @brief 绘制按钮行，选中项实心高亮、悬停项加亮描边
static void Draw(int active, Vector2 mouse) {
    for (size_t i = 0; i < kAlgorithms.size(); ++i) {
        Rectangle const bounds = Bounds(i);
        bool const hovered = CheckCollisionPointRec(mouse, bounds);
        bool const selected = static_cast<int>(i) == active;
        Color const border = hovered ? Color{190, 194, 198, 255} : Color{75, 80, 85, 255};
        Color const text = selected ? Color{18, 20, 23, 255} : Color{178, 182, 187, 255};

        if (selected) {
            DrawRectangleRec(bounds, Color{232, 233, 235, 255});
        }
        else {
            DrawRectangleLinesEx(bounds, 1.0f, border);
        }

        char const* name = kAlgorithms[i].name;
        int const text_width = MeasureText(name, kFontSize);
        DrawText(name, static_cast<int>(bounds.x + (bounds.width - static_cast<float>(text_width)) * 0.5f),
                 static_cast<int>(bounds.y + (bounds.height - static_cast<float>(kFontSize)) * 0.5f), kFontSize, text);
    }
}

/// @brief 把窗口标题同步为当前存储格式
static void SyncWindowTitle(int active) {
    SetWindowTitle(std::format("Naive Echo — {}", kAlgorithms[static_cast<size_t>(active)].name).c_str());
}
} // namespace algorithm_bar

// ------------------------------------------------------------
//  miniaudio 回调
// ------------------------------------------------------------
extern "C" void MaCallback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;

    auto* input = static_cast<const float*>(pInput);
    auto* output = static_cast<float*>(pOutput);

    if (input == nullptr || output == nullptr)
        return;

    constexpr size_t kChunk = 512;
    float buffer[kChunk];
    size_t remaining = static_cast<size_t>(frameCount);
    size_t offset = 0;

    while (remaining > 0) {
        size_t const nf = std::min(remaining, kChunk);

        // 处理前：测量输入峰值
        s_peak_meter.MeasureInput(input + offset, nf);

        // DSP 处理（当前选中的存储格式）
        int const active = s_active.load(std::memory_order_acquire);
        kAlgorithms[static_cast<size_t>(active)].Process(input + offset, buffer, nf);

        // 处理后：测量输出峰值
        s_peak_meter.MeasureOutput(buffer, nf);

        std::copy_n(buffer, nf, output + offset);

        offset += nf;
        remaining -= nf;
    }
}

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------
int main(void) {
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, "Naive Echo");
    SetTargetFPS(60);

    static constexpr float kSampleRate = 48000.0f;

    // ── 初始化全部格式（切换时不重建，各自保留参数与延迟线内容）──
    for (auto const& algo : kAlgorithms) {
        algo.Init(kSampleRate);
        algo.CreateKnobs();
    }
    algorithm_bar::SyncWindowTitle(0);

    // ── miniaudio full-duplex（单声道入 → 单声道出）──
    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.playback.format = ma_format_f32;
    config.playback.channels = 1;
    config.sampleRate = static_cast<ma_uint32>(kSampleRate);
    config.dataCallback = MaCallback;
    config.pUserData = nullptr;
    config.periodSizeInMilliseconds = 10;

    ma_device device;
    ma_result result = ma_device_init(nullptr, &config, &device);
    bool audio_ok = (result == MA_SUCCESS);
    if (audio_ok) {
        ma_device_start(&device);
    }

    // ── main loop ──
    while (!WindowShouldClose()) {
        Vector2 const mouse = GetMousePosition();

        // ── 顶部按钮：点击切换存储格式 ──
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            int const hit = algorithm_bar::HitTest(mouse);
            if (hit >= 0 && hit != s_active.load(std::memory_order_relaxed)) {
                s_active.store(hit, std::memory_order_relaxed);
                algorithm_bar::SyncWindowTitle(hit);
            }
        }
        int const active = s_active.load(std::memory_order_relaxed);

        BeginDrawing();
        ClearBackground(Color{10, 10, 12, 255});

        // ── 按钮行 ──
        algorithm_bar::Draw(active, mouse);

        // ── 当前格式的旋钮 ──
        kAlgorithms[static_cast<size_t>(active)].DisplayKnobs();

        // ── 当前格式说明 ──
        DrawText(kAlgorithms[static_cast<size_t>(active)].description, 14, 120, 10, Color{150, 155, 165, 255});

        // ── 峰值电平表 ──
        {
            constexpr float kDecay = 0.88f;
            constexpr float kBarW = 18.0f;
            constexpr float kGap = 12.0f;
            constexpr float kBarX = static_cast<float>(kWindowWidth) - 2.0f * kBarW - kGap - 14.0f;
            constexpr float kBarY = 14.0f;
            constexpr float kBarH = 180.0f;

            s_peak_meter.Update(kDecay);
            s_peak_meter.Draw(kBarX, kBarY, kBarH, kBarW, kGap);
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
