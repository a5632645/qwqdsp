#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>

#include "miniaudio.h"
#include "raylib.h"

#include "reverb/fdn_reverb_gui.hpp"
#include "reverb/gardner_reverb_gui.hpp"
#include "reverb/peak_meter.hpp"
#include "reverb/plate_reverb_gui.hpp"
#include "reverb/shimmer_reverb_gui.hpp"

// ------------------------------------------------------------
//  ReverbGuiConcept — 混响 GUI 包装器必须满足的接口
//
// 任何用作算法槽位的类型都需要提供以下成员：
//   - static constexpr int kWindowWidth
//   - static constexpr int kWindowHeight
//   - void Init(float sample_rate)
//   - void Process(const float* input, float* output_l,
//                  float* output_r, size_t frame_count) noexcept
//   - void CreateKnobs()
//   - void DisplayKnobs()
// ------------------------------------------------------------
template <typename T>
concept ReverbGuiConcept = requires(T& obj, const float* in, float* ol, float* or_, size_t n, float sr) {
    // ── 静态编译期常量 ──
    { T::kWindowWidth } -> std::convertible_to<int>;
    { T::kWindowHeight } -> std::convertible_to<int>;

    // ── 方法签名 ──
    { obj.Init(sr) } -> std::same_as<void>;
    { obj.Process(in, ol, or_, n) } noexcept -> std::same_as<void>;
    { obj.CreateKnobs() } -> std::same_as<void>;
    { obj.DisplayKnobs() } -> std::same_as<void>;
};

// ------------------------------------------------------------
//  算法槽位 — 新增一种混响：加一个全局实例，再往 kAlgorithms 里加一项
// ------------------------------------------------------------
static ShimmerReverbGui s_shimmer;
static PlateReverbGui s_plate;
static GardnerRoomGui s_gardner_room;
static GardnerHallGui s_gardner_hall;
static Fdn4SmallGui s_fdn_small;
static Fdn4MediumGui s_fdn_medium;
static Fdn4LargeGui s_fdn_large;

static_assert(ReverbGuiConcept<ShimmerReverbGui>,
              "ShimmerReverbGui 必须满足 ReverbGuiConcept：需要 kWindowWidth/kWindowHeight "
              "静态常量，以及 Init/Process/CreateKnobs/DisplayKnobs 方法");
static_assert(ReverbGuiConcept<PlateReverbGui>,
              "PlateReverbGui 必须满足 ReverbGuiConcept：需要 kWindowWidth/kWindowHeight "
              "静态常量，以及 Init/Process/CreateKnobs/DisplayKnobs 方法");
static_assert(ReverbGuiConcept<GardnerRoomGui>,
              "GardnerRoomGui 必须满足 ReverbGuiConcept：需要 kWindowWidth/kWindowHeight "
              "静态常量，以及 Init/Process/CreateKnobs/DisplayKnobs 方法");
static_assert(ReverbGuiConcept<GardnerHallGui>,
              "GardnerHallGui 必须满足 ReverbGuiConcept：需要 kWindowWidth/kWindowHeight "
              "静态常量，以及 Init/Process/CreateKnobs/DisplayKnobs 方法");
static_assert(ReverbGuiConcept<Fdn4SmallGui>,
              "Fdn4SmallGui 必须满足 ReverbGuiConcept：需要 kWindowWidth/kWindowHeight "
              "静态常量，以及 Init/Process/CreateKnobs/DisplayKnobs 方法");
static_assert(ReverbGuiConcept<Fdn4MediumGui>,
              "Fdn4MediumGui 必须满足 ReverbGuiConcept：需要 kWindowWidth/kWindowHeight "
              "静态常量，以及 Init/Process/CreateKnobs/DisplayKnobs 方法");
static_assert(ReverbGuiConcept<Fdn4LargeGui>,
              "Fdn4LargeGui 必须满足 ReverbGuiConcept：需要 kWindowWidth/kWindowHeight "
              "静态常量，以及 Init/Process/CreateKnobs/DisplayKnobs 方法");

struct ReverbAlgorithm {
    const char* name;
    void (*Init)(float sample_rate);
    void (*Process)(const float* input, float* output_l, float* output_r, size_t frame_count) noexcept;
    void (*CreateKnobs)();
    void (*DisplayKnobs)();
};

static const std::array<ReverbAlgorithm, 7> kAlgorithms{{
    {
        "Shimmer",
        [](float sample_rate) { s_shimmer.Init(sample_rate); },
        [](const float* in, float* l, float* r, size_t n) noexcept { s_shimmer.Process(in, l, r, n); },
        [] { s_shimmer.CreateKnobs(); },
        [] { s_shimmer.DisplayKnobs(); },
    },
    {
        "Plate",
        [](float sample_rate) { s_plate.Init(sample_rate); },
        [](const float* in, float* l, float* r, size_t n) noexcept { s_plate.Process(in, l, r, n); },
        [] { s_plate.CreateKnobs(); },
        [] { s_plate.DisplayKnobs(); },
    },
    {
        "Gardner Room",
        [](float sample_rate) { s_gardner_room.Init(sample_rate); },
        [](const float* in, float* l, float* r, size_t n) noexcept { s_gardner_room.Process(in, l, r, n); },
        [] { s_gardner_room.CreateKnobs(); },
        [] { s_gardner_room.DisplayKnobs(); },
    },
    {
        "Gardner Hall",
        [](float sample_rate) { s_gardner_hall.Init(sample_rate); },
        [](const float* in, float* l, float* r, size_t n) noexcept { s_gardner_hall.Process(in, l, r, n); },
        [] { s_gardner_hall.CreateKnobs(); },
        [] { s_gardner_hall.DisplayKnobs(); },
    },
    {
        "FDN Small",
        [](float sample_rate) { s_fdn_small.Init(sample_rate); },
        [](const float* in, float* l, float* r, size_t n) noexcept { s_fdn_small.Process(in, l, r, n); },
        [] { s_fdn_small.CreateKnobs(); },
        [] { s_fdn_small.DisplayKnobs(); },
    },
    {
        "FDN Medium",
        [](float sample_rate) { s_fdn_medium.Init(sample_rate); },
        [](const float* in, float* l, float* r, size_t n) noexcept { s_fdn_medium.Process(in, l, r, n); },
        [] { s_fdn_medium.CreateKnobs(); },
        [] { s_fdn_medium.DisplayKnobs(); },
    },
    {
        "FDN Large",
        [](float sample_rate) { s_fdn_large.Init(sample_rate); },
        [](const float* in, float* l, float* r, size_t n) noexcept { s_fdn_large.Process(in, l, r, n); },
        [] { s_fdn_large.CreateKnobs(); },
        [] { s_fdn_large.DisplayKnobs(); },
    },
}};

// 当前选中的算法下标，UI 线程写入、音频线程读取
static std::atomic<int> s_active{0};

// 窗口取各算法布局尺寸的最大值
static constexpr int kWindowWidth = std::max(ShimmerReverbGui::kWindowWidth, PlateReverbGui::kWindowWidth);
static constexpr int kWindowHeight = std::max(ShimmerReverbGui::kWindowHeight, PlateReverbGui::kWindowHeight);

static PeakMeter s_peak_meter;

// ------------------------------------------------------------
//  顶部算法按钮
// ------------------------------------------------------------
namespace algorithm_bar {
// 占用窗口顶部空白区：旋钮自 y=42 起、电平表自 x=514 起，按钮行均不与两者重叠。
// 算法名长短不一，按钮宽度按文字实测宽度自适应，整行从 x 起依次排开。
static constexpr float kX = 14.0f;
static constexpr float kY = 6.0f;
static constexpr float kHeight = 26.0f;
static constexpr float kGap = 6.0f;
static constexpr float kPadding = 6.0f;
static constexpr int kFontSize = 10;

/// @brief 第 index 个按钮的宽度（文字宽度 + 两侧留白）
static float WidthOf(size_t index) {
    return static_cast<float>(MeasureText(kAlgorithms[index].name, kFontSize)) + 2 * kPadding;
}

/// @brief 第 index 个按钮的矩形
static Rectangle Bounds(size_t index) {
    float x = kX;
    for (size_t i = 0; i < index; ++i) {
        x += WidthOf(i) + kGap;
    }
    return Rectangle{x, kY, WidthOf(index), kHeight};
}

/// @brief 返回鼠标位置命中的算法下标，未命中返回 -1
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

/// @brief 把窗口标题同步为当前算法名
static void SyncWindowTitle(int active) {
    SetWindowTitle(kAlgorithms[static_cast<size_t>(active)].name);
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
    float buf_l[kChunk];
    float buf_r[kChunk];
    size_t remaining = static_cast<size_t>(frameCount);
    size_t offset = 0;

    while (remaining > 0) {
        size_t nf = std::min(remaining, kChunk);

        // 处理前：测量输入峰值
        s_peak_meter.MeasureInput(input + offset, nf);

        // DSP 处理（当前选中的算法）
        int const active = s_active.load(std::memory_order_acquire);
        kAlgorithms[static_cast<size_t>(active)].Process(input + offset, buf_l, buf_r, nf);

        // 处理后：测量输出峰值
        s_peak_meter.MeasureOutput(buf_l, buf_r, nf);

        for (size_t i = 0; i < nf; ++i) {
            output[2 * (offset + i)] = buf_l[i];
            output[2 * (offset + i) + 1] = buf_r[i];
        }

        offset += nf;
        remaining -= nf;
    }
}

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------
int main(void) {
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, kAlgorithms[0].name);
    SetTargetFPS(60);

    static constexpr float kSampleRate = 48000.0f;

    // ── 初始化全部算法（切换时不重建，各自保留参数与内部状态）──
    for (auto const& algo : kAlgorithms) {
        algo.Init(kSampleRate);
        algo.CreateKnobs();
    }

    // ── miniaudio full-duplex（立体声入 → 立体声出）──
    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 2;
    config.playback.format = ma_format_f32;
    config.playback.channels = 2;
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

    // ── main loop ──
    while (!WindowShouldClose()) {
        Vector2 const mouse = GetMousePosition();

        // ── 顶部按钮：点击切换算法 ──
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

        // ── 当前算法的旋钮 ──
        kAlgorithms[static_cast<size_t>(active)].DisplayKnobs();

        // ── 峰值电平表（hold 衰减更新 + 绘制）──
        {
            constexpr float kDecay = 0.88f;
            constexpr float kBarW = 18.0f;
            constexpr float kGap = 4.0f;
            constexpr float kGroupGap = 12.0f;
            // 所有电平条放在右侧：4 条 ×18 + 3 个间距(4+12+4) = 92px 偏移
            constexpr float kBarX = static_cast<float>(kWindowWidth) - 92.0f - 14.0f;
            constexpr float kBarY = 14.0f;
            constexpr float kBarH = 252.0f;

            s_peak_meter.Update(kDecay);
            s_peak_meter.Draw(kBarX, kBarY, kBarH, kBarW, kGap, kGroupGap);
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
