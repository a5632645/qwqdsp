/// @brief 提供 GUI demo 的共享截图验证入口

#pragma once

#include "raylib.h"

namespace qwqdsp_gui_test {

/// @brief 截图模式；模板参数为 false 时所有操作都是内联空操作
template <bool kEnabled>
class ShotMode;

template <>
class ShotMode<false> {
public:
    explicit ShotMode(char const* = nullptr, int = 10) noexcept {}

    void ApplyWindowFlags() const noexcept {}
    [[nodiscard]] bool IsEnabled() const noexcept { return false; }
    [[nodiscard]] char const* Path() const noexcept { return nullptr; }
    bool Tick() noexcept { return false; }
};

template <>
class ShotMode<true> {
public:
    /// @param path PNG 输出路径；nullptr 表示不截图
    /// @param capture_frame 第几帧导出 PNG
    explicit ShotMode(char const* path = nullptr, int capture_frame = 10) noexcept
        : path_(path), capture_frame_(capture_frame) {}

    /// @brief 在 InitWindow 之前调用，截图模式下隐藏窗口
    void ApplyWindowFlags() const noexcept {
        if (path_ != nullptr) {
            SetConfigFlags(FLAG_WINDOW_HIDDEN);
        }
    }

    [[nodiscard]] bool IsEnabled() const noexcept {
        return path_ != nullptr;
    }

    [[nodiscard]] char const* Path() const noexcept {
        return path_;
    }

    /// @brief 每帧 EndDrawing 之后调用；到达截图帧时导出 PNG 并返回 true
    bool Tick() noexcept {
        if (path_ == nullptr) {
            return false;
        }
        if (++frame_ < capture_frame_) {
            return false;
        }

        Image image = LoadImageFromScreen();
        bool const ok = ExportImage(image, path_);
        UnloadImage(image);
        if (ok) {
            TraceLog(LOG_INFO, "GUI shot saved: %s", path_);
        }
        else {
            TraceLog(LOG_WARNING, "GUI shot failed: %s", path_);
        }
        return true;
    }

private:
    char const* path_;
    int capture_frame_;
    int frame_ = 0;
};

} // namespace qwqdsp_gui_test
