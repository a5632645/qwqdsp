#pragma once

#include <atomic>
#include <cstring>
#include <span>
#include <vector>

#include "raylib.h"

/**
 * @brief 从右侧接收新像素列，内部维护滚动缓冲区 + raylib 纹理
 *
 * PushColumn 由音频线程调用，Draw 由渲染线程调用。
 */
struct ScrollingImage {
    void Init(int width, int height) noexcept {
        width_ = width;
        height_ = height;
        buf_.resize(width_ * height_);
        std::fill(buf_.begin(), buf_.end(), Color{0, 0, 0, 255});

        Image img = GenImageColor(width_, height_, Color{0, 0, 0, 255});
        tex_ = LoadTextureFromImage(img);
        UnloadImage(img);
    }

    /** @brief 从右侧推入一列新像素，图像左移一列 */
    void PushColumn(std::span<const Color> column) noexcept {
        for (int y = 0; y < height_; ++y) {
            std::memmove(&buf_[y * width_], &buf_[y * width_ + 1], (width_ - 1) * sizeof(Color));
            buf_[y * width_ + width_ - 1] = column[y];
        }
        dirty_.store(true, std::memory_order_release);
    }

    /** @brief 渲染线程：必要时更新纹理并绘制到画布，自动缩放填充 targetSize */
    void Draw(int x, int y, int targetWidth, int targetHeight) noexcept {
        if (dirty_.load(std::memory_order_acquire)) {
            UpdateTexture(tex_, buf_.data());
            dirty_.store(false, std::memory_order_release);
        }
        DrawTexturePro(tex_, {0, 0, static_cast<float>(width_), static_cast<float>(height_)},
                       {static_cast<float>(x), static_cast<float>(y), static_cast<float>(targetWidth),
                        static_cast<float>(targetHeight)},
                       {0, 0}, 0.0f, WHITE);
    }

    void Unload() noexcept {
        UnloadTexture(tex_);
    }
private:
    int width_{}, height_{};
    std::vector<Color> buf_;
    Texture2D tex_{};
    std::atomic<bool> dirty_{false};
};
