#pragma once

#include <array>
#include <span>
#include <vector>

#include <qwqdsp/window/blackman.hpp>
#include <qwqdsp/window/blackman_harris.hpp>
#include <qwqdsp/window/blackman_harris_3term.hpp>
#include <qwqdsp/window/hamming.hpp>
#include <qwqdsp/window/hann.hpp>
#include <qwqdsp/window/helper.hpp>

#include "raylib.h"

/**
 * @brief 输入音频块，内部管理音频环回缓冲和加窗，每产出新列立即回调 push
 *
 * FrameProcessor 需提供:
 *   void Process(raw_frame, window, windowed_frame) noexcept;
 *   std::span<const Color> GetColumn() const noexcept;
 */
struct SpectrogramColumn {
    void Init(int height, int sampleRate, int fftSize, int hopSize) noexcept {
        height_ = height;
        fft_size_ = fftSize;
        hop_size_ = hopSize;

        window_.resize(fft_size_);
        qwqdsp_window::BlackmanHarrisThreeTerm::Window(window_, true);
        qwqdsp_window::Helper::Normalize(window_);

        in_buffer_.resize(fft_size_, 0.0f);
        fft_in_.resize(fft_size_);
    }

    template <typename FrameProcessor, typename PushFunc>
    void ProcessAudio(std::span<const float> audio, FrameProcessor&& proc, PushFunc&& push) noexcept {
        const float* src = audio.data();
        int num = static_cast<int>(audio.size());

        while (num != 0) {
            int need = fft_size_ - in_count_;
            int can_do = std::min(need, num);
            std::copy_n(src, can_do, in_buffer_.begin() + in_count_);

            in_count_ += can_do;
            num -= can_do;
            src += can_do;

            if (in_count_ == fft_size_) {
                // 加窗后委托给 FrameProcessor (传入原始帧 + 窗口 + 加窗帧)
                for (int i = 0; i < fft_size_; ++i)
                    fft_in_[i] = window_[i] * in_buffer_[i];
                proc.Process(in_buffer_, window_, fft_in_);

                // 环回平移
                for (int i = 0; i < (fft_size_ - hop_size_); ++i)
                    in_buffer_[i] = in_buffer_[i + hop_size_];
                in_count_ -= hop_size_;

                // 每产出一列立即推送
                push(proc.GetColumn());
            }
        }
    }

    int ColumnHeight() const noexcept {
        return height_;
    }
private:
    std::vector<float> window_;
    std::vector<float> in_buffer_;
    std::vector<float> fft_in_;
    int in_count_{};
    int fft_size_{}, hop_size_{}, height_{};
};
