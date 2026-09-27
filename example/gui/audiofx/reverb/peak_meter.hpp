#pragma once

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>

#include "raylib.h"

// ------------------------------------------------------------
//  PeakMeter — 独立峰值电平表
//
//  音频线程：
//    MeasureInput()  处理前调用，扫描输入峰值
//    MeasureOutput() 处理后调用，扫描输出峰值
//
//  UI 线程：
//    Update(decay)   每帧调用，带 hold 衰减
//    Draw(...)       绘制电平条
// ------------------------------------------------------------
struct PeakMeter {
    // ── 原子原始峰值（音频线程写入，UI 线程读取）──
    std::atomic<float> raw_in_l{0.0f};
    std::atomic<float> raw_in_r{0.0f};
    std::atomic<float> raw_out_l{0.0f};
    std::atomic<float> raw_out_r{0.0f};

    // --------------------------------------------------------
    //  MeasureInput — 处理前：扫描交错立体声输入缓冲的峰值
    // --------------------------------------------------------
    void MeasureInput(const float* interleaved, size_t n) noexcept {
        float p_l = 0.0f, p_r = 0.0f;
        for (size_t i = 0; i < n; ++i) {
            float in_l = interleaved[2 * i];
            float in_r = interleaved[2 * i + 1];
            p_l = std::max(p_l, std::abs(in_l));
            p_r = std::max(p_r, std::abs(in_r));
        }
        raw_in_l.store(p_l, std::memory_order_relaxed);
        raw_in_r.store(p_r, std::memory_order_relaxed);
    }

    // --------------------------------------------------------
    //  MeasureOutput — 处理后：扫描分离立体声输出缓冲的峰值
    // --------------------------------------------------------
    void MeasureOutput(const float* l, const float* r, size_t n) noexcept {
        float p_l = 0.0f, p_r = 0.0f;
        for (size_t i = 0; i < n; ++i) {
            p_l = std::max(p_l, std::abs(l[i]));
            p_r = std::max(p_r, std::abs(r[i]));
        }
        raw_out_l.store(p_l, std::memory_order_relaxed);
        raw_out_r.store(p_r, std::memory_order_relaxed);
    }

    // --------------------------------------------------------
    //  Update — UI 线程每帧调用，带 hold 衰减更新显示值
    // --------------------------------------------------------
    void Update(float decay) noexcept {
        auto upd = [decay](float& disp, std::atomic<float>& raw) {
            float r = raw.load(std::memory_order_relaxed);
            disp = (r > disp) ? r : disp * decay;
        };
        upd(disp_in_l_, raw_in_l);
        upd(disp_in_r_, raw_in_r);
        upd(disp_out_l_, raw_out_l);
        upd(disp_out_r_, raw_out_r);
    }

    // --------------------------------------------------------
    //  归一化 dB 值获取 [0, 1]（-60dB ~ 0dB）
    // --------------------------------------------------------
    float NormInL() const noexcept {
        return ToNorm(disp_in_l_);
    }
    float NormInR() const noexcept {
        return ToNorm(disp_in_r_);
    }
    float NormOutL() const noexcept {
        return ToNorm(disp_out_l_);
    }
    float NormOutR() const noexcept {
        return ToNorm(disp_out_r_);
    }

    // --------------------------------------------------------
    //  Draw — 绘制 4 条电平条
    //
    //  布局： L in | R in  ← group_gap →  L out | R out
    //  所有条都位于窗口右侧。
    //  每条的宽度 bar_w，组内间距 gap，in/out 组间距 group_gap。
    // --------------------------------------------------------
    void Draw(float x, float y, float total_h, float bar_w, float gap, float group_gap) const {
        // ── 颜色映射（绿→黄→红）──
        auto bar_color = [](float norm) -> Color {
            if (norm < 0.667f) {
                float t = norm / 0.667f;
                return Color{(unsigned char)(t * 255.0f), 255, 0, 255};
            }
            else {
                float t = (norm - 0.667f) / 0.333f;
                return Color{255, (unsigned char)((1.0f - t) * 255.0f), 0, 255};
            }
        };

        // ── 绘制单条 ──
        auto draw_one = [&](float bx, float norm, const char* label) {
            DrawRectangle(static_cast<int>(bx), static_cast<int>(y), static_cast<int>(bar_w), static_cast<int>(total_h),
                          Color{25, 25, 30, 255});
            DrawRectangleLines(static_cast<int>(bx), static_cast<int>(y), static_cast<int>(bar_w),
                               static_cast<int>(total_h), Color{60, 60, 70, 255});
            if (norm > 0.001f) {
                int fill_h = static_cast<int>(norm * total_h);
                DrawRectangle(static_cast<int>(bx) + 1, static_cast<int>(y + total_h - fill_h),
                              static_cast<int>(bar_w) - 2, fill_h, bar_color(norm));
            }
            DrawText(label, static_cast<int>(bx), static_cast<int>(y + total_h + 4), 9, Color{140, 140, 150, 255});
        };

        float norms[4] = {NormInL(), NormInR(), NormOutL(), NormOutR()};
        const char* labels[4] = {"L in", "R in", "L out", "R out"};

        // 组内间距序列
        float gaps[3] = {gap, group_gap, gap};
        float bx = x;
        for (int i = 0; i < 4; ++i) {
            draw_one(bx, norms[i], labels[i]);
            if (i < 3)
                bx += bar_w + gaps[i];
        }
    }
private:
    float disp_in_l_ = 0.0f;
    float disp_in_r_ = 0.0f;
    float disp_out_l_ = 0.0f;
    float disp_out_r_ = 0.0f;

    static float ToNorm(float amp) noexcept {
        float db = 20.0f * std::log10(std::max(amp, 1e-10f));
        return std::clamp((db + 60.0f) / 60.0f, 0.0f, 1.0f);
    }
};
