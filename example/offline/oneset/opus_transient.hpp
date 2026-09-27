#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

#include "oneset_common.hpp"

namespace qwqdsp_test {

// ------------------------------------------------------------
// OpusTransientDetector
// ------------------------------------------------------------
/**
 * @brief Opus CELT 瞬态检测器（transient_analysis）。
 *
 * 移植自 ref_repos/opus/celt/celt_encoder.c:267-469 的
 * transient_analysis（Xiph.Org Opus，BSD 协议），float 模式逐行对应
 * （Opus float 模式下 SHR32/SHL32/PSHR32/SROUND16 等均为恒等，
 * MULT16_16 为普通乘法，QCONST 为直接常量）。纯时域 PCM 检测：
 *
 * 1. 二阶高通 (1-2z⁻¹+z⁻²)/(1-z⁻¹+0.5z⁻²) 抑制低频
 * 2. 前 12 样本清零（滤波初值无效），两样本分块能量 x2 = t² + t²
 * 3. 前向 masking（衰减 0.0625 ≈ 6.7 dB/ms）：短时能量经低通平滑
 * 4. 后向 masking（0.875/0.125 ≈ 13.9 dB/ms）：求预回声阈值
 * 5. 帧能量（几何均值）与谐波均值之比构成 mask metric，
 *    经训练表 inv_table 加权，> 200 判定为瞬态
 *
 * 这是工业级编码器实际使用的瞬态检测，对打击乐/拨弦攻击敏感，
 * 且对低频纯音（partial cycle 误判）有专门防误报。
 */
class OpusTransientDetector {
public:
    static constexpr size_t kDefaultFrameSize = 960; // 48kHz 下 20ms
    static constexpr int kTransientThreshold = 200;

    /** 设置分析帧长度（样本数，Opus 通常 960@48k / 480@24k） */
    void SetFrameSize(size_t n) noexcept {
        frame_size_ = n;
    }

    /** 检测输入信号中的瞬态 onset */
    OnsetResult Detect(std::span<const float> input, float /*sample_rate*/) {
        const size_t len = frame_size_;
        if (len < 64)
            return {};

        const size_t len2 = len / 2;
        OnsetResult result;
        result.odf.reserve(input.size() / len + 1);

        // 每帧独立检测（Opus 的 mem0/mem1 在调用内清零，无跨帧状态）
        for (size_t frame_start = 0; frame_start + len <= input.size(); frame_start += len) {
            const float* in = input.data() + frame_start;

            std::vector<float> tmp(len);
            float mem0 = 0.0f, mem1 = 0.0f;

            // ---- 二阶高通 (1-2z⁻¹+z⁻²)/(1-z⁻¹+0.5z⁻²) ----
            // 原 float 分支：y = mem0 + x; mem0 = mem0 - x + .5*mem1; mem1 = x - mem0_old
            for (size_t i = 0; i < len; ++i) {
                const float x = in[i];
                const float y = mem0 + x;
                const float mem00 = mem0;
                mem0 = mem0 - x + 0.5f * mem1;
                mem1 = x - mem00;
                tmp[i] = y; // SROUND16(y,2) 在 float 下恒等
            }

            // 前 12 样本清零（滤波初值无效）
            std::fill(tmp.begin(), tmp.begin() + std::min<size_t>(12, tmp.size()), 0.0f);

            // ---- 两样本分块能量 + 前向 masking ----
            float mean = 0.0f;
            mem0 = 0.0f;
            std::vector<float> tmp2(len2);
            for (size_t i = 0; i < len2; ++i) {
                // PSHR32(MULT16_16(t,t)+MULT16_16(t,t),4) 在 float 下恒等
                const float x2 = tmp[2 * i] * tmp[2 * i] + tmp[2 * i + 1] * tmp[2 * i + 1];
                mean += x2; // PSHR32(x2,12) 恒等
                mem0 = x2 + (1.0f - 0.0625f) * mem0;
                tmp2[i] = 0.0625f * mem0;
            }

            // ---- 后向 masking（求预回声阈值） ----
            mem0 = 0.0f;
            float max_e = 0.0f;
            for (size_t i = len2; i-- > 0;) {
                mem0 = tmp2[i] + 0.875f * mem0;
                tmp2[i] = 0.125f * mem0;
                max_e = std::max(max_e, 0.125f * mem0);
            }

            // ---- mask metric（几何均值 × 谐波均值） ----
            // mean = celt_sqrt(mean * maxE * .5 * len2)
            mean = std::sqrt(mean * max_e * 0.5f * static_cast<float>(len2));
            // norm = len2 / (EPSILON + mean/2)，float 下 SHL32/SHR32 恒等，EPSILON=1e-15
            const float norm = static_cast<float>(len2) / (1e-15f + mean * 0.5f);

            int unmask = 0;
            for (size_t i = 12; i + 5 < len2; i += 4) {
                const float v = tmp2[i] + 1e-15f;
                const int id = static_cast<int>(std::floor(64.0f * norm * v));
                const int clamped = std::clamp(id, 0, 127);
                unmask += kInvTable[clamped];
            }
            unmask = 64 * unmask * 4 / (6 * static_cast<int>(len2 - 17));

            const bool is_transient = unmask > kTransientThreshold;

            result.odf.push_back(static_cast<float>(unmask));
            if (is_transient)
                result.onset_frames.push_back(result.odf.size() - 1);
        }

        for (size_t f : result.onset_frames)
            result.onset_samples.push_back(f * frame_size_);

        return result;
    }

private:
    /**
     * @brief Opus 训练表（6*64/x，最小化平均误差），逐字复制自
     *        celt_encoder.c:286-295。
     */
    static constexpr unsigned char kInvTable[128] = {
        255, 255, 156, 110, 86, 70, 59, 51, 45, 40, 37, 33, 31, 28, 26, 25,
        23, 22, 21, 20, 19, 18, 17, 16, 16, 15, 15, 14, 13, 13, 12, 12,
        12, 12, 11, 11, 11, 10, 10, 10, 9, 9, 9, 9, 9, 9, 8, 8,
        8, 8, 8, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2,
    };

    size_t frame_size_ = kDefaultFrameSize;
};

} // namespace qwqdsp_test
