#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

#include "oneset_common.hpp"

namespace qwqdsp_test {

// ------------------------------------------------------------
// FdkAttackDetector
// ------------------------------------------------------------
/**
 * @brief FDK-AAC 块切换攻击检测器（FDKaacEnc_BlockSwitching）。
 *
 * 移植自 ref_repos/fdk-aac/libAACenc/src/block_switch.cpp:229-409
 * （Fraunhofer FDK AAC，自定义许可）。纯时域检测，AAC 编码器用它决定
 * 长窗/短窗切换：
 *
 * 1. 每帧（granule，默认 1024）分 8 个子窗（每窗 128 样本）
 * 2. 每子窗：二阶 IIR 高通（hiPassCoeff={-0.5095,0.7548}）→
 *    未滤波平方和 windowNrg 与滤波平方和 windowNrgF
 * 3. 跨窗平滑历史能量：accWindowNrg = 0.7·accWindowNrg + 0.3·前一窗
 * 4. attack：当前滤波能量 × 0.1 > 历史能量（即能量比 > 10）
 *    且最大能量超过绝对门限 minAttackNrg
 * 5. 跨帧延续：上一帧 attack 且上一帧最后一窗仍远大于当前首窗
 *
 * 状态（IIR、accWindowNrg、上一帧子窗能量）跨帧保持，与编码器一致。
 */
class FdkAttackDetector {
public:
    static constexpr size_t kDefaultGranule = 1024; // AAC 帧长
    static constexpr size_t kNumWindows = 8;        // 子窗数
    /** 跨窗历史能量系数：accWindowNrg = 0.7·acc + 0.3·prev */
    static constexpr float kAccWindowNrgFac = 0.3f;
    static constexpr float kOneMinusAccWindowNrgFac = 0.7f;
    /** 攻击比值下限的倒数（attackRatio = 10） */
    static constexpr float kInvAttackRatio = 0.1f;
    /**
     * @brief 最小攻击能量门限（float 尺度）。
     *
     * FDK 定点：minAttackNrg = (1e6·NORM_PCM_ENERGY) >> 7，其中
     * NORM_PCM_ENERGY = (1/32768)²，windowNrgF 累加含 /32 定标。
     * 换算到本实现（window_nrg_f = Σ iir_state1²，iir_state1 由
     * x/2 高通得到）的等效门限 ≈ 1e6·2⁻³² ≈ 2.33e-4。
     */
    static constexpr float kMinAttackNrg = 1e6f / (32768.0f * 32768.0f * 128.0f) * 32.0f;

    /** 设置帧长（granule，AAC 通常 1024） */
    void SetGranuleSize(size_t n) noexcept {
        granule_ = n;
    }

    /** 检测输入信号中的瞬态 onset */
    OnsetResult Detect(std::span<const float> input, float /*sample_rate*/) {
        const size_t granule = granule_;
        const size_t window_len = granule / kNumWindows;

        if (granule < kNumWindows * 4)
            return {};

        // 状态（跨帧保持，FDK 的 BLOCK_SWITCHING_CONTROL）
        float acc_window_nrg = 0.0f;
        float last_attack = 0;          // 上一帧是否 attack
        size_t last_attack_index = 0;   // 上一帧 attack 的子窗序号
        std::vector<float> prev_filtered(kNumWindows, 0.0f); // 上一帧滤波能量
        float iir_state0 = 0.0f, iir_state1 = 0.0f;          // IIR 状态

        OnsetResult result;
        result.odf.reserve(input.size() / granule + 1);

        // 每 granule 处理一帧（与 FDK 编码器一致的 hop）
        for (size_t frame_start = 0; frame_start + granule <= input.size(); frame_start += granule) {
            const float* pcm = input.data() + frame_start;

            std::vector<float> cur_filtered(kNumWindows);
            float en_max = 0.0f;

            // ---- 计算每子窗能量 ----
            for (size_t w = 0; w < kNumWindows; ++w) {
                float window_nrg = 0.0f;
                float window_nrg_f = 0.0f;
                for (size_t i = 0; i < window_len; ++i) {
                    const float x = pcm[w * window_len + i] * 0.5f; // tempUnfiltered = x>>1
                    // 二阶 IIR 高通：t1 = .5*c1*(x-s0); t2 = .5*c0*s1; s0=x; s1=(t1-t2)*2
                    const float t1 = 0.5f * kHiPassCoeff1 * (x - iir_state0);
                    const float t2 = 0.5f * kHiPassCoeff0 * iir_state1;
                    iir_state0 = x;
                    iir_state1 = (t1 - t2) * 2.0f;

                    window_nrg += iir_state0 * iir_state0;     // 未滤波能量
                    window_nrg_f += iir_state1 * iir_state1;   // 滤波能量
                }
                cur_filtered[w] = window_nrg_f;
                en_max = std::max(en_max, window_nrg_f);
            }

            // ---- attack 检测 ----
            bool attack = false;
            size_t attack_index = 0;
            float en_m1 = prev_filtered[kNumWindows - 1]; // 上一帧最后一窗

            for (size_t i = 0; i < kNumWindows; ++i) {
                acc_window_nrg = kOneMinusAccWindowNrgFac * acc_window_nrg + kAccWindowNrgFac * en_m1;
                if (cur_filtered[i] * kInvAttackRatio > acc_window_nrg) {
                    attack = true;
                    attack_index = i;
                }
                en_m1 = cur_filtered[i];
            }

            if (en_max < kMinAttackNrg)
                attack = false;

            // 跨帧延续：上一帧 attack 且上一帧最后一窗仍远大于当前首窗
            // 原版：prev_last/16 > 10·cur_1/16（fMult 的定点换算），即 prev_last > 10·cur_1
            if (!attack && last_attack) {
                if (prev_filtered[kNumWindows - 1] > 10.0f * cur_filtered[1]
                    && last_attack_index == kNumWindows - 1) {
                    attack = true;
                    attack_index = 0;
                }
            }

            // ODF：最大滤波能量相对历史能量比（近似 FDK 的攻击强度）
            const float odf_value = (acc_window_nrg > 1e-12f) ? en_max / acc_window_nrg : 0.0f;

            // ---- 状态推进 ----
            last_attack = attack ? 1 : 0;
            last_attack_index = attack_index;
            prev_filtered = cur_filtered;

            result.odf.push_back(odf_value);
            if (attack) {
                const size_t onset_frame = result.odf.size() - 1;
                result.onset_frames.push_back(onset_frame);
                result.onset_samples.push_back(frame_start + attack_index * window_len);
            }
        }

        return result;
    }

private:
    /** FDK-AAC 二阶高通系数 */
    static constexpr float kHiPassCoeff0 = -0.5095f;
    static constexpr float kHiPassCoeff1 = 0.7548f;

    size_t granule_ = kDefaultGranule;
};

} // namespace qwqdsp_test
