#pragma once
#include <algorithm>
#include <cmath>
#include <qwqdsp/fx/phase_vocoder2.hpp>
#include <qwqdsp/segement/analyze_synthsis_offline2.hpp>
#include <span>
#include <vector>

#include "phase_locked_vocoder.hpp" // detail::LinearResample

namespace qwqdsp_test {

// ------------------------------------------------------------
// RunPhaseVocoder3
// ------------------------------------------------------------
/**
 * @brief 库版相位梯度声码器 + 离线分析/合成分帧器。
 *
 * 使用 qwqdsp_fx::PhaseGradientVocoder 处理每一帧，由
 * qwqdsp_segement::AnalyzeSynthsisOffline2 负责分帧与重叠相加，
 * 最后通过线性重采样还原音高移动部分。
 *
 * @param input 输入单声道信号
 * @param kt    时间拉伸比
 * @param kp    音高比例（>1 升高，<1 降低）
 * @return 处理后的输出信号
 */
static inline std::vector<float> RunPhaseVocoder3(std::span<const float> input, float kt, float kp) {
    qwqdsp_fx::PhaseGradientVocoder dsp;
    dsp.Init(4096, 8192, 1024);
    dsp.Reset();
    int hop_analyze = dsp.SetTimeStretch(kt * kp);

    qwqdsp_segement::AnalyzeSynthsisOffline2 as;
    as.SetInputHop(hop_analyze);
    as.SetOutputHop(1024);
    as.SetInputSize(4096);
    as.SetOutputSize(4096);

    std::vector<float> out;
    as.Process(input, out,
               [&dsp](std::span<const float> input, std::span<float> output) { dsp.Process(input, output); });

    out = detail::LinearResample(out, kp);
    return out;
}

} // namespace qwqdsp_test
