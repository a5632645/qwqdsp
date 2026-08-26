#pragma once
#include <algorithm>
#include <audio_ops.hpp>
#include <cmath>
#include <format>
#include <iostream>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

namespace qwqdsp_test {

// ------------------------------------------------------------
// PSOLA — Pitch Synchronous Overlap and Add
// ------------------------------------------------------------
/**
 * @brief 离线 TD-PSOLA 时间拉伸/音高移动，支持独立控制音高和共振峰。
 *
 * ## 算法概述
 *
 * 1. **音高标记检测**：基于自相关的基音周期估计 + 波形峰值细化。
 * 2. **颗粒提取**：以每个音高标记为中心，用汉宁窗截取 2 周期长度的颗粒。
 * 3. **重叠相加合成**：将颗粒按目标音高周期重定位到合成时间轴并叠加。
 *
 * ## 共振峰控制
 *
 * 先对输入信号做重采样（改变共振峰），再用 PSOLA 将音高修正到目标值：
 * - `formant_shift` 控制共振峰缩放（>1 升高，<1 降低）
 * - `pitch_shift` 独立控制音高偏移（半音程为单位）
 */
struct Psola {
    /** 采样率 (Hz) */
    float sample_rate_ = 44100.0f;
    /** 音高偏移，单位为半音（0 = 不变，+12 = 升八度） */
    float pitch_shift_semitones_ = 0.0f;
    /** 共振峰缩放因子（1.0 = 不变，>1 升高，<1 降低） */
    float formant_shift_ = 1.0f;
    /** 最小基频 (Hz) */
    float min_f0_ = 75.0f;
    /** 最大基频 (Hz) */
    float max_f0_ = 500.0f;

    /** 处理主函数 */
    std::vector<float> Process(std::span<const float> input) {
        const float pitch_factor = std::pow(2.0f, pitch_shift_semitones_ / 12.0f);

        // ----------------------------------------------------
        // Step 0: 若 formant_shift != 1，先重采样改变共振峰
        // ----------------------------------------------------
        // 重采样后：音高偏移 formant_shift，共振峰偏移 formant_shift，时长缩放 formant_shift
        std::vector<float> proc_input(input.begin(), input.end());
        if (std::abs(formant_shift_ - 1.0f) > 1e-6f) {
            proc_input = LinearResample(proc_input, formant_shift_);
        }

        // ----------------------------------------------------
        // Step 1: 检测音高标记
        // ----------------------------------------------------
        PitchMarks marks = DetectPitchMarks(proc_input);
        if (marks.positions.size() < 2) {
            std::cout << std::format("[PSOLA] 音高标记不足 ({}), 返回原始信号\n", marks.positions.size()) << std::flush;
            return proc_input;
        }

        // ----------------------------------------------------
        // Step 2: PSOLA 合成 — 同时修正音高和时长
        // ----------------------------------------------------
        // PSOLA 需要把 proc_input 的时长（= input 时长 x formant_shift）拉回原始时长
        // 同时把音高从 formant_shift 倍修正到 pitch_factor 倍
        const float psola_pitch = pitch_factor / formant_shift_; // 音高修正比
        const float psola_stretch = 1.0f / formant_shift_;       // 时长修正比

        std::vector<float> output = PsolaSynthesize(proc_input, marks, psola_pitch, psola_stretch);

        // 归一化
        qwqdsp_support::AudioOps::Normalize(output);
        return output;
    }

private:
    // ------------------------------------------------------------
    // 音高标记检测
    // ------------------------------------------------------------
    struct PitchMarks {
        std::vector<size_t> positions; // 采样点位置
        std::vector<float> periods;    // 各标记处的周期（采样点数）
        std::vector<bool> voiced;      // 浊音标记
    };

    PitchMarks DetectPitchMarks(const std::vector<float>& data) {
        const float kMinPeriod = sample_rate_ / max_f0_;
        const float kMaxPeriod = sample_rate_ / min_f0_;
        const size_t min_period = static_cast<size_t>(kMinPeriod);
        const size_t max_period = static_cast<size_t>(kMaxPeriod);
        const size_t def_period = (min_period + max_period) / 2;

        // --- 自相关音高轮廓检测 ---
        const size_t hop = std::max<size_t>(12, static_cast<size_t>(min_period * 0.75f));
        const size_t start = max_period * 2;
        const size_t end = (data.size() > max_period * 2) ? data.size() - max_period * 2 : 0;

        std::vector<float> periods;
        std::vector<float> scores;
        std::vector<bool> voiced;
        std::vector<size_t> positions;

        float prev_period = static_cast<float>(def_period);

        for (size_t center = start; center <= end; center += hop) {
            auto [period, score] =
                DetectPeriod(data, center - max_period, min_period, max_period, data.size(), prev_period);
            const bool is_voiced = (score >= 0.55f && period > 0.0f);
            const float p = is_voiced ? period : prev_period;
            periods.push_back(p);
            scores.push_back(score);
            voiced.push_back(is_voiced);
            positions.push_back(center);
            if (is_voiced) {
                prev_period = period;
            }
        }

        // 周期平滑
        SmoothPeriods(periods, voiced, static_cast<float>(def_period));

        // --- 从轮廓提取音高标记（峰值细化） ---
        return ExtractMarks(data, periods, voiced, scores, positions, start, hop, min_period, max_period);
    }

    /** 在指定位置检测基音周期 */
    static std::pair<float, float> DetectPeriod(const std::vector<float>& data, size_t pos, size_t min_period,
                                                size_t max_period, size_t end, float prev_period) {
        if (pos + max_period * 2 > end) {
            return {0.0f, 0.0f};
        }

        // 如果 prev_period > 0，先在小范围内搜索
        if (prev_period > 0.0f) {
            const size_t local_min = std::clamp(static_cast<size_t>(prev_period * 0.78f), min_period, max_period);
            const size_t local_max = std::clamp(static_cast<size_t>(prev_period * 1.28f), min_period, max_period);
            auto result = DetectPeriodRange(data, pos, local_min, local_max, prev_period);
            if (result.second >= 0.58f)
                return result;
        }

        return DetectPeriodRange(data, pos, min_period, max_period, prev_period);
    }

    /** 在指定滞后范围内检测周期（归一化自相关） */
    static std::pair<float, float> DetectPeriodRange(const std::vector<float>& data, size_t pos, size_t min_lag,
                                                     size_t max_lag, float prev_period) {
        if (min_lag > max_lag)
            return {0.0f, 0.0f};

        const size_t n = max_lag;
        // 计算 e1 = Σ data[pos+i]²
        float e1 = 0.0f;
        for (size_t i = 0; i < n && pos + i < data.size(); ++i) {
            e1 += data[pos + i] * data[pos + i];
        }

        struct LagScore {
            size_t lag;
            float score;
        };
        std::vector<float> corr(max_lag + 2, 0.0f);

        for (size_t lag = min_lag; lag <= max_lag; ++lag) {
            float sum = 0.0f, e2 = 0.0f;
            for (size_t i = 0; i < n && pos + i + lag < data.size(); ++i) {
                const float b = data[pos + i + lag];
                sum += data[pos + i] * b;
                e2 += b * b;
            }
            const float denom = std::sqrt(e1 * e2);
            corr[lag] = (denom > 1e-10f) ? sum / denom : 0.0f;
        }

        // 寻找局部峰值 + 倍周期抑制
        size_t best_lag = 0;
        float best_score = 0.0f;
        float best_metric = -std::numeric_limits<float>::infinity();

        for (size_t lag = min_lag + 1; lag < max_lag; ++lag) {
            if (corr[lag] < 0.3f)
                continue;
            if (corr[lag] < corr[lag - 1] || corr[lag] < corr[lag + 1])
                continue;

            float metric = corr[lag];
            if (prev_period > 0.0f) {
                metric += 0.18f * std::max(-1.0f, 1.0f - std::abs(std::log(static_cast<float>(lag) / prev_period)));
            }

            // 倍周期检测：如果 2*lag 处也有高相关，取之
            const size_t doubled = lag * 2;
            if (doubled < max_lag && corr[doubled] >= corr[lag] * 0.88f && corr[doubled] >= corr[doubled - 1]
                && corr[doubled] >= corr[doubled + 1]) {
                float dm = corr[doubled];
                if (prev_period > 0.0f) {
                    dm += 0.18f * std::max(-1.0f, 1.0f - std::abs(std::log(static_cast<float>(doubled) / prev_period)));
                }
                if (dm > metric) {
                    metric = dm;
                    corr[lag] = corr[doubled]; // 用倍周期的相关值
                }
            }

            if (metric > best_metric) {
                best_metric = metric;
                best_score = corr[lag];
                best_lag = lag;
            }
        }

        if (best_lag == 0) {
            // 退化为最大值查找
            float best_val = -std::numeric_limits<float>::infinity();
            for (size_t lag = min_lag; lag <= max_lag; ++lag) {
                if (corr[lag] > best_val) {
                    best_val = corr[lag];
                    best_lag = lag;
                    best_score = corr[lag];
                }
            }
            // 若极值在搜索边界上，认为周期不在范围内
            if ((best_lag == min_lag && corr[min_lag] > corr[min_lag + 1])
                || (best_lag == max_lag && corr[max_lag] > corr[max_lag - 1])) {
                return {0.0f, 0.0f};
            }
        }

        if (best_score <= 0.35f || best_lag == 0) {
            return {0.0f, std::max(0.0f, best_score)};
        }

        // 抛物线插值子采样精度
        float period = static_cast<float>(best_lag);
        if (best_lag > min_lag && best_lag < max_lag) {
            const float ym1 = corr[best_lag - 1];
            const float y0 = corr[best_lag];
            const float yp1 = corr[best_lag + 1];
            const float denom = ym1 - 2.0f * y0 + yp1;
            if (std::abs(denom) > 1e-8f) {
                period += 0.5f * (ym1 - yp1) / denom;
                period = std::clamp(period, static_cast<float>(best_lag) - 0.5f, static_cast<float>(best_lag) + 0.5f);
            }
        }

        return {period, best_score};
    }

    /** 中值滤波平滑周期轨迹 */
    static void SmoothPeriods(std::vector<float>& periods, const std::vector<bool>& voiced, float def_period) {
        if (periods.empty())
            return;

        // 两遍中值滤波
        for (int pass = 0; pass < 2; ++pass) {
            auto next = periods;
            for (size_t i = 0; i < periods.size(); ++i) {
                if (!voiced[i])
                    continue;
                std::vector<float> scratch;
                for (int k = -2; k <= 2; ++k) {
                    const int idx = static_cast<int>(i) + k;
                    if (idx >= 0 && idx < static_cast<int>(periods.size()) && voiced[idx]) {
                        scratch.push_back(periods[idx]);
                    }
                }
                if (scratch.size() < 3)
                    continue;
                std::sort(scratch.begin(), scratch.end());
                const float med = scratch[scratch.size() / 2];
                next[i] = std::clamp(0.6f * periods[i] + 0.4f * med, med * 0.75f, med * 1.35f);
            }
            periods = next;
        }

        // 前后向填充非浊音段
        size_t first_voiced = 0;
        while (first_voiced < voiced.size() && !voiced[first_voiced])
            ++first_voiced;
        float seed = (first_voiced < voiced.size()) ? periods[first_voiced] : def_period;

        float prev = seed;
        for (size_t i = 0; i < periods.size(); ++i) {
            if (voiced[i]) {
                periods[i] = std::clamp(periods[i], prev * 0.8f, prev * 1.25f);
                prev = periods[i];
            }
            else {
                periods[i] = prev;
            }
        }

        float next_val = prev;
        for (size_t i = periods.size(); i-- > 0;) {
            if (voiced[i])
                next_val = periods[i];
            else
                periods[i] = next_val;
        }
    }

    /** 在波形峰值附近细化位置 */
    static size_t PeakNear(const std::vector<float>& data, size_t center, size_t radius) {
        const size_t start = (center > radius) ? center - radius : 1;
        const size_t end = std::min(data.size() - 2, center + radius);

        // 寻找局部极大值（考虑绝对值）
        size_t best = std::clamp(center, start, end);
        float best_val = -std::numeric_limits<float>::infinity();
        for (size_t i = start; i <= end; ++i) {
            const float v = std::abs(data[i]);
            if (v >= std::abs(static_cast<long double>(data[i - 1]))
                && v >= std::abs(static_cast<long double>(data[i + 1]))) {
                const float score =
                    v
                    - 0.08f * static_cast<float>(std::abs(static_cast<long long>(i) - static_cast<long long>(center)))
                          / std::max(1.0f, static_cast<float>(radius));
                if (score > best_val) {
                    best = i;
                    best_val = score;
                }
            }
        }
        if (best_val > -std::numeric_limits<float>::infinity())
            return best;

        // 退化为距离加权最大值
        float bv = -std::numeric_limits<float>::infinity();
        for (size_t i = start; i <= end; ++i) {
            const float score =
                std::abs(data[i])
                - 0.08f * static_cast<float>(std::abs(static_cast<long long>(i) - static_cast<long long>(center)))
                      / std::max(1.0f, static_cast<float>(radius));
            if (score > bv) {
                bv = score;
                best = i;
            }
        }
        return best;
    }

    /** 从轮廓提取音高标记 */
    PitchMarks ExtractMarks(const std::vector<float>& data, const std::vector<float>& periods,
                            const std::vector<bool>& voiced, const std::vector<float>& scores,
                            const std::vector<size_t>& positions, size_t start, size_t hop, size_t min_period,
                            size_t max_period) {
        // 找一个稳定的锚点
        int anchor_idx = -1;
        for (int i = 1; i < static_cast<int>(voiced.size()) - 1; ++i) {
            if (voiced[i - 1] && voiced[i] && voiced[i + 1]) {
                anchor_idx = i;
                break;
            }
        }
        if (anchor_idx < 0) {
            float best_score = 0.0f;
            for (int i = 0; i < static_cast<int>(voiced.size()); ++i) {
                if (voiced[i] && scores[i] > best_score) {
                    best_score = scores[i];
                    anchor_idx = i;
                }
            }
        }
        if (anchor_idx < 0)
            return {};

        const size_t anchor_center = start + static_cast<size_t>(anchor_idx) * hop;
        const float anchor_period = periods[anchor_idx];
        const size_t anchor_radius = std::max<size_t>(4, static_cast<size_t>(anchor_period * 0.35f));
        const size_t anchor_mark = PeakNear(data, anchor_center, anchor_radius);

        // --- 反向扩展（向左） ---
        std::vector<size_t> head_pos;
        std::vector<float> head_periods;
        std::vector<bool> head_voiced;
        size_t pos = anchor_mark;

        while (pos > min_period) {
            const float period = PeriodAt(periods, positions, start, hop, static_cast<float>(pos));
            const float predicted_f = static_cast<float>(pos) - period;
            if (predicted_f <= 0.0f)
                break;
            const size_t predicted = static_cast<size_t>(predicted_f);
            const bool is_voiced = VoicedAt(voiced, positions, start, hop, static_cast<float>(predicted));

            const size_t radius = std::max<size_t>(
                4,
                static_cast<size_t>(PeriodAt(periods, positions, start, hop, static_cast<float>(predicted)) * 0.35f));
            size_t mark = is_voiced ? PeakNear(data, predicted, radius) : predicted;

            const size_t min_step = std::max<size_t>(1, static_cast<size_t>(period * 0.55f));
            const size_t max_step = std::max(min_step + 1, static_cast<size_t>(period * 1.8f));
            const size_t step = pos - mark;
            if (step < min_step)
                mark = pos - min_step;
            if (step > max_step)
                mark = pos - max_step;
            if (mark <= 0 || mark >= pos)
                break;

            head_pos.push_back(mark);
            head_periods.push_back(PeriodAt(periods, positions, start, hop, static_cast<float>(mark)));
            head_voiced.push_back(is_voiced);
            pos = mark;
        }

        // 反转头部（因为是从锚点向左扩展的）
        std::reverse(head_pos.begin(), head_pos.end());
        std::reverse(head_periods.begin(), head_periods.end());
        std::reverse(head_voiced.begin(), head_voiced.end());

        // --- 组装 ---
        PitchMarks result;
        result.positions = std::move(head_pos);
        result.periods = std::move(head_periods);
        result.voiced = std::move(head_voiced);

        result.positions.push_back(anchor_mark);
        result.periods.push_back(anchor_period);
        result.voiced.push_back(true);

        // --- 正向扩展（向右） ---
        pos = anchor_mark;
        while (pos + min_period < data.size()) {
            const float period = PeriodAt(periods, positions, start, hop, static_cast<float>(pos));
            const float predicted_f = static_cast<float>(pos) + period;
            if (predicted_f + static_cast<float>(min_period) >= static_cast<float>(data.size()))
                break;
            const size_t predicted = static_cast<size_t>(predicted_f);
            const bool is_voiced = VoicedAt(voiced, positions, start, hop, static_cast<float>(predicted));

            const size_t radius = std::max<size_t>(
                4,
                static_cast<size_t>(PeriodAt(periods, positions, start, hop, static_cast<float>(predicted)) * 0.35f));
            size_t mark = is_voiced ? PeakNear(data, predicted, radius) : predicted;

            const size_t min_step = std::max<size_t>(1, static_cast<size_t>(period * 0.55f));
            const size_t max_step = std::max(min_step + 1, static_cast<size_t>(period * 1.8f));
            const size_t step = mark - pos;
            if (step < min_step)
                mark = pos + min_step;
            if (step > max_step)
                mark = pos + max_step;
            if (mark <= pos || mark >= data.size())
                break;

            result.positions.push_back(mark);
            result.periods.push_back(PeriodAt(periods, positions, start, hop, static_cast<float>(mark)));
            result.voiced.push_back(is_voiced);
            pos = mark;
        }

        return result;
    }

    /** 在给定采样点位置插值周期值 */
    static float PeriodAt(const std::vector<float>& periods, const std::vector<size_t>& positions, size_t start,
                          size_t hop, float sample_pos) {
        const float x = (sample_pos - static_cast<float>(start)) / static_cast<float>(hop);
        if (x <= 0.0f)
            return periods.front();
        if (x >= static_cast<float>(periods.size()) - 1.0f)
            return periods.back();
        const size_t i = static_cast<size_t>(x);
        const float frac = x - static_cast<float>(i);
        return periods[i] * (1.0f - frac) + periods[i + 1] * frac;
    }

    /** 判断给定位置是否为浊音 */
    static bool VoicedAt(const std::vector<bool>& voiced, const std::vector<size_t>& positions, size_t start,
                         size_t hop, float sample_pos) {
        const int idx =
            static_cast<int>(std::round((sample_pos - static_cast<float>(start)) / static_cast<float>(hop)));
        const int clamped = std::clamp(idx, 0, static_cast<int>(voiced.size()) - 1);
        return voiced[clamped];
    }

    // ------------------------------------------------------------
    // PSOLA 合成
    // ------------------------------------------------------------

    /** 添加一个颗粒（双半汉宁窗）到输出缓冲区 */
    static void AddGrain(const std::vector<float>& data, size_t src_pos, float left, float right, std::span<float> out,
                         std::span<float> norm, size_t dst_pos) {
        const size_t l = std::max<size_t>(1, static_cast<size_t>(std::round(left)));
        const size_t r = std::max<size_t>(1, static_cast<size_t>(std::round(right)));

        for (size_t i = 0; i < l; ++i) {
            const int si = static_cast<int>(src_pos) - static_cast<int>(l) + static_cast<int>(i);
            const int di = static_cast<int>(dst_pos) - static_cast<int>(l) + static_cast<int>(i);
            if (si < 0 || si >= static_cast<int>(data.size()) || di < 0 || di >= static_cast<int>(out.size()))
                continue;
            // 左半汉宁：从 0 上升到 1
            const float w =
                0.5f * (1.0f - std::cos(std::numbers::pi_v<float> * static_cast<float>(i) / static_cast<float>(l)));
            out[di] += data[si] * w;
            norm[di] += w;
        }

        for (size_t i = 0; i < r; ++i) {
            const int si = static_cast<int>(src_pos) + static_cast<int>(i);
            const int di = static_cast<int>(dst_pos) + static_cast<int>(i);
            if (si < 0 || si >= static_cast<int>(data.size()) || di < 0 || di >= static_cast<int>(out.size()))
                continue;
            // 右半汉宁：从 1 下降到 0
            const float w =
                0.5f * (1.0f + std::cos(std::numbers::pi_v<float> * static_cast<float>(i) / static_cast<float>(r)));
            out[di] += data[si] * w;
            norm[di] += w;
        }
    }

    /** PSOLA 重合成
     *
     * @param data          输入信号
     * @param marks         音高标记
     * @param pitch_factor  音高偏移比（>1 升高，<1 降低）
     * @param stretch_factor 时长拉伸比（>1 变长，<1 变短，=1 不变）
     * @return 合成输出，长度为 data.size() * stretch_factor
     */
    std::vector<float> PsolaSynthesize(const std::vector<float>& data, const PitchMarks& marks, float pitch_factor,
                                       float stretch_factor) {
        const size_t out_len = static_cast<size_t>(std::round(static_cast<float>(data.size()) * stretch_factor));
        std::vector<float> out(out_len, 0.0f);
        std::vector<float> norm(out_len, 0.0f);

        if (marks.positions.size() < 2)
            return out;

        const float min_period_f = sample_rate_ / max_f0_;
        const float max_period_f = sample_rate_ / min_f0_;

        const size_t last = marks.positions.size() - 1;

        // 合成循环：在输出时间轴上以目标周期步进放置颗粒
        float syn_pos_f = 0.0f;
        size_t cursor = 0;

        while (syn_pos_f < static_cast<float>(out_len)) {
            // 当前合成位置对应的原始时间位置（考虑时长拉伸）
            const float src_time = syn_pos_f / stretch_factor;
            if (src_time > static_cast<float>(marks.positions[last]) + marks.periods[last])
                break;

            // 找到最近的原始标记
            while (cursor + 1 < marks.positions.size() && static_cast<float>(marks.positions[cursor + 1]) <= src_time) {
                ++cursor;
            }

            size_t best = cursor;
            if (cursor + 1 < marks.positions.size()) {
                const float dist_curr = std::abs(static_cast<float>(marks.positions[cursor]) - src_time);
                const float dist_next = std::abs(static_cast<float>(marks.positions[cursor + 1]) - src_time);
                if (dist_next < dist_curr)
                    best = cursor + 1;
            }

            // 左右窗口大小 = 到相邻标记的距离
            const float left = (best > 0) ? static_cast<float>(marks.positions[best] - marks.positions[best - 1])
                                          : marks.periods[best];
            const float right = (best < last) ? static_cast<float>(marks.positions[best + 1] - marks.positions[best])
                                              : marks.periods[best];

            const float cl = std::clamp(left, min_period_f, max_period_f * 2.0f);
            const float cr = std::clamp(right, min_period_f, max_period_f * 2.0f);

            AddGrain(data, marks.positions[best], cl, cr, out, norm, static_cast<size_t>(std::round(syn_pos_f)));

            // 合成步长 = 原始周期 / pitch_factor（音高越高周期越短 → 步长越小）
            const float raw_step = marks.periods[best] / pitch_factor;
            const float step = marks.voiced[best] ? raw_step : 0.5f * (cl + cr) / pitch_factor;
            syn_pos_f += std::clamp(step, min_period_f * 0.75f, max_period_f * 1.25f);
        }

        // 归一化（OLA 幅度补偿）
        for (size_t i = 0; i < out_len; ++i) {
            if (norm[i] > 1e-8f) {
                out[i] /= norm[i];
            }
        }

        return out;
    }

    // ------------------------------------------------------------
    // 线性重采样（用于共振峰控制）
    // ------------------------------------------------------------
    static std::vector<float> LinearResample(std::span<const float> input, float ratio) {
        if (std::abs(ratio - 1.0f) < 1e-6f) {
            return std::vector<float>(input.begin(), input.end());
        }
        const size_t out_len = static_cast<size_t>(std::round(static_cast<float>(input.size()) * ratio));
        if (out_len < 2)
            return {};
        std::vector<float> output(out_len);
        if (out_len == 1) {
            output[0] = input[0];
            return output;
        }
        const float step = static_cast<float>(input.size() - 1) / static_cast<float>(out_len - 1);
        for (size_t i = 0; i < out_len; ++i) {
            const float pos = step * static_cast<float>(i);
            const size_t idx = static_cast<size_t>(pos);
            const float frac = pos - static_cast<float>(idx);
            const size_t nxt = std::min(idx + 1, input.size() - 1);
            output[i] = input[idx] + frac * (input[nxt] - input[idx]);
        }
        return output;
    }
};

} // namespace qwqdsp_test
