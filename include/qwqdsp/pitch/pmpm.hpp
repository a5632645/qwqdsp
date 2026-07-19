#pragma once
#include "pmpm/pmpm_core.hpp"
#include "pyin/pyin_hmm.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

namespace qwqdsp_pitch {

// ------------------------------------------------------------
// Pmpm
// ------------------------------------------------------------
// 概率 MPM (pMPM) 基频估计算法。
// 多 cutoff 扫描 NSDF 峰值，HMM 维特比解码平滑轨迹。
//
// 用法:
//   Pmpm pmpm;
//   pmpm.Init(sample_rate, block_size, step_size);
//
//   for (each frame) {
//       pmpm.Process(audio_block);
//   }
//
//   auto pitch_track = pmpm.GetPitchTrack();  // HMM 平滑后的基频轨迹

class Pmpm {
public:
    void Init(float sample_rate, int block_size, int step_size = 0) {
        fs_ = sample_rate;
        block_size_ = block_size;
        step_size_ = (step_size > 0) ? step_size : block_size / 4;

        core_.Init(sample_rate, block_size);
        pitch_hmm_.Init();

        frame_count_ = 0;
        pitch_probs_.clear();
        result_.clear();
    }

    // --------------------------------------------------------
    // Process
    // --------------------------------------------------------
    // 处理一帧音频，提取基频候选概率分布并存储。
    // 需在所有帧处理完毕后调用 GetPitchTrack()。
    void Process(std::span<const float> block) {
        auto candidates = core_.Process(block, 10, 0.001f);

        // 将 Hz → MIDI pitch 并存储
        std::vector<std::pair<double, double>> frame_probs;
        for (auto const& c : candidates) {
            if (c.pitch_hz > 0.0f) {
                double midi = 12.0 * std::log(static_cast<double>(c.pitch_hz) / 440.0) / std::log(2.0) + 69.0;
                frame_probs.emplace_back(midi, static_cast<double>(c.probability));
            }
        }
        pitch_probs_.push_back(std::move(frame_probs));
        ++frame_count_;
    }

    // --------------------------------------------------------
    // GetPitchTrack
    // --------------------------------------------------------
    // 对所有已处理的帧运行 HMM 维特比解码，返回平滑后的基频轨迹。
    // 返回的 vector 长度为 frame_count_，元素为 Hz，无声帧为 0。
    std::vector<float> GetPitchTrack() {
        if (pitch_probs_.empty())
            return {};

        // 多候选: 用不同中心音高偏置运行多次 HMM
        size_t const n_candidate = 13;
        size_t const n_frame = pitch_probs_.size();

        // 预备正态分布权重
        auto normal_pdf = [](double x, double mean, double sd) -> double {
            double z = (x - mean) / sd;
            return std::exp(-0.5 * z * z) / (sd * std::sqrt(2.0 * std::numbers::pi));
        };
        double const max_norm = normal_pdf(0.0, 0.0, 8.0);

        std::vector<std::vector<float>> pitch_tracks(n_candidate, std::vector<float>(n_frame, 0.0f));
        std::vector<float> freq_sum(n_candidate, 0.0f);
        std::vector<int> freq_count(n_candidate, 0);

        for (size_t ic = 0; ic < n_candidate; ++ic) {
            double centre_pitch = 45.0 + 3.0 * static_cast<double>(ic);

            // 用 Gaussian 加权每个帧的候选概率
            std::vector<std::vector<double>> obs_prob(n_frame);
            for (size_t f = 0; f < n_frame; ++f) {
                auto const& frame_candidates = pitch_probs_[f];
                size_t n_c = frame_candidates.size();

                std::vector<std::pair<double, double>> weighted;
                double sum_prob = 0.0;
                for (size_t i = 0; i < n_c; ++i) {
                    double pitch = frame_candidates[i].first;
                    double prob =
                        frame_candidates[i].second * normal_pdf(pitch - centre_pitch, 0.0, 8.0) / max_norm * 2.0;
                    weighted.emplace_back(pitch, prob);
                    sum_prob += prob;
                }
                if (sum_prob > 0.0) {
                    for (auto& w : weighted) {
                        w.second /= sum_prob;
                    }
                }

                obs_prob[f] = pitch_hmm_.CalculateObsProb(weighted);
            }

            // 维特比解码
            auto path = pitch_hmm_.DecodeViterbi(obs_prob);

            // 将状态路径转换为频率
            for (size_t f = 0; f < n_frame && f < path.size(); ++f) {
                double freq = pitch_hmm_.GetFrequency(path[f]);
                if (freq > 0) {
                    pitch_tracks[ic][f] = static_cast<float>(freq);
                    freq_sum[ic] += static_cast<float>(freq);
                    freq_count[ic]++;
                }
            }
        }

        // 去重: 去除 80% 以上帧一致的重复轨迹
        std::vector<bool> is_duplicate(n_candidate, false);
        for (size_t i = 0; i < n_candidate; ++i) {
            if (freq_count[i] <= 0) {
                is_duplicate[i] = true;
                continue;
            }
            for (size_t j = i + 1; j < n_candidate; ++j) {
                if (freq_count[j] <= 0)
                    continue;
                size_t count_equal = 0;
                for (size_t f = 0; f < n_frame; ++f) {
                    if ((pitch_tracks[j][f] == 0.0f && pitch_tracks[i][f] == 0.0f)
                        || (pitch_tracks[i][f] > 0.0f && pitch_tracks[j][f] > 0.0f
                            && std::abs(pitch_tracks[i][f] / pitch_tracks[j][f] - 1.0f) < 0.01f)) {
                        ++count_equal;
                    }
                }
                if (static_cast<float>(count_equal) / static_cast<float>(n_frame) > 0.8f) {
                    if (freq_count[i] > freq_count[j]) {
                        is_duplicate[j] = true;
                    }
                    else if (i < j) {
                        is_duplicate[i] = true;
                    }
                }
            }
        }

        // 选择最佳轨迹（帧数最多的非重复轨迹）
        int best_idx = -1;
        int best_count = 0;
        for (size_t i = 0; i < n_candidate; ++i) {
            if (!is_duplicate[i] && freq_count[i] > best_count) {
                best_count = freq_count[i];
                best_idx = static_cast<int>(i);
            }
        }

        if (best_idx >= 0) {
            result_ = std::move(pitch_tracks[best_idx]);
        }
        else {
            result_.assign(n_frame, 0.0f);
        }

        return result_;
    }

    // 直接访问每帧的候选概率（用于自定义后处理）
    std::span<const std::vector<std::pair<double, double>>> PitchProbabilities() const noexcept {
        return pitch_probs_;
    }

    void SetMinPitch(float hz) noexcept {
        core_.SetMinPitch(hz);
    }

    void SetMaxPitch(float hz) noexcept {
        core_.SetMaxPitch(hz);
    }
private:
    float fs_{};
    int block_size_{};
    int step_size_{};

    PmpmCore core_;
    MonoPitchHmm pitch_hmm_;

    int frame_count_{0};
    std::vector<std::vector<std::pair<double, double>>> pitch_probs_;
    std::vector<float> result_;
};

} // namespace qwqdsp_pitch
