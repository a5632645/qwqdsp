#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace qwqdsp_pitch {

// ------------------------------------------------------------
// SparseHmm
// ------------------------------------------------------------
// 稀疏维特比解码器基类。
// 使用 from/to/transProb 稀疏表示转移矩阵。
// 参考: pYIN SparseHMM (Queen Mary, Centre for Digital Music)

class SparseHmm {
public:
    virtual ~SparseHmm() = default;

    // 子类必须实现：计算观测概率
    virtual std::vector<double> CalculateObsProb(std::span<const std::pair<double, double>> pitch_prob) const = 0;

    // 维特比解码
    // obs_prob: [nFrame][nState] 观测概率
    // 返回最优状态路径
    std::vector<int> DecodeViterbi(std::span<const std::vector<double>> obs_prob) const noexcept {
        if (obs_prob.empty())
            return {};

        size_t const n_state = init_.size();
        size_t const n_frame = obs_prob.size();
        size_t const n_trans = from_.size();

        std::vector<double> delta(n_state, 0.0);
        std::vector<double> old_delta(n_state);
        std::vector<std::vector<int>> psi(n_frame);
        std::vector<int> path(n_frame, static_cast<int>(n_state) - 1);

        // ----- 初始化第一帧 -----
        double delta_sum = 0.0;
        for (size_t s = 0; s < n_state; ++s) {
            old_delta[s] = init_[s] * obs_prob[0][s];
            delta_sum += old_delta[s];
        }
        if (delta_sum > 0.0) {
            for (size_t s = 0; s < n_state; ++s) {
                old_delta[s] /= delta_sum;
            }
        }
        psi[0].assign(n_state, 0);

        // ----- 前向递推 -----
        for (size_t f = 1; f < n_frame; ++f) {
            delta.assign(n_state, 0.0);
            psi[f].assign(n_state, 0);

            // 稀疏转移：只遍历非零转移
            for (size_t t = 0; t < n_trans; ++t) {
                size_t from = from_[t];
                size_t to = to_[t];
                double val = old_delta[from] * trans_prob_[t];
                if (val > delta[to]) {
                    delta[to] = val;
                    psi[f][to] = static_cast<int>(from);
                }
            }

            // 乘以观测概率并归一化
            delta_sum = 0.0;
            for (size_t s = 0; s < n_state; ++s) {
                delta[s] *= obs_prob[f][s];
                delta_sum += delta[s];
            }

            if (delta_sum > 0.0) {
                for (size_t s = 0; s < n_state; ++s) {
                    old_delta[s] = delta[s] / delta_sum;
                }
            }
            else {
                // 零概率保护
                for (size_t s = 0; s < n_state; ++s) {
                    old_delta[s] = 1.0 / static_cast<double>(n_state);
                }
            }
        }

        // ----- 回溯 -----
        double best_val = 0.0;
        for (size_t s = 0; s < n_state; ++s) {
            if (old_delta[s] > best_val) {
                best_val = old_delta[s];
                path[n_frame - 1] = static_cast<int>(s);
            }
        }
        for (int f = static_cast<int>(n_frame) - 2; f >= 0; --f) {
            path[f] = psi[f + 1][path[f + 1]];
        }

        return path;
    }
protected:
    std::vector<double> init_;
    std::vector<size_t> from_;
    std::vector<size_t> to_;
    std::vector<double> trans_prob_;
};

// ------------------------------------------------------------
// MonoPitchHmm
// ------------------------------------------------------------
// 音高跟踪 HMM: 每个音高有 voiced/unvoiced 两种状态。
// 状态 = 音高索引 + m_nPitch（无声），共 2*m_nPitch 个状态。
// 参考: pYIN MonoPitchHMM

class MonoPitchHmm : public SparseHmm {
public:
    MonoPitchHmm() = default;

    void Init(double min_freq = 61.735, int n_bps = 5, double self_trans = 0.99, double yin_trust = 0.5) noexcept {
        m_min_freq_ = min_freq;
        m_n_bps_ = n_bps;
        m_self_trans_ = self_trans;
        m_yin_trust_ = yin_trust;
        m_transition_width_ = 5 * (m_n_bps_ / 2) + 1;
        m_n_pitch_ = 69 * m_n_bps_; // 69 semitones range

        // 预计算频率
        m_freqs_.resize(2 * m_n_pitch_);
        for (int i = 0; i < m_n_pitch_; ++i) {
            m_freqs_[i] = m_min_freq_ * std::pow(2.0, i / (12.0 * m_n_bps_));
            m_freqs_[i + m_n_pitch_] = -m_freqs_[i]; // 负频率 = 无声
        }

        Build();
    }

    // 计算观测概率
    // pitch_prob:  (midi_pitch, probability) 对列表
    std::vector<double> CalculateObsProb(std::span<const std::pair<double, double>> pitch_prob) const override {
        size_t const n_state = 2 * static_cast<size_t>(m_n_pitch_) + 1;
        std::vector<double> out(n_state, 0.0);

        double prob_yin_pitched = 0.0;

        // 将基频候选分配到最近的 frequency bin
        for (auto const& pp : pitch_prob) {
            double freq = 440.0 * std::pow(2.0, (pp.first - 69.0) / 12.0);
            if (freq <= m_min_freq_)
                continue;

            double old_d = 1000.0;
            for (int i = 0; i < m_n_pitch_; ++i) {
                double d = std::abs(freq - m_freqs_[i]);
                if (old_d < d && i > 0) {
                    out[static_cast<size_t>(i) - 1] = pp.second;
                    prob_yin_pitched += pp.second;
                    break;
                }
                old_d = d;
            }
        }

        double prob_really_pitched = m_yin_trust_ * prob_yin_pitched;

        // 分配 voiced/unvoiced 概率
        for (int i = 0; i < m_n_pitch_; ++i) {
            auto idx = static_cast<size_t>(i);
            if (prob_yin_pitched > 0.0) {
                out[idx] *= (prob_really_pitched / prob_yin_pitched);
            }
            out[idx + static_cast<size_t>(m_n_pitch_)] = (1.0 - prob_really_pitched) / static_cast<double>(m_n_pitch_);
        }

        return out;
    }

    // 获取频率值 (Hz)，负值表示无声
    double GetFrequency(int state) const noexcept {
        if (state >= 0 && state < static_cast<int>(m_freqs_.size())) {
            return m_freqs_[state];
        }
        return 0.0;
    }

    int NumPitch() const noexcept {
        return m_n_pitch_;
    }
private:
    void Build() noexcept {
        size_t const n_state = 2 * static_cast<size_t>(m_n_pitch_);

        // 初始概率：均匀分布
        init_.assign(n_state, 1.0 / static_cast<double>(n_state));

        from_.clear();
        to_.clear();
        trans_prob_.clear();

        for (int i_pitch = 0; i_pitch < m_n_pitch_; ++i_pitch) {
            // 转移范围
            int half_width = static_cast<int>(m_transition_width_) / 2;
            int min_next = std::max(0, i_pitch - half_width);
            int max_next = std::min(m_n_pitch_ - 1, i_pitch + half_width);

            // 三角权重
            std::vector<double> weights;
            double weight_sum = 0.0;
            int theoretical_min = i_pitch - half_width;
            for (int j = min_next; j <= max_next; ++j) {
                double w;
                if (j <= i_pitch) {
                    w = static_cast<double>(j - theoretical_min + 1);
                }
                else {
                    w = static_cast<double>(i_pitch - theoretical_min + 1 - (j - i_pitch));
                }
                weights.push_back(w);
                weight_sum += w;
            }

            auto i_p = static_cast<size_t>(i_pitch);
            auto i_p_u = i_p + static_cast<size_t>(m_n_pitch_);

            for (int j = min_next; j <= max_next; ++j) {
                auto j_p = static_cast<size_t>(j);
                auto j_p_u = j_p + static_cast<size_t>(m_n_pitch_);
                double w = weights[static_cast<size_t>(j - min_next)] / weight_sum;

                // voiced → voiced
                from_.push_back(i_p);
                to_.push_back(j_p);
                trans_prob_.push_back(w * m_self_trans_);

                // voiced → unvoiced
                from_.push_back(i_p);
                to_.push_back(j_p_u);
                trans_prob_.push_back(w * (1.0 - m_self_trans_));

                // unvoiced → unvoiced
                from_.push_back(i_p_u);
                to_.push_back(j_p_u);
                trans_prob_.push_back(w * m_self_trans_);

                // unvoiced → voiced
                from_.push_back(i_p_u);
                to_.push_back(j_p);
                trans_prob_.push_back(w * (1.0 - m_self_trans_));
            }
        }
    }

    double m_min_freq_{61.735};
    int m_n_bps_{5};
    int m_n_pitch_{0};
    size_t m_transition_width_{0};
    double m_self_trans_{0.99};
    double m_yin_trust_{0.5};
    std::vector<double> m_freqs_;
};

} // namespace qwqdsp_pitch
