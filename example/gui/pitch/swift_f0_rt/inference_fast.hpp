#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>

#include <Eigen/Dense>

#include <qwqdsp/pitch/hide/swift_f0_model.hpp>

namespace swift_f0_rt {

/**
 * @brief 移位 GEMM 版 swift_f0 推理
 *
 * 与 qwqdsp_swift_f0::SwiftF0Inference 数值等价（浮点求和顺序不同，差异在 1e-5 量级）。
 *
 * 5 层 5x5 SAME 卷积被改写为：对每个 (kh, kw) 偏移，输入缓冲与输出缓冲之间
 * 只差一个常量行位移（因为时间维补零后行列式布局不变），因此一次连续行块的
 * GEMM 即可完成整层该偏移的累加，避免朴素六重循环与 im2col 的额外内存流量。
 *
 * 时间维在缓冲内左右各补 kPad 帧零；每层结束后把非数据列与守护行清零，
 * 使下一层看到的边界仍是零填充。
 */
class FastInference {
public:
    using Mat = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    /// @brief 预计算各层权重矩阵，必须在 Process 之前调用一次
    void Init() {
        using namespace qwqdsp_swift_f0;
        BuildConv(1, 8, ConvLayer1::kWeight.data(), cw1_);
        BuildConv(8, 16, ConvLayer2::kWeight.data(), cw2_);
        BuildConv(16, 32, ConvLayer3::kWeight.data(), cw3_);
        BuildConv(32, 64, ConvLayer4::kWeight.data(), cw4_);
        BuildConv(64, 1, ConvLayer5::kWeight.data(), cw5_);

        // FreqProjection 权重 [200, 132, 1] → (132 x 200)
        wfp_t_.resize(kNumMelBins, kNumPitchBins);
        for (int p = 0; p < kNumPitchBins; ++p) {
            for (int f = 0; f < kNumMelBins; ++f) {
                wfp_t_(f, p) = FreqProjection::kWeight[static_cast<size_t>(p) * kNumMelBins + f];
            }
        }
        fp_bias_.resize(kNumPitchBins);
        for (int p = 0; p < kNumPitchBins; ++p) {
            fp_bias_[p] = FreqProjection::kBias[p];
        }
        inited_ = true;
    }

    /**
     * @brief 对一段 log-magnitude 频谱推理
     * @param log_mag     [num_frames, 132] 行主序
     * @param num_frames  帧数
     * @param pitch_hz    输出 [num_frames] 基频 (Hz)
     * @param confidence  输出 [num_frames] 置信度
     */
    void Process(const float* log_mag, int num_frames, float* pitch_hz, float* confidence) {
        using namespace qwqdsp_swift_f0;
        if (!inited_ || num_frames <= 0) {
            return;
        }

        constexpr int kNumFreq = kNumMelBins; // 132
        const int tp = num_frames + 2 * kPad;
        const int r = kNumFreq * tp;
        const int rb = r + 2 * kPad;
        EnsureBuffers(rb);

        // ---- 输入填零 + 写入数据列 ----
        buf_[0].setZero(rb, 1);
        for (int t = 0; t < num_frames; ++t) {
            for (int f = 0; f < kNumFreq; ++f) {
                buf_[0](kPad + f * tp + kPad + t, 0) = log_mag[static_cast<size_t>(t) * kNumFreq + f];
            }
        }

        // ---- 5 层卷积 ----
        Conv(buf_[0], cw1_, ConvLayer1::kBias.data(), buf_[1], kNumFreq, num_frames, tp, r, rb);
        Conv(buf_[1], cw2_, ConvLayer2::kBias.data(), buf_[2], kNumFreq, num_frames, tp, r, rb);
        Conv(buf_[2], cw3_, ConvLayer3::kBias.data(), buf_[3], kNumFreq, num_frames, tp, r, rb);
        Conv(buf_[3], cw4_, ConvLayer4::kBias.data(), buf_[4], kNumFreq, num_frames, tp, r, rb);
        Conv(buf_[4], cw5_, ConvLayer5::kBias.data(), buf_[5], kNumFreq, num_frames, tp, r, rb);

        // ---- FreqProjection: (T x 132) * (132 x 200) ----
        x_.resize(num_frames, kNumFreq);
        for (int t = 0; t < num_frames; ++t) {
            for (int f = 0; f < kNumFreq; ++f) {
                x_(t, f) = buf_[5](kPad + f * tp + kPad + t, 0);
            }
        }
        proj_.noalias() = x_ * wfp_t_;
        for (int t = 0; t < num_frames; ++t) {
            for (int p = 0; p < kNumPitchBins; ++p) {
                proj_(t, p) += fp_bias_[p];
            }
        }

        // ---- Softmax → voicing mask → centroid ----
        sm_.resize(num_frames, kNumPitchBins);
        masked_.resize(num_frames, kNumPitchBins);
        constexpr int kVoicingHalfWidth = 9;
        for (int t = 0; t < num_frames; ++t) {
            float max_val = proj_(t, 0);
            int argmax = 0;
            for (int p = 1; p < kNumPitchBins; ++p) {
                if (proj_(t, p) > max_val) {
                    max_val = proj_(t, p);
                    argmax = p;
                }
            }

            float sum_exp = 0.0f;
            for (int p = 0; p < kNumPitchBins; ++p) {
                float e = std::exp(proj_(t, p) - max_val);
                sm_(t, p) = e;
                sum_exp += e;
            }
            for (int p = 0; p < kNumPitchBins; ++p) {
                sm_(t, p) /= sum_exp;
            }

            float sum_masked = 0.0f;
            for (int p = 0; p < kNumPitchBins; ++p) {
                if (std::abs(p - argmax) <= kVoicingHalfWidth) {
                    masked_(t, p) = sm_(t, p);
                    sum_masked += sm_(t, p);
                }
                else {
                    masked_(t, p) = 0.0f;
                }
            }

            float const norm = sum_masked + 1.0e-7f;
            float acc = 0.0f;
            for (int p = 0; p < kNumPitchBins; ++p) {
                acc += (masked_(t, p) / norm) * kPitchBinCenters[p];
            }
            pitch_hz[t] = acc;
            confidence[t] = sum_masked;
        }
    }

private:
    /// 时间维每侧补零帧数
    static constexpr int kPad = 2;

    struct ConvWeights {
        int in_c{};
        int out_c{};
        // w[kh][kw] 形状 (in_c x out_c)
        Mat w[5][5];
    };

    /// @brief 把 [out_c][in_c][5][5] 权重拆成 25 个 (in_c x out_c) 矩阵
    static void BuildConv(int in_c, int out_c, const float* w, ConvWeights& cw) {
        cw.in_c = in_c;
        cw.out_c = out_c;
        for (int kh = 0; kh < 5; ++kh) {
            for (int kw = 0; kw < 5; ++kw) {
                auto& m = cw.w[kh][kw];
                m.resize(in_c, out_c);
                for (int ic = 0; ic < in_c; ++ic) {
                    for (int oc = 0; oc < out_c; ++oc) {
                        size_t idx = (static_cast<size_t>(oc) * in_c + ic) * 25 + kh * 5 + kw;
                        m(ic, oc) = w[idx];
                    }
                }
            }
        }
    }

    void EnsureBuffers(int rb) {
        for (auto& b : buf_) {
            if (b.rows() != rb) {
                b.resize(rb, b.cols() > 0 ? b.cols() : 1);
            }
        }
    }

    /**
     * @brief 一层 5x5 SAME 卷积 + bias + ReLU
     *
     * 缓冲行号 = 数据行号 + kPad；输出数据列 [kPad, kPad+num_frames)。
     * 对固定 (kh, kw)，输入行号 = 输出行号 + (kh-2)*tp + (kw-2)，是常量位移。
     */
    static void Conv(const Mat& in, const ConvWeights& cw, const float* bias, Mat& out, int num_freq, int num_frames,
                     int tp, int r, int rb) {
        out.setZero(rb, cw.out_c);
        for (int kh = 0; kh < 5; ++kh) {
            int const dh = kh - 2;
            int const f_lo = std::max(0, -dh);
            int const f_hi = std::min(num_freq - 1, num_freq - 1 - dh);
            if (f_lo > f_hi) {
                continue;
            }
            int const n = (f_hi - f_lo + 1) * tp;
            int const row_lo = f_lo * tp + kPad;
            int const src_base = row_lo + dh * tp;
            for (int kw = 0; kw < 5; ++kw) {
                int const dw = kw - 2;
                out.middleRows(row_lo, n).noalias() += in.middleRows(src_base + dw, n) * cw.w[kh][kw];
            }
        }

        // ---- bias + ReLU（仅数据列），其余清零 ----
        for (int row = 0; row < rb; ++row) {
            int const data_row = row - kPad;
            int const col = (data_row >= 0 && data_row < r) ? (data_row % tp) : -1;
            bool const is_data = col >= kPad && col < kPad + num_frames;
            for (int oc = 0; oc < cw.out_c; ++oc) {
                if (is_data) {
                    float const v = out(row, oc) + bias[oc];
                    out(row, oc) = v > 0.0f ? v : 0.0f;
                }
                else {
                    out(row, oc) = 0.0f;
                }
            }
        }
    }

    ConvWeights cw1_, cw2_, cw3_, cw4_, cw5_;
    Mat buf_[6];
    Mat x_, proj_, sm_, masked_;
    Mat wfp_t_;
    Eigen::VectorXf fp_bias_;
    bool inited_{false};
};

} // namespace swift_f0_rt
