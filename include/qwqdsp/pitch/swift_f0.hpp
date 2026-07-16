#pragma once

#ifdef QWQDSP_HAVE_EIGEN

#include "qwqdsp/pitch/swift_f0_model.hpp"
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cstring>

namespace qwqdsp_swift_f0 {

// ------------------------------------------------------------
// SwiftF0Inference
// ------------------------------------------------------------
// 基于 Eigen 的 swift_f0 卷积网络推理。
// 输入：log-magnitude 频谱 [T, 132]（由外部 STFT 提供）
// 输出：pitch_hz [T], confidence [T]
//
// 前处理（外部实现）：
//   音频 16kHz → STFT(hop=256, n_fft=1024, hann) → mag → slice[3:135) → log(1e-8)
// ------------------------------------------------------------

class SwiftF0Inference {
public:
    using Tensor2f = Eigen::Tensor<float, 2>;
    using Tensor3f = Eigen::Tensor<float, 3>;
    using Tensor4f = Eigen::Tensor<float, 4>;

    SwiftF0Inference() = default;

    /**
     * @brief 运行完整推理
     * @param log_mag_spec   [T, 132] log-magnitude 频谱
     * @param pitch_hz       输出 [T] 基频估计 (Hz)
     * @param confidence     输出 [T] 置信度 (0~1)
     */
    void Process(const Tensor2f& log_mag_spec,
                 Eigen::Ref<Eigen::VectorXf> pitch_hz,
                 Eigen::Ref<Eigen::VectorXf> confidence) {
        const int T = static_cast<int>(log_mag_spec.dimension(0));
        const int F = static_cast<int>(log_mag_spec.dimension(1)); // 132

        // ONNX internal layout: [N=1, C=1, H=freq, W=time]
        Tensor4f x(1, 1, F, T);
        for (int t = 0; t < T; ++t)
            for (int f = 0; f < F; ++f)
                x(0, 0, f, t) = log_mag_spec(t, f);

        // ---- Conv layers (all SAME padding, stride=1) ----
        // Each keeps spatial shape [F, T]
        x = Conv2DSameReLU(x, ConvLayer1::kWeight.data(), ConvLayer1::kBias.data(),
                           8, 1, 5, 5);  // [1, 8,  F, T]
        x = Conv2DSameReLU(x, ConvLayer2::kWeight.data(), ConvLayer2::kBias.data(),
                           16, 8, 5, 5); // [1, 16, F, T]
        x = Conv2DSameReLU(x, ConvLayer3::kWeight.data(), ConvLayer3::kBias.data(),
                           32, 16, 5, 5); // [1, 32, F, T]
        x = Conv2DSameReLU(x, ConvLayer4::kWeight.data(), ConvLayer4::kBias.data(),
                           64, 32, 5, 5); // [1, 64, F, T]
        x = Conv2DSameReLU(x, ConvLayer5::kWeight.data(), ConvLayer5::kBias.data(),
                           1, 64, 5, 5);  // [1, 1,  F, T]

        // Squeeze dim 1 → [1, F, T]
        // Transpose → [T, F] for FreqProjection
        Tensor2f squeezed(T, F);
        for (int t = 0; t < T; ++t)
            for (int f = 0; f < F; ++f)
                squeezed(t, f) = x(0, 0, f, t);

        // ---- FreqProjection (1×1 conv) ----
        // weight [200, 132, 1], bias [200]
        // output[t, p] = sum_f squeezed[t, f] * w[p, f, 0] + b[p]
        Tensor2f proj(T, kNumPitchBins);
        const float* w_fp = FreqProjection::kWeight.data();
        const float* b_fp = FreqProjection::kBias.data();
        for (int t = 0; t < T; ++t) {
            for (int p = 0; p < kNumPitchBins; ++p) {
                float sum = b_fp[p];
                for (int f = 0; f < kNumMelBins; ++f) {
                    sum += squeezed(t, f) * w_fp[p * kNumMelBins + f];
                }
                proj(t, p) = sum;
            }
        }

        // ---- Softmax over pitch bins (dim=1) ----
        Tensor2f softmax_out(T, kNumPitchBins);
        Eigen::VectorXi argmax(T);
        for (int t = 0; t < T; ++t) {
            float max_val = proj(t, 0);
            argmax(t) = 0;
            for (int p = 1; p < kNumPitchBins; ++p) {
                if (proj(t, p) > max_val) {
                    max_val = proj(t, p);
                    argmax(t) = p;
                }
            }

            float sum_exp = 0.0f;
            for (int p = 0; p < kNumPitchBins; ++p) {
                float e = std::exp(proj(t, p) - max_val);
                softmax_out(t, p) = e;
                sum_exp += e;
            }
            for (int p = 0; p < kNumPitchBins; ++p) {
                softmax_out(t, p) /= sum_exp;
            }
        }

        // ---- Voicing mask: |p - argmax| <= 9 ----
        // mask out pitch bins far from the argmax peak
        constexpr int kVoicingHalfWidth = 9;
        Tensor2f masked(T, kNumPitchBins);
        for (int t = 0; t < T; ++t) {
            float sum_masked = 0.0f;
            for (int p = 0; p < kNumPitchBins; ++p) {
                if (std::abs(p - argmax(t)) <= kVoicingHalfWidth) {
                    masked(t, p) = softmax_out(t, p);
                    sum_masked += softmax_out(t, p);
                } else {
                    masked(t, p) = 0.0f;
                }
            }
            // Renormalize
            float norm = sum_masked + 1e-7f;
            for (int p = 0; p < kNumPitchBins; ++p) {
                masked(t, p) /= norm;
            }
            // Confidence = sum of masked softmax (before renormalization)
            confidence(t) = sum_masked;
        }

        // ---- Centroid (weighted sum): pitch_hz ----
        for (int t = 0; t < T; ++t) {
            float sum = 0.0f;
            for (int p = 0; p < kNumPitchBins; ++p) {
                sum += masked(t, p) * kPitchBinCenters[p];
            }
            pitch_hz(t) = sum;
        }
    }

private:
    // ------------------------------------------------------------
    // Conv2D (SAME padding) + Bias + ReLU
    // ------------------------------------------------------------
    // input  shape: [N, C, H, W]
    // output shape: [N, OutC, H, W]  (SAME padding, stride=1)
    // weight_data: OutC * C * kH * kW  flat float array
    // bias_data:   OutC flat float array
    // ------------------------------------------------------------
    static Tensor4f Conv2DSameReLU(
        const Tensor4f& input,
        const float* weight_data,
        const float* bias_data,
        int out_c, int in_c, int kH, int kW) {

        int N = static_cast<int>(input.dimension(0));
        int H = static_cast<int>(input.dimension(2));
        int W = static_cast<int>(input.dimension(3));

        int pad_h = kH / 2;  // SAME: pad = floor(k/2)
        int pad_w = kW / 2;

        Tensor4f output(N, out_c, H, W);

        for (int n = 0; n < N; ++n) {
            for (int oc = 0; oc < out_c; ++oc) {
                for (int h = 0; h < H; ++h) {
                    for (int w = 0; w < W; ++w) {
                        float sum = bias_data[oc];
                        for (int ic = 0; ic < in_c; ++ic) {
                            for (int kh = 0; kh < kH; ++kh) {
                                int ih = h + kh - pad_h;
                                if (ih < 0 || ih >= H) continue;
                                for (int kw = 0; kw < kW; ++kw) {
                                    int iw = w + kw - pad_w;
                                    if (iw < 0 || iw >= W) continue;

                                    float in_val = input(n, ic, ih, iw);
                                    float w_val = weight_data[
                                        oc * (in_c * kH * kW) +
                                        ic * (kH * kW) +
                                        kh * kW +
                                        kw
                                    ];
                                    sum += in_val * w_val;
                                }
                            }
                        }
                        // ReLU
                        output(n, oc, h, w) = (sum > 0.0f) ? sum : 0.0f;
                    }
                }
            }
        }

        return output;
    }
};

} // namespace qwqdsp_swift_f0

#endif // QWQDSP_HAVE_EIGEN
