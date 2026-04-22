#include <algorithm>
#include <array>
#include <limits>
#include <numbers>
#include "denoise_weights.h"
#include "mfcc.h"
#include "nnom/inc/nnom.h"

#define NUM_FILTER 20

#define NUM_ORDER 1

#define NUM_COEFF_PAIR 3

#define FILTER_COEFF_A                                                                         \
    {                                                                                          \
        1., -1.96123368, 0.96283383, 1., -1.95164666, 0.95824788, 1., -1.93714371, 0.95310919, \
        1., -1.91617215, 0.9473545,  1., -1.88677298, 0.94091414, 1., -1.84649779, 0.93371159, \
        1., -1.79232148, 0.9256631,  1., -1.72056016, 0.91667725, 1., -1.62681198, 0.90665466, \
        1., -1.50595143, 0.89548775, 1., -1.35222737, 0.88306056, 1., -1.15954321, 0.86924871, \
        1., -0.92203548, 0.85391956, 1., -0.63511305, 0.83693243, 1., -0.29716294, 0.81813906, \
        1., 0.08785538,  0.79738417, 1., 0.50677509,  0.77450622, 1., 0.93265768,  0.74933815, \
        1., 1.32004587,  0.72170809, 1., 1.60202043,  0.6914399}

#define FILTER_COEFF_B                                                                         \
    {                                                                                          \
        0.01858308, 0., -0.01858308, 0.02087606, 0., -0.02087606, 0.0234454,  0., -0.0234454,  \
        0.02632275, 0., -0.02632275, 0.02954293, 0., -0.02954293, 0.0331442,  0., -0.0331442,  \
        0.03716845, 0., -0.03716845, 0.04166138, 0., -0.04166138, 0.04667267, 0., -0.04667267, \
        0.05225612, 0., -0.05225612, 0.05846972, 0., -0.05846972, 0.06537564, 0., -0.06537564, \
        0.07304022, 0., -0.07304022, 0.08153378, 0., -0.08153378, 0.09093047, 0., -0.09093047, \
        0.10130791, 0., -0.10130791, 0.11274689, 0., -0.11274689, 0.12533093, 0., -0.12533093, \
        0.13914595, 0., -0.13914595, 0.15428005, 0., -0.15428005}

#define NUM_FEATURES NUM_FILTER

class RnnNoise {
public:
    static constexpr int kFrameSize = 512;

    RnnNoise() {
        model_ = nnom_model_create();
        mfcc_ = mfcc_create(NUM_FEATURES, 0, NUM_FEATURES, kFrameSize, 0, true);
    }

    void Reset() noexcept {
        in_count_ = kFrameSize / 2;
        out_count_ = 0;
    }

    void Process(float* ptr, int num_samples) noexcept {
        float* src_ptr = ptr;
        int samples_left = num_samples;

        while (samples_left > 0) {
            // 1. 填充输入缓冲区
            int need = kFrameSize - in_count_;
            int fill = std::min(need, samples_left);
            std::copy_n(src_ptr, fill, in_buffer_.begin() + in_count_);

            in_count_ += fill;
            src_ptr += fill;
            samples_left -= fill;

            // 2. 当缓冲区满时处理一帧 (注意：这里包含 50% 重叠逻辑)
            if (in_count_ == kFrameSize) {
                ProcessFrame();
                // 帧移：保留后 50% 数据供下次重叠
                std::copy_n(in_buffer_.begin() + kFrameSize / 2, kFrameSize / 2, in_buffer_.begin());
                in_count_ = kFrameSize / 2;
            }
            
            // 3. 输出处理后的数据
            // 注意：实时处理会有固定延迟（kRawFrameSize/2），确保 out_buffer 有足够数据
            int can_write = std::min(num_samples, out_count_);
            std::copy_n(out_buffer_.data(), can_write, ptr);
        
            // 移除已输出的数据
            out_count_ -= can_write;
            if (out_count_ > 0) {
                std::memmove(out_buffer_.data(), out_buffer_.data() + can_write, out_count_ * sizeof(float));
            }
        }

    }

    void ProcessFrame() noexcept {
        // 临时存储 16k 的音频块（长度为 512，包含 256 旧数据 + 256 新数据）
        std::array<float, kFrameSize> audio_16k_float;
        std::array<int16_t, kFrameSize> audio_16k_int16;

        // --- 步骤 1: 将当前 48k 缓冲区的全部数据降采样到 16k ---
        // 注意：in_buffer_ 此时包含 1536 个采样 (48k 下的 50% 重叠窗口)
        for (int i = 0; i < kFrameSize; ++i) {
            audio_16k_float[i] = in_buffer_[i];                                                           // 滤波用 float
            audio_16k_int16[i] = static_cast<int16_t>(std::clamp(in_buffer_[i], -1.0f, 1.0f) * 32767.0f); // MFCC 用 int16
        }

        // --- NN 推理部分 (保持不变) ---
        mfcc_compute(mfcc_, audio_16k_int16.data(), mfcc_feature);
        for (uint32_t i = 0; i < NUM_FEATURES; i++) {
            mfcc_feature_diff[i] = mfcc_feature[i] - mfcc_feature_prev[i];
            mfcc_feature_diff1[i] = mfcc_feature_diff[i] - mfcc_feature_diff_prev[i];
        }
        memcpy(mfcc_feature_prev, mfcc_feature, NUM_FEATURES * sizeof(float));
        memcpy(mfcc_feature_diff_prev, mfcc_feature_diff, NUM_FEATURES * sizeof(float));

        // combine MFCC with derivatives
        memcpy(nn_features, mfcc_feature, NUM_FEATURES * sizeof(float));
        memcpy(&nn_features[NUM_FEATURES], mfcc_feature_diff, 10 * sizeof(float));
        memcpy(&nn_features[NUM_FEATURES + 10], mfcc_feature_diff1, 10 * sizeof(float));

        quantize_data(nn_features, nn_features_q7, NUM_FEATURES + 20, 3);
        memcpy(nnom_input_data, nn_features_q7, sizeof(nnom_input_data));
        model_run(model_);

        for (int i = 0; i < NUM_FEATURES; i++) {
            band_gains[i] = (float)(nnom_output_data[i]) / 127.f;
        }

        // for (int i = 0; i < NUM_FEATURES; i++) {
        //     band_gains[i] = std::max(band_gains_prev[i] * 0.8f, (float)nnom_output_data[i] / 127.f);
        // }
        // memcpy(band_gains_prev, band_gains, NUM_FEATURES * sizeof(float));

        // --- 滤波应用部分 ---
        // 【关键修复 2】设置增益。注意：如果 coeff_b 是针对 16k 的，
        // 你不能直接在 48k 的 equalizer 上用，除非 set_gains 内部做了频率映射。
        set_gains((float*)coeff_b, (float*)b_, band_gains, NUM_FILTER, NUM_ORDER);

        std::array<float, kFrameSize / 2> filtered_16k;

        // 输入是 audio_16k_float 的后半段
        equalizer(audio_16k_float.data() + kFrameSize / 2, filtered_16k.data(), kFrameSize / 2, (float*)b_,
                  (float*)coeff_a, NUM_FILTER, NUM_ORDER);

        // --- 步骤 4: 升采样回 48k ---
        // 将 256 个 16k 采样还原为 768 个 48k 采样 (kRawFrameSize / 2)
        float* out_ptr = out_buffer_.data() + out_count_;
        for (int i = 0; i < kFrameSize / 2; ++i) {
            float sample = filtered_16k[i]; // 增益补偿

            // 简单的线性插值升采样（或使用专门的 Upsampler）
            out_ptr[i] = sample;
        }

        out_count_ += kFrameSize / 2;
    }

    float band_gains[NUM_FILTER] = {0};
private:
    void quantize_data(float* din, int8_t* dout, uint32_t size, uint32_t int_bit) {
        float limit = (1 << int_bit);
        for (uint32_t i = 0; i < size; i++) dout[i] = (int8_t)(std::max(std::min(din[i], limit), -limit) / limit * 127);
    }

    void set_gains(float* b_in, float* b_out, float* gains, uint32_t num_band, uint32_t num_order) {
        uint32_t num_coeff = num_order * 2 + 1;
        for (uint32_t i = 0; i < num_band; i++)
            for (uint32_t c = 0; c < num_coeff; c++)
                b_out[num_coeff * i + c] = b_in[num_coeff * i + c] * gains[i]; // only need to set b.
    }

    void equalizer(float* x, float* y, uint32_t signal_len, float* b, float* a, uint32_t num_band, uint32_t num_order) {
        // the y history for each band
        static float y_h[NUM_FILTER][NUM_COEFF_PAIR] = {0};
        static float x_h[NUM_COEFF_PAIR * 2] = {0};
        uint32_t num_coeff = num_order * 2 + 1;

        // i <= num_coeff (where historical x is involved in the first few points)
        // combine state and new data to get a continual x input.
        memcpy(x_h + num_coeff, x, num_coeff * sizeof(float));
        for (uint32_t i = 0; i < num_coeff; i++) {
            y[i] = 0;
            for (uint32_t n = 0; n < num_band; n++) {
                y_h_update(y_h[n], num_coeff);
                y_h[n][0] = b[n * num_coeff] * x_h[i + num_coeff];
                for (uint32_t c = 1; c < num_coeff; c++)
                    y_h[n][0] += b[n * num_coeff + c] * x_h[num_coeff + i - c] - a[n * num_coeff + c] * y_h[n][c];
                y[i] += y_h[n][0];
            }
        }
        // store the x for the state of next round
        memcpy(x_h, &x[signal_len - num_coeff], num_coeff * sizeof(float));

        // i > num_coeff; the rest data not involed the x history
        for (uint32_t i = num_coeff; i < signal_len; i++) {
            y[i] = 0;
            for (uint32_t n = 0; n < num_band; n++) {
                y_h_update(y_h[n], num_coeff);
                y_h[n][0] = b[n * num_coeff] * x[i];
                for (uint32_t c = 1; c < num_coeff; c++)
                    y_h[n][0] += b[n * num_coeff + c] * x[i - c] - a[n * num_coeff + c] * y_h[n][c];
                y[i] += y_h[n][0];
            }
        }
    }

    void y_h_update(float* y_h, uint32_t len) {
        for (uint32_t i = len - 1; i > 0; i--) y_h[i] = y_h[i - 1];
    }

    std::array<float, kFrameSize> in_buffer_{};
    int in_count_{kFrameSize / 2};
    std::array<float, 8192> out_buffer_{};
    int out_count_{};

    nnom_model_t* model_;
    mfcc_t* mfcc_;

    // mfcc features and their derivatives
    float mfcc_feature[NUM_FEATURES] = {0};
    float mfcc_feature_prev[NUM_FEATURES] = {0};
    float mfcc_feature_diff[NUM_FEATURES] = {0};
    float mfcc_feature_diff_prev[NUM_FEATURES] = {0};
    float mfcc_feature_diff1[NUM_FEATURES] = {0};
    // features for NN
    float nn_features[64] = {0};
    int8_t nn_features_q7[64] = {0};

    // NN results, which is the gains for each frequency band
    float band_gains_prev[NUM_FILTER] = {0};

    // 0db gains coefficient
    float coeff_b[NUM_FILTER][NUM_COEFF_PAIR] = FILTER_COEFF_B;
    float coeff_a[NUM_FILTER][NUM_COEFF_PAIR] = FILTER_COEFF_A;
    // dynamic gains coefficient
    float b_[NUM_FILTER][NUM_COEFF_PAIR] = {0};
};
