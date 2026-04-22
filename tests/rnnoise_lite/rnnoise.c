#include "rnnoise.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "weights.h"

#define RNNLITE_MIN(a, b)            ((a) > (b) ? (b) : (a))
#define RNNLITE_MAX(a, b)            ((a) > (b) ? (a) : (b))
#define RNNLITE_CLAMP(x, xmin, xmax) RNNLITE_MAX(RNNLITE_MIN(x, xmax), xmin)

// clang-format off
static const float kFilterCoeffA1A2[] = {
    -1.96123368f, 0.96283383f,
    -1.95164666f, 0.95824788f,
    -1.93714371f, 0.95310919f,
    -1.91617215f, 0.9473545f,
    -1.88677298f, 0.94091414f,
    -1.84649779f, 0.93371159f,
    -1.79232148f, 0.9256631f,
    -1.72056016f, 0.91667725f,
    -1.62681198f, 0.90665466f,
    -1.50595143f, 0.89548775f,
    -1.35222737f, 0.88306056f,
    -1.15954321f, 0.86924871f,
    -0.92203548f, 0.85391956f,
    -0.63511305f, 0.83693243f,
    -0.29716294f, 0.81813906f,
    0.08785538f,  0.79738417f,
    0.50677509f,  0.77450622f,
    0.93265768f,  0.74933815f,
    1.32004587f,  0.72170809f,
    1.60202043f,  0.6914399f
};

static const float kFilterCoeffB[] = {
    0.01858308f,
    0.02087606f,
    0.0234454f,
    0.02632275f,
    0.02954293f,
    0.0331442f,
    0.03716845f,
    0.04166138f,
    0.04667267f,
    0.05225612f,
    0.05846972f,
    0.06537564f,
    0.07304022f,
    0.08153378f,
    0.09093047f,
    0.10130791f,
    0.11274689f,
    0.12533093f,
    0.13914595f,
    0.15428005f
};
// clang-format on

static void RnnoiseState_ProcessFrame(struct RnnoiseState* state);
static void RnnoiseState_ShiftInBuffer(struct RnnoiseState* state);
static int RnnoiseState_ReadOutput(struct RnnoiseState* state, int16_t* ptr, int num_samples);
static void RnnoiseState_UpdateEq(struct RnnoiseState* state, float band_gains[NUM_FILTER]);
static void RnnoiseState_ProcessEq(struct RnnoiseState* state, int16_t const* in_ptr, int16_t* out_ptr,
                                   int num_samples);

#define Q24_FROM_FLOAT(f) ((int32_t)((f) * 16777216.0f))
#define Q24_FROM_Q15(x)   ((int32_t)((int32_t)(x) << 9))
#define Q15_FROM_Q24(x)   (RNNLITE_CLAMP((x) >> 9, INT16_MIN, INT16_MAX))

static void BandpassBiquad_Process(struct BandpassBiquad* biquad, int16_t const* in_ptr, int16_t* out_ptr,
                                   int num_samples);
static void BandpassBiquad_ProcessAdd(struct BandpassBiquad* biquad, int16_t const* in_ptr, int16_t* out_ptr,
                                      int num_samples);

// ----------------------------------------
// rnnoise private
// ----------------------------------------

static void quantize_data(float* din, int8_t* dout, uint32_t size, uint32_t int_bit) {
    float limit = (1 << int_bit);
    for (uint32_t i = 0; i < size; i++)
        dout[i] = (int8_t)(RNNLITE_MAX(RNNLITE_MIN(din[i], limit), -limit) / limit * 127);
}

void RnnoiseState_Init(struct RnnoiseState* state) {
    state->model_ = nnom_model_create();
    mfcc_create(&state->mfcc_, NUM_FEATURES, 0, NUM_FEATURES, kFrameSize, 0, true);
    for (int i = 0; i < NUM_FILTER; ++i) {
        state->eq[i].a1 = Q24_FROM_FLOAT(kFilterCoeffA1A2[2 * i]);
        state->eq[i].a2 = Q24_FROM_FLOAT(kFilterCoeffA1A2[2 * i + 1]);
    }
}

void RnnoiseState_Reset(struct RnnoiseState* state) {
    state->in_count = kFrameSize / 2;
    state->out_count = 0;
    memset(state->in_buffer, 0, sizeof(int16_t) * kFrameSize);
    memset(state->out_buffer, 0, sizeof(int16_t) * kFrameSize);
    memset(state->mfcc_feature_prev, 0, sizeof(state->mfcc_feature_prev));
    memset(state->mfcc_feature_diff_prev, 0, sizeof(state->mfcc_feature_diff_prev));
    memset(state->band_gains_prev, 0, sizeof(state->band_gains_prev));
    for (int i = 0; i < NUM_FILTER; ++i) {
        state->eq[i].s1 = 0;
        state->eq[i].s2 = 0;
    }
}

static void RnnoiseState_ProcessFrame(struct RnnoiseState* state) {
    float mfcc_feature[NUM_FEATURES] = {0};
    mfcc_compute(&state->mfcc_, state->in_buffer, mfcc_feature);

    float mfcc_feature_diff[NUM_FEATURES] = {0};
    float mfcc_feature_diff1[NUM_FEATURES] = {0};
    for (uint32_t i = 0; i < NUM_FEATURES; i++) {
        mfcc_feature_diff[i] = mfcc_feature[i] - state->mfcc_feature_prev[i];
        mfcc_feature_diff1[i] = mfcc_feature_diff[i] - state->mfcc_feature_diff_prev[i];
    }
    memcpy(state->mfcc_feature_prev, mfcc_feature, NUM_FEATURES * sizeof(float));
    memcpy(state->mfcc_feature_diff_prev, mfcc_feature_diff, NUM_FEATURES * sizeof(float));

    float nn_features[64] = {0};
    memcpy(nn_features, mfcc_feature, NUM_FEATURES * sizeof(float));
    memcpy(&nn_features[NUM_FEATURES], mfcc_feature_diff, 10 * sizeof(float));
    memcpy(&nn_features[NUM_FEATURES + 10], mfcc_feature_diff1, 10 * sizeof(float));

    int8_t nn_features_q7[64] = {0};
    quantize_data(nn_features, nn_features_q7, NUM_FEATURES + 20, 3);
    memcpy(nnom_input_data, nn_features_q7, sizeof(nnom_input_data));
    model_run(state->model_);

    float band_gains[NUM_FILTER] = {0};
    for (int i = 0; i < NUM_FEATURES; i++) {
        band_gains[i] = (float)(nnom_output_data[i]) / 127.f;
    }
    // for (int i = 0; i < NUM_FEATURES; i++) {
    //     band_gains[i] = RNNLITE_MAX(state->band_gains_prev[i] * 0.8f, band_gains[i]);
    // }
    // memcpy(state->band_gains_prev, band_gains, NUM_FEATURES * sizeof(float));

    RnnoiseState_UpdateEq(state, band_gains);
    RnnoiseState_ProcessEq(state, state->in_buffer + kFrameSize / 2, state->out_buffer + state->out_count,
                           kFrameSize / 2);
    state->out_count += kFrameSize / 2;
}

static void RnnoiseState_ShiftInBuffer(struct RnnoiseState* state) {
    memmove(state->in_buffer, state->in_buffer + kFrameSize / 2, (kFrameSize / 2) * sizeof(int16_t));
    state->in_count -= kFrameSize / 2;
}

static int RnnoiseState_ReadOutput(struct RnnoiseState* state, int16_t* ptr, int num_samples) {
    if (num_samples > state->out_count) {
        num_samples = state->out_count;
    }
    if (num_samples == 0) return 0;

    memcpy(ptr, state->out_buffer, num_samples * sizeof(int16_t));
    memmove(state->out_buffer, state->out_buffer + num_samples, (state->out_count - num_samples) * sizeof(int16_t));
    state->out_count -= num_samples;
    return num_samples;
}

static void RnnoiseState_UpdateEq(struct RnnoiseState* state, float band_gains[NUM_FILTER]) {
    for (int i = 0; i < NUM_FILTER; ++i) {
        state->eq[i].b = Q24_FROM_FLOAT(kFilterCoeffB[i] * band_gains[i]);
    }
}

static void RnnoiseState_ProcessEq(struct RnnoiseState* state, int16_t const* in_ptr, int16_t* out_ptr,
                                   int num_samples) {
    BandpassBiquad_Process(&state->eq[0], in_ptr, out_ptr, num_samples);
    for (int i = 1; i < NUM_FILTER; ++i) {
        BandpassBiquad_ProcessAdd(&state->eq[i], in_ptr, out_ptr, num_samples);
    }
}

void RnnoiseState_Process(struct RnnoiseState* state, int16_t const* src, int16_t* dst, int num_samples) {
    int num_write_samples_need = num_samples;
    while (num_samples != 0) {
        int need = kFrameSize - state->in_count;
        int process = num_samples > need ? need : num_samples;
        memcpy(state->in_buffer + state->in_count, src, sizeof(int16_t) * process);
        state->in_count += process;
        num_samples -= process;
        src += process;

        if (state->in_count == kFrameSize) {
            RnnoiseState_ProcessFrame(state);
            RnnoiseState_ShiftInBuffer(state);
        }

        int readed = RnnoiseState_ReadOutput(state, dst, RNNLITE_MIN(num_write_samples_need, state->out_count));
        num_write_samples_need -= readed;
        dst += readed;
    }
    memset(dst, 0, sizeof(int16_t) * num_write_samples_need);
}

// ----------------------------------------
// bandpass biquad private
// ----------------------------------------

static void BandpassBiquad_Process(struct BandpassBiquad* biquad, int16_t const* in_ptr, int16_t* out_ptr,
                                   int num_samples) {
    int32_t s1 = biquad->s1;
    int32_t s2 = biquad->s2;
    int32_t b0 = biquad->b;
    int32_t a1 = biquad->a1;
    int32_t a2 = biquad->a2;
    for (int i = 0; i < num_samples; ++i) {
        int32_t x = Q24_FROM_Q15(in_ptr[i]);

        int64_t acc = (int64_t)x * b0 + ((int64_t)s1 << 24);
        int32_t output = acc >> 24;

        acc = ((int64_t)s2 << 24) - (int64_t)output * a1;
        s1 = acc >> 24;

        acc = -(int64_t)x * b0 - (int64_t)output * a2;
        s2 = acc >> 24;

        out_ptr[i] = Q15_FROM_Q24(output);
    }
    biquad->s1 = s1;
    biquad->s2 = s2;
}

static void BandpassBiquad_ProcessAdd(struct BandpassBiquad* biquad, int16_t const* in_ptr, int16_t* out_ptr,
                                      int num_samples) {
    int32_t s1 = biquad->s1;
    int32_t s2 = biquad->s2;
    int32_t b0 = biquad->b;
    int32_t a1 = biquad->a1;
    int32_t a2 = biquad->a2;
    for (int i = 0; i < num_samples; ++i) {
        int32_t x = Q24_FROM_Q15(in_ptr[i]);

        int64_t acc = (int64_t)x * b0 + ((int64_t)s1 << 24);
        int32_t output = acc >> 24;

        acc = ((int64_t)s2 << 24) - (int64_t)output * a1;
        s1 = acc >> 24;

        acc = -(int64_t)x * b0 - (int64_t)output * a2;
        s2 = acc >> 24;

        out_ptr[i] += Q15_FROM_Q24(output);
    }
    biquad->s1 = s1;
    biquad->s2 = s2;
}
