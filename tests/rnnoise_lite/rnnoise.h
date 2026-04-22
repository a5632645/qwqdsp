#ifndef __RNNOISE_LITE_H__
#define __RNNOISE_LITE_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "mfcc.h"
#include "nnom/inc/nnom.h"

#define NUM_FILTER     20
#define NUM_ORDER      1
#define NUM_COEFF_PAIR 3
#define NUM_FEATURES   NUM_FILTER

enum {
    kFrameSize = 512,
};

struct BandpassBiquad {
    int32_t s1;
    int32_t s2;
    int32_t b;
    int32_t a1;
    int32_t a2;
};

struct RnnoiseState {
    int16_t in_buffer[512];
    int in_count;
    int16_t out_buffer[1024];
    int out_count;

    nnom_model_t* model_;
    mfcc_t mfcc_;

    float mfcc_feature_prev[NUM_FEATURES];
    float mfcc_feature_diff_prev[NUM_FEATURES];
    float band_gains_prev[NUM_FILTER];
    struct BandpassBiquad eq[NUM_FILTER];
};

void RnnoiseState_Init(struct RnnoiseState* state);
void RnnoiseState_Reset(struct RnnoiseState* state);
void RnnoiseState_Process(struct RnnoiseState* state, int16_t const* src, int16_t* dst, int num_samples);

#ifdef __cplusplus
}
#endif

#endif // ifndef __RNNOISE_LITE_H__
