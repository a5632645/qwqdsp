#ifndef __MFCC_H__
#define __MFCC_H__

#include <math.h>
#include <stdint.h>
#include "string.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SAMP_FREQ     16000
#define MEL_LOW_FREQ  20
#define MEL_HIGH_FREQ 8000

#define M_2PI 6.283185307179586476925286766559005

#define NUM_FILTER 20
#define NUM_FRAME  512

typedef struct _mfcc_t {
    int num_mfcc_features;
    int num_features_offset;
    int num_fbank;
    int frame_len;
    int frame_len_padded;
    int is_append_energy;
    float preempha;
    float frame[NUM_FRAME];
    float buffer[NUM_FRAME];
    float mel_energies[NUM_FILTER];
    float window_func[NUM_FRAME];
    int fbank_filter_first[NUM_FILTER];
    int fbank_filter_last[NUM_FILTER];
    float* mel_fbank[NUM_FILTER];
    float dct_matrix[NUM_FILTER * NUM_FILTER];
    float fft_buffer[NUM_FRAME * 2];
} mfcc_t;

void mfcc_create(mfcc_t* mfcc, int num_mfcc_features, int feature_offset, int num_fbank, int frame_len, float preempha,
                 int is_append_energy);

void mfcc_compute(mfcc_t* mfcc, const int16_t* audio_data, float* mfcc_out);

#ifdef __cplusplus
}
#endif

#endif
