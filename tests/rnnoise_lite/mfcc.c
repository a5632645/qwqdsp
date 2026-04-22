#include <stdlib.h>
#include <string.h>

#include "float.h"
#include "mfcc.h"

#define M_PI 3.14159265358979323846264338327950288

#ifndef MFCC_PLATFORM_ARM
// FFT code from arduino_fft: https://github.com/lloydroc/arduino_fft
// change to float data£¬ modify to fit within this file
// see the above link for license( MIT license).
#include <math.h>
#include <stdio.h>

static void rearrange(float data_re[], float data_im[], const unsigned int N) {
    unsigned int target = 0;
    for (unsigned int position = 0; position < N; position++) {
        if (target > position) {
            const float temp_re = data_re[target];
            const float temp_im = data_im[target];
            data_re[target] = data_re[position];
            data_im[target] = data_im[position];
            data_re[position] = temp_re;
            data_im[position] = temp_im;
        }
        unsigned int mask = N;
        while (target & (mask >>= 1)) target &= ~mask;
        target |= mask;
    }
}

static void compute(float data_re[], float data_im[], const unsigned int N) {
    const float pi = -3.14159265358979323846f;
    for (unsigned int step = 1; step < N; step <<= 1) {
        const unsigned int jump = step << 1;
        const float step_d = (float)step;
        float twiddle_re = 1.0;
        float twiddle_im = 0.0;
        for (unsigned int group = 0; group < step; group++) {
            for (unsigned int pair = group; pair < N; pair += jump) {
                const unsigned int match = pair + step;
                const float product_re = twiddle_re * data_re[match] - twiddle_im * data_im[match];
                const float product_im = twiddle_im * data_re[match] + twiddle_re * data_im[match];
                data_re[match] = data_re[pair] - product_re;
                data_im[match] = data_im[pair] - product_im;
                data_re[pair] += product_re;
                data_im[pair] += product_im;
            }
            // we need the factors below for the next iteration
            // if we don't iterate then don't compute
            if (group + 1 == step) {
                continue;
            }
            float angle = pi * ((float)group + 1) / step_d;
            twiddle_re = cosf(angle);
            twiddle_im = sinf(angle);
        }
    }
}

static void fft(float data_re[], float data_im[], const int N) {
    rearrange(data_re, data_im, N);
    compute(data_re, data_im, N);
}

#endif /* end of FFT implmentation*/

// ----------------------------------------
//
// ----------------------------------------

static void create_mel_fbank(mfcc_t* mfcc);
static void create_dct_matrix(mfcc_t* mfcc, int32_t input_length, int32_t coefficient_count);

// ----------------------------------------
//
// ----------------------------------------

static const float kMelFBank0[] = {
    0.130661f, 0.483331f, 0.821832f, 0.852739f, 0.539412f, 0.237317f,
};
static const float kMelFBank1[] = {
    0.147261f, 0.460588f, 0.762683f, 0.945677f, 0.663792f, 0.391032f, 0.126825f,
};
static const float kMelFBank2[] = {
    0.054323f, 0.336208f, 0.608968f, 0.873176f, 0.870648f, 0.622030f, 0.380537f, 0.145770f,
};
static const float kMelFBank3[] = {
    0.129352f, 0.377970f, 0.619463f, 0.854230f, 0.917367f, 0.694991f, 0.478333f, 0.267105f, 0.061042f,
};
static const float kMelFBank4[] = {
    0.082633f, 0.305009f, 0.521667f, 0.732895f, 0.938958f, 0.859898f, 0.663443f, 0.471463f, 0.283758f, 0.100146f,
};
static const float kMelFBank5[] = {
    0.140102f, 0.336557f, 0.528537f, 0.716242f, 0.899854f, 0.920446f,
    0.744500f, 0.572151f, 0.403257f, 0.237680f, 0.075294f,
};
static const float kMelFBank6[] = {
    0.079554f, 0.255500f, 0.427849f, 0.596743f, 0.762320f, 0.924706f, 0.915979f,
    0.759619f, 0.606108f, 0.455343f, 0.307227f, 0.161670f, 0.018584f,
};
static const float kMelFBank7[] = {
    0.084021f, 0.240381f, 0.393892f, 0.544657f, 0.692773f, 0.838330f, 0.981416f,
    0.877887f, 0.739500f, 0.603350f, 0.469364f, 0.337474f, 0.207618f, 0.079730f,
};
static const float kMelFBank8[] = {
    0.122113f, 0.260500f, 0.396650f, 0.530636f, 0.662526f, 0.792382f, 0.920270f, 0.953756f,
    0.829637f, 0.707318f, 0.586751f, 0.467884f, 0.350669f, 0.235064f, 0.121022f, 0.008503f,
};
static const float kMelFBank9[] = {
    0.046244f, 0.170363f, 0.292682f, 0.413249f, 0.532116f, 0.649331f, 0.764936f, 0.878978f, 0.991497f,
    0.897466f, 0.787873f, 0.679689f, 0.572873f, 0.467397f, 0.363223f, 0.260321f, 0.158661f, 0.058213f,
};
static const float kMelFBank10[] = {
    0.102534f, 0.212127f, 0.320311f, 0.427127f, 0.532603f, 0.636777f, 0.739679f, 0.841339f, 0.941787f, 0.958947f,
    0.860836f, 0.763855f, 0.667978f, 0.573178f, 0.479434f, 0.386720f, 0.295017f, 0.204298f, 0.114548f, 0.025743f,
};
static const float kMelFBank11[] = {
    0.041053f, 0.139164f, 0.236145f, 0.332022f, 0.426822f, 0.520566f, 0.613280f, 0.704983f,
    0.795702f, 0.885452f, 0.974257f, 0.937862f, 0.850891f, 0.764807f, 0.679593f, 0.595233f,
    0.511708f, 0.429004f, 0.347101f, 0.265989f, 0.185650f, 0.106069f, 0.027232f,
};
static const float kMelFBank12[] = {
    0.062138f, 0.149109f, 0.235193f, 0.320407f, 0.404767f, 0.488292f, 0.570996f, 0.652899f, 0.734011f,
    0.814350f, 0.893931f, 0.972768f, 0.949126f, 0.871739f, 0.795054f, 0.719061f, 0.643747f, 0.569101f,
    0.495110f, 0.421762f, 0.349048f, 0.276956f, 0.205475f, 0.134595f, 0.064308f,
};
static const float kMelFBank13[] = {
    0.050874f, 0.128261f, 0.204946f, 0.280939f, 0.356253f, 0.430899f, 0.504890f, 0.578238f, 0.650952f, 0.723044f,
    0.794525f, 0.865405f, 0.935692f, 0.994599f, 0.925463f, 0.856891f, 0.788873f, 0.721398f, 0.654458f, 0.588047f,
    0.522155f, 0.456774f, 0.391896f, 0.327515f, 0.263621f, 0.200208f, 0.137268f, 0.074796f, 0.012782f,
};
static const float kMelFBank14[] = {
    0.005401f, 0.074537f, 0.143109f, 0.211127f, 0.278602f, 0.345542f, 0.411953f, 0.477845f, 0.543226f,
    0.608104f, 0.672485f, 0.736379f, 0.799792f, 0.862732f, 0.925204f, 0.987218f, 0.951221f, 0.890109f,
    0.829434f, 0.769195f, 0.709381f, 0.649991f, 0.591015f, 0.532450f, 0.474288f, 0.416524f, 0.359155f,
    0.302173f, 0.245576f, 0.189353f, 0.133504f, 0.078023f, 0.022906f,
};
static const float kMelFBank15[] = {
    0.048779f, 0.109891f, 0.170566f, 0.230805f, 0.290619f, 0.350009f, 0.408985f, 0.467550f, 0.525712f,
    0.583476f, 0.640845f, 0.697827f, 0.754424f, 0.810647f, 0.866496f, 0.921977f, 0.977094f, 0.968144f,
    0.913738f, 0.859681f, 0.805968f, 0.752595f, 0.699558f, 0.646854f, 0.594476f, 0.542421f, 0.490686f,
    0.439266f, 0.388159f, 0.337362f, 0.286866f, 0.236672f, 0.186773f, 0.137170f, 0.087859f, 0.038832f,
};
static const float kMelFBank16[] = {
    0.031856f, 0.086262f, 0.140319f, 0.194032f, 0.247405f, 0.300442f, 0.353147f, 0.405524f, 0.457579f,
    0.509314f, 0.560734f, 0.611841f, 0.662638f, 0.713134f, 0.763328f, 0.813227f, 0.862830f, 0.912141f,
    0.961168f, 0.990090f, 0.941628f, 0.893442f, 0.845532f, 0.797891f, 0.750518f, 0.703412f, 0.656566f,
    0.609979f, 0.563648f, 0.517570f, 0.471741f, 0.426163f, 0.380829f, 0.335739f, 0.290888f, 0.246273f,
    0.201895f, 0.157747f, 0.113831f, 0.070143f, 0.026678f,
};
static const float kMelFBank17[] = {
    0.009910f, 0.058372f, 0.106558f, 0.154468f, 0.202109f, 0.249482f, 0.296588f, 0.343434f, 0.390021f, 0.436352f,
    0.482430f, 0.528259f, 0.573837f, 0.619171f, 0.664261f, 0.709112f, 0.753727f, 0.798105f, 0.842253f, 0.886169f,
    0.929857f, 0.973322f, 0.983440f, 0.940420f, 0.897619f, 0.855032f, 0.812663f, 0.770503f, 0.728553f, 0.686812f,
    0.645276f, 0.603945f, 0.562815f, 0.521884f, 0.481152f, 0.440615f, 0.400272f, 0.360122f, 0.320162f, 0.280391f,
    0.240805f, 0.201408f, 0.162191f, 0.123155f, 0.084299f, 0.045625f, 0.007126f,
};
static const float kMelFBank18[] = {
    0.016560f, 0.059580f, 0.102381f, 0.144968f, 0.187337f, 0.229497f, 0.271447f, 0.313188f, 0.354724f,
    0.396055f, 0.437185f, 0.478116f, 0.518848f, 0.559385f, 0.599728f, 0.639878f, 0.679838f, 0.719609f,
    0.759195f, 0.798592f, 0.837809f, 0.876845f, 0.915701f, 0.954375f, 0.992875f, 0.968800f, 0.930651f,
    0.892670f, 0.854863f, 0.817221f, 0.779748f, 0.742443f, 0.705300f, 0.668321f, 0.631503f, 0.594846f,
    0.558348f, 0.522006f, 0.485820f, 0.449789f, 0.413911f, 0.378187f, 0.342613f, 0.307189f, 0.271912f,
    0.236783f, 0.201798f, 0.166959f, 0.132265f, 0.097711f, 0.063297f, 0.029026f,
};
static const float kMelFBank19[] = {
    0.031200f, 0.069349f, 0.107330f, 0.145137f, 0.182779f, 0.220252f, 0.257557f, 0.294700f, 0.331679f, 0.368497f,
    0.405154f, 0.441652f, 0.477994f, 0.514180f, 0.550211f, 0.586089f, 0.621813f, 0.657387f, 0.692811f, 0.728088f,
    0.763217f, 0.798202f, 0.833041f, 0.867735f, 0.902289f, 0.936703f, 0.970974f, 0.994894f, 0.960896f, 0.927038f,
    0.893315f, 0.859725f, 0.826270f, 0.792944f, 0.759754f, 0.726691f, 0.693757f, 0.660954f, 0.628276f, 0.595723f,
    0.563298f, 0.530997f, 0.498818f, 0.466761f, 0.434827f, 0.403012f, 0.371317f, 0.339744f, 0.308284f, 0.276943f,
    0.245719f, 0.214610f, 0.183614f, 0.152734f, 0.121966f, 0.091308f, 0.060763f, 0.030326f,
};

static const int kMelFBankStart[] = {
    1, 4, 7, 11, 15, 20, 25, 31, 38, 45, 54, 63, 74, 86, 99, 115, 132, 151, 173, 198
};
static const int kMelFBankEnd[] = {
    6, 10, 14, 19, 24, 30, 37, 44, 53, 62, 73, 85, 98, 114, 131, 150, 172, 197, 224, 255
};

static inline float InverseMelScale(float mel_freq) {
    return 700.0f * (expf(mel_freq / 1127.0f) - 1.0f);
}

static inline float MelScale(float freq) {
    return 1127.0f * logf(1.0f + freq / 700.0f);
}

static void create_mel_fbank(mfcc_t* mfcc) {
    mfcc->mel_fbank[0] = kMelFBank0;
    mfcc->mel_fbank[1] = kMelFBank1;
    mfcc->mel_fbank[2] = kMelFBank2;
    mfcc->mel_fbank[3] = kMelFBank3;
    mfcc->mel_fbank[4] = kMelFBank4;
    mfcc->mel_fbank[5] = kMelFBank5;
    mfcc->mel_fbank[6] = kMelFBank6;
    mfcc->mel_fbank[7] = kMelFBank7;
    mfcc->mel_fbank[8] = kMelFBank8;
    mfcc->mel_fbank[9] = kMelFBank9;
    mfcc->mel_fbank[10] = kMelFBank10;
    mfcc->mel_fbank[11] = kMelFBank11;
    mfcc->mel_fbank[12] = kMelFBank12;
    mfcc->mel_fbank[13] = kMelFBank13;
    mfcc->mel_fbank[14] = kMelFBank14;
    mfcc->mel_fbank[15] = kMelFBank15;
    mfcc->mel_fbank[16] = kMelFBank16;
    mfcc->mel_fbank[17] = kMelFBank17;
    mfcc->mel_fbank[18] = kMelFBank18;
    mfcc->mel_fbank[19] = kMelFBank19;
}

static void create_dct_matrix(mfcc_t* mfcc, int32_t input_length, int32_t coefficient_count) {
    int32_t k, n;
    float normalizer;
    normalizer = sqrtf(2.0f / (float)input_length);
    for (k = 0; k < coefficient_count; k++) {
        for (n = 0; n < input_length; n++) {
            mfcc->dct_matrix[k * input_length + n] = normalizer * cosf(((float)M_PI) / input_length * (n + 0.5f) * k);
        }
    }
}

// ----------------------------------------
//
// ----------------------------------------

void mfcc_create(mfcc_t* mfcc, int num_mfcc_features, int feature_offset, int num_fbank, int frame_len, float preempha,
                 int is_append_energy) {
    mfcc->num_mfcc_features = num_mfcc_features;
    mfcc->num_features_offset = feature_offset;
    mfcc->num_fbank = num_fbank;
    mfcc->frame_len = frame_len;
    mfcc->preempha = preempha;
    mfcc->is_append_energy = is_append_energy;

    // Round-up to nearest power of 2.
    mfcc->frame_len_padded = (int)powf(2, ceilf((logf(frame_len) / logf(2))));

    memset(mfcc->frame, 0, sizeof(mfcc->frame));
    memset(mfcc->buffer, 0, sizeof(mfcc->buffer));
    memset(mfcc->mel_energies, 0, sizeof(mfcc->mel_energies));

    // create window function, hanning
    for (int i = 0; i < frame_len; i++) {
        mfcc->window_func[i] = 0.5f - 0.5f * cosf((float)M_2PI * ((float)i) / (float)(frame_len));
    }

    // create mel filterbank
    memcpy(mfcc->fbank_filter_first, kMelFBankStart, sizeof(mfcc->fbank_filter_first));
    memcpy(mfcc->fbank_filter_last, kMelFBankEnd, sizeof(mfcc->fbank_filter_last));
    create_mel_fbank(mfcc);

    // create DCT matrix
    create_dct_matrix(mfcc, mfcc->num_fbank, num_mfcc_features);

    memset(mfcc->fft_buffer, 0, sizeof(mfcc->fft_buffer));
}

void mfcc_compute(mfcc_t* mfcc, const int16_t* audio_data, float* mfcc_out) {
    int32_t i, j, bin;

    // 1. TensorFlow way of normalizing .wav data to (-1,1) and 2. do pre-emphasis.
    float last = (float)audio_data[0] / (1 << 15);
    mfcc->frame[0] = last;
    for (i = 1; i < mfcc->frame_len; i++) {
        mfcc->frame[i] = ((float)audio_data[i] - last * mfcc->preempha) / (1 << 15);
        last = (float)audio_data[i];
    }

    // Fill up remaining with zeros
    if (mfcc->frame_len_padded - mfcc->frame_len)
        memset(&mfcc->frame[mfcc->frame_len], 0, sizeof(float) * (mfcc->frame_len_padded - mfcc->frame_len));

    // windows filter
    for (i = 0; i < mfcc->frame_len; i++) {
        mfcc->frame[i] *= mfcc->window_func[i];
    }

#ifdef MFCC_PLATFORM_ARM // ToDo add other fft implementation
    // Compute FFT
    arm_rfft_fast_f32(mfcc->rfft, mfcc->frame, mfcc->buffer, 0);

    // Convert to power spectrum
    // frame is stored as [real0, realN/2-1, real1, im1, real2, im2, ...]
    int32_t half_dim = mfcc->frame_len_padded / 2;
    float first_energy = mfcc->buffer[0] * mfcc->buffer[0];
    float last_energy = mfcc->buffer[1] * mfcc->buffer[1]; // handle this special case
    for (i = 1; i < half_dim; i++) {
        float real = mfcc->buffer[i * 2];
        float im = mfcc->buffer[i * 2 + 1];
        mfcc->buffer[i] = real * real + im * im;
    }
    mfcc->buffer[0] = first_energy;
    mfcc->buffer[half_dim] = last_energy;

#else // end of ARM_fft
    // not yet optimized for memory
    float* data_re = mfcc->fft_buffer;
    float* data_im = &mfcc->fft_buffer[mfcc->frame_len_padded];

    memcpy(data_re, mfcc->frame, mfcc->frame_len_padded * sizeof(float));
    memset(data_im, 0, mfcc->frame_len_padded * sizeof(float));

    fft(data_re, data_im, mfcc->frame_len_padded);
    // only need half (N/2+1)
    for (int i = 0; i <= mfcc->frame_len_padded / 2; i++) {
        mfcc->buffer[i] = data_re[i] * data_re[i] + data_im[i] * data_im[i];
    }
#endif

    float sqrt_data;
    // Apply mel filterbanks
    for (bin = 0; bin < mfcc->num_fbank; bin++) {
        j = 0;
        float mel_energy = 0;
        int32_t first_index = mfcc->fbank_filter_first[bin];
        int32_t last_index = mfcc->fbank_filter_last[bin];
        for (i = first_index; i <= last_index; i++) {
            mel_energy += mfcc->buffer[i] * mfcc->mel_fbank[bin][j++];
        }
        mfcc->mel_energies[bin] = mel_energy / mfcc->frame_len_padded;

        // avoid log of zero
        if (mel_energy == 0.0f) mfcc->mel_energies[bin] = FLT_MIN;
    }

    // Take log
    float total_energy = 0;
    for (bin = 0; bin < mfcc->num_fbank; bin++) {
        total_energy += mfcc->mel_energies[bin];
        mfcc->mel_energies[bin] = logf(mfcc->mel_energies[bin]);
    }

    // Take DCT. Uses matrix mul.
    int out_index = 0;
    for (i = mfcc->num_features_offset; i < mfcc->num_mfcc_features; i++) {
        float sum = 0.0;
        for (j = 0; j < mfcc->num_fbank; j++) {
            sum += mfcc->dct_matrix[i * mfcc->num_fbank + j] * mfcc->mel_energies[j];
        }
        mfcc_out[out_index] = sum;
        out_index++;
    }

    // whether replace the first energy by log of total energy
    if (mfcc->is_append_energy) {
        mfcc_out[0] = logf(total_energy);
    }
}
