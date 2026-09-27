// #include "librkstudio.h"
// #if defined(ENABLE_HAL)
// #include "hal_bsp.h"
// #include "hal_base.h"
// #endif

#include "AudioFile.h"
#include "work_dir.hpp"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define M_PI 3.14159265358979323846f

#define REVERB_MIN(a, b)     ((a) > (b) ? (b) : (a))
#define REVERB_MAX(a, b)     ((a) < (b) ? (b) : (a))
#define REVERB_LERP(a, b, t) ((a) + ((b) - (a)) * (t))

// ----------------------------------------
// Delayline
// ----------------------------------------

struct DelayLine {
    float* buffer;
};

// ----------------------------------------
// OnePoleFilter
// ----------------------------------------

struct OnePoleFilter {
    float lag_;
};

static void OnePoleFilter_Reset(struct OnePoleFilter* f) {
    f->lag_ = 0;
}

static float OnePoleFilter_ComputeCoeff(float w) {
    float kMaxOmega = M_PI - 1e-5f;

    if (w < 0.0f) {
        return 0.0f;
    }
    else if (w > kMaxOmega) {
        return 1.0f;
    }
    else {
        float k = tanf(w / 2);
        return k / (1 + k);
    }
}

static float OnePoleFilter_TickLowpass(struct OnePoleFilter* f, float x, float coeff) {
    float lag = f->lag_;
    float delta = coeff * (x - lag);
    lag += delta;
    float y = lag;
    lag += delta;
    f->lag_ = lag;
    return y;
}

// ----------------------------------------
// Reverb
// ----------------------------------------

struct ReverbParam {
    float chorus_amount;   // [0, 1]
    float chorus_freq;     // [0.003, 8.0]
    float wet;             // [0, 1]
    float pre_lowpass;     // [0, 130]
    float pre_highpass;    // [0, 130]
    float low_damp_pitch;  // [0, 130]
    float high_damp_pitch; // [0, 130]
    float low_damp_db;     // [-6, 0]
    float high_damp_db;    // [-6, 0]
    float size;            // [0, 1]
    float decay_ms;        // [15ms, 64s]
    float pre_delay;       // [0, 300ms]
};

static void ReverbParam_SetDefault(struct ReverbParam* self) {
    self->chorus_amount = 0.05f;
    self->chorus_freq = 0.25f;
    self->wet = 0.25f;
    self->pre_lowpass = 0.0f;
    self->pre_highpass = 110.0f;
    self->low_damp_pitch = 0.0f;
    self->high_damp_pitch = 90.0f;
    self->low_damp_db = 0.0f;
    self->high_damp_db = -1.0f;
    self->size = 0.5f;
    self->decay_ms = 1000.0f;
    self->pre_delay = 0.0f;
}

struct Reverb {
    struct ReverbParam param;

    // allpass network
    float* allpass_lookups_;
    int poly_allpass_mask_;
    int allpass_write_pos_;
    int max_allpass_size_;

    // delay network
    float* feedback_memorie_;
    float* feedback_lookup_[16];
    int write_index_;
    int feedback_mask_;
    int max_feedback_size_;
    float decays_[16];
    float feedback_offsets_[16];

    struct OnePoleFilter low_shelf_filters_[16];
    struct OnePoleFilter high_shelf_filters_[16];

    float low_coefficient_;
    float low_amplitude_;
    float high_coefficient_;
    float high_amplitude_;

    // pre filter
    struct OnePoleFilter low_pre_filter_[2];
    struct OnePoleFilter high_pre_filter_[2];

    float low_pre_coefficient_;
    float high_pre_coefficient_;

    // chorus
    float chorus_phase_;
    float chorus_amount_[4];

    // predelay
    float sample_delay_;
    float sample_delay_increment_;

    // mix
    float dry_;
    float wet_;

    // other
    float fs_;
    float fs_ratio_;
    float buffer_scale_ratio_;
};

static const float kAllpassDelays[] = {1001, 799,  933, 876, 895, 807, 907, 853,
                                       957,  1019, 711, 567, 833, 779, 663, 997};
static const float kFeedbackDelays[] = {6753.2f,  9278.4f, 7704.5f,  11328.5f, 9701.12f, 5512.5f,  8480.45f, 5638.65f,
                                        3120.73f, 3429.5f, 3626.37f, 7713.52f, 4521.54f, 6518.97f, 5265.56f, 5630.25};

static const float kT60Amplitude = 0.001f;
static const float kAllpassFeedback = 0.6f;
static const float kMinDelay = 3.0f;

static const int kBaseSampleRate = 44100;
static const int kDefaultSampleRate = 88200;
static const int kNetworkSize = 16;
static const int kBaseFeedbackBits = 14;
static const int kExtraLookupSample = 4;
static const int kBaseAllpassBits = 10;
static const int kMinSizePower = -3;
static const int kMaxSizePower = 1;
static const float kSizePowerRange = kMaxSizePower - kMinSizePower;

static const float kMaxChorusDrift = 2500.0f;
static const float kMinDecayTime = 0.1f;
static const float kMaxDecayTime = 100.0f;
static const float kMaxChorusFrequency = 16.0f;
static const float kChorusShiftAmount = 0.9f;
static const float kSampleDelayMultiplier = 0.05f;
static const float kSampleIncrementMultiplier = 0.05f;

// -------------------- reverb private --------------------

static float _Pitch2Freq(float pitch) {
    return 440.0f * exp2f((pitch - 69.0f) / 12.0f);
}

static float _Freq2W(float f, float fs) {
    return f * M_PI * 2 / fs;
}

static float _Db2Gain(float db) {
    return powf(10.0f, db / 20.0f);
}

static void Reverb_WriteFeedback(struct Reverb* self, float const src[16]) {
    int const idx = self->write_index_;
    for (int i = 0; i < 16; ++i) {
        float* ptr = self->feedback_lookup_[i];
        ptr[idx] = src[i];
    }
}

static void Reverb_ReadFeedback(struct Reverb* self, float dst[16], float offset[16]) {
    float const warp_add = (float)(self->write_index_ + self->feedback_mask_);
    for (int i = 0; i < 16; ++i) {
        float rpos = warp_add - offset[i];
        int idx = ((int)(rpos)-1) & self->feedback_mask_;
        float t = rpos - (float)(int)rpos;

        float const* ptr = self->feedback_lookup_[i];
        float yn1 = ptr[idx];
        float y0 = ptr[idx + 1];
        float y1 = ptr[idx + 2];
        float y2 = ptr[idx + 3];

        float d0 = (y1 - yn1) * (0.5f);
        float d1 = (y2 - y0) * (0.5f);
        float d = y1 - y0;
        float m0 = (3.0f) * d - (2.0f) * d0 - d1;
        float m1 = d0 - (2.0f) * d + d1;
        dst[i] = y0 + t * (d0 + t * (m0 + t * m1));
    }
}

static void Reverb_ReadAllpass(struct Reverb* self, float dst[16], int const offset[16]) {
    int const warp_add = self->allpass_write_pos_ + self->poly_allpass_mask_;
    for (int i = 0; i < 16; ++i) {
        int rpos = (warp_add - offset[i]) & self->poly_allpass_mask_;
        int idx = rpos * 16 + i;
        dst[i] = self->allpass_lookups_[idx];
    }
}

static void Reverb_WriteAllpass(struct Reverb* self, float const src[16]) {
    float* dst = self->allpass_lookups_ + self->allpass_write_pos_ * 16;
    for (int i = 0; i < 16; ++i) {
        dst[i] = src[i];
    }
}

static void Reverb_Scatter(float const src[16], float dst[16]) {
    float t[4];
    t[0] = 0.5f * (src[0] + src[4] + src[8] + src[12]);
    t[1] = 0.5f * (src[1] + src[5] + src[9] + src[13]);
    t[2] = 0.5f * (src[2] + src[6] + src[10] + src[14]);
    t[3] = 0.5f * (src[3] + src[7] + src[11] + src[15]);

    float s[4];
    s[0] = 0.5f * (src[0] + src[1] + src[2] + src[3]);
    s[1] = 0.5f * (src[4] + src[5] + src[6] + src[7]);
    s[2] = 0.5f * (src[8] + src[9] + src[10] + src[11]);
    s[3] = 0.5f * (src[12] + src[13] + src[14] + src[15]);

    float sum_all = 0.5f * (s[0] + s[1] + s[2] + s[3]);

    dst[0] = sum_all + src[0] - s[0] - t[0];
    dst[1] = sum_all + src[1] - s[0] - t[1];
    dst[2] = sum_all + src[2] - s[0] - t[2];
    dst[3] = sum_all + src[3] - s[0] - t[3];

    dst[4] = sum_all + src[4] - s[1] - t[0];
    dst[5] = sum_all + src[5] - s[1] - t[1];
    dst[6] = sum_all + src[6] - s[1] - t[2];
    dst[7] = sum_all + src[7] - s[1] - t[3];

    dst[8] = sum_all + src[8] - s[2] - t[0];
    dst[9] = sum_all + src[9] - s[2] - t[1];
    dst[10] = sum_all + src[10] - s[2] - t[2];
    dst[11] = sum_all + src[11] - s[2] - t[3];

    dst[12] = sum_all + src[12] - s[3] - t[0];
    dst[13] = sum_all + src[13] - s[3] - t[1];
    dst[14] = sum_all + src[14] - s[3] - t[2];
    dst[15] = sum_all + src[15] - s[3] - t[3];
}

static void Reverb_WarpBuffer(struct Reverb* self) {
    int const where = self->max_feedback_size_;
    for (int i = 0; i < 16; ++i) {
        float* ptr = self->feedback_lookup_[i];
        ptr[where] = ptr[0];
        ptr[where + 1] = ptr[1];
        ptr[where + 2] = ptr[2];
        ptr[where + 3] = ptr[3];
    }
}

// -------------------- reverb public --------------------

static void Reverb_Init(struct Reverb* self, float fs) {
    self->fs_ = fs;
    self->fs_ratio_ = fs / kBaseSampleRate;

    // self.predelay_.Init(fs, 300.0f);

    self->low_pre_coefficient_ = 0.1f;
    self->high_pre_coefficient_ = 0.1f;
    self->low_coefficient_ = 0.1f;
    self->high_coefficient_ = 0.1f;
    self->low_amplitude_ = 0.0f;
    self->high_amplitude_ = 0.0f;
    self->sample_delay_ = kMinDelay;

    float const buffer_scale_ratio = fs / kBaseSampleRate;
    self->buffer_scale_ratio_ = buffer_scale_ratio;
    int const base_feedback_size = ceilf(buffer_scale_ratio * (1 << (kBaseFeedbackBits + kMaxSizePower)));
    int max_feedback_size = 1;
    while (max_feedback_size < base_feedback_size) {
        max_feedback_size *= 2;
    }
    self->max_feedback_size_ = (int)(max_feedback_size);
    self->feedback_mask_ = self->max_feedback_size_ - 1;

    int each_size = max_feedback_size + kExtraLookupSample;
    self->feedback_memorie_ = (float*)malloc(each_size * kNetworkSize * sizeof(float));
    for (int i = 0; i < kNetworkSize; ++i) {
        self->feedback_lookup_[i] = &self->feedback_memorie_[i * each_size];
    }

    int const base_allpass_size = ceilf(buffer_scale_ratio * (1 << kBaseAllpassBits));
    int max_allpass_size = 1;
    while (max_allpass_size < base_allpass_size) {
        max_allpass_size *= 2;
    }
    self->max_allpass_size_ = (int)(max_allpass_size);
    self->poly_allpass_mask_ = self->max_allpass_size_ - 1;
    self->allpass_lookups_ = (float*)malloc(max_allpass_size * kNetworkSize * sizeof(float));

    self->write_index_ = 0;
    self->allpass_write_pos_ = 0;
}

static void Reverb_FreeMeme(struct Reverb* self) {
    free(self->feedback_memorie_);
    self->feedback_memorie_ = NULL;

    free(self->allpass_lookups_);
    self->allpass_lookups_ = NULL;

    for (int i = 0; i < 16; ++i) {
        self->feedback_lookup_[i] = NULL;
    }
}

static void Reverb_Reset(struct Reverb* self) {
    self->wet_ = 0;
    self->dry_ = 0;

    for (int i = 0; i < 4; ++i) {
        self->chorus_amount_[i] = self->param.chorus_amount * kMaxChorusDrift;
    }

    for (int i = 0; i < 16; ++i) {
        OnePoleFilter_Reset(&self->low_shelf_filters_[i]);
    }
    for (int i = 0; i < 16; ++i) {
        OnePoleFilter_Reset(&self->high_shelf_filters_[i]);
    }

    for (int i = 0; i < 2; ++i) {
        OnePoleFilter_Reset(&self->low_pre_filter_[i]);
    }
    for (int i = 0; i < 2; ++i) {
        OnePoleFilter_Reset(&self->high_pre_filter_[i]);
    }

    for (int i = 0; i < 16; ++i) {
        self->decays_[i] = 0;
    }

    for (size_t i = 0; i < 16; ++i) {
        self->feedback_offsets_[i] = kFeedbackDelays[i];
    }

    memset(self->allpass_lookups_, 0, sizeof(float) * 16 * self->max_allpass_size_);
    memset(self->feedback_memorie_, 0, sizeof(float) * 16 * (self->max_feedback_size_ + kExtraLookupSample));
}

static void Reverb_Update(struct Reverb* self) {}

static void Reverb_ProcessStereo(struct Reverb* self, float* left, float* right, int num_samples) {
    Reverb_WarpBuffer(self);

    struct ReverbParam* param = &self->param;
    float const tick_increment = 1.0f / (float)(num_samples);

    float current_dry = self->dry_;
    float current_wet = self->wet_;
    float current_low_pre_coefficient = self->low_pre_coefficient_;
    float current_high_pre_coefficient = self->high_pre_coefficient_;
    float current_low_coefficient = self->low_coefficient_;
    float current_low_amplitude = self->low_amplitude_;
    float current_high_coefficient = self->high_coefficient_;
    float current_high_amplitude = self->high_amplitude_;

    self->wet_ = sinf(param->wet * M_PI / 2);
    self->dry_ = cosf(param->wet * M_PI / 2);
    float delta_wet = (self->wet_ - current_wet) * tick_increment;
    float delta_dry = (self->dry_ - current_dry) * tick_increment;

    float const low_pre_cutoff_frequency = _Pitch2Freq(param->pre_lowpass);
    self->low_pre_coefficient_ = OnePoleFilter_ComputeCoeff(_Freq2W(low_pre_cutoff_frequency, self->fs_));
    float const high_pre_cutoff_frequency = _Pitch2Freq(param->pre_highpass);
    self->high_pre_coefficient_ = OnePoleFilter_ComputeCoeff(_Freq2W(high_pre_cutoff_frequency, self->fs_));
    float delta_low_pre_coefficient = (self->low_pre_coefficient_ - current_low_pre_coefficient) * tick_increment;
    float delta_high_pre_coefficient = (self->high_pre_coefficient_ - current_high_pre_coefficient) * tick_increment;

    float const low_cutoff_frequency = _Pitch2Freq(param->low_damp_pitch);
    self->low_coefficient_ = OnePoleFilter_ComputeCoeff(_Freq2W(low_cutoff_frequency, self->fs_));
    float const high_cutoff_frequency = _Pitch2Freq(param->high_damp_pitch);
    self->high_coefficient_ = OnePoleFilter_ComputeCoeff(_Freq2W(high_cutoff_frequency, self->fs_));
    float delta_low_coefficient = (self->low_coefficient_ - current_low_coefficient) * tick_increment;
    float delta_high_coefficient = (self->high_coefficient_ - current_high_coefficient) * tick_increment;

    self->low_amplitude_ = 1.0f - _Db2Gain(param->low_damp_db);
    self->high_amplitude_ = _Db2Gain(param->high_damp_db);
    float delta_low_amplitude = (self->low_amplitude_ - current_low_amplitude) * tick_increment;
    float delta_high_amplitude = (self->high_amplitude_ - current_high_amplitude) * tick_increment;

    float const size_mult = exp2f(param->size * kSizePowerRange + kMinSizePower);

    float const decay_samples = param->decay_ms * kBaseSampleRate / 1000.0f;
    float const decay_period = size_mult / decay_samples;

    float current_decay[16];
    float delta_decay[16];
    for (int i = 0; i < 16; ++i) {
        current_decay[i] = self->decays_[i];
    }
    for (int i = 0; i < 16; ++i) {
        self->decays_[i] = powf(kT60Amplitude, kFeedbackDelays[i] * decay_period);
    }
    for (int i = 0; i < 16; ++i) {
        delta_decay[i] = (self->decays_[i] - current_decay[i]) * tick_increment;
    }

    int allpass_offset[16];
    for (int i = 0; i < 16; ++i) {
        allpass_offset[i] = (int)(kAllpassDelays[i] * self->buffer_scale_ratio_);
    }

    float const chorus_phase_increment = param->chorus_freq / self->fs_;

    float const network_offset = 2.0f * M_PI / kNetworkSize;
    float phase_offset[4] = {network_offset * 0, network_offset * 1, network_offset * 2, network_offset * 3};

    float container_phase[4];
    for (int i = 0; i < 4; ++i) {
        container_phase[i] = phase_offset[i] + self->chorus_phase_ * 2.0f * M_PI;
    }
    self->chorus_phase_ += (float)(num_samples)*chorus_phase_increment;
    self->chorus_phase_ -= floorf(self->chorus_phase_);

    float chorus_increment_real = cosf(chorus_phase_increment * (2.0f * M_PI));
    float chorus_increment_imaginary = sinf(chorus_phase_increment * (2.0f * M_PI));

    float current_chorus_real[4];
    float current_chorus_imaginary[4];
    for (size_t i = 0; i < 4; ++i) {
        current_chorus_real[i] = cosf(container_phase[i]);
    }
    for (size_t i = 0; i < 4; ++i) {
        current_chorus_imaginary[i] = sinf(container_phase[i]);
    }

    float delay14[16];
    for (int i = 0; i < 16; ++i) {
        delay14[i] = kFeedbackDelays[i] * (self->fs_ / kBaseSampleRate * size_mult);
    }

    float current_chorus_amount[4];
    for (int i = 0; i < 4; ++i) {
        current_chorus_amount[i] = self->chorus_amount_[i];
    }

    for (int i = 0; i < 4; ++i) {
        self->chorus_amount_[i] = param->chorus_amount * kMaxChorusDrift * self->fs_ratio_;
    }
    for (int i = 0; i < 4; ++i) {
        self->chorus_amount_[i] = REVERB_MIN(self->chorus_amount_[i], delay14[0 + i] - 8 * 4);
    }
    for (int i = 0; i < 4; ++i) {
        self->chorus_amount_[i] = REVERB_MIN(self->chorus_amount_[i], delay14[4 + i] - 8 * 4);
    }
    for (int i = 0; i < 4; ++i) {
        self->chorus_amount_[i] = REVERB_MIN(self->chorus_amount_[i], delay14[8 + i] - 8 * 4);
    }
    for (int i = 0; i < 4; ++i) {
        self->chorus_amount_[i] = REVERB_MIN(self->chorus_amount_[i], delay14[12 + i] - 8 * 4);
    }

    float delta_chorus_amount[4];
    for (int i = 0; i < 4; ++i) {
        delta_chorus_amount[i] = (self->chorus_amount_[i] - current_chorus_amount[i]) * tick_increment;
    }
    for (int i = 0; i < 4; ++i) {
        current_chorus_amount[i] *= size_mult;
    }

    float current_sample_delay = self->sample_delay_;
    float current_delay_increment = self->sample_delay_increment_;
    float end_target = current_sample_delay + current_delay_increment * (float)(num_samples);
    float target_delay = REVERB_MAX(kMinDelay, param->pre_delay * self->fs_ / 1000.0f);
    target_delay = REVERB_LERP(self->sample_delay_, target_delay, kSampleDelayMultiplier);
    float makeup_delay = target_delay - end_target;
    float delta_delay_increment =
        makeup_delay / (0.5f * (float)(num_samples * num_samples)) * kSampleIncrementMultiplier;

    float feedback_offset[16];
    for (int i = 0; i < 16; ++i) {
        feedback_offset[i] = self->feedback_offsets_[i];
    }

    for (int i = 0; i < num_samples; ++i) {
        // paralle chorus delaylines
        for (int j = 0; j < 4; ++j) {
            current_chorus_amount[j] += delta_chorus_amount[j];
        }
        for (int j = 0; j < 4; ++j) {
            current_chorus_real[j] = current_chorus_real[j] * chorus_increment_real
                                   - current_chorus_imaginary[j] * chorus_increment_imaginary;
        }
        for (int j = 0; j < 4; ++j) {
            current_chorus_imaginary[j] = current_chorus_imaginary[j] * chorus_increment_real
                                        + current_chorus_real[j] * chorus_increment_imaginary;
        }

        for (int j = 0; j < 4; ++j) {
            feedback_offset[j + 0] = delay14[0 + j] + current_chorus_real[j] * current_chorus_amount[j];
        }
        for (int j = 0; j < 4; ++j) {
            feedback_offset[j + 4] = delay14[4 + j] - current_chorus_real[j] * current_chorus_amount[j];
        }
        for (int j = 0; j < 4; ++j) {
            feedback_offset[j + 8] = delay14[8 + j] + current_chorus_imaginary[j] * current_chorus_amount[j];
        }
        for (int j = 0; j < 4; ++j) {
            feedback_offset[j + 12] = delay14[12 + j] - current_chorus_imaginary[j] * current_chorus_amount[j];
        }

        float feedback_read[16];
        Reverb_ReadFeedback(self, feedback_read, feedback_offset);

        float leftx = left[i];
        float rightx = right[i];
        leftx = OnePoleFilter_TickLowpass(&self->high_pre_filter_[0], leftx, current_high_pre_coefficient);
        rightx = OnePoleFilter_TickLowpass(&self->high_pre_filter_[1], rightx, current_high_pre_coefficient);
        leftx = OnePoleFilter_TickLowpass(&self->low_pre_filter_[0], leftx, current_low_pre_coefficient) - leftx;
        rightx = OnePoleFilter_TickLowpass(&self->low_pre_filter_[1], rightx, current_low_pre_coefficient) - rightx;

        // paralle polyphase allpass
        float allpass_read[16];
        Reverb_ReadAllpass(self, allpass_read, allpass_offset);

        float allpass_delay_input[16];
        for (int j = 0; j < 16; ++j) {
            allpass_delay_input[j] = feedback_read[j] - allpass_read[j] * kAllpassFeedback;
        }

        float allpass_write[16];
        allpass_write[0] = allpass_delay_input[0] + leftx * 0.5f;
        allpass_write[1] = allpass_delay_input[1] + rightx * 0.5f;
        allpass_write[2] = allpass_delay_input[2] + leftx * 0.5f;
        allpass_write[3] = allpass_delay_input[3] + rightx * 0.5f;
        allpass_write[4] = allpass_delay_input[4] + leftx * 0.5f;
        allpass_write[5] = allpass_delay_input[5] + rightx * 0.5f;
        allpass_write[6] = allpass_delay_input[6] + leftx * 0.5f;
        allpass_write[7] = allpass_delay_input[7] + rightx * 0.5f;
        allpass_write[8] = allpass_delay_input[8] + leftx * 0.5f;
        allpass_write[9] = allpass_delay_input[9] + rightx * 0.5f;
        allpass_write[10] = allpass_delay_input[10] + leftx * 0.5f;
        allpass_write[11] = allpass_delay_input[11] + rightx * 0.5f;
        allpass_write[12] = allpass_delay_input[12] + leftx * 0.5f;
        allpass_write[13] = allpass_delay_input[13] + rightx * 0.5f;
        allpass_write[14] = allpass_delay_input[14] + leftx * 0.5f;
        allpass_write[15] = allpass_delay_input[15] + rightx * 0.5f;
        Reverb_WriteAllpass(self, allpass_write);
        self->allpass_write_pos_ = (self->allpass_write_pos_ + 1) & self->poly_allpass_mask_;

        float allpass_output[16];
        for (int j = 0; j < 16; ++j) {
            allpass_output[j] = allpass_read[j] + allpass_delay_input[j] * kAllpassFeedback;
        }

        // scatter matrix
        // write1 = 0.25 * Sum(All) + Ao0 - 0.5 * Sum(Ao0) - 0.5 * (Ao0 + Ao1 + Ao2 + Ao3)
        // write2 = 0.25 * Sum(All) + Ao1 - 0.5 * Sum(Ao1) - 0.5 * (Ao0 + Ao1 + Ao2 + Ao3)
        // write3 = 0.25 * Sum(All) + Ao2 - 0.5 * Sum(Ao2) - 0.5 * (Ao0 + Ao1 + Ao2 + Ao3)
        // write4 = 0.25 * Sum(All) + Ao3 - 0.5 * Sum(Ao3) - 0.5 * (Ao0 + Ao1 + Ao2 + Ao3)
        float write[16];
        Reverb_Scatter(allpass_output, write);

        // damp filter
        float high_filter[16];
        for (int j = 0; j < 16; ++j) {
            high_filter[j] =
                OnePoleFilter_TickLowpass(&self->high_shelf_filters_[j], write[j], current_high_coefficient);
        }
        for (int j = 0; j < 16; ++j) {
            write[j] = high_filter[j] + (current_high_amplitude) * (write[j] - high_filter[j]);
        }

        float low_filter[16];
        for (int j = 0; j < 16; ++j) {
            low_filter[j] = OnePoleFilter_TickLowpass(&self->low_shelf_filters_[j], write[j], current_low_coefficient);
        }
        for (int j = 0; j < 16; ++j) {
            write[j] -= low_filter[j] * current_low_amplitude;
        }

        // decay block
        for (int j = 0; j < 16; ++j) {
            current_decay[j] += delta_decay[j];
        }

        float store[16];
        for (int j = 0; j < 16; ++j) {
            store[j] = current_decay[j] * write[j];
        }
        self->write_index_ = (self->write_index_ + 1) & self->feedback_mask_;
        Reverb_WriteFeedback(self, store);

        // scatter matrix
        float feed_forward[16];
        Reverb_Scatter(store, feed_forward);

        // predelay
        float lefty = 0;
        float righty = 0;
        for (int j = 0; j < 8; ++j) {
            lefty += write[2 * j];
            lefty += current_decay[2 * j] * feed_forward[2 * j] * 0.125f;
        }
        for (int j = 0; j < 8; ++j) {
            righty += write[2 * j + 1];
            righty += current_decay[2 * j + 1] * feed_forward[2 * j + 1] * 0.125f;
        }

        // self.predelay_.Push(output);
        // auto audio_out = current_wet * self.predelay_.GetAfterPush(current_sample_delay) + current_dry * input;
        // left[i] = audio_out[0];
        // if constexpr (!kMono) {
        // right[i] = audio_out[1];
        // }

        left[i] = current_wet * lefty + current_dry * leftx;
        right[i] = current_wet * righty + current_dry * rightx;

        current_delay_increment += delta_delay_increment;
        current_sample_delay += current_delay_increment;
        current_sample_delay = REVERB_MAX(current_sample_delay, kMinDelay);
        current_dry += delta_dry;
        current_wet += delta_wet;
        current_high_coefficient += delta_high_coefficient;
        current_high_amplitude += delta_high_amplitude;
        current_low_pre_coefficient += delta_low_pre_coefficient;
        current_high_pre_coefficient += delta_high_pre_coefficient;
        current_low_coefficient += delta_low_coefficient;
        current_low_amplitude += delta_low_amplitude;
    }

    self->sample_delay_increment_ = current_delay_increment;
    self->sample_delay_ = current_sample_delay;
    for (int i = 0; i < 16; ++i) {
        self->feedback_offsets_[i] = feedback_offset[i];
    }
}

// int Custom_Algo_SetParam(BlockAlgo* algoSt, void* param, uint32_t nParamOff, uint32_t nParamSize)
// {
//     int chanNum = algoSt->nOutputs;
//     float *pParam = (float *)algoSt->pParam;
//     CustomAlgoSt *st = (CustomAlgoSt *)algoSt->pAlgoCore;

//     if (param)
//     {
//         memcpy(pParam + nParamOff, param, nParamSize);
//     }
//     st->gain = pParam[1];
//     // todo use the pParam
//     return 0;
// }

// int Custom_Algo_GetAlgoSize(BlockAlgo *algoSt)
// {
//     return sizeof(CustomAlgoSt);
// }

// int Custom_Algo_Init(BlockAlgo *algoSt)
// {

//     CustomAlgoSt *st = (CustomAlgoSt *)algoSt->pAlgoCore;

//     // algoSt->pParam ready now, don't need care
//     return Custom_Algo_SetParam(algoSt, NULL, 0, 0);
// }

// int Custom_Algo_Process(BlockAlgo *algoSt, float **dataIn, float **dataOut)
// {
//     CustomAlgoSt *st = (CustomAlgoSt *)algoSt->pAlgoCore;
//     int chanNum = algoSt->nInputs;
//     int blockSize = algoSt->pBlockPro->nBlockSize;
//     for (int chan = 0; chan < chanNum; chan++)
//     {
//         for (int i = 0; i < blockSize; i++)
//             // just test for scale input data
//             dataOut[chan][i] = dataIn[chan][i] * st->gain;
//     }
//     return 0;
// }

// int Custom_Algo_Destroy(BlockAlgo *algoSt)
// {
//     return 0;
// }

// int Custom_Algo_GetParam(BlockAlgo* algoSt, void* pParam, uint32_t nParamOff, uint32_t *nParamSize)
// {
//     memcpy(pParam, algoSt->pParam + nParamOff, *nParamSize);
//     for (int i = nParamOff; i < *nParamSize / 4; i++)
//         LOGD("param_%d:%f\n", i, algoSt->pParam[i]);
//     return 0;
// }

// int librkstudio_register(void *param)
// {
//     int ret;
//     RkStudioFunc func =
//     {
//         .algoType        = 0x1101,                       // must same with algoType in custom_algo.json
//         .algoInit        = Custom_Algo_Init,
//         .algoSetParam    = Custom_Algo_SetParam,
//         .algoGetParam    = Custom_Algo_GetParam,
//         .algoGetAlgoSize = Custom_Algo_GetAlgoSize,      // You can also set this function as NULL
//         .algoProcess     = Custom_Algo_Process,
//         .algoDestroy     = Custom_Algo_Destroy,
//     };

//     ret = RKST_ALGO_Register(&func);
//     LOGD("register algoType = %d, ret = %d\n", func.algoType, ret);
//     return ret;
// }

int main() {
    AudioFile<float> file;
    file.load(qwqdsp_support::WormholeWav());
    auto leftx = file.samples.front();
    auto rightx = leftx;

    Reverb dsp;
    Reverb_Init(&dsp, file.getSampleRate());
    Reverb_Reset(&dsp);
    ReverbParam_SetDefault(&dsp.param);
    dsp.param.wet = 1.0f;
    dsp.param.decay_ms = 8000.0f;

    constexpr int block_size = 256;
    int wpos = 0;
    int size = leftx.size();
    while (size != 0) {
        Reverb_WarpBuffer(&dsp);
        int can_do = REVERB_MIN(size, block_size);
        Reverb_ProcessStereo(&dsp, leftx.data() + wpos, rightx.data() + wpos, can_do);
        wpos += can_do;
        size -= can_do;
    }

    file.setNumChannels(2);
    file.samples.clear();
    file.samples.reserve(2);
    file.samples.emplace_back(std::move(leftx));
    file.samples.emplace_back(std::move(rightx));

    file.save(qwqdsp_support::OutputFile("reverb.wav"));
}
