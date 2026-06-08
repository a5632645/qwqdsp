#include <string.h>
#include <stdio.h>
#include <math.h>
#include "spectrum_driver.h"

int rfft_api(int32_t *buf, int n, int inverse);

/*################################### 宏定义 ###################################*/
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef DBG
// #define DBG(...) printf(__VA_ARGS__)
#define DBG(...)
#endif

/*################################### 变量(静态声明) ###################################*/

/*################################### 变量(全局声明) ###################################*/

/*################################### 变量(静态定义) ###################################*/
static const uint16_t SpectrumUserBinTab[SPECTRUM_USER_BIN_NUM] = { //数组元素范围[0, SPECTRUM_LEN/2]，其中2个不可用：0、SPECTRUM_LEN/2
     1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
    51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
    61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
    71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
    81,
};
static float window_compensation;

/*################################### 变量(全局定义) ###################################*/
SpectrumContext gsc;

/*################################### 函数(静态声明) ###################################*/

/*################################### 函数(全局声明) ###################################*/

/*################################### 函数(静态定义) ###################################*/
static inline float rfft_bin_2_amp(int32_t *fft_out, uint16_t bin)
{
    float re, im, amp;

    switch(bin)
    {
        case 0:
            re = (float)fft_out[0];
            im = 0.0f;
            amp = sqrtf(re * re + im * im) / SPECTRUM_LEN;
            break;
        case SPECTRUM_LEN / 2:
            re = (float)fft_out[1];
            im = 0.0f;
            amp = sqrtf(re * re + im * im) / SPECTRUM_LEN;
            break;
        default:
            re = (float)fft_out[2 * bin + 0];
            im = (float)fft_out[2 * bin + 1];
            amp = 2.0f * sqrtf(re * re + im * im) / SPECTRUM_LEN;
            break;
    }

    return amp * window_compensation;
}

static inline float rfft_amp_2_dbfs_raw(float amp)
{
    if(amp <= 0.0f)
    {
        return -120.0f;
    }
    return 20.0f * log10f(amp / FFT_FULL_SCALE);
}
static inline float rfft_amp_2_dbfs(float amp)
{
    float dbfs = rfft_amp_2_dbfs_raw(amp);
    if(dbfs <= FFT_DBFS_FLOOR)
    {
        return FFT_DBFS_FLOOR;
    }
    return dbfs;
}

#if SPECTRUM_ENABLE_DEBUG_MEL
static inline float hz_2_mel(float hz)
{
    return 2595.0f * log10f(1.0f + hz / 700.0f);
}
static inline float mel_2_hz(float mel)
{
    return 700.0f * (powf(10.0f, mel / 2595.0f) - 1.0f);
}
static inline float mel_filter_weight(float hz, float left_hz, float center_hz, float right_hz)
{
    if((hz <= left_hz) || (hz >= right_hz))
    {
        return 0.0f;
    }

    if(hz < center_hz)
    {
        return (hz - left_hz) / (center_hz - left_hz);
    }

    return (right_hz - hz) / (right_hz - center_hz);
}
#endif

/*################################### 函数(全局定义) ###################################*/
void spectrum_init(SpectrumContext *ct)
{
    memset(ct, 0x00, sizeof(SpectrumContext));
    ct->frame_flag = -(2000 / (256/(SPECTRUM_SAMPLE_RATE/1000))); //256数据帧大概2秒
    
    /*** 窗函数 + 窗补偿 ***/
    float sum = 0.0f;
    for(int i=0; i<SPECTRUM_LEN; i++)
    {
        ct->win[i] = (float)(0.5 - 0.5 * cos(2.0 * M_PI * i / (SPECTRUM_LEN - 1)));
        sum += ct->win[i];
    }
    window_compensation = (float)SPECTRUM_LEN / sum;
#if RFFT_SCALE_DOWN
	window_compensation *= (float)(SPECTRUM_LEN >> 1);
#endif

    /*** Mel滤波器组 ***/
#if SPECTRUM_ENABLE_DEBUG_MEL
    float mel_min = hz_2_mel(SPECTRUM_MEL_MIN_HZ);
    float mel_max = hz_2_mel(SPECTRUM_MEL_MAX_HZ);
    float mel_step = (mel_max - mel_min) / (SPECTRUM_MEL_FILTER_NUM + 1);

    for(int i=0; i<SPECTRUM_MEL_FILTER_NUM; i++)
    {
        ct->mel_left_hz[i] = mel_2_hz(mel_min + (float)i * mel_step);
        ct->mel_center_hz[i] = mel_2_hz(mel_min + (float)(i + 1) * mel_step);
        ct->mel_right_hz[i] = mel_2_hz(mel_min + (float)(i + 2) * mel_step);
    }
#endif

    /*** 打印频率 ***/
    #if 0
    for(int i=0; i<SPECTRUM_USER_BIN_NUM; i++)
    {
        uint32_t freq_hz = (uint32_t)SpectrumUserBinTab[i] * SPECTRUM_SAMPLE_RATE / SPECTRUM_LEN;

        DBG("spectrum bin[%02d] fft[%03u] = %luHz\n",
               i,
               (unsigned int)SpectrumUserBinTab[i],
               (unsigned long)freq_hz);
    }
    #endif
}

uint8_t spectrum_process(SpectrumContext *ct, fft_pcm_type_t *pcm)
{
    int32_t *buf = ct->fft_buf;
    float *win   = ct->win;
    float *dB    = ct->dbfs;

    /*** 加窗 ***/
	for(int i=0; i<SPECTRUM_LEN; i++)
	{
		buf[i] = (int32_t)((float)pcm[i] * win[i]);
	}

    /*** FFT ***/
    if(rfft_api(buf, SPECTRUM_LEN, 1) == 1)
    {
        for(int i=0; i<SPECTRUM_USER_BIN_NUM; i++)
        {
            uint16_t user_bin = SpectrumUserBinTab[i];
            dB[i] = rfft_amp_2_dbfs(rfft_bin_2_amp(buf, user_bin));
        }
        spectrum_update_level(ct);
        return TRUE;
    }
    else
    {
        return FALSE;
    }
}

void spectrum_update_level(SpectrumContext *ct)
{
    for(int i=0; i<SPECTRUM_USER_BIN_NUM; i++)
    {
        float dbfs = ct->dbfs[i];

        if(dbfs <= FFT_DBFS_FLOOR)
        {
            ct->level_cur[i] = 0;
        }
        else if(dbfs >= 0.0f)
        {
            ct->level_cur[i] = SPECTRUM_MAX_LEVEL;
        }
        else
        {
            float level = (dbfs - FFT_DBFS_FLOOR) * (float)SPECTRUM_MAX_LEVEL / (0.0f - FFT_DBFS_FLOOR);
            ct->level_cur[i] = (uint8_t)(level + 0.5f);
        }
        DBG("bin(%2d) = %3d\n", SpectrumUserBinTab[i], ct->level_cur[i]);
    }
}
















//[XH]打印bin + hz + db
#if SPECTRUM_ENABLE_DEBUG_BIN
uint8_t spectrum_process_debug(SpectrumContext *ct, fft_pcm_type_t *pcm)
{
    int32_t *buf = ct->fft_buf;
    float *win   = ct->win;

    /*** 加窗 ***/
	for(int i=0; i<SPECTRUM_LEN; i++)
	{
		buf[i] = (int32_t)((float)pcm[i] * win[i]);
	}

    /*** FFT ***/
    if(rfft_api(buf, SPECTRUM_LEN, 1) == 1)
    {
        for(int i=SPECTRUM_FFT_BIN_MIN; i<=SPECTRUM_FFT_BIN_MAX; i++)
        {
            float dB = rfft_amp_2_dbfs(rfft_bin_2_amp(buf, i));
            int Hz = ((i * SPECTRUM_SAMPLE_RATE + SPECTRUM_LEN / 2) / SPECTRUM_LEN);
            DBG("%3d, %5d, %.2f\n", i, Hz, dB);
        }
        return TRUE;
    }
    else
    {
        return FALSE;
    }
}
#endif

//[XH]峰值估算
#if SPECTRUM_ENABLE_DEBUG_PEAK
uint8_t spectrum_peak_estimate(int32_t *fft_out, SpectrumPeakEstimate *peak)
{
	uint16_t max_bin = 1;
	float max_amplitude = 0.0f;

	for(uint16_t bin = 1; bin < SPECTRUM_FFT_BIN_MAX; bin++)
	{
		float amplitude = rfft_bin_2_amp(fft_out, bin);
		if(amplitude > max_amplitude)
		{
			max_amplitude = amplitude;
			max_bin = bin;
		}
	}

    if(rfft_amp_2_dbfs_raw(max_amplitude) <= FFT_DBFS_FLOOR) //静音
    {
        return 0;
    }

    float alpha = rfft_amp_2_dbfs_raw(rfft_bin_2_amp(fft_out, max_bin - 1));
    float beta  = rfft_amp_2_dbfs_raw(rfft_bin_2_amp(fft_out, max_bin));
    float gamma = rfft_amp_2_dbfs_raw(rfft_bin_2_amp(fft_out, max_bin + 1));
	float denominator = alpha - 2.0f * beta + gamma;
	float delta = 0.0f;

	if(denominator != 0.0f)
	{
		delta = 0.5f * (alpha - gamma) / denominator;
		if(delta > 0.5f)
		{
			delta = 0.5f;
		}
		else if(delta < -0.5f)
		{
			delta = -0.5f;
		}
	}

	peak->bin  = (float)max_bin + delta;
	peak->freq = peak->bin * SPECTRUM_SAMPLE_RATE / SPECTRUM_LEN;
	peak->dbfs = beta - 0.25f * (alpha - gamma) * delta;

    DBG("\nEstimated peak: bin=%.2f, freq=%.2fHz, db=%.2f\n\n", peak->bin, peak->freq, peak->dbfs);

	return 1;
}
#endif

//[XH]打印Mel滤波器组频谱
#if SPECTRUM_ENABLE_DEBUG_MEL
void spectrum_update_mel_dbfs(SpectrumContext *ct)
{
    int32_t *buf = ct->fft_buf;

    for(int band=0; band<SPECTRUM_MEL_FILTER_NUM; band++)
    {
        float power_sum = 0.0f;
        float weight_sum = 0.0f;

        for(int bin=0; bin<=SPECTRUM_LEN/2; bin++)
        {
            float hz = (float)bin * SPECTRUM_SAMPLE_RATE / SPECTRUM_LEN;
            float weight = mel_filter_weight(hz,
                                             ct->mel_left_hz[band],
                                             ct->mel_center_hz[band],
                                             ct->mel_right_hz[band]);

            if(weight > 0.0f)
            {
                float amp = rfft_bin_2_amp(buf, bin);
                power_sum += amp * amp * weight;
                weight_sum += weight;
            }
        }

        if(weight_sum > 0.0f)
        {
            ct->mel_dbfs[band] = rfft_amp_2_dbfs(sqrtf(power_sum / weight_sum));
        }
        else
        {
            ct->mel_dbfs[band] = FFT_DBFS_FLOOR;
        }
    }
}
uint8_t spectrum_process_mel(SpectrumContext *ct, fft_pcm_type_t *pcm)
{
    int32_t *buf = ct->fft_buf;
    float *win   = ct->win;

    /*** 加窗 ***/
	for(int i=0; i<SPECTRUM_LEN; i++)
	{
		buf[i] = (int32_t)((float)pcm[i] * win[i]);
	}

    /*** FFT ***/
    if(rfft_api(buf, SPECTRUM_LEN, 1) == 1)
    {
        spectrum_update_mel_dbfs(ct);

        DBG("spectrum_process_mel()\n");
        for(int i=0; i<SPECTRUM_MEL_FILTER_NUM; i++)
        {
            DBG("%d, %d, %.2f\n", i, (int)(ct->mel_center_hz[i] + 0.5f), ct->mel_dbfs[i]);
        }
        return TRUE;
    }
    else
    {
        return FALSE;
    }
}
#endif

//测试方法
/*

#define FFT_TEST_TONE_HZ	(1750)
static fft_pcm_type_t source[SPECTRUM_LEN];
static void generate_sine_source(void)
{
	for(int i=0; i<SPECTRUM_LEN; i++)
	{
		float sample = sinf(2.0f * M_PI * FFT_TEST_TONE_HZ * i / SPECTRUM_SAMPLE_RATE);
		source[i] = (fft_pcm_type_t)roundf(sample * (FFT_FULL_SCALE - 1.0f));
	}
}



	generate_sine_source();

	{
		spectrum_init(&gsc);

#if SPECTRUM_ENABLE_DEBUG_BIN
		spectrum_process_debug(&gsc, source);
#else
		spectrum_process(&gsc, source);
#endif

#if SPECTRUM_ENABLE_DEBUG_PEAK
		SpectrumPeakEstimate peak;
		spectrum_peak_estimate(gsc.fft_buf, &peak);
#endif

#if SPECTRUM_ENABLE_DEBUG_MEL
		spectrum_process_mel(&gsc, source);
#endif
	}

*/
