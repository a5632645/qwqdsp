#ifndef _SPECTRUM_DRIVER_H_
#define _SPECTRUM_DRIVER_H_

#include <cstdint>

/*################################### define/enum ###################################*/
#define SPECTRUM_SAMPLE_RATE    (48000)
#define SPECTRUM_LEN            (256)

#define SPECTRUM_FFT_BIN_MIN    (1) //0是直流分量
#define SPECTRUM_FFT_BIN_MAX    ((SPECTRUM_LEN>>1) - 1)

#define SPECTRUM_USER_BIN_NUM   (81)
#define SPECTRUM_MAX_LEVEL      (116)

//梅尔滤波器组
#define SPECTRUM_MEL_FILTER_NUM (32)
#define SPECTRUM_MEL_MIN_HZ     (0.0f)
#define SPECTRUM_MEL_MAX_HZ     (15000.0f)

//硬件FFT
#define RFFT_PCM_BIT            (16)
#if RFFT_PCM_BIT == 16
    #define RFFT_SCALE_DOWN     (0)
    typedef int16_t fft_pcm_type_t;
    #define FFT_FULL_SCALE		(32768.0f)
#else
    #define RFFT_SCALE_DOWN     (1)
    typedef int32_t fft_pcm_type_t;
    #define FFT_FULL_SCALE		(8388608.0f)
#endif
#define FFT_DBFS_FLOOR			(-60.0f)

//调试功能开关
#define SPECTRUM_ENABLE_DEBUG_BIN       (1)
#define SPECTRUM_ENABLE_DEBUG_PEAK      (1)
#define SPECTRUM_ENABLE_DEBUG_MEL       (1)

/*################################### typedef ###################################*/
typedef struct _SpectrumContext {
    int32_t fft_buf[SPECTRUM_LEN];
    float win[SPECTRUM_LEN];
    float dbfs[SPECTRUM_USER_BIN_NUM]; //单位dB

#if SPECTRUM_ENABLE_DEBUG_MEL
    float mel_dbfs[SPECTRUM_MEL_FILTER_NUM]; //单位dB
    float mel_left_hz[SPECTRUM_MEL_FILTER_NUM];
    float mel_center_hz[SPECTRUM_MEL_FILTER_NUM];
    float mel_right_hz[SPECTRUM_MEL_FILTER_NUM];
#endif

    int32_t frame_flag;
    uint8_t level_cur[SPECTRUM_USER_BIN_NUM];
    uint8_t level_tar[SPECTRUM_USER_BIN_NUM];
} SpectrumContext;

typedef struct {
	float bin;
	float freq;
	float dbfs;
} SpectrumPeakEstimate;

/*################################### 全局变量 ###################################*/
extern SpectrumContext gsc;

/*################################### 全局函数 ###################################*/
void spectrum_init(SpectrumContext *ct);
uint8_t spectrum_process(SpectrumContext *ct, fft_pcm_type_t *pcm);
void spectrum_update_level(SpectrumContext *ct);

//调试用
#if SPECTRUM_ENABLE_DEBUG_BIN
uint8_t spectrum_process_debug(SpectrumContext *ct, fft_pcm_type_t *pcm);
#endif
#if SPECTRUM_ENABLE_DEBUG_PEAK
uint8_t spectrum_peak_estimate(int32_t *fft_out, SpectrumPeakEstimate *peak);
#endif

//梅尔滤波器组
#if SPECTRUM_ENABLE_DEBUG_MEL
void spectrum_update_mel_dbfs(SpectrumContext *ct);
uint8_t spectrum_process_mel(SpectrumContext *ct, fft_pcm_type_t *pcm);
#endif

#endif
