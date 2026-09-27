#ifndef _SPECTRUM_DRIVER_H_
#define _SPECTRUM_DRIVER_H_

#include <stdint.h>

#define SPECTRUM_SAMPLE_RATE    (44100)
#define SPECTRUM_SIZE           (2048)
#define SPECTRUM_BIN_SIZE       (25)
#define SPECTRUM_LEVEL_NUM      (126)
#define SMOOTH_TAU_MS           (100)

/************************ 频谱显示参数 ************************/
typedef struct _SpectrumContext_ {
    uint16_t sample_rate;
    uint16_t spt_len;
    uint8_t bin_num;
    uint8_t level_num;
    uint8_t is_init;

    uint16_t bin_frq[SPECTRUM_BIN_SIZE];//the frequency band for every bin
    uint8_t step_table[SPECTRUM_BIN_SIZE];
    int32_t fft_buf[SPECTRUM_SIZE];
    int32_t hm[SPECTRUM_SIZE];

    //[XH]
    uint32_t spectrum_out[SPECTRUM_BIN_SIZE];
    int16_t sample_flag;
    
    uint8_t level_cur[SPECTRUM_BIN_SIZE];
    uint8_t level_dst[SPECTRUM_BIN_SIZE];
    uint8_t blocks_bottom_level[SPECTRUM_BIN_SIZE]; //blocks底部在几楼

    int16_t window_gain;
    uint32_t smooth_alpha;      // Q16.16 平滑因子
    uint32_t smooth_invalpha;   // Q16.16 1−α
    uint64_t smooth_buf[SPECTRUM_BIN_SIZE];
} SpectrumContext;
extern SpectrumContext gsc;

/**
 * @brief Initialization for Spectrum
 * 
 * @param sptct A piece of memory that user provided as context for process Spectrum
 * @param sample_rate pcm sample rate 
 * @param spt_len pcm samples for every Spectrum process, can be 128,256,512..., not larger than SPECTRUM_SIZE
 * @param bin_num bin numbers for output, supported 16..., not larger than SPECTRUM_BIN_SIZE
 * @return true 
 * @return false 
 */
uint8_t spectrum_init(SpectrumContext *sptct, uint16_t sample_rate, uint16_t spt_len, uint8_t bin_num);

/**
 * @brief Sprctrum processing routine
 * 
 * @param sptct Spectrum Context
 * @param pcm input pcm data, len is spt_len*2 bytes
 * @param bin_out output bin data, len is bin_num*2 bytes. sum of (real)^2+(img)^2 in the current bin
 * @return true 
 * @return false 
 */
uint8_t spectrum_process(SpectrumContext *sptct, int16_t *pcm, uint32_t *bin_out);

/**
 * @brief Show Frequency Interval for every bin
 *
 * @param sptct Spectrum Context
 * @return uint16_t * an arry which length is bin_num. the unit is Hz
  */
uint16_t *spectrum_get_freq_interval(SpectrumContext *sptct);

/**
 * @brief quantized 16 levels for spectrum_process's output
 *
 * @param sptct Spectrum Context
 * @param *bin_value output by spectrum_process's bin_out array.
 * @param *level_out quantized level value. shows 0~15. Please check LevelTableDB for every level's dB
 * @return void
  */
void spectrum_quantized_level(SpectrumContext *sptct, uint32_t *bin_value);

#endif
