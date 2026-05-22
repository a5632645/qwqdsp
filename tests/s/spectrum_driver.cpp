#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spectrum_driver.h"
#include <math.h>
#include <stdbool.h>

/* rfft_api 函数声明 (由 spectrum.cpp 使用 qwqdsp 浮点 FFT 实现) */
int rfft_api(int32_t *buf, int n, int inverse);

/*################################### 宏定义 ###################################*/
#define M_2PI 6.283185307179586476925286766559005

/*################################### 变量(静态声明) ###################################*/

/*################################### 变量(全局声明) ###################################*/

/*################################### 变量(静态定义) ###################################*/

static constexpr int kMul = SPECTRUM_SIZE / 256;

static const uint8_t FreqStepTable[] = {
    2 * kMul, //步进累计( 2-1)，频率172Hz
    1 * kMul, //步进累计( 3-1)，频率344Hz
    1 * kMul, //步进累计( 4-1)，频率516Hz
    1 * kMul, //步进累计( 5-1)，频率689Hz
    1 * kMul, //步进累计( 6-1)，频率861Hz
    1 * kMul, //步进累计( 7-1)，频率1033Hz
    1 * kMul, //步进累计( 8-1)，频率1205Hz
    1 * kMul, //步进累计( 9-1)，频率1378Hz
    1 * kMul, //步进累计(10-1)，频率1550Hz
    1 * kMul, //步进累计(11-1)，频率1722Hz
    1 * kMul, //步进累计(12-1)，频率1894Hz     //低频 end
    2 * kMul, //步进累计(14-1)，频率2239Hz     //高频(中间idx=11) start
    2 * kMul, //步进累计(16-1)，频率2583Hz     
    2 * kMul, //步进累计(18-1)，频率2928Hz
    2 * kMul, //步进累计(20-1)，频率3273Hz
    3 * kMul, //步进累计(23-1)，频率3789Hz
    4 * kMul, //步进累计(27-1)，频率4479Hz
    4 * kMul, //步进累计(31-1)，频率5168Hz
    4 * kMul, //步进累计(35-1)，频率5857Hz
    4 * kMul, //步进累计(39-1)，频率6546Hz
    4 * kMul, //步进累计(43-1)，频率7235Hz
    4 * kMul, //步进累计(47-1)，频率7924Hz
    4 * kMul, //步进累计(51-1)，频率8613Hz
    5 * kMul, //步进累计(56-1)，频率9475Hz
    5 * kMul, //步进累计(61-1)，频率10336Hz
};

/*################################### 变量(全局定义) ###################################*/
SpectrumContext gsc;

/*################################### 函数(静态声明) ###################################*/

/*################################### 函数(全局声明) ###################################*/

/*################################### 函数(静态定义) ###################################*/

/*################################### 函数(全局定义) ###################################*/
uint8_t spectrum_init(SpectrumContext *sptct, uint16_t sample_rate, uint16_t spt_len, uint8_t bin_num)
{
    /*** 局部变量 ***/
	int i, sum = 0;
	float frq_one_unit;

    /*** 初始化全局变量 ***/
    memset(sptct, 0, sizeof(SpectrumContext));
    sptct->sample_rate = sample_rate;
    sptct->spt_len     = spt_len;
    sptct->bin_num     = bin_num;
    sptct->is_init     = true;

    /*** 不明数据：FFT计算中用到 ***/
    float fwindow_sum = 0.0f;
    for(i=0; i<sptct->spt_len; i++)
    {
        float fwin = 0.5 - 0.5 * cos(M_2PI * i / (sptct->spt_len-1));
    	sptct->hm[i] = (int32_t)((fwin) * 32767);
        fwindow_sum += fwin;
    }
    sptct->window_gain = fwindow_sum / 2;

    /*** 频率步进表 ***/
	for(i=0; i<sptct->bin_num; i++)
    {
    	sptct->step_table[i] = FreqStepTable[i];
    }

    /*** 预计算平滑系数 (Q16.16 定点) ***/
    {
        int tau_ms = SMOOTH_TAU_MS;
        int64_t x = (int64_t)sptct->spt_len * 655360000 / tau_ms / sptct->sample_rate;
        int64_t x2 = x * x / 65536;
        int64_t x3 = x2 * x / 65536;
        int64_t a = x - x2 / 2 + x3 / 6;
        if (a < 0)     a = 0;
        if (a > 65536) a = 65536;
        sptct->smooth_alpha    = (uint32_t)a;
        sptct->smooth_invalpha = 65536 - (uint32_t)a;
    }

    /*** 计算得到：每段频率倍数表(基波频率*此倍数) ***/
    frq_one_unit = (float)sample_rate / spt_len;
    for(i=0; i<bin_num; i++)
    {
    	sptct->bin_frq[i] = frq_one_unit * (sptct->step_table[i] + sum);
        sum += sptct->step_table[i];
    }

    /*** 打印：各个频段对应的频率值 ***/
    #ifdef CFG_FUNC_DEBUG_EN
    DBG("\n---------------------------------\n");
	for (i=0; i < SPECTRUM_BIN_SIZE; i++)
	{
		DBG("%02d: %05dHz\n", i, sptct->bin_frq[i]);
	}
	DBG("---------------------------------\n");
    #endif

    return true;
}

uint8_t spectrum_process(SpectrumContext *sptct, int16_t *pcm, uint32_t *bin_out)
{
    int i, j, k;

    memset(bin_out, 0, sizeof(uint32_t) * sptct->bin_num);
    for (i = 0; i < sptct->spt_len; i++)
    {
    	sptct->fft_buf[i] = (pcm[i] * sptct->hm[i]) >> 15;
    }

    if(rfft_api(sptct->fft_buf, sptct->spt_len, 1) == 1)
    {
        sptct->fft_buf[0] = 0;
        sptct->fft_buf[1] = 0;

        int64_t psd_max = 0;
        for (i = 0, j = 0, k = 0; i < sptct->spt_len; i+=2, j++)
        {
            int64_t re = sptct->fft_buf[i];
            int64_t im = sptct->fft_buf[i+1];
            int64_t psd = (re * re) + (im * im);

            if (psd > psd_max) {
                psd_max = psd;
            }

            if(j >= sptct->step_table[k])
            {
            	j -= sptct->step_table[k];
                bin_out[k] = psd_max / sptct->window_gain / sptct->window_gain;
                psd_max = 0;
            	k++;
                if(k >= SPECTRUM_BIN_SIZE)
                {
                    break;
                }
            }
        }

        // 指数平滑 (单独循环, 系数已在 init 中预计算)
        {
            uint32_t a    = sptct->smooth_alpha;
            uint32_t ia   = sptct->smooth_invalpha;
            for (int _i = 0; _i < SPECTRUM_BIN_SIZE; ++_i) {
                uint64_t sm = (uint64_t)a * bin_out[_i]
                            + (uint64_t)ia * sptct->smooth_buf[_i];
                sm >>= 16;
                sptct->smooth_buf[_i] = sm;
                bin_out[_i] = (uint32_t)sm;
            }
        }

        return true;
    }
    else
    {
        return false;
    }
}

void spectrum_quantized_level(SpectrumContext *sptct, uint32_t *bin_value)
{
    for (int i = 0; i < SPECTRUM_BIN_SIZE; ++i) {
        float psd = (float)bin_value[i];
        if (psd <= 0.0f) {
            sptct->level_dst[i] = 0;
            continue;
        }

        // PSD → 幅度 dB: 10·log10(psd) − 20·log10(32768)  [−∞, 0] dB
        float db = 10.0f * log10f(psd) - 20.0f * log10f(32768.0f);
        if (db < -80.0f) db = -80.0f;
        if (db >   0.0f) db =   0.0f;

        // 线性映射 [−80, 0] dB → [0, 127]
        sptct->level_dst[i] = (uint8_t)((db + 80.0f) * 127.0f / 80.0f + 0.5f);
    }
}

uint16_t *spectrum_get_freq_interval(SpectrumContext *sptct)
{
	return sptct->bin_frq;
}
