/*
 * Copyright (c) 2018-2020
 * Jianjia Ma
 * majianjia@live.com
 *
 * SPDX-License-Identifier: Apache-2.0
 * 
 * Notice: 
 * Code in this file inlcudes derivative works from CMSIS, which is released under alternative license.
 * Please check the LICENSE file for detial.
 *
 * Change Logs:
 * Date           Author       Notes
 * 2019-02-05     Jianjia Ma   The first version
 * 2019-03-19     Jianjia Ma   Local C implementation partly from CMSIS-NN
 */

#ifndef __NNOM_LOCAL_H__
#define __NNOM_LOCAL_H__

#ifdef __cplusplus
extern "C" {
#endif


#include "stdint.h"
#include "../port/nnom_port.h"
#include "nnom.h"

#ifdef ARM_NN_TRUNCATE
#define NNOM_TRUNCATE
#endif

// SSAT implementation with C code
#ifndef __NNOM_SSAT
static inline int __NNOM_SSAT(int32_t value, int32_t bit) {
    int32_t min = -(1<<(bit-1));
    int32_t max = (1<<(bit-1)) - 1;
    if (value < min)
        return min;
    else if (value > max)
        return max;
    else
        return value;
}
#endif

// USAT implementation with C code
#ifndef __NNOM_USAT
static inline int __NNOM_USAT(int32_t value, int32_t bit) {
    int32_t max = (1<<(bit-1)) - 1;
    if (value < 0)
        return 0;
    else if (value > max)
        return max;
    else
        return value;
}
#endif

#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define MIN(A, B) ((A) < (B) ? (A) : (B))


// Those functions/tables below are partially modifed from CMSIS-NN lib
// https://github.com/ARM-software/CMSIS_5
//
void local_avepool_q7_HWC(const q7_t *Im_in,           // input image
	const uint16_t dim_im_in_x,  // input image dimension x or W
	const uint16_t dim_im_in_y,  // input image dimension y or H
	const uint16_t ch_im_in,     // number of input image channels
	const uint16_t dim_kernel_x, // window kernel size
	const uint16_t dim_kernel_y, // window kernel size
	const uint16_t padding_x,    // padding sizes
	const uint16_t padding_y,    // padding sizes
	const uint16_t stride_x,     // stride
	const uint16_t stride_y,     // stride
	const uint16_t dim_im_out_x, // output image dimension x or W
	const uint16_t dim_im_out_y, // output image dimension y or H
	const uint16_t output_shift, // output right shift
	q7_t *bufferA,               // a buffer for local storage, NULL by now
	q7_t *Im_out);

void local_avepool_q7_CHW(const q7_t *Im_in,           // input image
	const uint16_t dim_im_in_x,  // input image dimension x or W
	const uint16_t dim_im_in_y,  // input image dimension y or H
	const uint16_t ch_im_in,     // number of input image channels
	const uint16_t dim_kernel_x, // window kernel size
	const uint16_t dim_kernel_y, // window kernel size
	const uint16_t padding_x,    // padding sizes
	const uint16_t padding_y,    // padding sizes
	const uint16_t stride_x,     // stride
	const uint16_t stride_y,     // stride
	const uint16_t dim_im_out_x, // output image dimension x or W
	const uint16_t dim_im_out_y, // output image dimension y or H
	const uint16_t output_shift, // output right shift
	q7_t *bufferA,               // a buffer for local storage, NULL by now
	q7_t *Im_out);

// modified from CMSIS-NN test_ref                            
void local_maxpool_q7_HWC(const q7_t * Im_in, 				// input image
	const uint16_t dim_im_in_x,   	// input image dimension x or W
	const uint16_t dim_im_in_y,   	// input image dimension y or H
	const uint16_t ch_im_in,    	// number of input image channels
	const uint16_t dim_kernel_x,  	// window kernel size
	const uint16_t dim_kernel_y,  	// window kernel size
	const uint16_t padding_x, 		// padding sizes
	const uint16_t padding_y, 		// padding sizes
	const uint16_t stride_x,  		// stride
	const uint16_t stride_y,  		// stride
	const uint16_t dim_im_out_x,  	// output image dimension x or W
	const uint16_t dim_im_out_y,  	// output image dimension y or H
	q7_t * bufferA, 				// a buffer for local storage, NULL by now
	q7_t * Im_out);

void local_maxpool_q7_CHW(const q7_t * Im_in, 				// input image
	const uint16_t dim_im_in_x,   	// input image dimension x or W
	const uint16_t dim_im_in_y,   	// input image dimension y or H
	const uint16_t ch_im_in,    	// number of input image channels
	const uint16_t dim_kernel_x,  	// window kernel size
	const uint16_t dim_kernel_y,  	// window kernel size
	const uint16_t padding_x, 		// padding sizes
	const uint16_t padding_y, 		// padding sizes
	const uint16_t stride_x,  		// stride
	const uint16_t stride_y,  		// stride
	const uint16_t dim_im_out_x,  	// output image dimension x or W
	const uint16_t dim_im_out_y,  	// output image dimension y or H
	q7_t * bufferA, 				// a buffer for local storage, NULL by now
	q7_t * Im_out);
							
void local_sumpool_q7_HWC(const q7_t * Im_in, // input image
	const uint16_t dim_im_in_x,   	// input image dimension x or W
	const uint16_t dim_im_in_y,   	// input image dimension y or H
	const uint16_t ch_im_in,    	// number of input image channels
	const uint16_t dim_kernel_x,  	// window kernel size
	const uint16_t dim_kernel_y,  	// window kernel size
	const uint16_t padding_x, 		// padding sizes
	const uint16_t padding_y, 		// padding sizes
	const uint16_t stride_x,  		// stride
	const uint16_t stride_y,  		// stride
	const uint16_t dim_im_out_x,  	// output image dimension x or W
	const uint16_t dim_im_out_y,  	// output image dimension y or H
	q7_t * bufferA, 				// a buffer for local storage, size = 4*output_size
	q7_t * Im_out);
							
void local_sumpool_q7_CHW(const q7_t * Im_in, // input image
	const uint16_t dim_im_in_x,   	// input image dimension x or W
	const uint16_t dim_im_in_y,   	// input image dimension y or H
	const uint16_t ch_im_in,    	// number of input image channels
	const uint16_t dim_kernel_x,  	// window kernel size
	const uint16_t dim_kernel_y,  	// window kernel size
	const uint16_t padding_x, 		// padding sizes
	const uint16_t padding_y, 		// padding sizes
	const uint16_t stride_x,  		// stride
	const uint16_t stride_y,  		// stride
	const uint16_t dim_im_out_x,  	// output image dimension x or W
	const uint16_t dim_im_out_y,  	// output image dimension y or H
	q7_t * bufferA, 				// a buffer for local storage, size = 4*output_size
	q7_t * Im_out);

// customised up sample pooling
void local_up_sampling_q7_HWC(const q7_t *Im_in,       // input image
	const uint16_t dim_im_in_x,  // input image dimension x or W
	const uint16_t dim_im_in_y,  // input image dimension y or H
	const uint16_t ch_im_in,     // number of input image channels
	const uint16_t dim_kernel_x, // window kernel size
	const uint16_t dim_kernel_y, // window kernel size
	const uint16_t dim_im_out_x, // output image dimension x or W
	const uint16_t dim_im_out_y, // output image dimension y or H
	q7_t *bufferA,               // NULL
	q7_t *Im_out);
						  
void local_up_sampling_q7_CHW(const q7_t *Im_in,       // input image
	const uint16_t dim_im_in_x,  // input image dimension x or W
	const uint16_t dim_im_in_y,  // input image dimension y or H
	const uint16_t ch_im_in,     // number of input image channels
	const uint16_t dim_kernel_x, // window kernel size
	const uint16_t dim_kernel_y, // window kernel size
	const uint16_t dim_im_out_x, // output image dimension x or W
	const uint16_t dim_im_out_y, // output image dimension y or H
	q7_t *bufferA,               // NULL
	q7_t *Im_out);

void local_convolve_HWC_q7_nonsquare(const q7_t *Im_in,                // input image
	const uint16_t dim_im_in_x,     // input image dimention x
	const uint16_t dim_im_in_y,     // input image dimention y
	const uint16_t ch_im_in,        // number of input image channels
	const q7_t *wt,                 // kernel weights
	const uint16_t ch_im_out,       // number of filters, i.e., output image channels
	const uint16_t dim_kernel_x,    // filter kernel size x
	const uint16_t dim_kernel_y,    // filter kernel size y
	const uint16_t padding_x,       // padding sizes x
	const uint16_t padding_y,       // padding sizes y
	const uint16_t stride_x,        // stride x
	const uint16_t stride_y,        // stride y
    const uint16_t dilation_x,      // dilation x
	const uint16_t dilation_y,      // dilation y
	const q7_t *bias,               // bias
	const nnom_qformat_param_t *bias_shift,                            // bias shifts
    const nnom_qformat_param_t *out_shift,                             // output shift
    const nnom_qtype_t q_type,      // per channel or per tensor
    q7_t *Im_out,                   // output image
	const uint16_t dim_im_out_x,    // output image dimension x
	const uint16_t dim_im_out_y,    // output image dimension y
	q15_t *bufferA,                 //buffer space for input
	q7_t *bufferB                   //buffer space for output
);
					   
void local_convolve_CHW_q7_nonsquare(const q7_t *Im_in,                // input image
	const uint16_t dim_im_in_x,     // input image dimention x
	const uint16_t dim_im_in_y,     // input image dimention y
	const uint16_t ch_im_in,        // number of input image channels
	const q7_t *wt,                 // kernel weights
	const uint16_t ch_im_out,       // number of filters, i.e., output image channels
	const uint16_t dim_kernel_x,    // filter kernel size x
	const uint16_t dim_kernel_y,    // filter kernel size y
	const uint16_t padding_x,       // padding sizes x
	const uint16_t padding_y,       // padding sizes y
	const uint16_t stride_x,        // stride x
	const uint16_t stride_y,        // stride y
    const uint16_t dilation_x,      // dilation x
	const uint16_t dilation_y,      // dilation y
	const q7_t *bias,               // bias
	const nnom_qformat_param_t *bias_shift,                            // bias shifts
    const nnom_qformat_param_t *out_shift,                             // output shift
    const nnom_qtype_t q_type,      // per channel or per tensor
    q7_t *Im_out,                   // output image
	const uint16_t dim_im_out_x,    // output image dimension x
	const uint16_t dim_im_out_y,    // output image dimension y
	q15_t *bufferA,                 //buffer space for input
	q7_t *bufferB                   //buffer space for output
);

void local_conv_trans_HWC_q7_nonsquare(const int8_t * Im_in,
	const uint16_t dim_im_in_x,     // input image dimention x
	const uint16_t dim_im_in_y,     // input image dimention y
	const uint16_t ch_im_in,        // number of input image channels
	const q7_t *wt,                 // kernel weights
	const uint16_t ch_im_out,       // number of filters, i.e., output image channels
	const uint16_t dim_kernel_x,    // filter kernel size x
	const uint16_t dim_kernel_y,    // filter kernel size y
	const uint16_t padding_x,       // padding sizes x
	const uint16_t padding_y,       // padding sizes y
	const uint16_t stride_x,        // stride x
	const uint16_t stride_y,        // stride y
    const uint16_t dilation_x,      // dilation x
	const uint16_t dilation_y,      // dilation y
	const q7_t *bias,               // bias
	const uint16_t bias_shift, const uint16_t out_shift, q7_t *Im_out, // output image
	const uint16_t dim_im_out_x,    // output image dimension x
	const uint16_t dim_im_out_y,    // output image dimension y
	q15_t *bufferA,                 //buffer space for input
	q7_t *bufferB                   //buffer space for output
);

void local_depthwise_separable_conv_HWC_q7_nonsquare(const q7_t *Im_in,// input image
	const uint16_t dim_im_in_x,     // input image dimention x
	const uint16_t dim_im_in_y,     // input image dimention y
	const uint16_t ch_im_in,        // number of input image channels
	const q7_t *wt,                 // kernel weights
	const uint16_t ch_im_out,       // number of filters, i.e., output image channels
	const uint16_t dim_kernel_x,    // filter kernel size x
	const uint16_t dim_kernel_y,    // filter kernel size y
	const uint16_t padding_x,       // padding sizes x
	const uint16_t padding_y,       // padding sizes y
	const uint16_t stride_x,        // stride x
	const uint16_t stride_y,        // stride y
    const uint16_t dilation_x,      // dilation x
	const uint16_t dilation_y,      // dilation y
	const q7_t *bias,               // bias
	const nnom_qformat_param_t *bias_shift,                            // bias shifts
    const nnom_qformat_param_t *out_shift,                             // output shift
    const nnom_qtype_t q_type,      // per channel or per tensor
    q7_t *Im_out,                   // output image
	const uint16_t dim_im_out_x,    // output image dimension x
	const uint16_t dim_im_out_y,    // output image dimension y
	q15_t *bufferA,                 //buffer space for input
	q7_t *bufferB                   //buffer space for output
);

void local_depthwise_separable_conv_CHW_q7_nonsquare(const q7_t *Im_in,// input image
	const uint16_t dim_im_in_x,     // input image dimention x
	const uint16_t dim_im_in_y,     // input image dimention y
	const uint16_t ch_im_in,        // number of input image channels
	const q7_t *wt,                 // kernel weights
	const uint16_t ch_im_out,       // number of filters, i.e., output image channels
	const uint16_t dim_kernel_x,    // filter kernel size x
	const uint16_t dim_kernel_y,    // filter kernel size y
	const uint16_t padding_x,       // padding sizes x
	const uint16_t padding_y,       // padding sizes y
	const uint16_t stride_x,        // stride x
	const uint16_t stride_y,        // stride y
    const uint16_t dilation_x,      // dilation x
	const uint16_t dilation_y,      // dilation y
	const q7_t *bias,               // bias
	const nnom_qformat_param_t *bias_shift,                            // bias shifts
    const nnom_qformat_param_t *out_shift,                             // output shift
    const nnom_qtype_t q_type,      // per channel or per tensor
    q7_t *Im_out,                   // output image
	const uint16_t dim_im_out_x,    // output image dimension x
	const uint16_t dim_im_out_y,    // output image dimension y
	q15_t *bufferA,                 //buffer space for input
	q7_t *bufferB                   //buffer space for output
);

void local_zero_padding_HWC_q7(const q7_t *Im_in,           // input image
	const uint16_t dim_im_in_x,    // input image dimention x
	const uint16_t dim_im_in_y,    // input image dimention y
	const uint16_t ch_im_in,       // number of input image channels
	const uint16_t padding_top,    // padding sizes y
	const uint16_t padding_bottom, // padding sizes y
	const uint16_t padding_left,   // padding sizes x
	const uint16_t padding_right,  // padding sizes x
	q7_t *Im_out,                  // output image
	const uint16_t dim_im_out_x,   // output image dimension x
	const uint16_t dim_im_out_y);  // output image dimension y 
						 
void local_zero_padding_CHW_q7(const q7_t *Im_in,           // input image
	const uint16_t dim_im_in_x,    // input image dimention x
	const uint16_t dim_im_in_y,    // input image dimention y
	const uint16_t ch_im_in,       // number of input image channels
	const uint16_t padding_top,    // padding sizes y
	const uint16_t padding_bottom, // padding sizes y
	const uint16_t padding_left,   // padding sizes x
	const uint16_t padding_right,  // padding sizes x
	q7_t *Im_out,                  // output image
	const uint16_t dim_im_out_x,   // output image dimension x
	const uint16_t dim_im_out_y);  // output image dimension y 
						 
void local_cropping_HWC_q7(const q7_t *Im_in,           // input image
	const uint16_t dim_im_in_x,    // input image dimention x
	const uint16_t dim_im_in_y,    // input image dimention y
	const uint16_t ch_im_in,       // number of input image channels
	const uint16_t padding_top,    // padding sizes y
	const uint16_t padding_bottom, // padding sizes y
	const uint16_t padding_left,   // padding sizes x
	const uint16_t padding_right,  // padding sizes x
	q7_t *Im_out,                  // output image
	const uint16_t dim_im_out_x,   // output image dimension x
	const uint16_t dim_im_out_y);  // output image dimension y 
						 
void local_cropping_CHW_q7(const q7_t *Im_in,           // input image
	const uint16_t dim_im_in_x,    // input image dimention x
	const uint16_t dim_im_in_y,    // input image dimention y
	const uint16_t ch_im_in,       // number of input image channels
	const uint16_t padding_top,    // padding sizes y
	const uint16_t padding_bottom, // padding sizes y
	const uint16_t padding_left,   // padding sizes x
	const uint16_t padding_right,  // padding sizes x
	q7_t *Im_out,                  // output image
	const uint16_t dim_im_out_x,   // output image dimension x
	const uint16_t dim_im_out_y);  // output image dimension y 

void local_fully_connected_q7_opt(const q7_t * pV,    // pointer to vector
	const q7_t * pM,    // pointer to matrix
	const uint16_t dim_vec, // length of the vector
	const uint16_t num_of_rows, // numCol of A
	const uint16_t bias_shift,  // amount of left-shift for bias
	const uint16_t out_shift,   // amount of right-shift for output
	const q7_t * bias, q7_t * pOut, // output operand
	q15_t * vec_buffer);


void local_fully_connected_q7(const q7_t * pV,    // pointer to vector
	const q7_t * pM,    		// pointer to matrix
	const uint16_t dim_vec, 	// length of the vector
	const uint16_t num_of_rows, // numCol of A
	const uint16_t bias_shift,  // amount of left-shift for bias
	const uint16_t out_shift,   // amount of right-shift for output
	const q7_t * bias, q7_t * pOut, // output operand
	q15_t * vec_buffer);

// matrix dot, 
// it takes reorderd weight as input, (see dense layer for detail. this is basiclly a dense opt without bias)
void local_dot_q7_opt(const q7_t *pV, // pointer to vector
	const q7_t *pM,               // pointer to matrix
	const uint16_t dim_vec,       // length of the vector
	const uint16_t num_of_rows,   // numCol of A
	const uint16_t out_shift,   // amount of right-shift for output
	q7_t *pOut);				// result buffer  

void local_dot_q7(const q7_t *pV, // pointer to vector
	const q7_t *pM,               // pointer to matrix
	const uint16_t dim_vec,       // length of the vector
	const uint16_t num_of_rows,   // numCol of A
	const uint16_t out_shift,     // amount of right-shift for output
	 q7_t *pOut);                   // output operand)



// softmax
void local_softmax_q7(const q7_t * vec_in, const uint32_t dim_vec, q7_t * p_out);

// sigmoid
void local_sigmoid_q7(q7_t * data, uint32_t size, int16_t int_width);

// tanh
void local_tanh_q7(q7_t * data, uint32_t size, int16_t int_width);

// relu
void local_relu_q7(q7_t * data, uint32_t size);

// leaky relu
void local_leaky_relu_q7(q7_t *data, q7_t alpha, uint32_t size);

// alpha in q7 format with dec_bit=7
// max and threshold has the same Q format with the activation
void local_adv_relu_q7(q7_t *data, q7_t alpha, q7_t max, q7_t threshold, uint32_t size);

// hard sigmoid, 
// y=-1 if x < -2.5
// y=1  if x > 2.5
// otherwise y = 0.2 * x + 0.5 (y=0.20315 * x + 0.5)
void local_hard_sigmoid_q7(q7_t *data, uint32_t size, int16_t dec_bit);

// hard tanh
// y=-1 if x < -1
// y=1  if x > 1
// otherwise y = x
void local_hard_tanh_q7(q7_t *data, uint32_t size, int16_t dec_bit);

// matrix ops
void local_mult_q7(q7_t * pSrcA, q7_t * pSrcB, q7_t * pDst, const uint16_t out_shift, uint32_t blockSize);

// add 
void local_add_q7(q7_t * pSrcA, q7_t * pSrcB, q7_t * pDst, const uint16_t out_shift,  uint32_t blockSize);

// sub 
void local_sub_q7(q7_t * pSrcA, q7_t * pSrcB, q7_t * pDst, const uint16_t out_shift, uint32_t blockSize);

// take multiple blocks (>2) as input
void local_multiple_add_q7( q7_t *p_dst,
                  const int16_t out_shift,
                  uint32_t block_size,
                  uint32_t num_block,
                  q7_t **p_src);

void local_multiple_mult_q7( q7_t *p_dst,
                  const int16_t out_shift,
                  uint32_t block_size,
                  uint32_t num_block,
                  q7_t **p_src);

void local_multiple_sub_q7( q7_t *p_dst,
                  const int16_t out_shift,
                  uint32_t block_size,
                  uint32_t num_block,
                  q7_t **p_src);


// Below tables credit to CMSIS
// For more info. check CMSIS-NN lib
// https://github.com/ARM-software/CMSIS_5/blob/develop/CMSIS/NN/Source/NNSupportFunctions/arm_nntables.c
static const q7_t nnom_sigmoid_table_q7[256] = {
    0x40, 0x42, 0x44, 0x46, 0x48, 0x4a, 0x4c, 0x4e,
    0x50, 0x52, 0x53, 0x55, 0x57, 0x59, 0x5a, 0x5c,
    0x5e, 0x5f, 0x61, 0x62, 0x63, 0x65, 0x66, 0x67,
    0x69, 0x6a, 0x6b, 0x6c, 0x6d, 0x6e, 0x6f, 0x70,
    0x71, 0x72, 0x72, 0x73, 0x74, 0x74, 0x75, 0x76,
    0x76, 0x77, 0x77, 0x78, 0x78, 0x79, 0x79, 0x7a,
    0x7a, 0x7a, 0x7b, 0x7b, 0x7b, 0x7c, 0x7c, 0x7c,
    0x7c, 0x7c, 0x7d, 0x7d, 0x7d, 0x7d, 0x7d, 0x7e,
    0x7e, 0x7e, 0x7e, 0x7e, 0x7e, 0x7e, 0x7e, 0x7f,
    0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f,
    0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f,
    0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f,
    0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f,
    0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f,
    0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f,
    0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
    0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
    0x01, 0x01, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02,
    0x02, 0x02, 0x03, 0x03, 0x03, 0x03, 0x03, 0x04,
    0x04, 0x04, 0x04, 0x04, 0x05, 0x05, 0x05, 0x06,
    0x06, 0x06, 0x07, 0x07, 0x08, 0x08, 0x09, 0x09,
    0x0a, 0x0a, 0x0b, 0x0c, 0x0c, 0x0d, 0x0e, 0x0e,
    0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16,
    0x17, 0x19, 0x1a, 0x1b, 0x1d, 0x1e, 0x1f, 0x21,
    0x22, 0x24, 0x26, 0x27, 0x29, 0x2b, 0x2d, 0x2e,
    0x30, 0x32, 0x34, 0x36, 0x38, 0x3a, 0x3c, 0x3e,
};


static const q7_t nnom_tanh_table_q7[256] = {
    (q7_t)0x00, (q7_t)0x08, (q7_t)0x10, (q7_t)0x18, (q7_t)0x1f, (q7_t)0x27, (q7_t)0x2e, (q7_t)0x35,
    (q7_t)0x3b, (q7_t)0x41, (q7_t)0x47, (q7_t)0x4c, (q7_t)0x51, (q7_t)0x56, (q7_t)0x5a, (q7_t)0x5e,
    (q7_t)0x61, (q7_t)0x65, (q7_t)0x68, (q7_t)0x6a, (q7_t)0x6d, (q7_t)0x6f, (q7_t)0x71, (q7_t)0x72,
    (q7_t)0x74, (q7_t)0x75, (q7_t)0x76, (q7_t)0x78, (q7_t)0x78, (q7_t)0x79, (q7_t)0x7a, (q7_t)0x7b,
    (q7_t)0x7b, (q7_t)0x7c, (q7_t)0x7c, (q7_t)0x7d, (q7_t)0x7d, (q7_t)0x7e, (q7_t)0x7e, (q7_t)0x7e,
    (q7_t)0x7e, (q7_t)0x7e, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f,
    (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f,
    (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f,
    (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f,
    (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f,
    (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f,
    (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f,
    (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f,
    (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f,
    (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f,
    (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f, (q7_t)0x7f,
    (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80,
    (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80,
    (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80,
    (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80,
    (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80,
    (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80,
    (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80,
    (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80,
    (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80,
    (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x80, (q7_t)0x81,
    (q7_t)0x81, (q7_t)0x81, (q7_t)0x81, (q7_t)0x81, (q7_t)0x81, (q7_t)0x81, (q7_t)0x81, (q7_t)0x82,
    (q7_t)0x82, (q7_t)0x82, (q7_t)0x82, (q7_t)0x82, (q7_t)0x83, (q7_t)0x83, (q7_t)0x84, (q7_t)0x84,
    (q7_t)0x85, (q7_t)0x85, (q7_t)0x86, (q7_t)0x87, (q7_t)0x88, (q7_t)0x88, (q7_t)0x8a, (q7_t)0x8b,
    (q7_t)0x8c, (q7_t)0x8e, (q7_t)0x8f, (q7_t)0x91, (q7_t)0x93, (q7_t)0x96, (q7_t)0x98, (q7_t)0x9b,
    (q7_t)0x9f, (q7_t)0xa2, (q7_t)0xa6, (q7_t)0xaa, (q7_t)0xaf, (q7_t)0xb4, (q7_t)0xb9, (q7_t)0xbf,
    (q7_t)0xc5, (q7_t)0xcb, (q7_t)0xd2, (q7_t)0xd9, (q7_t)0xe1, (q7_t)0xe8, (q7_t)0xf0, (q7_t)0xf8,
};


// ------------ 16bit ops --------------------

void local_avepool_q15_HWC(const q15_t *Im_in,           // input image
	const uint16_t dim_im_in_x,  // input image dimension x or W
	const uint16_t dim_im_in_y,  // input image dimension y or H
	const uint16_t ch_im_in,     // number of input image channels
	const uint16_t dim_kernel_x, // window kernel size
	const uint16_t dim_kernel_y, // window kernel size
	const uint16_t padding_x,    // padding sizes
	const uint16_t padding_y,    // padding sizes
	const uint16_t stride_x,     // stride
	const uint16_t stride_y,     // stride
	const uint16_t dim_im_out_x, // output image dimension x or W
	const uint16_t dim_im_out_y, // output image dimension y or H
	const uint16_t output_shift, // output right shift
	q7_t *bufferA,               // a buffer for local storage, NULL by now
	q15_t *Im_out);

void local_avepool_q15_CHW(const q15_t *Im_in,           // input image
	const uint16_t dim_im_in_x,  // input image dimension x or W
	const uint16_t dim_im_in_y,  // input image dimension y or H
	const uint16_t ch_im_in,     // number of input image channels
	const uint16_t dim_kernel_x, // window kernel size
	const uint16_t dim_kernel_y, // window kernel size
	const uint16_t padding_x,    // padding sizes
	const uint16_t padding_y,    // padding sizes
	const uint16_t stride_x,     // stride
	const uint16_t stride_y,     // stride
	const uint16_t dim_im_out_x, // output image dimension x or W
	const uint16_t dim_im_out_y, // output image dimension y or H
	const uint16_t output_shift, // output right shift
	q7_t *bufferA,               // a buffer for local storage, NULL by now
	q15_t *Im_out);

void local_maxpool_q15_HWC(const q15_t *Im_in,           // input image
	const uint16_t dim_im_in_x,  // input image dimension x or W
	const uint16_t dim_im_in_y,  // input image dimension y or H
	const uint16_t ch_im_in,     // number of input image channels
	const uint16_t dim_kernel_x, // window kernel size
	const uint16_t dim_kernel_y, // window kernel size
	const uint16_t padding_x,    // padding sizes
	const uint16_t padding_y,    // padding sizes
	const uint16_t stride_x,     // stride
	const uint16_t stride_y,     // stride
	const uint16_t dim_im_out_x, // output image dimension x or W
	const uint16_t dim_im_out_y, // output image dimension y or H
	q7_t *bufferA,               // a buffer for local storage, NULL by now
	q15_t *Im_out);

void local_maxpool_q15_CHW(const q15_t *Im_in,           // input image
	const uint16_t dim_im_in_x,  // input image dimension x or W
	const uint16_t dim_im_in_y,  // input image dimension y or H
	const uint16_t ch_im_in,     // number of input image channels
	const uint16_t dim_kernel_x, // window kernel size
	const uint16_t dim_kernel_y, // window kernel size
	const uint16_t padding_x,    // padding sizes
	const uint16_t padding_y,    // padding sizes
	const uint16_t stride_x,     // stride
	const uint16_t stride_y,     // stride
	const uint16_t dim_im_out_x, // output image dimension x or W
	const uint16_t dim_im_out_y, // output image dimension y or H
	q7_t *bufferA,               // a buffer for local storage, NULL by now
	q15_t *Im_out);

void local_sumpool_q15_HWC(const q15_t *Im_in,           // input image
	const uint16_t dim_im_in_x,  // input image dimension x or W
	const uint16_t dim_im_in_y,  // input image dimension y or H
	const uint16_t ch_im_in,     // number of input image channels
	const uint16_t dim_kernel_x, // window kernel size
	const uint16_t dim_kernel_y, // window kernel size
	const uint16_t padding_x,    // padding sizes
	const uint16_t padding_y,    // padding sizes
	const uint16_t stride_x,     // stride
	const uint16_t stride_y,     // stride
	const uint16_t dim_im_out_x, // output image dimension x or W
	const uint16_t dim_im_out_y, // output image dimension y or H
    const uint16_t output_shift, // output right shift
	q7_t *bufferA,               // a buffer for local storage, size = 4*output_size
	q15_t *Im_out);

void local_sumpool_q15_CHW(const q15_t *Im_in,           // input image
	const uint16_t dim_im_in_x,  // input image dimension x or W
	const uint16_t dim_im_in_y,  // input image dimension y or H
	const uint16_t ch_im_in,     // number of input image channels
	const uint16_t dim_kernel_x, // window kernel size
	const uint16_t dim_kernel_y, // window kernel size
	const uint16_t padding_x,    // padding sizes
	const uint16_t padding_y,    // padding sizes
	const uint16_t stride_x,     // stride
	const uint16_t stride_y,     // stride
	const uint16_t dim_im_out_x, // output image dimension x or W
	const uint16_t dim_im_out_y, // output image dimension y or H
	const uint16_t output_shift, // output right shift
	q7_t *bufferA,               // a buffer for local storage, size = 4*output_size
	q15_t *Im_out);  

void local_up_sampling_q15_HWC(const q15_t *Im_in,       // input image
	const uint16_t dim_im_in_x,  // input image dimension x or W
	const uint16_t dim_im_in_y,  // input image dimension y or H
	const uint16_t ch_im_in,     // number of input image channels
	const uint16_t dim_kernel_x, // window kernel size
	const uint16_t dim_kernel_y, // window kernel size
	const uint16_t dim_im_out_x, // output image dimension x or W
	const uint16_t dim_im_out_y, // output image dimension y or H
	q7_t *bufferA,               // a buffer for local storage, NULL by now
	q15_t *Im_out);

 void local_up_sampling_q15_CHW(const q15_t *Im_in,       // input image
	const uint16_t dim_im_in_x,  // input image dimension x or W
	const uint16_t dim_im_in_y,  // input image dimension y or H
	const uint16_t ch_im_in,     // number of input image channels
	const uint16_t dim_kernel_x, // window kernel size
	const uint16_t dim_kernel_y, // window kernel size
	const uint16_t dim_im_out_x, // output image dimension x or W
	const uint16_t dim_im_out_y, // output image dimension y or H
	q7_t *bufferA,               // a buffer for local storage, NULL by now
	q15_t *Im_out);

void local_convolve_HWC_q15_nonsquare(const q15_t *Im_in,   // input image
	const uint16_t dim_im_in_x,   // input image dimention x
	const uint16_t dim_im_in_y,   // input image dimention y
	const uint16_t ch_im_in,      // number of input image channels
	const q7_t *wt,               // kernel weights
	const uint16_t ch_im_out,     // number of filters, i.e., output image channels
	const uint16_t dim_kernel_x,  // filter kernel size x
	const uint16_t dim_kernel_y,  // filter kernel size y
	const uint16_t padding_x,     // padding sizes x
	const uint16_t padding_y,     // padding sizes y
	const uint16_t stride_x,      // stride x
	const uint16_t stride_y,      // stride y
    const uint16_t dilation_x,    // dilation x
	const uint16_t dilation_y,    // dilation y
	const q7_t *bias,             // bias
	const nnom_qformat_param_t *bias_shift, // bias shifts
    const nnom_qformat_param_t *out_shift, // output shift
    const nnom_qtype_t q_type,      // per channel or per tensor
    q15_t *Im_out,                  // output image
	const uint16_t dim_im_out_x,    // output image dimension x
	const uint16_t dim_im_out_y,    // output image dimension y
	q15_t *bufferA,                 //buffer space for input
	q7_t *bufferB                   //buffer space for output
);
void local_convolve_CHW_q15_nonsquare(const q15_t *Im_in,                // input image
	const uint16_t dim_im_in_x,     // input image dimention x
	const uint16_t dim_im_in_y,     // input image dimention y
	const uint16_t ch_im_in,        // number of input image channels
	const q7_t *wt,                 // kernel weights
	const uint16_t ch_im_out,       // number of filters, i.e., output image channels
	const uint16_t dim_kernel_x,    // filter kernel size x
	const uint16_t dim_kernel_y,    // filter kernel size y
	const uint16_t padding_x,       // padding sizes x
	const uint16_t padding_y,       // padding sizes y
	const uint16_t stride_x,        // stride x
	const uint16_t stride_y,        // stride y
    const uint16_t dilation_x,      // dilation x
	const uint16_t dilation_y,      // dilation y
	const q7_t *bias,               // bias
	const nnom_qformat_param_t *bias_shift,     // bias shifts
    const nnom_qformat_param_t *out_shift,      // output shift
    const nnom_qtype_t q_type,      // per channel or per tensor
    q15_t *Im_out,                   // output image
	const uint16_t dim_im_out_x,    // output image dimension x
	const uint16_t dim_im_out_y,    // output image dimension y
	q15_t *bufferA,                 //buffer space for input
	q7_t *bufferB                   //buffer space for output
);

void local_conv_trans_HWC_q15_nonsquare(const int8_t * Im_in,
	const uint16_t dim_im_in_x,     // input image dimention x
	const uint16_t dim_im_in_y,     // input image dimention y
	const uint16_t ch_im_in,        // number of input image channels
	const q7_t *wt,                 // kernel weights
	const uint16_t ch_im_out,       // number of filters, i.e., output image channels
	const uint16_t dim_kernel_x,    // filter kernel size x
	const uint16_t dim_kernel_y,    // filter kernel size y
	const uint16_t padding_x,       // padding sizes x
	const uint16_t padding_y,       // padding sizes y
	const uint16_t stride_x,        // stride x
	const uint16_t stride_y,        // stride y
    const uint16_t dilation_x,      // dilation x
	const uint16_t dilation_y,      // dilation y
	const q7_t *bias,               // bias
	const uint16_t bias_shift, const uint16_t out_shift, q15_t *Im_out, // output image
	const uint16_t dim_im_out_x,    // output image dimension x
	const uint16_t dim_im_out_y,    // output image dimension y
	q15_t *bufferA,                 //buffer space for input
	q7_t *bufferB                   //buffer space for output
);

void local_depthwise_separable_conv_HWC_q15_nonsquare(const q15_t *Im_in,// input image
	const uint16_t dim_im_in_x,     // input image dimention x
	const uint16_t dim_im_in_y,     // input image dimention y
	const uint16_t ch_im_in,        // number of input image channels
	const q7_t *wt,                 // kernel weights
	const uint16_t ch_im_out,       // number of filters, i.e., output image channels
	const uint16_t dim_kernel_x,    // filter kernel size x
	const uint16_t dim_kernel_y,    // filter kernel size y
	const uint16_t padding_x,       // padding sizes x
	const uint16_t padding_y,       // padding sizes y
	const uint16_t stride_x,        // stride x
	const uint16_t stride_y,        // stride y
    const uint16_t dilation_x,      // dilation x
	const uint16_t dilation_y,      // dilation y
	const q7_t *bias,               // bias
	const nnom_qformat_param_t *bias_shift,                            // bias shifts
    const nnom_qformat_param_t *out_shift,                             // output shift
    const nnom_qtype_t q_type,      // per channel or per tensor
    q15_t *Im_out,                   // output image
	const uint16_t dim_im_out_x,    // output image dimension x
	const uint16_t dim_im_out_y,    // output image dimension y
	q15_t *bufferA,                 //buffer space for input
	q7_t *bufferB                   //buffer space for output
);

void local_depthwise_separable_conv_CHW_q15_nonsquare(const q15_t *Im_in,// input image
	const uint16_t dim_im_in_x,     // input image dimention x
	const uint16_t dim_im_in_y,     // input image dimention y
	const uint16_t ch_im_in,        // number of input image channels
	const q7_t *wt,                 // kernel weights
	const uint16_t ch_im_out,       // number of filters, i.e., output image channels
	const uint16_t dim_kernel_x,    // filter kernel size x
	const uint16_t dim_kernel_y,    // filter kernel size y
	const uint16_t padding_x,       // padding sizes x
	const uint16_t padding_y,       // padding sizes y
	const uint16_t stride_x,        // stride x
	const uint16_t stride_y,        // stride y
    const uint16_t dilation_x,      // dilation x
	const uint16_t dilation_y,      // dilation y
	const q7_t *bias,               // bias
	const nnom_qformat_param_t *bias_shift,                            // bias shifts
    const nnom_qformat_param_t *out_shift,                             // output shift
    const nnom_qtype_t q_type,      // per channel or per tensor
    q15_t *Im_out,                   // output image
	const uint16_t dim_im_out_x,    // output image dimension x
	const uint16_t dim_im_out_y,    // output image dimension y
	q15_t *bufferA,                 //buffer space for input
	q7_t *bufferB                   //buffer space for output
);

void local_zero_padding_HWC_q15(const q15_t *Im_in,           // input image
	const uint16_t dim_im_in_x,    // input image dimention x
	const uint16_t dim_im_in_y,    // input image dimention y
	const uint16_t ch_im_in,       // number of input image channels
	const uint16_t padding_top,    // padding sizes y
	const uint16_t padding_bottom, // padding sizes y
	const uint16_t padding_left,   // padding sizes x
	const uint16_t padding_right,  // padding sizes x
	q15_t *Im_out,                  // output image
	const uint16_t dim_im_out_x,   // output image dimension x
	const uint16_t dim_im_out_y);   // output image dimension y 

void local_zero_padding_CHW_q15(const q15_t *Im_in,           // input image
	const uint16_t dim_im_in_x,    // input image dimention x
	const uint16_t dim_im_in_y,    // input image dimention y
	const uint16_t ch_im_in,       // number of input image channels
	const uint16_t padding_top,    // padding sizes y
	const uint16_t padding_bottom, // padding sizes y
	const uint16_t padding_left,   // padding sizes x
	const uint16_t padding_right,  // padding sizes x
	q15_t *Im_out,                  // output image
	const uint16_t dim_im_out_x,   // output image dimension x
	const uint16_t dim_im_out_y);   // output image dimension y 

void local_cropping_HWC_q15(const q15_t *Im_in,           // input image
	const uint16_t dim_im_in_x,    // input image dimention x
	const uint16_t dim_im_in_y,    // input image dimention y
	const uint16_t ch_im_in,       // number of input image channels
	const uint16_t padding_top,    // padding sizes y
	const uint16_t padding_bottom, // padding sizes y
	const uint16_t padding_left,   // padding sizes x
	const uint16_t padding_right,  // padding sizes x
	q15_t *Im_out,                  // output image
	const uint16_t dim_im_out_x,   // output image dimension x
	const uint16_t dim_im_out_y);   // output image dimension y 

void local_cropping_CHW_q15(const q15_t *Im_in,           // input image
	const uint16_t dim_im_in_x,    // input image dimention x
	const uint16_t dim_im_in_y,    // input image dimention y
	const uint16_t ch_im_in,       // number of input image channels
	const uint16_t padding_top,    // padding sizes y
	const uint16_t padding_bottom, // padding sizes y
	const uint16_t padding_left,   // padding sizes x
	const uint16_t padding_right,  // padding sizes x
	q15_t *Im_out,                  // output image
	const uint16_t dim_im_out_x,   // output image dimension x
	const uint16_t dim_im_out_y);   // output image dimension y 


void local_dot_q15(const q15_t *pV, // pointer to vector
	const q15_t *pM,               // pointer to matrix
	const uint16_t dim_vec,       // length of the vector
	const uint16_t num_of_rows,   // numCol of A
	const uint16_t out_shift,     // amount of right-shift for output
	 q15_t *pOut);                   // output operand)

void local_dot_q15_opt(const q15_t * pV,
	const q15_t * pM,
	const uint16_t dim_vec,
	const uint16_t num_of_rows,
	const uint16_t out_shift, 
	q15_t * pOut);

// original implementation
// this support none bias, the it will perform like a dot. 
// set the `bias=NULL` to work
void local_fully_connected_mat_q7_vec_q15(const q15_t * pV, // pointer to vector
	const q7_t * pM,			// pointer to matrix
	const uint16_t dim_vec,  	// length of the vector
	const uint16_t num_of_rows, // numCol of A
	const uint16_t bias_shift, 	// amount of left-shift for bias
	const uint16_t out_shift, 	// amount of right-shift for output
	const q7_t * bias,			// bias
	q15_t * pOut,				// output
	q15_t * vec_buffer);     	// not used but to keep the interface same as the ARM's version

// work on recorder matrix
// this support none bias, set the bias=NULL to work								   
void local_fully_connected_mat_q7_vec_q15_opt(const q15_t * pV,
	const q7_t * pM,
	const uint16_t dim_vec,
	const uint16_t num_of_rows,
	const uint16_t bias_shift,
	const uint16_t out_shift, 
	const q7_t * bias, 
	q15_t * pOut, 
	q15_t * vec_buffer);

// matrix operation Q15
void local_multiple_add_q15( q15_t *p_dst,
                  const int16_t out_shift,
                  uint32_t block_size,
                  uint32_t num_block,
                  q15_t **p_src);
                  
void local_multiple_mult_q15( q15_t *p_dst,
                  const int16_t out_shift,
                  uint32_t block_size,
                  uint32_t num_block,
                  q15_t **p_src);

void local_multiple_sub_q15( q15_t *p_dst,
                  const int16_t out_shift,
                  uint32_t block_size,
                  uint32_t num_block,
                  q15_t **p_src);

void local_mult_q15(q15_t * pSrcA, q15_t * pSrcB, q15_t * pDst, const uint16_t out_shift, uint32_t blockSize);

// add 
void local_add_q15(q15_t * pSrcA, q15_t * pSrcB, q15_t * pDst, const uint16_t out_shift,  uint32_t blockSize);

// sub 
void local_sub_q15(q15_t * pSrcA, q15_t * pSrcB, q15_t * pDst, const uint16_t out_shift, uint32_t blockSize);

// Convert Q7 to Q15
void local_q7_to_q15_no_shift(const q7_t *src, q15_t *des, uint32_t size);
void local_q7_to_q15(const q7_t *src, q15_t *des, uint32_t size);

// q15 shift to q7
void local_q15_to_q7(const q15_t *src, q7_t *des,  uint32_t shift, uint32_t size);

// y = 1 - x
void local_1_minor_z_q15(q15_t *src, q15_t *des, uint16_t dec_bit, uint32_t size);

void local_softmax_q15(const q15_t * vec_in, const uint16_t dim_vec, q15_t * p_out);
void local_hard_sigmoid_q15(q15_t *data, uint32_t size, int16_t dec_bit);
void local_hard_tanh_q15(q15_t *data, uint32_t size, int16_t dec_bit);
void local_relu_q15(q15_t *data, uint32_t size);
void local_leaky_relu_q15(q15_t *data, q7_t alpha, uint32_t size);
void local_adv_relu_q15(q15_t *data, q7_t negative_slope, q15_t max, q15_t threshold, uint32_t size);
void local_sigmoid_q15(q15_t * data, uint32_t size, uint16_t int_width);
void local_tanh_q15(q15_t * data, uint32_t size, uint16_t int_width);


static const q15_t nnom_sigmoid_table_q15[256] = {
    0x4000, 0x4200, 0x43ff, 0x45fc, 0x47f5, 0x49eb, 0x4bdc, 0x4dc8,
    0x4fad, 0x518a, 0x5360, 0x552c, 0x56ef, 0x58a8, 0x5a57, 0x5bfb,
    0x5d93, 0x5f20, 0x60a1, 0x6216, 0x637f, 0x64db, 0x662b, 0x676f,
    0x68a6, 0x69d2, 0x6af1, 0x6c05, 0x6d0d, 0x6e09, 0x6efb, 0x6fe2,
    0x70be, 0x7190, 0x7258, 0x7316, 0x73cc, 0x7478, 0x751b, 0x75b7,
    0x764a, 0x76d6, 0x775b, 0x77d8, 0x784f, 0x78c0, 0x792a, 0x798f,
    0x79ee, 0x7a48, 0x7a9d, 0x7aed, 0x7b39, 0x7b80, 0x7bc4, 0x7c03,
    0x7c3f, 0x7c78, 0x7cad, 0x7ce0, 0x7d0f, 0x7d3c, 0x7d66, 0x7d8d,
    0x7db3, 0x7dd6, 0x7df7, 0x7e16, 0x7e33, 0x7e4f, 0x7e69, 0x7e81,
    0x7e98, 0x7eae, 0x7ec2, 0x7ed5, 0x7ee7, 0x7ef8, 0x7f08, 0x7f17,
    0x7f25, 0x7f32, 0x7f3e, 0x7f4a, 0x7f55, 0x7f5f, 0x7f69, 0x7f72,
    0x7f7b, 0x7f83, 0x7f8a, 0x7f91, 0x7f98, 0x7f9e, 0x7fa4, 0x7faa,
    0x7faf, 0x7fb4, 0x7fb8, 0x7fbd, 0x7fc1, 0x7fc5, 0x7fc8, 0x7fcc,
    0x7fcf, 0x7fd2, 0x7fd5, 0x7fd7, 0x7fda, 0x7fdc, 0x7fde, 0x7fe0,
    0x7fe2, 0x7fe4, 0x7fe6, 0x7fe7, 0x7fe9, 0x7fea, 0x7feb, 0x7fed,
    0x7fee, 0x7fef, 0x7ff0, 0x7ff1, 0x7ff2, 0x7ff3, 0x7ff4, 0x7ff4,
    0x000b, 0x000c, 0x000c, 0x000d, 0x000e, 0x000f, 0x0010, 0x0011,
    0x0012, 0x0013, 0x0015, 0x0016, 0x0017, 0x0019, 0x001a, 0x001c,
    0x001e, 0x0020, 0x0022, 0x0024, 0x0026, 0x0029, 0x002b, 0x002e,
    0x0031, 0x0034, 0x0038, 0x003b, 0x003f, 0x0043, 0x0048, 0x004c,
    0x0051, 0x0056, 0x005c, 0x0062, 0x0068, 0x006f, 0x0076, 0x007d,
    0x0085, 0x008e, 0x0097, 0x00a1, 0x00ab, 0x00b6, 0x00c2, 0x00ce,
    0x00db, 0x00e9, 0x00f8, 0x0108, 0x0119, 0x012b, 0x013e, 0x0152,
    0x0168, 0x017f, 0x0197, 0x01b1, 0x01cd, 0x01ea, 0x0209, 0x022a,
    0x024d, 0x0273, 0x029a, 0x02c4, 0x02f1, 0x0320, 0x0353, 0x0388,
    0x03c1, 0x03fd, 0x043c, 0x0480, 0x04c7, 0x0513, 0x0563, 0x05b8,
    0x0612, 0x0671, 0x06d6, 0x0740, 0x07b1, 0x0828, 0x08a5, 0x092a,
    0x09b6, 0x0a49, 0x0ae5, 0x0b88, 0x0c34, 0x0cea, 0x0da8, 0x0e70,
    0x0f42, 0x101e, 0x1105, 0x11f7, 0x12f3, 0x13fb, 0x150f, 0x162e,
    0x175a, 0x1891, 0x19d5, 0x1b25, 0x1c81, 0x1dea, 0x1f5f, 0x20e0,
    0x226d, 0x2405, 0x25a9, 0x2758, 0x2911, 0x2ad4, 0x2ca0, 0x2e76,
    0x3053, 0x3238, 0x3424, 0x3615, 0x380b, 0x3a04, 0x3c01, 0x3e00,
};


static const q15_t nnom_tanh_table_q15[256] = {
    (q15_t)0x0000, (q15_t)0x07fd, (q15_t)0x0feb, (q15_t)0x17b9, (q15_t)0x1f59, (q15_t)0x26bf, (q15_t)0x2ddf, (q15_t)0x34ae,
    (q15_t)0x3b27, (q15_t)0x4142, (q15_t)0x46fd, (q15_t)0x4c56, (q15_t)0x514d, (q15_t)0x55e2, (q15_t)0x5a1a, (q15_t)0x5df6,
    (q15_t)0x617c, (q15_t)0x64b0, (q15_t)0x6797, (q15_t)0x6a37, (q15_t)0x6c95, (q15_t)0x6eb5, (q15_t)0x709e, (q15_t)0x7254,
    (q15_t)0x73dc, (q15_t)0x753a, (q15_t)0x7672, (q15_t)0x7788, (q15_t)0x787f, (q15_t)0x795b, (q15_t)0x7a1e, (q15_t)0x7acb,
    (q15_t)0x7b65, (q15_t)0x7bee, (q15_t)0x7c66, (q15_t)0x7cd1, (q15_t)0x7d30, (q15_t)0x7d84, (q15_t)0x7dce, (q15_t)0x7e0f,
    (q15_t)0x7e49, (q15_t)0x7e7d, (q15_t)0x7eaa, (q15_t)0x7ed2, (q15_t)0x7ef5, (q15_t)0x7f14, (q15_t)0x7f30, (q15_t)0x7f48,
    (q15_t)0x7f5e, (q15_t)0x7f71, (q15_t)0x7f82, (q15_t)0x7f91, (q15_t)0x7f9e, (q15_t)0x7fa9, (q15_t)0x7fb3, (q15_t)0x7fbc,
    (q15_t)0x7fc4, (q15_t)0x7fcb, (q15_t)0x7fd1, (q15_t)0x7fd7, (q15_t)0x7fdc, (q15_t)0x7fe0, (q15_t)0x7fe4, (q15_t)0x7fe7,
    (q15_t)0x7fea, (q15_t)0x7fed, (q15_t)0x7fef, (q15_t)0x7ff1, (q15_t)0x7ff3, (q15_t)0x7ff4, (q15_t)0x7ff6, (q15_t)0x7ff7,
    (q15_t)0x7ff8, (q15_t)0x7ff9, (q15_t)0x7ffa, (q15_t)0x7ffa, (q15_t)0x7ffb, (q15_t)0x7ffc, (q15_t)0x7ffc, (q15_t)0x7ffd,
    (q15_t)0x7ffd, (q15_t)0x7ffd, (q15_t)0x7ffe, (q15_t)0x7ffe, (q15_t)0x7ffe, (q15_t)0x7ffe, (q15_t)0x7fff, (q15_t)0x7fff,
    (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff,
    (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff,
    (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff,
    (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff,
    (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff, (q15_t)0x7fff,
    (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000,
    (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000,
    (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000,
    (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8000,
    (q15_t)0x8000, (q15_t)0x8000, (q15_t)0x8001, (q15_t)0x8001, (q15_t)0x8001, (q15_t)0x8001, (q15_t)0x8001, (q15_t)0x8001,
    (q15_t)0x8001, (q15_t)0x8001, (q15_t)0x8001, (q15_t)0x8002, (q15_t)0x8002, (q15_t)0x8002, (q15_t)0x8002, (q15_t)0x8003,
    (q15_t)0x8003, (q15_t)0x8003, (q15_t)0x8004, (q15_t)0x8004, (q15_t)0x8005, (q15_t)0x8006, (q15_t)0x8006, (q15_t)0x8007,
    (q15_t)0x8008, (q15_t)0x8009, (q15_t)0x800a, (q15_t)0x800c, (q15_t)0x800d, (q15_t)0x800f, (q15_t)0x8011, (q15_t)0x8013,
    (q15_t)0x8016, (q15_t)0x8019, (q15_t)0x801c, (q15_t)0x8020, (q15_t)0x8024, (q15_t)0x8029, (q15_t)0x802f, (q15_t)0x8035,
    (q15_t)0x803c, (q15_t)0x8044, (q15_t)0x804d, (q15_t)0x8057, (q15_t)0x8062, (q15_t)0x806f, (q15_t)0x807e, (q15_t)0x808f,
    (q15_t)0x80a2, (q15_t)0x80b8, (q15_t)0x80d0, (q15_t)0x80ec, (q15_t)0x810b, (q15_t)0x812e, (q15_t)0x8156, (q15_t)0x8183,
    (q15_t)0x81b7, (q15_t)0x81f1, (q15_t)0x8232, (q15_t)0x827c, (q15_t)0x82d0, (q15_t)0x832f, (q15_t)0x839a, (q15_t)0x8412,
    (q15_t)0x849b, (q15_t)0x8535, (q15_t)0x85e2, (q15_t)0x86a5, (q15_t)0x8781, (q15_t)0x8878, (q15_t)0x898e, (q15_t)0x8ac6,
    (q15_t)0x8c24, (q15_t)0x8dac, (q15_t)0x8f62, (q15_t)0x914b, (q15_t)0x936b, (q15_t)0x95c9, (q15_t)0x9869, (q15_t)0x9b50,
    (q15_t)0x9e84, (q15_t)0xa20a, (q15_t)0xa5e6, (q15_t)0xaa1e, (q15_t)0xaeb3, (q15_t)0xb3aa, (q15_t)0xb903, (q15_t)0xbebe,
    (q15_t)0xc4d9, (q15_t)0xcb52, (q15_t)0xd221, (q15_t)0xd941, (q15_t)0xe0a7, (q15_t)0xe847, (q15_t)0xf015, (q15_t)0xf803,
};

#ifdef __cplusplus
}
#endif

#endif  /* __NNOM_LOCAL_H__ */
