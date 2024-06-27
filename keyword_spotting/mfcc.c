/*
 * Copyright (C) 2018 Arm Limited or its affiliates. All rights reserved.
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * Licensed under the Apache License, Version 2.0 (the License); you may
 * not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Modified by Jianjia Ma for C implementation.
 *
 */

/*
 * Description: MFCC feature extraction to match with TensorFlow MFCC Op
 */

#include "mfcc.h"

#include <stdlib.h>
#include <string.h>

#include "float.h"
#include "nnom_port.h"

#define M_PI 3.14159265358979323846264338327950288

// const int16_t
// cos_lookup[]={9999,9996,9993,9987,9981,9972,9963,9951,9939,9924,9909,9891,9873,9852,9831,9807,9783,9757,9729,9700,9669,9637,9604,9569,9533,9495,9456,9415,9373,9329,9285,9238,9191,9142,9091,9039,8986,8932,8876,8819,8760,8700,8639,8577,8513,8448,8382,8314,8245,8175,8104,8032,7958,7883,7807,7730,7651,7572,7491,7409,7326,7242,7157,7070,6983,6895,6806,6715,6624,6531,6438,6343,6248,6152,6055,5956,5857,5758,5657,5555,5453,5349,5245,5141,5035,4928,4821,4713,4605,4496,4386,4275,4164,4052,3939,3826,3713,3598,3484,3368,3253,3136,3020,2902,2785,2667,2548,2429,2310,2191,2071,1950,1830,1709,1588,1467,1345,1224,1102,9801,8579,7356,6132,4906,3680,2454,1227,0000,-1227,-2454,-3680,-4906,-6132,-7356,-8579,-9801,-1102,-1224,-1345,-1467,-1588,-1709,-1830,-1950,-2071,-2191,-2310,-2429,-2548,-2667,-2785,-2902,-3020,-3136,-3253,-3368,-3484,-3598,-3713,-3826,-3939,-4052,-4164,-4275,-4386,-4496,-4605,-4713,-4821,-4928,-5035,-5141,-5245,-5349,-5453,-5555,-5657,-5758,-5857,-5956,-6055,-6152,-6248,-6343,-6438,-6531,-6624,-6715,-6806,-6895,-6983,-7070,-7157,-7242,-7326,-7409,-7491,-7572,-7651,-7730,-7807,-7883,-7958,-8032,-8104,-8175,-8245,-8314,-8382,-8448,-8513,-8577,-8639,-8700,-8760,-8819,-8876,-8932,-8986,-9039,-9091,-9142,-9191,-9238,-9285,-9329,-9373,-9415,-9456,-9495,-9533,-9569,-9604,-9637,-9669,-9700,-9729,-9757,-9783,-9807,-9831,-9852,-9873,-9891,-9909,-9924,-9939,-9951,-9963,-9972,-9981,-9987,-9993,-9996,-9999};

const int16_t sin_lookup[] = {
    0,    122,  245,  368,  490,  613,  735,  857,  980,   1102, 1224, 1345,
    1467, 1588, 1709, 1830, 1950, 2071, 2191, 2310, 2429,  2548, 2667, 2785,
    2902, 3020, 3136, 3253, 3368, 3484, 3598, 3713, 3826,  3939, 4052, 4164,
    4275, 4386, 4496, 4605, 4713, 4821, 4928, 5035, 5141,  5245, 5349, 5453,
    5555, 5657, 5758, 5857, 5956, 6055, 6152, 6248, 6343,  6438, 6531, 6624,
    6715, 6806, 6895, 6983, 7071, 7157, 7242, 7326, 7409,  7491, 7572, 7651,
    7730, 7807, 7883, 7958, 8032, 8104, 8175, 8245, 8314,  8382, 8448, 8513,
    8577, 8639, 8700, 8760, 8819, 8876, 8932, 8986, 9039,  9091, 9142, 9191,
    9238, 9285, 9329, 9373, 9415, 9456, 9495, 9533, 9569,  9604, 9637, 9669,
    9700, 9729, 9757, 9783, 9807, 9831, 9852, 9873, 9891,  9909, 9924, 9939,
    9951, 9963, 9972, 9981, 9987, 9993, 9996, 9999, 10000, 9999, 9996, 9993,
    9987, 9981, 9972, 9963, 9951, 9939, 9924, 9909, 9891,  9873, 9852, 9831,
    9807, 9783, 9757, 9729, 9700, 9669, 9637, 9604, 9569,  9533, 9495, 9456,
    9415, 9373, 9329, 9285, 9238, 9191, 9142, 9091, 9039,  8986, 8932, 8876,
    8819, 8760, 8700, 8639, 8577, 8513, 8448, 8382, 8314,  8245, 8175, 8104,
    8032, 7958, 7883, 7807, 7730, 7651, 7572, 7491, 7409,  7326, 7242, 7157,
    7071, 6983, 6895, 6806, 6715, 6624, 6531, 6438, 6343,  6248, 6152, 6055,
    5956, 5857, 5758, 5657, 5555, 5453, 5349, 5245, 5141,  5035, 4928, 4821,
    4713, 4605, 4496, 4386, 4275, 4164, 4052, 3939, 3826,  3713, 3598, 3484,
    3368, 3253, 3136, 3020, 2902, 2785, 2667, 2548, 2429,  2310, 2191, 2071,
    1950, 1830, 1709, 1588, 1467, 1345, 1224, 1102, 980,   857,  735,  613,
    490,  368,  245,  122};
const int16_t cos_lookup[] = {
    10000, 9999,  9996,  9993,  9987,  9981,  9972,  9963,  9951,  9939,  9924,
    9909,  9891,  9873,  9852,  9831,  9807,  9783,  9757,  9729,  9700,  9669,
    9637,  9604,  9569,  9533,  9495,  9456,  9415,  9373,  9329,  9285,  9238,
    9191,  9142,  9091,  9039,  8986,  8932,  8876,  8819,  8760,  8700,  8639,
    8577,  8513,  8448,  8382,  8314,  8245,  8175,  8104,  8032,  7958,  7883,
    7807,  7730,  7651,  7572,  7491,  7409,  7326,  7242,  7157,  7071,  6983,
    6895,  6806,  6715,  6624,  6531,  6438,  6343,  6248,  6152,  6055,  5956,
    5857,  5758,  5657,  5555,  5453,  5349,  5245,  5141,  5035,  4928,  4821,
    4713,  4605,  4496,  4386,  4275,  4164,  4052,  3939,  3826,  3713,  3598,
    3484,  3368,  3253,  3136,  3020,  2902,  2785,  2667,  2548,  2429,  2310,
    2191,  2071,  1950,  1830,  1709,  1588,  1467,  1345,  1224,  1102,  980,
    857,   735,   613,   490,   368,   245,   122,   0,     -122,  -245,  -368,
    -490,  -613,  -735,  -857,  -980,  -1102, -1224, -1345, -1467, -1588, -1709,
    -1830, -1950, -2071, -2191, -2310, -2429, -2548, -2667, -2785, -2902, -3020,
    -3136, -3253, -3368, -3484, -3598, -3713, -3826, -3939, -4052, -4164, -4275,
    -4386, -4496, -4605, -4713, -4821, -4928, -5035, -5141, -5245, -5349, -5453,
    -5555, -5657, -5758, -5857, -5956, -6055, -6152, -6248, -6343, -6438, -6531,
    -6624, -6715, -6806, -6895, -6983, -7071, -7157, -7242, -7326, -7409, -7491,
    -7572, -7651, -7730, -7807, -7883, -7958, -8032, -8104, -8175, -8245, -8314,
    -8382, -8448, -8513, -8577, -8639, -8700, -8760, -8819, -8876, -8932, -8986,
    -9039, -9091, -9142, -9191, -9238, -9285, -9329, -9373, -9415, -9456, -9495,
    -9533, -9569, -9604, -9637, -9669, -9700, -9729, -9757, -9783, -9807, -9831,
    -9852, -9873, -9891, -9909, -9924, -9939, -9951, -9963, -9972, -9981, -9987,
    -9993, -9996, -9999};
#ifndef MFCC_PLATFORM_ARM
// FFT code from arduino_fft: https://github.com/lloydroc/arduino_fft
// change to float data£¬ modify to fit within this file
// see the above link for license( MIT license).
#include <math.h>
#include <stdio.h>

void rearrange(float data_re[], float data_im[], const unsigned int N) {
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
    while (target & (mask >>= 1))
      target &= ~mask;
    target |= mask;
  }
  return;
}

void compute(float data_re[], float data_im[], const unsigned int N) {
  // puts("in compute\n");
  const float pi = -3.14159265358979323846;
  for (unsigned int step = 1; step < N; step <<= 1) {
    const unsigned int jump = step << 1;
    const int step_d = step;
    const float step_d1=(float)step;
    float twiddle_re = 1.0;
    float twiddle_im = 0.0;
    for (unsigned int group = 0; group < step; group++) {
      for (unsigned int pair = group; pair < N; pair += jump) {
        const unsigned int match = pair + step;
        const float product_re =
            twiddle_re * data_re[match] - twiddle_im * data_im[match];
        const float product_im =
            twiddle_im * data_re[match] + twiddle_re * data_im[match];
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
 float angle1 = ((float)group + 1) /step_d1;
 angle1*=pi;
      uint16_t angle = ((group + 1) << 8) / step_d;

      twiddle_re = (float)cos_lookup[angle] / 10000;
    //   twiddle_im = (float)sin_lookup[angle] / -10000;



// twiddle_re=cosf(angle1);
  twiddle_im = sinf(angle1);

    }

  }
  return;
}



// printf("%f,", angle);
// printf("group: %d, step_d: %f,angle: %f",group,step_d,angle);

// twiddle_re = cosf((float)angle *pi/256);

// printf("a:%f, ",twiddle_re1);

// print("%f  ",twiddle_re);

// print("%f",twiddle_re*10000);
// printf(",cosf: %f",twiddle_re);
// printf("%f,", twiddle_re);
// printf("\n");

// printf(",sinf: %f\n",twiddle_im);

int fft(float data_re[], float data_im[], const int N) {
  rearrange(data_re, data_im, N);
  compute(data_re, data_im, N);
  // printf("stopping");
  // asm volatile ("wfi");
  // printf("exiting wfi");
  return 0;
}

#endif /* end of FFT implmentation*/

static void *mfcc_malloc(size_t size) {
  void *p = nnom_malloc(size);
  if (p == NULL)
    return NULL;
  memset(p, 0, size);
  return p;
}

static void mfcc_free(void *p) {
  if (p != NULL)
    nnom_free(p);
  return;
}

mfcc_t *mfcc_create(int num_mfcc_features, int feature_offset, int num_fbank,
                    int frame_len, float preempha, int is_append_energy) {
  mfcc_t *mfcc;
  mfcc = mfcc_malloc(sizeof(mfcc_t));

  mfcc->num_mfcc_features = num_mfcc_features;
  mfcc->num_features_offset = feature_offset;
  mfcc->num_fbank = num_fbank;
  mfcc->frame_len = frame_len;
  mfcc->preempha = preempha;
  mfcc->is_append_energy = is_append_energy;

  // Round-up to nearest power of 2.
  mfcc->frame_len_padded = (int)powf(2, ceilf((logf(frame_len) / logf(2))));

  mfcc->frame = mfcc_malloc(sizeof(float) * mfcc->frame_len_padded);
  mfcc->buffer = mfcc_malloc(sizeof(float) * mfcc->frame_len_padded);
  mfcc->mel_energies = mfcc_malloc(sizeof(float) * mfcc->num_fbank);

  // create window function, hanning
  mfcc->window_func = mfcc_malloc(sizeof(float) * frame_len);
  for (int i = 0; i < frame_len; i++) {
    // printf("%f", )
    // printf()
    // int16_t a=cos_lookup[i];

    int16_t ov = i & (1 << 8);
    if (!ov) {
      mfcc->window_func[i] = 0.5f - 0.5f * (float)cos_lookup[i] / 10000;

    } else {
      int index = index & ((1 << 8) - 1);
      mfcc->window_func[i] = 0.5f - 0.5f * (float)cos_lookup[index] / -10000;
    }

	//  mfcc->window_func[i] = 0.5f - 0.5f*cosf((float)M_2PI * ((float)i) / (frame_len));

    // printf(" cosf: %f\n",cosf((float)M_2PI * ((float)i) / (frame_len)));
  }

  // create mel filterbank
  mfcc->fbank_filter_first = mfcc_malloc(sizeof(int32_t) * mfcc->num_fbank);
  mfcc->fbank_filter_last = mfcc_malloc(sizeof(int32_t) * mfcc->num_fbank);
  mfcc->mel_fbank = create_mel_fbank(mfcc);

  // create DCT matrix
  mfcc->dct_matrix = create_dct_matrix(mfcc->num_fbank, num_mfcc_features);

#ifdef MFCC_PLATFORM_ARM
  // initialize FFT
  mfcc->rfft = mfcc_malloc(sizeof(arm_rfft_fast_instance_f32));
  arm_rfft_fast_init_f32(mfcc->rfft, mfcc->frame_len_padded);
#else
  mfcc->fft_buffer = mfcc_malloc(sizeof(float) * mfcc->frame_len_padded * 2);
#endif

  return mfcc;
}

int mfcc_delete(mfcc_t *mfcc) {
  mfcc_free(mfcc->frame);
  mfcc_free(mfcc->buffer);
  mfcc_free(mfcc->mel_energies);
  mfcc_free(mfcc->window_func);
  mfcc_free(mfcc->fbank_filter_first);
  mfcc_free(mfcc->fbank_filter_last);
  mfcc_free(mfcc->dct_matrix);
  for (int i = 0; i < mfcc->num_fbank; i++)
    mfcc_free(mfcc->mel_fbank[i]);
  mfcc_free(mfcc->mel_fbank);

#ifdef MFCC_PLATFORM_ARM
  mfcc_free(mfcc->rfft);
#else
  mfcc_free(mfcc->fft_buffer);
#endif
  mfcc_free(mfcc);

  return 0;
}

float *create_dct_matrix(int32_t input_length, int32_t coefficient_count) {
  int32_t k, n;
  float *M = mfcc_malloc(sizeof(float) * input_length * coefficient_count);
  float normalizer;
#ifdef MFCC_PLATFORM_ARM
  arm_sqrt_f32(2.0f / (float)input_length, &normalizer);
#else
  normalizer = sqrtf(2.0f / (float)input_length);
#endif
  for (k = 0; k < coefficient_count; k++) {
    for (n = 0; n < input_length; n++) {
      M[k * input_length + n] =
          normalizer * cosf(((float)M_PI) / input_length * (n + 0.5f) * k);
    }
  }

  return M;
}

float **create_mel_fbank(mfcc_t *mfcc) {
  // compute points evenly spaced in mels
  float mel_low_freq = MelScale(MEL_LOW_FREQ);
  float mel_high_freq = MelScale(MEL_HIGH_FREQ);
  float mel_freq_delta = (mel_high_freq - mel_low_freq) / (mfcc->num_fbank + 1);

  float *bin = mfcc_malloc(sizeof(float) * (mfcc->num_fbank + 2));
  for (int i = 0; i < mfcc->num_fbank + 2; i++) {
    bin[i] = mel_low_freq + mel_freq_delta * i;
    bin[i] = floor((mfcc->frame_len_padded + 1) * InverseMelScale(bin[i]) /
                   SAMP_FREQ);
  }

  float **mel_fbank = mfcc_malloc(sizeof(float *) * mfcc->num_fbank);

  for (int j = 0; j < mfcc->num_fbank; j++) {
    mel_fbank[j] =
        mfcc_malloc(sizeof(float) * (mfcc->frame_len_padded / 2 + 1));
    for (int i = (int)bin[j]; i < (int)bin[j + 1]; i++)
      mel_fbank[j][i] = (i - bin[j]) / (bin[j + 1] - bin[j]);
    for (int i = (int)bin[j + 1]; i < (int)bin[j + 2]; i++)
      mel_fbank[j][i] = (bin[j + 2] - i) / (bin[j + 2] - bin[j + 1]);
  }

  mfcc_free(bin);
  return mel_fbank;
}

void mfcc_compute(mfcc_t *mfcc, const int16_t *audio_data, float *mfcc_out) {
  int32_t i, j, bin;

  // 1. TensorFlow way of normalizing .wav data to (-1,1) and 2. do
  // pre-emphasis.
  float last = (float)audio_data[0];
  mfcc->frame[0] = last / (1 << 15);
  for (i = 1; i < mfcc->frame_len; i++) {
    mfcc->frame[i] = ((float)audio_data[i] - last * mfcc->preempha) / (1 << 15);
    last = (float)audio_data[i];
  }
  // Fill up remaining with zeros
  if (mfcc->frame_len_padded - mfcc->frame_len)
    memset(&mfcc->frame[mfcc->frame_len], 0,
           sizeof(float) * (mfcc->frame_len_padded - mfcc->frame_len));

  // windows filter
  for (i = 0; i < mfcc->frame_len; i++) {
    mfcc->frame[i] *= mfcc->window_func[i];
  }

#ifdef MFCC_PLATFORM_ARM  // ToDo add other fft implementation
  // Compute FFT
  arm_rfft_fast_f32(mfcc->rfft, mfcc->frame, mfcc->buffer, 0);

  // Convert to power spectrum
  // frame is stored as [real0, realN/2-1, real1, im1, real2, im2, ...]
  int32_t half_dim = mfcc->frame_len_padded / 2;
  float first_energy = mfcc->buffer[0] * mfcc->buffer[0];
  float last_energy =
      mfcc->buffer[1] * mfcc->buffer[1];  // handle this special case
  for (i = 1; i < half_dim; i++) {
    float real = mfcc->buffer[i * 2];
    float im = mfcc->buffer[i * 2 + 1];
    mfcc->buffer[i] = real * real + im * im;
  }
  mfcc->buffer[0] = first_energy;
  mfcc->buffer[half_dim] = last_energy;

#else  // end of ARM_fft
       // not yet optimized for memory
  float *data_re = mfcc->fft_buffer;
  float *data_im = &mfcc->fft_buffer[mfcc->frame_len_padded];

  memcpy(data_re, mfcc->frame, mfcc->frame_len_padded * sizeof(float));
  memset(data_im, 0, mfcc->frame_len_padded * sizeof(float));

  fft(data_re, data_im, mfcc->frame_len_padded);
  // only need half (N/2+1)
  for (int i = 0; i <= mfcc->frame_len_padded / 2; i++) {
    mfcc->buffer[i] = (data_re[i] * data_re[i] + data_im[i] * data_im[i]) /
                      mfcc->frame_len_padded;
  }
#endif

  float sqrt_data;
  // Apply mel filterbanks
  for (bin = 0; bin < mfcc->num_fbank; bin++) {
    float mel_energy = 0;
    for (i = 0; i < mfcc->frame_len_padded / 2 + 1; i++) {
      mel_energy += mfcc->buffer[i] * mfcc->mel_fbank[bin][i];
    }
    mfcc->mel_energies[bin] = mel_energy;

    // avoid log of zero
    if (mel_energy == 0.0f)
      mfcc->mel_energies[bin] = FLT_MIN;
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
  return;
}
