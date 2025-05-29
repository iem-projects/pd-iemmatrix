/*
Uniformly Partitioned, Time-Variant, 
Multichannel-Input-Mulichannel-Output Block Convolution
(and because signal processing folks like incomprehensible
 abbreviations: UPTVMIMOBC, yeah!)

useful for all kinds of real-time processes
virtual acoustic reality,
6dof spatial audio rendering/auralization,
convolution reverb,
Ambisonic decoding to headphones (binaural),
Ambisonic encoding of compact spherical microphone arrays
Ambisonic decoding to compact spherical loudspeaker arrays
Wave-Field Synthesis processing,
dynamic Binaural rendering,
Beamforming with loudspeaker or micrphone arrays,
adaptive signal processing,
frequency response equalization,
...

Franz Zotter
Hannes Pescoller
Sourena Mosleh

Email-address:
zotter@iem.at
h.pescoller@kug.ac.at
sourena.mosleh@student.kug.ac.at

Institute of Electronic Music and Acoustics (IEM)
University of Music and Performing Arts Graz
2024, 2025
*/

#include "../include/convolver.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define TRUE 1
#define FALSE 0

/* crossfade functions */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void crossfade(float *y, float *y_new, float *w, int len) {
  for (int n = 0; n < len; n++)
    y[n] = y[n] * w[n] + y_new[n] * (1 - w[n]);
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void setNewIR(conv_data *conv, _Bool status) {
  conv->convolver_switch = status;
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
_Bool getNewIR(conv_data *conv) {
  _Bool result;
  result = conv->convolver_switch;
  return result;
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
/* PARTITIONED CONVOLUTION CORE */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void conv_process(conv_data *conv, float **in, float **out) {

  for (int in_ch = 0; in_ch < conv->num_inputs; in_ch++) {
    copyArray(conv->x_old[in_ch], conv->xtemp, conv->L); // save old block
    copyArray(in[in_ch], &conv->xtemp[conv->L],
              conv->L);                 // combine with new block
    fftwf_execute(conv->fftplan_xtemp); // perform FFT
    copyComplexArray(conv->xftemp, conv->xf[in_ch][conv->current_rb],
                     conv->L + 1); // write into ring buffer
    copyArray(in[in_ch], conv->x_old[in_ch],
              conv->L); // store previous for next call
  }

  for (int out_ch = 0; out_ch < conv->num_outputs; out_ch++) {
    resetComplexArray(conv->yftemp,
                      conv->L + 1); // zero output of partition accumulator
    for (int in_ch = 0; in_ch < conv->num_inputs; in_ch++) {
      for (int p = 0, px = conv->current_rb; p < conv->P; p++, px++) {
        freq_mul_acc(conv->xf[in_ch][px % conv->P],
                     conv->hf[conv->current_cf][out_ch][in_ch][p], conv->yftemp,
                     conv->L + 1);
      }
    }

    fftwf_execute(conv->ifftplan_y); // perform accumulated partitions iFFT
    if (getNewIR(conv) == TRUE) {
      int current_cf = (conv->current_cf + 1) % NUM_CF;
      resetComplexArray(conv->yftemp, conv->L + 1);
      for (int in_ch = 0; in_ch < conv->num_inputs; in_ch++) {
        for (int p = 0, px = conv->current_rb; p < conv->P; p++) {
          freq_mul_acc(conv->xf[in_ch][px % conv->P],
                       conv->hf[current_cf][out_ch][in_ch][p], conv->yftemp,
                       conv->L + 1);
          px++;
        }
      }
      fftwf_execute(conv->ifftplan_y_cf);
      crossfade(conv->y + conv->L, conv->y_cf + conv->L, conv->w, conv->L);
    }
    copyArrayWithGain(&conv->y[conv->L], out[out_ch], conv->L,
                      2 * conv->L); // second half is result
  }
  if (getNewIR(conv) == TRUE) {
    conv->current_cf = (conv->current_cf + 1) % NUM_CF;
    setNewIR(conv, FALSE); // update status
  }
  conv->current_rb = (conv->current_rb + conv->P - 1) %
                     conv->P; // decrease ring buffer position
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
conv_data *initConvolution(int L, int P, int Hann_len, int num_inputs, int num_outputs) {
  conv_data *conv = (conv_data *)malloc(sizeof(conv_data));
  conv->L = L;
  conv->P = P;
  conv->num_inputs = num_inputs;
  conv->num_outputs = num_outputs;
  conv->xtemp = new1DArray(conv->L * 2);
  conv->htemp = new1DArray(conv->L * 2);
  conv->y = new1DArray(conv->L * 2);
  conv->y_cf = new1DArray(conv->L * 2);
  conv->w = new1DArray(conv->L);
  conv->current_rb = 0;
  conv->current_cf = 1;
  conv->Hann_len = Hann_len;
  resetArray(conv->xtemp, conv->L * 2);
  resetArray(conv->htemp, conv->L * 2);
  resetArray(&conv->htemp[conv->L], conv->L); // zero pad
  resetArray(conv->w, conv->L);
  hann(conv->w, conv->Hann_len);
  conv->convolver_switch = 0;

  conv->x_old = new2DArray(conv->num_inputs, conv->L);
  reset2DArray(conv->x_old, conv->num_inputs, conv->L);

  conv->xf = new3DComplexArray(conv->num_inputs, conv->P, conv->L + 1);
  reset3DComplexArray(conv->xf, conv->num_inputs, conv->P, conv->L + 1);
  conv->hf = new5DComplexArray(NUM_CF, conv->num_outputs, conv->num_inputs,
                               conv->P, conv->L + 1);
  reset5DComplexArray(conv->hf, NUM_CF, conv->num_outputs, conv->num_inputs,
                      conv->P, conv->L + 1);

  conv->xftemp = new1DComplexArray(conv->L + 1);
  conv->hftemp = new1DComplexArray(conv->L + 1);
  conv->yftemp = new1DComplexArray(conv->L + 1);

  conv->fftplan_xtemp = fftwf_plan_dft_r2c_1d(conv->L * 2, conv->xtemp,
                                              conv->xftemp, FFTW_ESTIMATE);
  conv->fftplan_htemp = fftwf_plan_dft_r2c_1d(conv->L * 2, conv->htemp,
                                              conv->hftemp, FFTW_ESTIMATE);
  conv->ifftplan_y =
      fftwf_plan_dft_c2r_1d(conv->L * 2, conv->yftemp, conv->y, FFTW_ESTIMATE);
  conv->ifftplan_y_cf = fftwf_plan_dft_c2r_1d(conv->L * 2, conv->yftemp,
                                              conv->y_cf, FFTW_ESTIMATE);
  return conv;
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void freeConvolution(conv_data *conv) {
  fftwf_destroy_plan(conv->fftplan_xtemp);
  fftwf_destroy_plan(conv->fftplan_htemp);
  fftwf_destroy_plan(conv->ifftplan_y);
  fftwf_destroy_plan(conv->ifftplan_y_cf);

  free3DComplexArray(conv->xf, conv->num_inputs, conv->P);
  free5DComplexArray(conv->hf, NUM_CF, conv->num_outputs, conv->num_inputs,
                     conv->P);
  free2DArray(conv->x_old, conv->num_inputs);

  fftwf_free(conv->yftemp);
  fftwf_free(conv->xtemp);
  fftwf_free(conv->y);
  fftwf_free(conv->w);
  fftwf_free(conv->y_cf);
  fftwf_free(conv->htemp);
  fftwf_free(conv->hftemp);
  fftwf_free(conv->xftemp);
  free(conv);
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void setImpulseResponse(conv_data *conv, float ***inh) {
  setNewIR(conv, TRUE);
  conv->convolver_switch = 1;

  for (int out_ch = 0; out_ch < conv->num_outputs; out_ch++) {
    for (int in_ch = 0; in_ch < conv->num_inputs; in_ch++) {
      for (int partition = 0; partition < conv->P; partition++) {
        copyArray(&inh[out_ch][in_ch][partition * conv->L], conv->htemp,
                  conv->L);
        fftwf_execute(conv->fftplan_htemp);
        copyComplexArray(
            conv->hftemp,
            conv->hf[(conv->current_cf + 1) % NUM_CF][out_ch][in_ch][partition],
            conv->L + 1);
      }
    }
  }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void setImpulseResponseZeroPad(conv_data *conv, float ***inh, int num_samples) {
  int offset;
  int copy_length;
  setNewIR(conv, TRUE);
  conv->convolver_switch = 1;
  for (int out_ch = 0; out_ch < conv->num_outputs; out_ch++) {
    for (int in_ch = 0; in_ch < conv->num_inputs; in_ch++) {
      offset = 0;
      for (int partition = 0; partition < conv->P; partition++) {
        copy_length = (offset+conv->L > num_samples) ? num_samples-offset : conv->L;
        copyArray(&inh[out_ch][in_ch][offset], conv->htemp,
                  copy_length);
        resetArray(&conv->htemp[copy_length], conv->L - copy_length);
        fftwf_execute(conv->fftplan_htemp);
        copyComplexArray(
            conv->hftemp,
            conv->hf[(conv->current_cf + 1) % NUM_CF][out_ch][in_ch][partition],
            conv->L + 1);
        offset+=conv->L;
      }
    }
  }
}