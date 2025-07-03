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

For information on usage and redistribution, and for a DISCLAIMER OF ALL
WARRANTIES, see the file, "LICENSE.txt," in this distribution.

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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif /* HAVE_CONFIG_H */

#include "convolver.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#if HAVE_FFTWF
//# define USE_FFTWF 1
#else
# if defined(__GNUC__)
#  warning "Building without FFTW3"
# endif
#endif


#ifdef USE_FFTWF
# if defined(__GNUC__)
#  warning USE_FFTWF
# endif
static const t_fftwf_malloc my_malloc = fftwf_malloc;
static const t_fftwf_free my_free = fftwf_free;
static const t_fftwf_execute my_execute = fftwf_execute;
static const t_fftwf_destroy_plan my_destroy_plan = fftwf_destroy_plan;
static const t_fftwf_plan_dft_r2c_1d my_plan_dft_r2c_1d;
static const t_fftwf_plan_dft_c2r_1d my_plan_dft_c2r_1d;
static const int have_fftwf = 1;
#else
static t_fftwf_malloc my_malloc = malloc;
static t_fftwf_free my_free = free;
static t_fftwf_execute my_execute = 0;
static t_fftwf_destroy_plan my_destroy_plan = 0;
static t_fftwf_plan_dft_r2c_1d my_plan_dft_r2c_1d = 0;
static t_fftwf_plan_dft_c2r_1d my_plan_dft_c2r_1d = 0;
static int have_fftwf = 0;
#endif
#include "m_pd.h"
int IEMCONVOLVE(convolver_set_fftwf_functions) (const t_fftwf_functions*funs) {
#ifndef USE_FFTWF
# if defined(__GNUC__)
#  warning overwrite fftwf functions
#endif
  if(funs && funs->malloc)
    my_malloc = funs->malloc;
  if(funs && funs->free)
    my_free = funs->free;
  if(funs && funs->execute)
    my_execute = funs->execute;
  if(funs && funs->destroy_plan)
    my_destroy_plan = funs->destroy_plan;
  if(funs && funs->plan_dft_r2c_1d)
    my_plan_dft_r2c_1d = funs->plan_dft_r2c_1d;
  if(funs && funs->plan_dft_c2r_1d)
    my_plan_dft_c2r_1d = funs->plan_dft_c2r_1d;
  if(1
     && my_malloc
     && my_free
     && my_execute
     && my_destroy_plan
     && my_plan_dft_r2c_1d
     && my_plan_dft_c2r_1d
     )
    have_fftwf = 1;
#endif
  return have_fftwf;
}



/* crossfade functions */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(crossFade) (float *y, float *y_new, float *w_old, float *w_new, unsigned int len) {
  for (unsigned int n = 0; n < len; n++)
    y[n] = y[n] * w_old[n] + y_new[n] * w_new[n];
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(registerCrossFade) (conv_data *conv) {
  conv->register_crossfade = 1;
}
void IEMCONVOLVE(registerCrossFadeComplete) (conv_data *conv) {
  conv->register_crossfade = 0;
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
_Bool IEMCONVOLVE(wasCrossFadeRegistered)(conv_data *conv) {
  return conv->register_crossfade;
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
/* PARTITIONED CONVOLUTION CORE */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(convProcess) (conv_data *conv, float **in, float **out) {
  if(have_fftwf) {
  for (unsigned int in_ch = 0; in_ch < conv->num_inputs; in_ch++) {
    IEMCONVOLVE(copyArray) (conv->x_old[in_ch], conv->xtemp, conv->blocksize); // copy old signal block
    IEMCONVOLVE(copyArray) (in[in_ch], &conv->xtemp[conv->blocksize],conv->blocksize);// append new signal block
    my_execute(conv->fftplan_xtemp); // perform 2L FFT to xftemp
    IEMCONVOLVE(copyComplexArray) (conv->xftemp, conv->xf[in_ch][conv->current_rb],
                     conv->blocksize + 1); // write FFT into ring buffer
    IEMCONVOLVE(copyArray) (in[in_ch], conv->x_old[in_ch],conv->blocksize); // store current signal block
  }
  for (unsigned int out_ch = 0; out_ch < conv->num_outputs; out_ch++) {
    IEMCONVOLVE(resetComplexArray) (conv->yftemp,conv->blocksize + 1); // zero output's partition accumulator
    for (unsigned int in_ch = 0; in_ch < conv->num_inputs; in_ch++) { // MIMO convolutions (main)
      for (unsigned int p = 0, px = conv->current_rb; p < conv->num_partitions; p++, px++) {
        IEMCONVOLVE(freq_mul_acc) (conv->xf[in_ch][px % conv->num_partitions],
                     conv->hf[conv->current_cf][out_ch][in_ch][p], conv->yftemp,
                     conv->blocksize + 1);
      }
    }
    my_execute(conv->ifftplan_y); // perform iFFT of the main-IR output partition
    /* IF IR WAS UPDATED: ALSO COMPUTE NEW OUTPUT FOR CROSSFADE AT THE OUT CHANNEL */
    if (IEMCONVOLVE(wasCrossFadeRegistered(conv))) {
      unsigned int current_cf = (conv->current_cf + 1) % NUM_CF;
      IEMCONVOLVE(resetComplexArray) (conv->yftemp, conv->blocksize + 1); // zero output's partitions accumulator
      for (unsigned int in_ch = 0; in_ch < conv->num_inputs; in_ch++) { // MIMO convolutions (update)
        for (unsigned int p = 0, px = conv->current_rb; p < conv->num_partitions; p++) {
          IEMCONVOLVE(freq_mul_acc) (conv->xf[in_ch][px % conv->num_partitions],
                       conv->hf[current_cf][out_ch][in_ch][p], conv->yftemp,
                       conv->blocksize + 1);
          px++;
        }
      }
      my_execute(conv->ifftplan_y_cf); // perform iFFT of the updated-IR output partition
      IEMCONVOLVE(crossFade) (conv->y + conv->blocksize, conv->y_cf + conv->blocksize, conv->w_old, conv->w_new, conv->blocksize); // xfade main->updated
    }
    /* IF IR WAS UPDATED: CROSSFADE TO NEW OUTPUT OF THE OUT CHANNEL COMPLETE */
    IEMCONVOLVE(copyArrayWithGain) (&conv->y[conv->blocksize], out[out_ch], conv->blocksize,
                      2 * conv->blocksize); // second half times N is the output's resulting signal block
  }
  if (IEMCONVOLVE(wasCrossFadeRegistered(conv))) {
    conv->current_cf = (conv->current_cf + 1) % NUM_CF;
    IEMCONVOLVE(registerCrossFadeComplete(conv)); // update status
  }
  conv->current_rb = (conv->current_rb + conv->num_partitions - 1) %
                     conv->num_partitions; // decrease ring for IR partitions
  }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
conv_data *IEMCONVOLVE(initConvolution) (unsigned int blocksize, unsigned int num_partitions, unsigned int xfade_length, unsigned int num_inputs, unsigned int num_outputs, _Bool coherent_xfade) {
  if (have_fftwf) {
  conv_data *conv = (conv_data *)malloc(sizeof(conv_data));
  conv->blocksize = blocksize; // block length (2L=FFT length)
  conv->num_partitions = num_partitions; // number of partitions
  conv->num_inputs = num_inputs;
  conv->num_outputs = num_outputs;
  conv->xtemp = IEMCONVOLVE(new1DArray) (conv->blocksize * 2); // input signal FFT buffer
  conv->htemp = IEMCONVOLVE(new1DArray) (conv->blocksize * 2); // IR signal partition FFT buffer
  conv->y = IEMCONVOLVE(new1DArray) (conv->blocksize * 2);     // output signal FFT buffer main/old
  conv->y_cf = IEMCONVOLVE(new1DArray) (conv->blocksize * 2);  // output signal FFT buffer new for crossfade
  conv->w_old = IEMCONVOLVE(new1DArray) (conv->blocksize);      // fade-out window signal buffer
  conv->w_new = IEMCONVOLVE(new1DArray) (conv->blocksize);      // fade-in window signal buffer
  conv->current_rb = 0;         // partition ringbuffer index [0;P-1]
  conv->current_cf = 0;         // crossfade-from index {0, 1}
  conv->xfade_length = (xfade_length<blocksize)? xfade_length : blocksize;    // crossfade length (<=L) for Hann window
  IEMCONVOLVE(resetArray) (conv->xtemp, conv->blocksize * 2);
  IEMCONVOLVE(resetArray) (conv->htemp, conv->blocksize * 2);
  IEMCONVOLVE(resetArray) (&conv->htemp[conv->blocksize], conv->blocksize); // zero pad
  IEMCONVOLVE(resetArray) (conv->w_old, conv->blocksize);
  IEMCONVOLVE(resetArray) (conv->w_new, conv->blocksize);
  if (coherent_xfade) { // amplitude-complementary crossfade for coherent signals
    IEMCONVOLVE(cos2win)(conv->w_old, conv->xfade_length);
    IEMCONVOLVE(sin2win)(conv->w_new, conv->xfade_length);
  } else { // power complementary crossfade for incoherent signals
    IEMCONVOLVE(coswin)(conv->w_old, conv->xfade_length);
    IEMCONVOLVE(sinwin)(conv->w_new, conv->xfade_length);
  }
  conv->register_crossfade = 0;

  conv->x_old = IEMCONVOLVE(new2DArray) (conv->num_inputs, conv->blocksize); // old input signal buffer
  IEMCONVOLVE(reset2DArray) (conv->x_old, conv->num_inputs, conv->blocksize);
  // frequency-domain buffers holding all input-signal partitions
  conv->xf = IEMCONVOLVE(new3DComplexArray) (conv->num_inputs, conv->num_partitions, conv->blocksize + 1);
  IEMCONVOLVE(reset3DComplexArray) (conv->xf, conv->num_inputs, conv->num_partitions, conv->blocksize + 1);
  // frequency-domain buffers holding all MIMO ir partitions
  conv->hf = IEMCONVOLVE(new5DComplexArray) (NUM_CF, conv->num_outputs, conv->num_inputs,
                               conv->num_partitions, conv->blocksize + 1);
  IEMCONVOLVE(reset5DComplexArray) (conv->hf, NUM_CF, conv->num_outputs, conv->num_inputs,
                      conv->num_partitions, conv->blocksize + 1);
  // frequency-domain single-channel temporary buffers
  conv->xftemp = IEMCONVOLVE(new1DComplexArray) (conv->blocksize + 1);
  conv->hftemp = IEMCONVOLVE(new1DComplexArray) (conv->blocksize + 1);
  conv->yftemp = IEMCONVOLVE(new1DComplexArray) (conv->blocksize + 1);
  // single-channel FFT plan for temporary signal and ir blocks
  conv->fftplan_xtemp = my_plan_dft_r2c_1d(conv->blocksize * 2, conv->xtemp,
                                              conv->xftemp, FFTW_ESTIMATE);
  conv->fftplan_htemp = my_plan_dft_r2c_1d(conv->blocksize * 2, conv->htemp,
                                              conv->hftemp, FFTW_ESTIMATE);
  // single-channel IFFT plan for temporary main/old and new output signal
  conv->ifftplan_y =
      my_plan_dft_c2r_1d(conv->blocksize * 2, conv->yftemp, conv->y, FFTW_ESTIMATE);
  conv->ifftplan_y_cf = my_plan_dft_c2r_1d(conv->blocksize * 2, conv->yftemp,
                                              conv->y_cf, FFTW_ESTIMATE);
  return conv;
  }
  return 0;
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(freeConvolution) (conv_data *conv) {
  if(have_fftwf) {
  // single-channel FFT and IFFT plans of temporary signals
  my_destroy_plan(conv->fftplan_xtemp);
  my_destroy_plan(conv->fftplan_htemp);
  my_destroy_plan(conv->ifftplan_y);
  my_destroy_plan(conv->ifftplan_y_cf);
  // frequency-domain array holding all input signal partitions
  IEMCONVOLVE(free3DComplexArray) (conv->xf, conv->num_inputs, conv->num_partitions);
  // frequency-domain array with all MIMO ir partitions
  IEMCONVOLVE(free5DComplexArray)(conv->hf, NUM_CF, conv->num_outputs, conv->num_inputs,
                     conv->num_partitions);
  // old input signal-block buffers
  IEMCONVOLVE(free2DArray) (conv->x_old, conv->num_inputs);
  // time-domain temporary FFT signal blocks
  }
  if(my_free) {
  my_free(conv->yftemp);
  my_free(conv->xtemp);
  my_free(conv->y);
  my_free(conv->w_old);
  my_free(conv->w_new);
  my_free(conv->y_cf);
  my_free(conv->htemp);
  // frequency-domain temporary FFT blocks
  my_free(conv->hftemp);
  my_free(conv->xftemp);
  }
  free(conv);
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(setImpulseResponseZeroPad) (conv_data *conv, float ***inh, unsigned int num_samples, _Bool no_xfade_init) {
  unsigned int offset;
  unsigned int copy_length;
  int hot_or_cold_stream;
  if(!have_fftwf)
    return;
  if (!conv)
    return;
  if (no_xfade_init) { // update hot stream directly, without crossfading, for initialization
    hot_or_cold_stream = conv->current_cf;
    IEMCONVOLVE(registerCrossFadeComplete) (conv);
  } else { // update cold stream and register a cross-fade from hot -> cold stream next dsp cycle
    hot_or_cold_stream = (conv->current_cf + 1) % NUM_CF;
    IEMCONVOLVE(registerCrossFade) (conv);
  }
  for (unsigned int out_ch = 0; out_ch < conv->num_outputs; out_ch++) {
    for (unsigned int in_ch = 0; in_ch < conv->num_inputs; in_ch++) {
      offset = 0;
      for (unsigned int partition = 0; partition < conv->num_partitions; partition++) {
        copy_length = (offset+conv->blocksize > num_samples) ? num_samples-offset : conv->blocksize;
        IEMCONVOLVE(copyArray) (&inh[out_ch][in_ch][offset], conv->htemp,
                  copy_length);
        IEMCONVOLVE(resetArray) (&conv->htemp[copy_length], conv->blocksize - copy_length);
        my_execute(conv->fftplan_htemp);
        IEMCONVOLVE(copyComplexArray) (
            conv->hftemp,
            conv->hf[hot_or_cold_stream][out_ch][in_ch][partition],
            conv->blocksize + 1);
        offset+=conv->blocksize;
      }
    }
  }
}
