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
#else
# warning "Building without FFTW3"
#endif


/* crossfade functions */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(crossFade) (float *y, float *y_new, float *w_old, float *w_new, int len) {
  for (int n = 0; n < len; n++)
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
#if HAVE_FFTWF
  for (int in_ch = 0; in_ch < conv->num_inputs; in_ch++) {
    IEMCONVOLVE(copyArray) (conv->x_old[in_ch], conv->xtemp, conv->blocksize); // copy old signal block
    IEMCONVOLVE(copyArray) (in[in_ch], &conv->xtemp[conv->blocksize],conv->blocksize);// append new signal block
    fftwf_execute(conv->fftplan_xtemp); // perform 2L FFT to xftemp
    IEMCONVOLVE(copyComplexArray) (conv->xftemp, conv->xf[in_ch][conv->current_rb],
                     conv->blocksize + 1); // write FFT into ring buffer
    IEMCONVOLVE(copyArray) (in[in_ch], conv->x_old[in_ch],conv->blocksize); // store current signal block
  }
  for (int out_ch = 0; out_ch < conv->num_outputs; out_ch++) {
    IEMCONVOLVE(resetComplexArray) (conv->yftemp,conv->blocksize + 1); // zero output's partition accumulator
    for (int in_ch = 0; in_ch < conv->num_inputs; in_ch++) { // MIMO convolutions (main)
      for (int p = 0, px = conv->current_rb; p < conv->num_partitions; p++, px++) {
        IEMCONVOLVE(freq_mul_acc) (conv->xf[in_ch][px % conv->num_partitions],
                     conv->hf[conv->current_cf][out_ch][in_ch][p], conv->yftemp,
                     conv->blocksize + 1);
      }
    }
    fftwf_execute(conv->ifftplan_y); // perform iFFT of the main-IR output partition
    /* IF IR WAS UPDATED: ALSO COMPUTE NEW OUTPUT FOR CROSSFADE AT THE OUT CHANNEL */
    if (IEMCONVOLVE(wasCrossFadeRegistered(conv))) {
      int current_cf = (conv->current_cf + 1) % NUM_CF;
      IEMCONVOLVE(resetComplexArray) (conv->yftemp, conv->blocksize + 1); // zero output's partitions accumulator
      for (int in_ch = 0; in_ch < conv->num_inputs; in_ch++) { // MIMO convolutions (update)
        for (int p = 0, px = conv->current_rb; p < conv->num_partitions; p++) {
          IEMCONVOLVE(freq_mul_acc) (conv->xf[in_ch][px % conv->num_partitions],
                       conv->hf[current_cf][out_ch][in_ch][p], conv->yftemp,
                       conv->blocksize + 1);
          px++;
        }
      }
      fftwf_execute(conv->ifftplan_y_cf); // perform iFFT of the updated-IR output partition
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
#endif
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
conv_data *IEMCONVOLVE(initConvolution) (int blocksize, int num_partitions, int xfade_length, int num_inputs, int num_outputs, _Bool coherent_xfade) {
#if HAVE_FFTWF
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
  // frequency-domain buffers holding all input-signal paritions
  conv->xf = IEMCONVOLVE(new3DComplexArray) (conv->num_inputs, conv->num_partitions, conv->blocksize + 1);
  IEMCONVOLVE(reset3DComplexArray) (conv->xf, conv->num_inputs, conv->num_partitions, conv->blocksize + 1);
  // frequency-domain buffers holding all MIMO ir paritions
  conv->hf = IEMCONVOLVE(new5DComplexArray) (NUM_CF, conv->num_outputs, conv->num_inputs,
                               conv->num_partitions, conv->blocksize + 1);
  IEMCONVOLVE(reset5DComplexArray) (conv->hf, NUM_CF, conv->num_outputs, conv->num_inputs,
                      conv->num_partitions, conv->blocksize + 1);
  // frequency-domain single-channel temporary buffers
  conv->xftemp = IEMCONVOLVE(new1DComplexArray) (conv->blocksize + 1);
  conv->hftemp = IEMCONVOLVE(new1DComplexArray) (conv->blocksize + 1);
  conv->yftemp = IEMCONVOLVE(new1DComplexArray) (conv->blocksize + 1);
  // single-channel FFT plan for temporary signal and ir blocks
  conv->fftplan_xtemp = fftwf_plan_dft_r2c_1d(conv->blocksize * 2, conv->xtemp,
                                              conv->xftemp, FFTW_ESTIMATE);
  conv->fftplan_htemp = fftwf_plan_dft_r2c_1d(conv->blocksize * 2, conv->htemp,
                                              conv->hftemp, FFTW_ESTIMATE);
  // single-channel IFFT plan for temporary main/old and new output signal
  conv->ifftplan_y =
      fftwf_plan_dft_c2r_1d(conv->blocksize * 2, conv->yftemp, conv->y, FFTW_ESTIMATE);
  conv->ifftplan_y_cf = fftwf_plan_dft_c2r_1d(conv->blocksize * 2, conv->yftemp,
                                              conv->y_cf, FFTW_ESTIMATE);
  return conv;
#endif
  return 0;
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(freeConvolution) (conv_data *conv) {
#if HAVE_FFTWF
  // single-channel FFT and IFFT plans of temporary signals
  fftwf_destroy_plan(conv->fftplan_xtemp);
  fftwf_destroy_plan(conv->fftplan_htemp);
  fftwf_destroy_plan(conv->ifftplan_y);
  fftwf_destroy_plan(conv->ifftplan_y_cf);
  // frequency-domain array holding all input signal partitions
  IEMCONVOLVE(free3DComplexArray) (conv->xf, conv->num_inputs, conv->num_partitions);
  // frequency-domain array with all MIMO ir partitions
  IEMCONVOLVE(free5DComplexArray)(conv->hf, NUM_CF, conv->num_outputs, conv->num_inputs,
                     conv->num_partitions);
  // old input signal-block buffers
  IEMCONVOLVE(free2DArray) (conv->x_old, conv->num_inputs);
  // time-domain temporary FFT signal blocks
  fftwf_free(conv->yftemp);
  fftwf_free(conv->xtemp);
  fftwf_free(conv->y);
  fftwf_free(conv->w_old);
  fftwf_free(conv->w_new);
  fftwf_free(conv->y_cf);
  fftwf_free(conv->htemp);
  // frequency-domain temporary FFT blocks
  fftwf_free(conv->hftemp);
  fftwf_free(conv->xftemp);
  free(conv);
#endif
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(setImpulseResponseZeroPad) (conv_data *conv, float ***inh, int num_samples, _Bool no_xfade_init) {
  int offset;
  int copy_length;
  int hot_or_cold_stream;
#if HAVE_FFTWF
  if (!conv)
    return;
  if (no_xfade_init) { // update hot stream directly, without crossfading, for initialization
    hot_or_cold_stream = conv->current_cf;
    IEMCONVOLVE(registerCrossFadeComplete) (conv);
  } else { // update cold stream and register a cross-fade from hot -> cold stream next dsp cycle
    hot_or_cold_stream = (conv->current_cf + 1) % NUM_CF;
    IEMCONVOLVE(registerCrossFade) (conv);
  }
  for (int out_ch = 0; out_ch < conv->num_outputs; out_ch++) {
    for (int in_ch = 0; in_ch < conv->num_inputs; in_ch++) {
      offset = 0;
      for (int partition = 0; partition < conv->num_partitions; partition++) {
        copy_length = (offset+conv->blocksize > num_samples) ? num_samples-offset : conv->blocksize;
        IEMCONVOLVE(copyArray) (&inh[out_ch][in_ch][offset], conv->htemp,
                  copy_length);
        IEMCONVOLVE(resetArray) (&conv->htemp[copy_length], conv->blocksize - copy_length);
        fftwf_execute(conv->fftplan_htemp);
        IEMCONVOLVE(copyComplexArray) (
            conv->hftemp,
            conv->hf[hot_or_cold_stream][out_ch][in_ch][partition],
            conv->blocksize + 1);
        offset+=conv->blocksize;
      }
    }
  }
#endif
}
