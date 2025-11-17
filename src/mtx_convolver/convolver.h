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
#ifndef _mtx_convolver_convolver_h
#define _mtx_convolver_convolver_h

#include "array.h"

#include <iemmatrix_fft.h>

int IEMCONVOLVE(convolver_set_fftwf_functions) (void);

#define NUM_CF 2 // there are 2 crossfase buffers (re-occurring array dimension)
typedef struct Conv_data {
  unsigned int blocksize;          // signal block length (fft length = 2L)
  unsigned int num_partitions;          // number of convolution partitions (length of h)
  unsigned int num_inputs; // number of input channels
  unsigned int num_outputs;
  unsigned int xfade_length;
  _Bool register_crossfade; // parameter for switching between different buffers

  t_float *xtemp;   // 2L fft-input time-domain signal x=[xprev, xcurr]
  t_float *htemp;   // 2L time-domain impulse response input h=[h, 0]
  t_float *y;       // 2L time-domain result y=[ycircular, ylinear]
  t_float *y_cf;    // 2L time-domain result y=[ycircular, ylinear]
  t_float **x_old;  // previous input data
  t_float *w_old;   // fade-out window
  t_float *w_new;   // fade-in window
  unsigned int current_rb; // current_rb ring buffer position
  unsigned int current_cf;

  t_complex ***xf;   // L+1 positive-half DFT, partition input ring buffer
  t_complex *****hf; // L+1 positive-half DFT, partition stack of h
  t_complex *yftemp; // L+1 positive-half DFT output buffer
  t_complex *xftemp; // L+1 FFTW buffer for previous+current block
  t_complex *hftemp; // L+1 FFTW buffer for response partition
  t_iemfft_plan fftplan_xtemp; // FFTW DFT plan for previous+current block
  t_iemfft_plan fftplan_htemp; // FFTW DFT plan for impulse response
  t_iemfft_plan ifftplan_y;    // FFTW iDFT plan for current output signal
  t_iemfft_plan ifftplan_y_cf;
} conv_data;
/* crossfade functions */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(crossFade) (t_float *y, t_float *y_new, t_float *w_old, t_float *w_new, unsigned int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(registerCrossFade) (conv_data *conv);
/*-----------------------------------------------------------------------------------------------------------------------------*/
_Bool IEMCONVOLVE(wasCrossFadeRegistered)(conv_data *conv);
/*-----------------------------------------------------------------------------------------------------------------------------*/
/* PARTITIONED CONVOLUTION CORE */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(convProcess) (conv_data *conv, t_float **in, t_float **out);
/*-----------------------------------------------------------------------------------------------------------------------------*/
conv_data *IEMCONVOLVE(initConvolution) (unsigned int blocksize, unsigned int num_partitions, unsigned int xfade_length, unsigned int num_inputs, unsigned int num_outputs, _Bool coherent_xfade);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(freeConvolution) (conv_data *conv);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(setImpulseResponseZeroPad) (conv_data *conv, t_float ***inh, unsigned int num_samples, _Bool no_xfade_init);

#endif /* _mtx_convolver_convolver_h */
