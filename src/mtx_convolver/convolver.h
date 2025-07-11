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
#if HAVE_FFTWF
# include <fftw3.h>
#else
# include "stub/fftwf.h"
#endif

int IEMCONVOLVE(convolver_set_fftwf_functions) (const t_fftwf_functions*funs);


#define NUM_CF 2 // there are 2 crossfase buffers (re-occurring array dimension)
typedef struct Conv_data {
  unsigned int blocksize;          // signal block length (fft length = 2L)
  unsigned int num_partitions;          // number of convolution partitions (length of h)
  unsigned int num_inputs; // number of input channels
  unsigned int num_outputs;
  unsigned int xfade_length;
  _Bool register_crossfade; // parameter for switching between different buffers

  float *xtemp;   // 2L fft-input time-domain signal x=[xprev, xcurr]
  float *htemp;   // 2L time-domain impulse response input h=[h, 0]
  float *y;       // 2L time-domain result y=[ycircular, ylinear]
  float *y_cf;    // 2L time-domain result y=[ycircular, ylinear]
  float **x_old;  // previous input data
  float *w_old;   // fade-out window
  float *w_new;   // fade-in window
  unsigned int current_rb; // current_rb ring buffer position
  unsigned int current_cf;

  fftwf_complex ***xf;   // L+1 positive-half DFT, partition input ring buffer
  fftwf_complex *****hf; // L+1 positive-half DFT, partition stack of h
  fftwf_complex *yftemp; // L+1 positive-half DFT output buffer
  fftwf_complex *xftemp; // L+1 FFTW buffer for previous+current block
  fftwf_complex *hftemp; // L+1 FFTW buffer for response partition
  fftwf_plan fftplan_xtemp; // FFTW DFT plan for previous+current block
  fftwf_plan fftplan_htemp; // FFTW DFT plan for impulse response
  fftwf_plan ifftplan_y;    // FFTW iDFT plan for current output signal
  fftwf_plan ifftplan_y_cf;
} conv_data;
/* crossfade functions */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(crossFade) (float *y, float *y_new, float *w_old, float *w_new, unsigned int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(registerCrossFade) (conv_data *conv);
/*-----------------------------------------------------------------------------------------------------------------------------*/
_Bool IEMCONVOLVE(wasCrossFadeRegistered)(conv_data *conv);
/*-----------------------------------------------------------------------------------------------------------------------------*/
/* PARTITIONED CONVOLUTION CORE */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(convProcess) (conv_data *conv, float **in, float **out);
/*-----------------------------------------------------------------------------------------------------------------------------*/
conv_data *IEMCONVOLVE(initConvolution) (unsigned int blocksize, unsigned int num_partitions, unsigned int xfade_length, unsigned int num_inputs, unsigned int num_outputs, _Bool coherent_xfade);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(freeConvolution) (conv_data *conv);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(setImpulseResponseZeroPad) (conv_data *conv, float ***inh, unsigned int num_samples, _Bool no_xfade_init);

#endif /* _mtx_convolver_convolver_h */
