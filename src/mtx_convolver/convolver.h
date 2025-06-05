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

#include "array.h"
#include <fftw3.h>
#define NUM_CF 2 // there are 2 crossfase buffers (re-occurring array dimension)
typedef struct Conv_data {
  int L;          // signal block length (fft length = 2L)
  int P;          // number of convolution partitions (length of h)
  int num_inputs; // number of input channels
  int num_outputs;
  int Hann_len;
  _Bool convolver_switch; // parameter for swiching bitween different buffers

  float *xtemp;   // 2L fft-input time-domain signal x=[xprev, xcurr]
  float *htemp;   // 2L time-domain impulse response input h=[h, 0]
  float *y;       // 2L time-domain result y=[ycircular, ylinear]
  float *y_cf;    // 2L time-domain result y=[ycircular, ylinear]
  float **x_old;  // previous input data
  float *w;       // hann-function
  int current_rb; // current_rb ring buffer position
  int current_cf;
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
void crossfade(float *y, float *y_new, float *w, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void setNewIR(conv_data *conv, _Bool status);
/*-----------------------------------------------------------------------------------------------------------------------------*/
_Bool getNewIR(conv_data *conv);
/*-----------------------------------------------------------------------------------------------------------------------------*/
/* PARTITIONED CONVOLUTION CORE */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void conv_process(conv_data *conv, float **in, float **out);
/*-----------------------------------------------------------------------------------------------------------------------------*/
conv_data *initConvolution(int L, int P, int Hann_len, int num_inputs, int num_outputs);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void freeConvolution(conv_data *conv);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void setImpulseResponse(conv_data *conv, float ***inh);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void setImpulseResponseZeroPad(conv_data *conv, float ***inh, int num_samples);
