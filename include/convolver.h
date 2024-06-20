#include "fftw3.h"
#include "array.h"
#define NUM_CF 2 // there are 2 crossfase buffers (re-occurring array dimension)
typedef struct Conv_data{
    int L;  // signal block length (fft length = 2L)
    int P;  // number of convolution partitions (length of h)
    int INPUT_Channel_Number; //number of input channels
    int OUTPUT_Channel_Number;
    int Hann_len;
    _Bool convolver_switch; //parameter for swiching bitween different buffers

    float *xtemp; // 2L fft-input time-domain signal x=[xprev, xcurr]  
    float *htemp; // 2L time-domain impulse response input h=[h, 0]
    float *y; // 2L time-domain result y=[ycircular, ylinear]
    float *y_cf; // 2L time-domain result y=[ycircular, ylinear] 
    float **x_old; //previous input data
    float *w; // hann-function
    int current_rb; // current_rb ring buffer position
    int current_cf;
    fftwf_complex ***xf; // L+1 positive-half DFT, partition input ring buffer
    fftwf_complex *****hf; // L+1 positive-half DFT, partition stack of h
    fftwf_complex *yftemp;  // L+1 positive-half DFT output buffer
    fftwf_complex *xftemp; // L+1 FFTW buffer for previous+current block
    fftwf_complex *hftemp; // L+1 FFTW buffer for response partition
    fftwf_plan fftplan_xtemp; // FFTW DFT plan for previous+current block
    fftwf_plan fftplan_htemp; // FFTW DFT plan for impulse response
    fftwf_plan ifftplan_y;    // FFTW iDFT plan for current output signal
    fftwf_plan ifftplan_y_cf;
    
} conv_data;
/* crossfade functions */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void crossfade(float *y,float* y_new, float* w, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void setNewIR(conv_data *conv, _Bool status);
/*-----------------------------------------------------------------------------------------------------------------------------*/
_Bool getNewIR(conv_data *conv);
/*-----------------------------------------------------------------------------------------------------------------------------*/
/* PARTITIONED CONVOLUTION CORE */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void conv_process(conv_data *conv, float **in, float **out);
/*-----------------------------------------------------------------------------------------------------------------------------*/
conv_data* initConvolution(int L, int P, int Hann_len, int in_ch, int out_ch);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void freeConvolution(conv_data *conv);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void setImpulseResponse(conv_data *conv, float***inh);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void setImpulseResponse2DZeropad(conv_data *conv, float**inh, int num_samples);