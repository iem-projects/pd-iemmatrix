/*
 *  iemmatrix
 *
 *  objects for manipulating simple matrices
 *  mostly referring to matlab/octave matrix functions
 *
 * Copyright (c) 2025, IOhannes m zmölnig
 * IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */

/* iemmatrix_mayer: fftw(f) compatible interface to Pd's mayer_fft */

#include "iemmatrix_fft.h"

#include <stdlib.h>

#if PD_FLOATSIZE == 64
# define x_in d_in
# define x_out d_out
# define c_in cd_in
# define c_out cd_out
#else
# define x_in f_in
# define x_out f_out
# define c_in cf_in
# define c_out cf_out
#endif


/* helpers */
void mayerfloat2complex(unsigned int N, const t_float *in, t_complex *out) {
  /* mayer complex to interleaved complex */
}
void complex2mayerfloat(unsigned int N, const t_complex *in, t_float *out) {
  /* interleaved complex to mayer complex */
}
void complex2mayercomplex(unsigned int N, const t_complex *in, t_float* re, t_float *im) {
  /* interleaved complex to separate complex */
}
void mayercomplex2complex(unsigned int N, const t_float *im, t_float *re, t_complex* out) {
  /* separate complex to interleaved complex */
}

typedef struct iemmatrix_fft_plan {
  unsigned int n0;

  t_complexfloat*cf_in, *cf_out;
  t_complexdouble*cd_in, *cd_out;

  float*f_in, *f_out;
  double*d_in, *d_out;

  t_float *_re, *_im;

  int inverse;
  void*fftw_plan;
} t_iemmatrix_fft_plan;

void*iemmatrix_fft_malloc(size_t size) {
  return malloc(size);
}
void iemmatrix_fft_free(void*data) {
  free(data);
}

void iemmatrix_fft_destroy_plan(t_iemmatrix_fft_plan*plan) {
  freebytes(plan, sizeof(*plan));
}

t_iemmatrix_fft_plan*iemmatrix_rifft_plan_1d(int n0, t_complex*in, t_float*out) {
  t_iemmatrix_fft_plan*plan = calloc(1, sizeof(t_iemmatrix_fft_plan));
  plan->n0 = n0;
  plan->c_in = in;
  plan->x_out = out;
  plan->inverse = 1;
  return plan;
}
t_iemmatrix_fft_plan*iemmatrix_rfft_plan_1d(int n0, t_float*in, t_complex*out) {
  t_iemmatrix_fft_plan*plan = calloc(1, sizeof(t_iemmatrix_fft_plan));
  plan->n0 = n0;
  plan->x_in = in;
  plan->c_out = out;
  plan->inverse = 0;
  return plan;
}

t_iemmatrix_fft_plan*iemmatrix_fft_plan_1d(int n0, t_complex*in, t_complex* out) {
  t_iemmatrix_fft_plan*plan = calloc(1, sizeof(t_iemmatrix_fft_plan));
  plan->n0 = n0;
  plan->c_in = in;
  plan->c_out = out;
  plan->inverse = 0;
  return plan;
}

t_iemmatrix_fft_plan*iemmatrix_ifft_plan_1d(int n0, t_complex*in, t_complex* out) {
  t_iemmatrix_fft_plan*plan = iemmatrix_fft_plan_1d(n0, in, out);
  plan->inverse = 1;
  return plan;
}



void mayer_execute(const t_iemmatrix_fft_plan*plan) {
#warning total nonsense
  if(plan->fftw_plan) {
#if 0
    /* FFTW */
    /*
      IN.real : f_in [0..N]
      OUT.real: c_out[0..(N>>1 +1)][0]
      OUT.imag: c_out[0..(N>>1 +1)][1]
    */
    my_plan_dft_r2c_1d(N, f_in, c_out);

    /*
      IN.real : c_in [0..(N>>1 +1)][0] ?
      IN.imag : c_in [0..(N>>1 +1)][1] ?
      OUT.real: f_out[0..N]            ?
    */
    my_plan_dft_c2r_1d(N, c_in, f_out);

    /*
      IN.real : c_in [0..N][0] ?
      IN.imag : c_in [0..N][1] ?
      OUT.real: c_out[0..N][0] ?
      OUT.imag: c_out[0..N][1] ?
    */
    fftw_plan_dft_1d(N, c_in, c_out, direction);
#endif
  } else {
    if (0) {
      // real-valued FFTs
    } else if(plan->x_in && plan->c_out && !plan->inverse) {
      /* native precision rFFT */

      /* inplace rFFT
	 IN.real :  x_in[    0..(N    -1)]
	 OUT.real:  x_in[    0..(N>>1   )]
	 OUT.imag: -x_in[(N-1)..(N>>1 +1)], y[0]=y[N>>1+1]=0
      */
      mayer_realfft(plan->n0, plan->x_in);
      mayerfloat2complex(plan->n0, plan->x_in, plan->c_out);
    } else if(plan->c_in && plan->x_out && plan->inverse) {
      /* native precision rIFFT */
      complex2mayerfloat(plan->n0, plan->c_in, plan->x_out);
	/* inplace rIFFT
	   IN.real :  x_in[    0..(N>>1   )]
	   IN.imag : -x_in[(N-1)..(N>>1 +1)], y[0]=y[N>>1+1]=0
	   OUT.real:  x_in[    0..(N    -1)]
	*/
      mayer_realifft(plan->n0, plan->x_out);

      // complex-valued FFTs
    } else if(plan->c_in && plan->c_out && !plan->inverse) {
      /* native precision FFT */
      complex2mayercomplex(plan->n0, plan->c_in, plan->_re, plan->_im);
	/* inplace FFT
	   IN.real :  x_re[0..(N-1)]
	   IN.imag :  x_im[0..(N-1)]
	   OUT.real:  x_re[0..(N-1)]
	   OUT.imag:  x_im[0..(N-1)]
	*/
      mayer_fft(plan->n0, plan->_re, plan->_im);
      mayercomplex2complex(plan->n0, plan->_re, plan->_im, plan->c_out);
    } else if(plan->c_in && plan->c_out && plan->inverse) {
      /* native precision IFFT */
      complex2mayercomplex(plan->n0, plan->c_in, plan->_re, plan->_im);
      /*
	   IN.real :  x_re[0..(N-1)]
	   IN.imag :  x_im[0..(N-1)]
	   OUT.real:  x_re[0..(N-1)]
	   OUT.imag:  x_im[0..(N-1)]
      */
      mayer_ifft(plan->n0, plan->_re, plan->_im);
      mayercomplex2complex(plan->n0, plan->_re, plan->_im, plan->c_out);
    } else {
      pd_error(0, "iemmatrix_fft: no valid plan!");
    }
  }
}


int iemmatrix_fft_init(t_class*c) {
  /* initialize stubs */
  return 0;
}
