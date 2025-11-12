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


/* use cases
   mtx_convolve~
   - FFTWF:
     - IN:  real-valued pd_in->in => convProcess => real-valued out->pd_out
     - OUT:
   - mayer
     - IN: ---
     - OUT: ---

   mtx_fft
   - FFTW:
     - IN: ---
     - OUT: ---
   - mayer
     - IN: separate re/im atomslists are converted to separate re/im floatlists
     - OUT: re/im floatlists are scaled & then converted into separate re/im atomlists

   mtx_ifft
   - FFTW:
     - IN: ---
     - OUT: ---
   - mayer
     - IN: separate re/im atomslists are converted to separate re/im floatlists
     - OUT: re/im floatlists are scaled & then converted into separate re/im atomlists

   mtx_rifft
   - FFTW:
     - IN: re/im atomlists is written directly into interleaved data vector
     - OUT: scale output vector & convert to atomlist
   - mayer:
     - IN: atomlist is converted to t_float & then rearranged for mayer blocked data
     - OUT: scale output vector & convert to atomlist

   mtx_rfft
   - FFTW:
     - IN: atomlist is written directly into real-valued data vector
     - OUT: interleaved output vector is converted to separate re/im atomlists
   - mayer:
     - IN: atomlist is written directly into real-valued data vector
     - OUT: imag-values are extracted from mayer blocked data & then written to separate re/im atomlists

 */

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
static void mayerfloat2complex(const unsigned int N, const t_float *in, t_complex *out) {
  /* mayer complex to interleaved complex

    IN.real :  in[    0..(N>>1   )]
    IN.imag : -in[(N-1)..(N>>1 +1)], y[0]=y[N>>1+1]=0
    OUT.real: out[0..(N>>1 +1)].re
    OUT.imag: out[0..(N>>1 +1)].im

    OUT.re: +in[0], +in[1], ..., +in[(N>>2]
  */
  const unsigned int N2 = (N>>1);

  for(unsigned int n=0; n<N2; n++) {
    out[   n].re =  in[n];
    out[N2-n].im = -in[N2+n];
  }
  out[N2].re = in[N2];
  out[0].im = out[N2].im = 0.;
}
static void complex2mayerfloat(const unsigned int N, const t_complex *in, t_float *out) {
  /* interleaved complex to mayer complex

  OUT: +in[0].re, +in[1].re, ..., +in[(N>>1)].re, -in[(N>>1)-1].im, -in[(N>>1)-2].im, ..., -in[0].im
  */
  const unsigned int N2 = (N>>1) + 1;

  for(unsigned int n=0; n+2<N2; n++) {
    out[n] = in[n].re;
    out[N2+n] = -in[N2-(n+2)].im;
  }
  out[N2-2] = in[N2-2].re;
  out[N2-1] = in[N2-1].re;
}
static void complex_deinterleave(unsigned int N, const t_complex *in, t_float* re, t_float *im) {
  /* interleaved complex to separate complex */
  for(unsigned int n=0; n<N; n++) {
    re[n] = in[n].re;
    im[n] = in[n].im;
  }
}
static void complex_interleave(unsigned int N, const t_float *re, t_float *im, t_complex* out) {
  /* separate complex to interleaved complex */
  for(unsigned int n=0; n<N; n++) {
    out[n].re = re[n];
    out[n].im = im[n];
  }
}

struct iemmatrix_fft_plan {
  unsigned int n0;

  t_complexfloat*cf_in, *cf_out;
  t_complexdouble*cd_in, *cd_out;

  float*f_in, *f_out;
  double*d_in, *d_out;

  t_float *_re, *_im;

  int inverse;
  void*fftw_plan;
};

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



void iemmatrix_fft_execute(const t_iemmatrix_fft_plan*plan) {
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
      mayer_realfft(plan->n0, plan->x_in);
      mayerfloat2complex(plan->n0, plan->x_in, plan->c_out);
    } else if(plan->c_in && plan->x_out && plan->inverse) {
      /* native precision rIFFT */
      complex2mayerfloat(plan->n0, plan->c_in, plan->x_out);
      mayer_realifft(plan->n0, plan->x_out);

      // complex-valued FFTs
    } else if(plan->c_in && plan->c_out && !plan->inverse) {
      /* native precision FFT */
      complex_deinterleave(plan->n0, plan->c_in, plan->_re, plan->_im);
	/* inplace FFT
	   IN.real :  x_re[0..(N-1)]
	   IN.imag :  x_im[0..(N-1)]
	   OUT.real:  x_re[0..(N-1)]
	   OUT.imag:  x_im[0..(N-1)]
	*/
      mayer_fft(plan->n0, plan->_re, plan->_im);
      complex_interleave(plan->n0, plan->_re, plan->_im, plan->c_out);
    } else if(plan->c_in && plan->c_out && plan->inverse) {
      /* native precision IFFT */
      complex_deinterleave(plan->n0, plan->c_in, plan->_re, plan->_im);
      /*
	   IN.real :  x_re[0..(N-1)]
	   IN.imag :  x_im[0..(N-1)]
	   OUT.real:  x_re[0..(N-1)]
	   OUT.imag:  x_im[0..(N-1)]
      */
      mayer_ifft(plan->n0, plan->_re, plan->_im);
      complex_interleave(plan->n0, plan->_re, plan->_im, plan->c_out);
    } else {
      pd_error(0, "iemmatrix_fft: no valid plan!");
    }
  }
}


int iemmatrix_fft_init(t_class*c) {
  /* initialize stubs */
  return 0;
}
