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
#include "iemmatrix_stub.h"

#define stringify(s) str(s)
#define str(s) #s


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


#ifdef HAVE_FFTW
# include <fftw3.h>
#else
# include "stub/fftwf.h"
#endif

#if PD_FLOATSIZE == 64
#define fftwx(a) fftw_##a
#else
#define fftwx(a) fftwf_##a
#endif




typedef fftwx(complex)  t_fftw_complex;
typedef fftwx(plan)     t_fftw_plan;

typedef void(*t_fftw_destroy_plan)(t_fftw_plan);
typedef void(*t_fftw_execute)(const t_fftw_plan);
typedef t_fftw_plan(*t_fftw_plan_dft_1d)(int, t_fftw_complex*, t_fftw_complex*, signed, unsigned);
typedef t_fftw_plan(*t_fftw_plan_dft_c2r_1d)(int, t_fftw_complex*, t_float*, unsigned);
typedef t_fftw_plan(*t_fftw_plan_dft_r2c_1d)(int, t_float*, t_fftw_complex*, unsigned);

static t_fftw_destroy_plan my_fftw_destroy_plan = 0;
static t_fftw_execute my_fftw_execute = 0;
static t_fftw_plan_dft_1d my_fftw_plan_dft_1d = 0;
static t_fftw_plan_dft_c2r_1d my_fftw_plan_dft_c2r_1d = 0;
static t_fftw_plan_dft_r2c_1d my_fftw_plan_dft_r2c_1d = 0;

static int have_fftw = 0;


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

/* the actual FFT wrapper implementation */

struct iemfft_plan {
  /* data owned by the caller of iemmatrix_..._plan() */
  unsigned int n0;

  t_complexfloat*cf_in, *cf_out;
  t_complexdouble*cd_in, *cd_out;

  float*f_in, *f_out;
  double*d_in, *d_out;


  /* the rest is owned by us */

  t_float *_re, *_im; /* de-interleave complex data for mayer-fft */
  void *_in; /* in case the input data must not be modified */

  t_iemfft_flag flags;
  int inverse;
  t_fftw_plan fftw_plan;
};

void*iemfft_malloc(size_t size) {
  return malloc(size);
}
void iemfft_free(void*data) {
  free(data);
}

void iemfft_destroy_plan(t_iemfft_plan plan) {
  if(plan->fftw_plan) {
    my_fftw_destroy_plan(plan->fftw_plan);
  }
  iemfft_free(plan->_re);
  iemfft_free(plan->_im);
  iemfft_free(plan->_in);
  freebytes(plan, sizeof(*plan));
}

t_iemfft_plan iemfft_plan_rifft_1d(int n0, t_complex*in, t_float*out, t_iemfft_flag flags) {
  t_iemfft_plan plan = calloc(1, sizeof(*plan));
  plan->flags = flags;
  plan->n0 = n0;
  plan->c_in = in;
  plan->x_out = out;

  plan->inverse = 1;
  if(have_fftw)
    plan->fftw_plan = my_fftw_plan_dft_c2r_1d(n0, (t_fftw_complex*)in, out, FFTW_ESTIMATE);

  if(!plan->fftw_plan && (flags & PRESERVE_INPUT)) {
    plan->_in = iemfft_malloc(n0 * sizeof(*in));
  }

  return plan;
}
t_iemfft_plan iemfft_plan_rfft_1d(int n0, t_float*in, t_complex*out, t_iemfft_flag flags) {
  t_iemfft_plan plan = calloc(1, sizeof(*plan));
  plan->flags = flags;
  plan->n0 = n0;
  plan->x_in = in;
  plan->c_out = out;
  plan->inverse = 0;
  if(have_fftw)
    plan->fftw_plan = my_fftw_plan_dft_r2c_1d(n0, in, (t_fftw_complex*)out, FFTW_ESTIMATE);

  if(!plan->fftw_plan && (flags & PRESERVE_INPUT)) {
    plan->_in = iemfft_malloc(n0 * sizeof(*in));
  }
  return plan;
}

static t_iemfft_plan _fft_plan_1d(int n0, t_complex*in, t_complex* out, int inverse, t_iemfft_flag flags) {
  t_iemfft_plan plan = calloc(1, sizeof(*plan));
  plan->flags = flags;
  plan->n0 = n0;
  plan->c_in = in;
  plan->c_out = out;
  plan->inverse = inverse;
  if(have_fftw)
    plan->fftw_plan = my_fftw_plan_dft_1d(n0, (t_fftw_complex*)in, (t_fftw_complex*)out, inverse?+1:-1, FFTW_ESTIMATE);

  if(!plan->fftw_plan) {
    plan->_re = iemfft_malloc(n0 * sizeof(*plan->_re));
    plan->_im = iemfft_malloc(n0 * sizeof(*plan->_im));
  }
  return plan;
}

t_iemfft_plan iemfft_plan_fft_1d(int n0, t_complex*in, t_complex* out, t_iemfft_flag flags) {
  return _fft_plan_1d(n0, in, out, 0, flags);
}

t_iemfft_plan iemfft_plan_ifft_1d(int n0, t_complex*in, t_complex* out, t_iemfft_flag flags) {
  return _fft_plan_1d(n0, in, out, 1, flags);
}



void iemfft_execute(const t_iemfft_plan plan) {
  if(plan->fftw_plan) {
    my_fftw_execute(plan->fftw_plan);
  } else {
    if (0) {
      // real-valued FFTs
    } else if(plan->x_in && plan->c_out && !plan->inverse) {
      t_float*in = plan->x_in;
      if(plan->_in) {
        memcpy(plan->_in, plan->x_in, sizeof(*in)*plan->n0);
        in = plan->_in;
      }
      /* native precision rFFT */
      mayer_realfft(plan->n0, in);
      mayerfloat2complex(plan->n0, in, plan->c_out);
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


t_iemfft_backend iemfft_init(t_class*c) {
  static int tried_fftw = 0;
  if(!tried_fftw) {
    /* initialize stubs */
    tried_fftw = 1;


    /*
#define register_my_fftw(prefix, name)                           \
    my_fftw_##name = iemmatrix_get_stub(stringify(prefix) "_" #name, c)
#define show_addr(x) post(#x "\t%p", x)
    */
#define register_my_fftw(name)                           \
    my_fftw_##name = iemmatrix_get_stub(stringify(fftwx(name)), c)


    register_my_fftw(destroy_plan);

    register_my_fftw(destroy_plan);
    register_my_fftw(execute);
    register_my_fftw(plan_dft_1d);
    register_my_fftw(plan_dft_c2r_1d);
    register_my_fftw(plan_dft_r2c_1d);

    have_fftw = (1
                 && my_fftw_destroy_plan
                 && my_fftw_execute
                 && my_fftw_plan_dft_1d
                 && my_fftw_plan_dft_c2r_1d
                 && my_fftw_plan_dft_r2c_1d
                 );

    if(have_fftw) {
      post("iemmatrix: using FFTW for Fourier transforms");
    } else {
      post("iemmatrix: using built in Mayer/Ooura for Fourier transforms");
    }
  }

  if(have_fftw) {
    return FFTW;
  }

  return MAYER;
}
