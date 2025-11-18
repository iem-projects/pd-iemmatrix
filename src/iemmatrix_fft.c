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
# include "stub/fftw.h"
#endif

typedef void*(*t_my_malloc)(size_t);
typedef void (*t_my_free)(void*);

/* fftwf */
typedef void(*t_fftw32_destroy_plan)(fftwf_plan);
typedef void(*t_fftw32_execute)(const fftwf_plan);
typedef fftwf_plan(*t_fftw32_plan_dft_1d)(int, t_complexfloat*, t_complexfloat*, signed, unsigned);
typedef fftwf_plan(*t_fftw32_plan_dft_c2r_1d)(int, t_complexfloat*, float*, unsigned);
typedef fftwf_plan(*t_fftw32_plan_dft_r2c_1d)(int, float*, t_complexfloat*, unsigned);

static t_my_malloc _fftw32_malloc = 0;
static t_my_free _fftw32_free = 0;
static t_fftw32_destroy_plan _fftw32_destroy_plan = 0;
static t_fftw32_execute _fftw32_execute = 0;
static t_fftw32_plan_dft_1d _fftw32_plan_dft_1d = 0;
static t_fftw32_plan_dft_c2r_1d _fftw32_plan_dft_c2r_1d = 0;
static t_fftw32_plan_dft_r2c_1d _fftw32_plan_dft_r2c_1d = 0;

static int have_fftw32 = 0;


/* fftw */
typedef void(*t_fftw64_destroy_plan)(fftw_plan);
typedef void(*t_fftw64_execute)(const fftw_plan);
typedef fftw_plan(*t_fftw64_plan_dft_1d)(int, t_complexdouble*, t_complexdouble*, signed, unsigned);
typedef fftw_plan(*t_fftw64_plan_dft_c2r_1d)(int, t_complexdouble*, double*, unsigned);
typedef fftw_plan(*t_fftw64_plan_dft_r2c_1d)(int, double*, t_complexdouble*, unsigned);

static t_my_malloc _fftw64_malloc = 0;
static t_my_free _fftw64_free = 0;
static t_fftw64_destroy_plan _fftw64_destroy_plan = 0;
static t_fftw64_execute _fftw64_execute = 0;
static t_fftw64_plan_dft_1d _fftw64_plan_dft_1d = 0;
static t_fftw64_plan_dft_c2r_1d _fftw64_plan_dft_c2r_1d = 0;
static t_fftw64_plan_dft_r2c_1d _fftw64_plan_dft_r2c_1d = 0;

static int have_fftw64 = 0;



/* helpers */
static t_iemfft_flag flags2precision(t_iemfft_flag flags) {
  /* returns: DEFAULT(mayer), SINGLE(fftwf), or DOUBLE(fftw) */
  int want_precision = 0;
  if(flags & PREFER_SINGLE)
    want_precision = 32;
  if(flags & PREFER_DOUBLE)
    want_precision = 64;

#if PD_FLOATSIZE == 32
  if(!(have_fftw32 || have_fftw64))
    return DEFAULT;
  if(!have_fftw64)
    return PREFER_SINGLE;
  if((64 == want_precision) || !have_fftw32)
    return PREFER_DOUBLE;
  return PREFER_SINGLE;
#else
  if((32 == want_precision) && have_fftw32)
    return PREFER_SINGLE;
  if(have_fftw64)
    return PREFER_DOUBLE;
#endif
  return DEFAULT;
}

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

static void single2double(double*dst, const float*src, size_t n) {
  while(n--)
    *dst++ = *src++;
}
static void double2single(float*dst, const double*src, size_t n) {
  while(n--)
    *dst++ = *src++;
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
  void *_in; /* in case the input data must not be modified; or we need to type-convert */
  void *_out; /* for type-converting the output data */

  /* these are needed for conversions */
  void*inptr,*outptr; /* just pointers to the actually used input/output pointers */
  size_t insize, outsize; /* the actual number of elements in the input/output poitners */

  t_iemfft_flag flags;
  int inverse;

  fftwf_plan fftw32_plan;
  fftw_plan  fftw64_plan;
};

void*iemfft_malloc(size_t size) {
  if(_fftw32_malloc)
    return  _fftw32_malloc(size);
  if(_fftw64_malloc)
    return  _fftw64_malloc(size);
  return malloc(size);
}
void iemfft_free(void*data) {
  if(_fftw32_free)
    return  _fftw32_free(data);
  if(_fftw64_free)
    return  _fftw64_free(data);
  free(data);
}

void iemfft_destroy_plan(t_iemfft_plan plan) {
  if(plan->fftw32_plan) {
    _fftw32_destroy_plan(plan->fftw32_plan);
  }
  if(plan->fftw64_plan) {
    _fftw64_destroy_plan(plan->fftw64_plan);
  }
  iemfft_free(plan->_re);
  iemfft_free(plan->_im);
  iemfft_free(plan->_in);
  iemfft_free(plan->_out);
  freebytes(plan, sizeof(*plan));
}

t_iemfft_plan iemfft_plan_rifft_1d(int n0, t_complex*in, t_float*out, t_iemfft_flag flags) {
  t_iemfft_plan plan = calloc(1, sizeof(*plan));
  const t_iemfft_flag precision = flags2precision(flags);
  plan->flags = flags;
  plan->n0 = n0;

  plan->inptr = plan->c_in = in;
  plan->outptr = plan->x_out = out;
  plan->insize = ((n0>>1)+1)<<1;
  plan->outsize = n0;

  plan->inverse = 1;

#if PD_FLOATSIZE == 32
  t_complexdouble*_in = 0;
  double*_out = 0;
#else
  t_complexfloat*_in = 0;
  float*_out = 0;
#endif

  switch(precision) {
  case PREFER_SINGLE:
#if PD_FLOATSIZE != 32
    /* need to do type-conversion */
    _in = plan->_in = iemfft_malloc(plan->insize * sizeof(*_in));
    _out = plan->_out = iemfft_malloc(plan->outsize * sizeof(*_out));
    plan->fftw32_plan = _fftw32_plan_dft_c2r_1d(n0, _in, _out, FFTW_ESTIMATE);
#else
    plan->fftw32_plan = _fftw32_plan_dft_c2r_1d(n0, in, out, FFTW_ESTIMATE);
#endif
    break;
  case PREFER_DOUBLE:
#if PD_FLOATSIZE != 64
    /* need to do type-conversion */
    _in = plan->_in = iemfft_malloc(plan->insize * sizeof(*_in));
    _out = plan->_out = iemfft_malloc(plan->outsize * sizeof(*_out));
    plan->fftw64_plan = _fftw64_plan_dft_c2r_1d(n0, _in, _out, FFTW_ESTIMATE);
#else
    plan->fftw64_plan = _fftw64_plan_dft_c2r_1d(n0, in, out, FFTW_ESTIMATE);
#endif
    break;
  default:
    break;
    if(flags & PRESERVE_INPUT) {
      plan->_in = iemfft_malloc(n0 * sizeof(*in));
    }
  }

  return plan;
}
t_iemfft_plan iemfft_plan_rfft_1d(int n0, t_float*in, t_complex*out, t_iemfft_flag flags) {
  t_iemfft_plan plan = calloc(1, sizeof(*plan));
  const t_iemfft_flag precision = flags2precision(flags);
  plan->flags = flags;
  plan->n0 = n0;
  plan->inptr = plan->x_in = in;
  plan->outptr = plan->c_out = out;

  plan->insize = n0;
  plan->outsize = ((n0>>1)+1)<<1;

  plan->inverse = 0;

#if PD_FLOATSIZE == 32
  double*_in = 0;
  t_complexdouble*_out = 0;
#else
  float*_in = 0;
  t_complexfloat*_out = 0;
#endif

  switch(precision) {
  case PREFER_SINGLE:
#if PD_FLOATSIZE != 32
    /* need to do type-conversion */
    _in = plan->_in = iemfft_malloc(n0 * sizeof(*_in));
    _out = plan->_out = iemfft_malloc(n0 * sizeof(*_out));
    plan->fftw32_plan = _fftw32_plan_dft_r2c_1d(n0, _in, _out, FFTW_ESTIMATE);
#else
    plan->fftw32_plan = _fftw32_plan_dft_r2c_1d(n0, in, out, FFTW_ESTIMATE);
#endif
    break;
  case PREFER_DOUBLE:
#if PD_FLOATSIZE != 64
    /* need to do type-conversion */
    _in = plan->_in = iemfft_malloc(n0 * sizeof(*_in));
    _out = plan->_out = iemfft_malloc(n0 * sizeof(*_out));
    plan->fftw64_plan = _fftw64_plan_dft_r2c_1d(n0, _in, _out, FFTW_ESTIMATE);
#else
    plan->fftw64_plan = _fftw64_plan_dft_r2c_1d(n0, in, out, FFTW_ESTIMATE);
#endif
    break;
  default:
    if(flags & PRESERVE_INPUT) {
      plan->_in = iemfft_malloc(n0 * sizeof(*in));
    }
    break;
  }
  return plan;
}

static t_iemfft_plan _fft_plan_1d(int n0, t_complex*in, t_complex* out, int inverse, t_iemfft_flag flags) {
  t_iemfft_plan plan = calloc(1, sizeof(*plan));
  const t_iemfft_flag precision = flags2precision(flags);
  plan->flags = flags;
  plan->n0 = n0;
  plan->inptr = plan->c_in = in;
  plan->outptr = plan->c_out = out;
  plan->insize = plan->outsize = n0<<1;

  plan->inverse = inverse;

#if PD_FLOATSIZE == 32
  t_complexdouble*_in = 0, *_out = 0;
#else
  t_complexfloat*_in = 0, *_out = 0;
#endif

  switch(precision) {
  case PREFER_SINGLE:
#if PD_FLOATSIZE != 32
    /* need to do type-conversion */
    _in = plan->_in = iemfft_malloc(n0 * sizeof(*_in));
    _out = plan->_out = iemfft_malloc(n0 * sizeof(*_out));
    plan->fftw32_plan = _fftw32_plan_dft_1d(n0, _in, _out, inverse?+1:-1, FFTW_ESTIMATE);
#else
    plan->fftw32_plan = _fftw32_plan_dft_1d(n0, in, out, inverse?+1:-1, FFTW_ESTIMATE);
#endif
    break;
  case PREFER_DOUBLE:
#if PD_FLOATSIZE != 64
    /* need to do type-conversion */
    _in = plan->_in = iemfft_malloc(n0 * sizeof(*_in));
    _out = plan->_out = iemfft_malloc(n0 * sizeof(*_out));
    plan->fftw64_plan = _fftw64_plan_dft_1d(n0, _in, _out, inverse?+1:-1, FFTW_ESTIMATE);
#else
    plan->fftw64_plan = _fftw64_plan_dft_1d(n0, in, out, inverse?+1:-1, FFTW_ESTIMATE);
#endif
    break;
  default:
    plan->_re = iemfft_malloc(n0 * sizeof(*plan->_re));
    plan->_im = iemfft_malloc(n0 * sizeof(*plan->_im));
    break;
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
  if(0) {
  } else if(plan->fftw32_plan) {
#if PD_FLOATSIZE != 32
    if(plan->_in) {
      double2single(plan->_in, plan->inptr, plan->insize);
    }
#endif
    _fftw32_execute(plan->fftw32_plan);
#if PD_FLOATSIZE != 32
    if(plan->_out) {
      single2double(plan->outptr, plan->_out, plan->outsize);
    }
#endif
  } else if(plan->fftw64_plan) {
#if PD_FLOATSIZE != 64
    if(plan->_in) {
      single2double(plan->_in, plan->inptr, plan->insize);
    }
#endif
    _fftw64_execute(plan->fftw64_plan);
#if PD_FLOATSIZE != 64
    if(plan->_out) {
      double2single(plan->outptr, plan->_out, plan->outsize);
    }
#endif
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
  static int have_fftw = 0;
  static int tried_fftw = 0;
  if(!tried_fftw) {
    /* initialize stubs */
    tried_fftw = 1;


    /*
#define register_my_fftw(prefix, name)                           \
    my_fftw_##name = iemmatrix_get_stub(stringify(prefix) "_" #name, c)
#define show_addr(x) post(#x "\t%p", x)
    */


#define register_fftw32(name) \
    _fftw32_##name = iemmatrix_get_stub(stringify(fftwf_##name), c)
#define register_fftw64(name) \
    _fftw64_##name = iemmatrix_get_stub(stringify(fftw_##name), c)

    register_fftw32(malloc);
    register_fftw32(free);
    register_fftw32(destroy_plan);
    register_fftw32(execute);
    register_fftw32(plan_dft_1d);
    register_fftw32(plan_dft_c2r_1d);
    register_fftw32(plan_dft_r2c_1d);

    have_fftw32 = (1
                   && _fftw32_malloc
                   && _fftw32_free
                   && _fftw32_destroy_plan
                   && _fftw32_execute
                   && _fftw32_plan_dft_1d
                   && _fftw32_plan_dft_c2r_1d
                   && _fftw32_plan_dft_r2c_1d
                   );

    register_fftw64(malloc);
    register_fftw64(free);
    register_fftw64(destroy_plan);
    register_fftw64(execute);
    register_fftw64(plan_dft_1d);
    register_fftw64(plan_dft_c2r_1d);
    register_fftw64(plan_dft_r2c_1d);

    have_fftw64 = (1
                   && _fftw64_malloc
                   && _fftw64_free
                   && _fftw64_destroy_plan
                   && _fftw64_execute
                   && _fftw64_plan_dft_1d
                   && _fftw64_plan_dft_c2r_1d
                   && _fftw64_plan_dft_r2c_1d
                   );

    have_fftw = have_fftw32 || have_fftw64;

    if(have_fftw32 && have_fftw64) {
      post("iemmatrix: using FFTW/FFTWF/Mayer for Fourier transforms");
    } else if (have_fftw64) {
      post("iemmatrix: using FFTW/Mayer for Fourier transforms");
    } else if (have_fftw32) {
      post("iemmatrix: using FFTWF/Mayer for Fourier transforms");
    } else {
      post("iemmatrix: using built in Mayer/Ooura for Fourier transforms");
    }
  }

  if(have_fftw) {
    return FFTW;
  }

  return MAYER;
}
