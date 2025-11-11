/* ************************************* */
/* iemmatrix                             */
/* ************************************* */
/*  objects for simple matrix operations */
/* ************************************* */

/*
 * Copyright (c) IOhannes m zmölnig (forum::für::umläute), IEM KUG, Graz, Austria; 2025
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 * there are ABSOLUTELY NO WARRANTIES for anything
 */
#ifndef _iemmatrix_fft_h_
#define _iemmatrix_fft_h_

#include "iemmatrix.h"

typedef struct _complexfloat {
  float re;
  float im;
} t_complexfloat;
typedef struct _complexdouble {
  double re;
  double im;
} t_complexdouble;
#if PD_FLOATSIZE == 64
typedef struct _complexdouble t_complex;
#else
typedef struct _complexfloat t_complex;
#endif


typedef struct iemfft_plan *t_iemfft_plan;


typedef enum {
  NONE = 0,
  MAYER,
  FFTW,
  INVALID
} t_iemfft_backend;

typedef enum {
  DEFAULT = 0,
  PREFER_SINGLE =  1 << 0, /* prefer single precision processing */
  PREFER_DOUBLE =  1 << 1, /* prefer double precision processing */
  PRESERVE_INPUT = 1 << 2, /* 'in' MUST not be modified during execute */
} t_iemfft_flag;

/* initialize FFT backends */
t_iemfft_backend iemfft_init(t_class*c);


void*iemfft_malloc(size_t size);
void iemfft_free(void*data);

t_iemfft_plan iemfft_plan_fft_1d(int n0, t_complex*in, t_complex* out, t_iemfft_flag flags);
t_iemfft_plan iemfft_plan_ifft_1d(int n0, t_complex*in, t_complex* out, t_iemfft_flag flags);
/* real-valued FFTs */
t_iemfft_plan iemfft_plan_rfft_1d(int n0, t_float*in, t_complex*out, t_iemfft_flag flags);
t_iemfft_plan iemfft_plan_rifft_1d(int n0, t_complex*in, t_float*out, t_iemfft_flag flags);

void iemfft_execute(const t_iemfft_plan plan);
void iemfft_destroy_plan(t_iemfft_plan plan);

#endif /* _iemmatrix_fft_h_ */
