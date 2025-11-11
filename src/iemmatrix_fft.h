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
void*iemmatrix_fft_malloc(size_t size);
void iemmatrix_fft_free(void*data);

struct iemmatrix_fft_plan*iemmatrix_fft_plan_1d(int n0, t_complex*in, t_complex* out);
struct iemmatrix_fft_plan*iemmatrix_ifft_plan_1d(int n0, t_complex*in, t_complex* out);

/* real-valued FFTs */
struct iemmatrix_fft_plan*iemmatrix_rfft_plan_1d(int n0, t_float*in, t_complex*out);
struct iemmatrix_fft_plan*iemmatrix_rifft_plan_1d(int n0, t_complex*in, t_float*out);

void iemmatrix_fft_execute(const struct iemmatrix_fft_plan*plan);
void iemmatrix_fft_destroy_plan(struct iemmatrix_fft_plan*plan);

#endif /* _iemmatrix_fft_h_ */
