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
#else
# define x_in f_in
# define x_out f_out
#endif

typedef struct iemmatrix_fft_plan {
  int n0;
  t_complex*c_in, *c_out;
  float*f_in, *f_out;
  double*d_in, *d_out;
  int inverse;
} t_iemmatrix_fft_plan;

void*mayer_malloc(size_t size) {
  return malloc(size);
}
void mayer_free(void*data) {
  free(data);
}

void mayer_destroy_plan(t_iemmatrix_fft_plan*plan) {
  freebytes(plan, sizeof(*plan));
}

t_iemmatrix_fft_plan*iemmatrix_rifft_plan_1d(int n0, t_complex*in, t_float*out, unsigned flags) {
  t_iemmatrix_fft_plan*plan = calloc(1, sizeof(t_iemmatrix_fft_plan));
  plan->n0 = n0;
  plan->c_in = in;
  plan->x_out = out;
  plan->inverse = 1;
  (void)flags;
  return plan;
}
t_iemmatrix_fft_plan*iemmatrix_rfft_plan_1d(int n0, t_float*in, t_complex*out, unsigned flags) {
  t_iemmatrix_fft_plan*plan = calloc(1, sizeof(t_iemmatrix_fft_plan));
  plan->n0 = n0;
  plan->x_in = in;
  plan->c_out = out;
  plan->inverse = 0;
  (void)flags;
  return plan;
}

t_iemmatrix_fft_plan*iemmatrix_fft_plan_1d(int n0, t_complex*in, t_complex* out, int direction, int flags) {
  t_iemmatrix_fft_plan*plan = calloc(1, sizeof(t_iemmatrix_fft_plan));
  plan->n0 = n0;
  plan->c_in = in;
  plan->c_out = out;
  plan->inverse = (1 == direction);
  (void)flags;
  return plan;
}


void mayer_execute(const t_iemmatrix_fft_plan*plan) {
  #warning total nonsense
  if(plan->x_in) {
    if(plan->inverse) {
      /* ifft */
      mayer_ifft(plan->n0, (t_float*)0, (t_float*)0);
    } else {
      /* fft */
      mayer_fft(plan->n0, (t_float*)0, (t_float*)0);
    }
  } else {
    if(plan->inverse) {
      /* rifft */
      mayer_realifft(plan->n0, plan->x_in);
    } else {
      /* rfft */
      mayer_realfft(plan->n0, plan->x_in);
    }
  }
}
