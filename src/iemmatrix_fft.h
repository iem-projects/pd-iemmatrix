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

typedef struct _complex {
  t_float re;
  t_float im;
} t_complex;


void*iemmatrix_fft_malloc(size_t size);
void iemmatrix_fft_free(void*data);

struct iemmatrix_fft_plan*iemmatrix_fft_plan_1d(int n0, t_complex*in, t_complex* out, int direction, int flags);

/* real-valued FFTs */
struct iemmatrix_fft_plan*iemmatrix_rfft_plan_1d(int n0, t_float*in, t_complex*out, unsigned flags);
struct iemmatrix_fft_plan*iemmatrix_rifft_plan_1d(int n0, t_complex*in, t_float*out, unsigned flags);

void iemmatrix_fft_execute(const struct iemmatrix_fft_plan*plan);
void iemmatrix_fft_destroy_plan(struct iemmatrix_fft_plan*plan);

#endif /* _iemmatrix_fft_h_ */
