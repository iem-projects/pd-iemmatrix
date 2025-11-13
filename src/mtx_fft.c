/*
 *  iemmatrix
 *
 *  objects for manipulating simple matrices
 *  mostly referring to matlab/octave matrix functions
 *
 * Copyright (c) 2005, Franz Zotter
 * IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */

#include "iemmatrix.h"
#include <stdlib.h>

static t_class *mtx_fft_class;

typedef struct _MtxFFT_ MtxFFT;
struct _MtxFFT_ {
  t_object x_obj;

  t_float *f_re;
  t_float *f_im;

  t_outlet *msg_re_out;
  t_outlet *msg_im_out;

  t_atom *msg_re;
  t_atom *msg_im;
  int msg_size;
};

static void deleteMtxFFT (MtxFFT *x)
{
  if (x->f_re)   free (x->f_re);
  if (x->f_im)   free (x->f_im);
  if (x->msg_re) free (x->msg_re);
  if (x->msg_im) free (x->msg_im);
}

static void *newMtxFFT ()
{
  MtxFFT *x = (MtxFFT *) pd_new (mtx_fft_class);
  inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("matrix"),gensym(""));
  x->msg_re_out = outlet_new (&x->x_obj, gensym("matrix"));
  x->msg_im_out = outlet_new (&x->x_obj, gensym("matrix"));

  x->msg_size=0;
  x->msg_re=x->msg_im=0;

  x->f_re=x->f_im=0;

  return ((void *) x);
}

static void mtxFFTBang (MtxFFT *x)
{
  if (x->msg_im) {
    outlet_anything(x->msg_im_out, gensym("matrix"), x->msg_size, x->msg_im);
    outlet_anything(x->msg_re_out, gensym("matrix"), x->msg_size, x->msg_re);
  }
}

static void mtxFFTMatrixCold (MtxFFT *x, t_symbol *s,
                              int argc, t_atom *argv)
{
  int rows, columns, size;
  t_atom *msg_re = x->msg_re;
  t_atom *msg_im = x->msg_im;
  t_float *f_re = x->f_re;
  t_float *f_im = x->f_im;
  (void)s; /* unused */

  /* fftsize check */
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  rows = atom_getint (argv++);
  columns = atom_getint (argv++);
  size = rows * columns;

  if (columns < 4) {
    pd_error(x, "[mtx_fft]: matrix must have at least 4 columns");
    return;
  } else if (columns != (1 << ilog2(columns))) {
    pd_error(x, "[mtx_fft]: rowvector size no power of 2!");
    return;
  }

  /* ok, prepare real-part of FFT! */


  /* memory things */
  f_re=(t_float*)realloc(f_re, sizeof(t_float)*size);
  f_im=(t_float*)realloc(f_im, sizeof(t_float)*size);
  msg_re=(t_atom*)realloc(msg_re, sizeof(t_atom)*(size+2));
  msg_im=(t_atom*)realloc(msg_im, sizeof(t_atom)*(size+2));

  x->msg_size = size;
  x->msg_im = msg_im;
  x->msg_re = msg_re;
  x->f_re = f_re;
  x->f_im = f_im;

  /* main part */
  iemmatrix_list2floats(f_im, argv, size);
}


static void mtxFFTMatrixHot (MtxFFT *x, t_symbol *s,
                             int argc, t_atom *argv)
{
  int rows, columns, msg_size;
  int fft_count;
  t_atom *msg_re = x->msg_re;
  t_atom *msg_im = x->msg_im;
  t_float *f_re = x->f_re;
  t_float *f_im = x->f_im;
  (void)s; /* unused */

  /* fftsize check */
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  rows = atom_getint (argv++);
  columns = atom_getint (argv++);
  msg_size = rows * columns;

  if (msg_size != x->msg_size) {
    pd_error(x, "[mtx_fft]: left matrix has other dimensions than right matrix");
    return;
  } else if (columns < 4) {
    pd_error(x, "[mtx_fft]: matrix must have at least 4 columns");
    return;
  } else if (columns != (1 << ilog2(columns))) {
    pd_error(x, "[mtx_fft]: rowvector size no power of 2!");
    return;
  }
  /* ok, do the FFT! */

  /* main part */
  iemmatrix_list2floats(f_re, argv, msg_size);

  fft_count = rows;
  msg_re += 2;
  msg_im += 2;
  while (fft_count--) {
    mayer_fft (columns, f_re, f_im);
    iemmatrix_floats2list(msg_re, f_re, columns);
    iemmatrix_floats2list(msg_im, f_im, columns);
    f_im += columns;
    f_re += columns;
    msg_re += columns;
    msg_im += columns;
  }

  msg_re = x->msg_re;
  msg_im = x->msg_im;

  SETSYMBOL(msg_re, gensym("matrix"));
  SETSYMBOL(msg_im, gensym("matrix"));
  SETFLOAT(msg_re, rows);
  SETFLOAT(msg_im, rows);
  SETFLOAT(msg_re+1, columns);
  SETFLOAT(msg_im+1, columns);
  outlet_anything(x->msg_im_out, gensym("matrix"),
                  x->msg_size+2, msg_im);
  outlet_anything(x->msg_re_out, gensym("matrix"),
                  x->msg_size+2, msg_re);
}

void mtx_fft_setup (void)
{
  mtx_fft_class = class_new
                  (gensym("mtx_fft"),
                   (t_newmethod) newMtxFFT,
                   (t_method) deleteMtxFFT,
                   sizeof (MtxFFT),
                   CLASS_DEFAULT, 0);
  class_addbang (mtx_fft_class, (t_method) mtxFFTBang);
  class_addmethod (mtx_fft_class, (t_method) mtxFFTMatrixHot,
                   gensym("matrix"), A_GIMME,0);
  class_addmethod (mtx_fft_class, (t_method) mtxFFTMatrixCold, gensym(""),
                   A_GIMME,0);
}

void iemtx_fft_setup(void)
{
  mtx_fft_setup();
}
