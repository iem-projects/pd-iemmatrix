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
  int size;

  t_float *f_re;
  t_float *f_im;

  t_outlet *list_re_out;
  t_outlet *list_im_out;

  t_atom *list_re;
  t_atom *list_im;
};

static void deleteMtxFFT (MtxFFT *x)
{
  if (x->f_re) {
    free (x->f_re);
  }
  if (x->f_im) {
    free (x->f_im);
  }
  if (x->list_re) {
    free (x->list_re);
  }
  if (x->list_im) {
    free (x->list_im);
  }
}

static void *newMtxFFT ()
{
  MtxFFT *x = (MtxFFT *) pd_new (mtx_fft_class);
  inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("matrix"),gensym(""));
  x->list_re_out = outlet_new (&x->x_obj, gensym("matrix"));
  x->list_im_out = outlet_new (&x->x_obj, gensym("matrix"));

  x->size=0;
  x->f_re=x->f_im=0;
  x->list_re=x->list_im=0;

  return ((void *) x);
}

static void mtxFFTBang (MtxFFT *x)
{
  if (x->list_im) {
    outlet_anything(x->list_im_out, gensym("matrix"), x->size, x->list_im);
    outlet_anything(x->list_re_out, gensym("matrix"), x->size, x->list_re);
  }
}

static void writeFloatIntoList (int n, t_atom *l, t_float *f)
{
  for (; n--; f++, l++) {
    SETFLOAT (l, *f);
  }
}
static void readFloatFromList (int n, t_atom *l, t_float *f)
{
  while (n--) {
    *f++ = atom_getfloat (l++);
  }
}

static void mtxFFTMatrixCold (MtxFFT *x, t_symbol *s,
                              int argc, t_atom *argv)
{
  int rows, columns, size;
  t_atom *list_re = x->list_re;
  t_atom *list_im = x->list_im;
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
  } else if (columns == (1 << ilog2(columns))) {
    /* ok, prepare real-part of FFT! */

    /* memory things */
    f_re=(t_float*)realloc(f_re, sizeof(t_float)*size);
    f_im=(t_float*)realloc(f_im, sizeof(t_float)*size);
    list_re=(t_atom*)realloc(list_re, sizeof(t_atom)*(size+2));
    list_im=(t_atom*)realloc(list_im, sizeof(t_atom)*(size+2));

    x->size = size;
    x->list_im = list_im;
    x->list_re = list_re;
    x->f_re = f_re;
    x->f_im = f_im;

    /* main part */
    iemmatrix_list2floats(f_im, argv, size);

  } else {
    pd_error(x, "[mtx_fft]: rowvector size no power of 2!");
  }
}


static void mtxFFTMatrixHot (MtxFFT *x, t_symbol *s,
                             int argc, t_atom *argv)
{
  int rows, columns, size;
  int fft_count;
  t_atom *list_re = x->list_re;
  t_atom *list_im = x->list_im;
  t_float *f_re = x->f_re;
  t_float *f_im = x->f_im;
  (void)s; /* unused */

  /* fftsize check */
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  rows = atom_getint (argv++);
  columns = atom_getint (argv++);
  size = rows * columns;

  if (size != x->size) {
    pd_error(x, "[mtx_fft]: left matrix has other dimensions than right matrix");
  } else if (columns < 4) {
    pd_error(x, "[mtx_fft]: matrix must have at least 4 columns");
  } else if (columns == (1 << ilog2(columns))) {
    /* ok, do the FFT! */

    /* main part */
    iemmatrix_list2floats(f_re, argv, size);

    fft_count = rows;
    list_re += 2;
    list_im += 2;
    while (fft_count--) {
      mayer_fft (columns, f_re, f_im);
      iemmatrix_floats2list(list_re, f_re, columns);
      iemmatrix_floats2list(list_im, f_im, columns);
      f_im += columns;
      f_re += columns;
      list_re += columns;
      list_im += columns;
    }

    list_re = x->list_re;
    list_im = x->list_im;

    SETSYMBOL(list_re, gensym("matrix"));
    SETSYMBOL(list_im, gensym("matrix"));
    SETFLOAT(list_re, rows);
    SETFLOAT(list_im, rows);
    SETFLOAT(list_re+1, columns);
    SETFLOAT(list_im+1, columns);
    outlet_anything(x->list_im_out, gensym("matrix"),
                    x->size+2, list_im);
    outlet_anything(x->list_re_out, gensym("matrix"),
                    x->size+2, list_re);
  } else {
    pd_error(x, "[mtx_fft]: rowvector size no power of 2!");
  }

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
