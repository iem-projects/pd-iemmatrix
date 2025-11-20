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
#include "iemmatrix_stub.h"

#include "iemmatrix_fft.h"

#include <stdlib.h>

static t_class *mtx_rfft_class;

typedef struct _MTXRfft_ MTXRfft;
struct _MTXRfft_ {
  t_object x_obj;

  int columns, rows;

  t_iemfft_plan plan;
  t_float *f_in;
  t_complex *c_out;

  t_outlet *list_re_out;
  t_outlet *list_im_out;

  t_atom *list_re;
  t_atom *list_im;
  int list_size;
};

static void deleteMTXRfft (MTXRfft *x)
{
  if(x->plan) {
      iemfft_destroy_plan(x->plan);
  }
  if (x->f_in)   free(x->f_in);
  if (x->c_out)  free(x->c_out);
  if (x->list_re)free(x->list_re);
  if (x->list_im)free(x->list_im);
}

static void *newMTXRfft (t_symbol *s, int argc, t_atom *argv)
{
  MTXRfft *x = (MTXRfft *) pd_new (mtx_rfft_class);
  (void)s; (void)argc; (void)argv; /* unused */

  x->list_re_out = outlet_new (&x->x_obj, gensym("matrix"));
  x->list_im_out = outlet_new (&x->x_obj, gensym("matrix"));
  return ((void *) x);
}

static void mTXRfftBang (MTXRfft *x)
{
  if (x->list_im) {
    outlet_anything(x->list_im_out, gensym("matrix"), x->list_size, x->list_im);
    outlet_anything(x->list_re_out, gensym("matrix"), x->list_size, x->list_re);
  }
}

static void complex2list (int n, t_complex*c,
                          t_atom*re, t_atom*im)
{
  while (n--) {
    SETFLOAT (re, c->re);
    SETFLOAT (im, c->im);
    c++;
    re++; im++;
  }
}

static void mTXRfftMatrix (MTXRfft *x, t_symbol *s,
                           int argc, t_atom *argv)
{
  int rows = atom_getint (argv++);
  int columns = atom_getint (argv++);
  int columns_re = (columns>>1)
                   +1; /* N/2+1 samples needed for real part of realfft */
  int size = rows * columns;
  int in_size = argc-2;
  int list_size = columns_re * rows +
              2; /* +2 since the list also contains matrix row+col */

  t_float *f_in = x->f_in;
  t_complex *c_out = x->c_out;

  t_atom *list_re = x->list_re;
  t_atom *list_im = x->list_im;

  (void)s; /* unused */



  /* fftsize check */
  if (!size) {
    pd_error(x, "[mtx_rfft]: invalid dimensions");
    return;
  } else if (in_size<size) {
    pd_error(x, "[mtx_rfft]: sparse matrix not yet supported: use \"mtx_check\"");
    return;
  } else if (columns < 4) {
    pd_error(x, "[mtx_rfft]: matrix must have at least 4 columns");
    return;
  } else if (columns == (1 << ilog2(columns))) {
    /* ok, do the FFT! */
  } else {
    pd_error(x, "[mtx_rfft]: rowvector size no power of 2!");
    return;
  }

  /* memory things */
  if ((x->rows != rows)||(x->columns != columns)) {
    /* size changed, so re-allocate */
    x->f_in = f_in  = (t_float*  )realloc(f_in , sizeof(*f_in ) * columns);
    x->c_out= c_out = (t_complex*)realloc(c_out, sizeof(*c_out) * columns_re);

    if(x->plan)
      iemfft_destroy_plan(x->plan);

    x->rows = rows;
    x->columns = columns;

    x->plan = iemfft_plan_rfft_1d(columns, f_in, c_out, PREFER_DOUBLE);
  }

  list_re=(t_atom*)realloc(list_re, sizeof(t_atom)*list_size);
  list_im=(t_atom*)realloc(list_im, sizeof(t_atom)*list_size);

  x->list_size = list_size;
  x->list_im = list_im;
  x->list_re = list_re;

  for(int r=0; r<rows; r++) {
    int offset = columns_re * r;
    iemmatrix_list2floats(f_in, argv+r*columns, columns);
    iemfft_execute(x->plan);
    complex2list(columns_re, c_out, list_re + offset + 2, list_im + offset + 2);
  }

  SETSYMBOL(list_re, gensym("matrix"));
  SETSYMBOL(list_im, gensym("matrix"));
  SETFLOAT(list_re, rows);
  SETFLOAT(list_im, rows);
  SETFLOAT(list_re+1, columns_re);
  SETFLOAT(list_im+1, columns_re);
  outlet_anything(x->list_im_out, gensym("matrix"),
                  x->list_size, list_im);
  outlet_anything(x->list_re_out, gensym("matrix"),
                  x->list_size, list_re);
}

void mtx_rfft_setup (void)
{
  mtx_rfft_class = class_new
                   (gensym("mtx_rfft"),
                    (t_newmethod) newMTXRfft,
                    (t_method) deleteMTXRfft,
                    sizeof (MTXRfft),
                    CLASS_DEFAULT, A_GIMME, 0);
  class_addbang (mtx_rfft_class, (t_method) mTXRfftBang);
  class_addmethod (mtx_rfft_class, (t_method) mTXRfftMatrix,
                   gensym("matrix"), A_GIMME,0);

  iemfft_init(mtx_rfft_class);
}

void iemtx_rfft_setup(void)
{
  mtx_rfft_setup();
}
