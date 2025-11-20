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

static t_class *mtx_rifft_class;

typedef struct _MTXRifft_ {
  t_object x_obj;
  t_outlet *list_re_out;

  int columns, rows;
  int size; // ???

  t_iemfft_plan*plan;
  t_complex *c_in;
  t_float *f_out;

  int list_size;
  t_atom *list_re;
} MTXRifft;


/*--------------inverse real fft */

static void multiplyVector (int n, t_float *f, t_float fac)
{
  while (n--) {
    *f++ *= fac;
  }
}

static void list2real(int n, t_atom*re, t_complex*c) {
  while (n--) {
    c->re = atom_getfloat(re++);
    c++;
  }
}
static void list2imag(int n, t_atom*im, t_complex*c) {
  while (n--) {
    c->im = atom_getfloat(im++);
    c++;
  }
}

static void *newMTXRifft ()
{
  MTXRifft *x = (MTXRifft *) pd_new (mtx_rifft_class);
  inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("matrix"),gensym(""));
  x->list_re_out = outlet_new (&x->x_obj, gensym("matrix"));
  return ((void *) x);
}


static void mTXRifftMatrixCold (MTXRifft *x, t_symbol *s,
                                int argc, t_atom *argv)
{
  int rows = atom_getint (argv++);
  int columns_re = atom_getint (argv++);
  int in_size = argc-2;
  int columns = (columns_re-1)<<1;
  int list_size = columns_re * rows;
  int size = rows * columns;
  t_atom *list_re = x->list_re;

  (void)s; /* unused */

  /* ifftsize check */
  if (columns_re < 3) {
    pd_error(x, "[mtx_rifft]: matrix must have at least 3 columns");
    return;
  } else if (!size) {
    pd_error(x, "[mtx_rifft]: invalid dimensions");
    return;
  } else if (in_size < list_size) {
    pd_error(x, "[mtx_rifft]: sparse matrix not yet supported: use [mtx_check]");
    return;
  } else if (columns<4) {
    pd_error(x, "[mtx_rifft]: too small matrices");
    return;
  } else if (columns == (1 << ilog2(columns))) {
    /* OK */
  } else {
    pd_error(x, "[mtx_rifft]: rowvector 2*(size+1) no power of 2!");
    return;
  }

  /* memory things */
  if ((x->rows!=rows)||(columns!=x->columns)) {
    /* size changed, so re-allocate */
    t_complex* c_in = x->c_in  = (t_complex*)realloc(x->c_in ,sizeof(*x->c_in )*list_size);
    t_float*  f_out = x->f_out = (t_float*  )realloc(x->f_out,sizeof(*x->f_out)*size);

    for (int r=0; r<x->rows; r++) {
      iemfft_destroy_plan(x->plan[r]);
    }
    x->rows = rows;
    x->columns = columns;
    x->plan=(t_iemfft_plan*)realloc(x->plan,sizeof(*x->plan)*x->rows);
    for (int r=0; r<x->rows; r++) {
      x->plan[r] = iemfft_plan_rifft_1d(columns, c_in + r * columns_re, f_out + r * columns, PREFER_DOUBLE);
    }

    x->list_re = list_re =(t_atom*)realloc(list_re, sizeof(*list_re)*(size+2));
  }

  x->list_size = list_size;
  x->rows = rows;
  x->columns = columns;

  /* main part: reading imaginary part */
  t_complex*c_in = x->c_in;
  for (int r=0; r<rows; r++) {
    list2imag(columns_re, argv, c_in);
    c_in += columns_re;
    argv += columns_re;
  }
}

static void mTXRifftMatrixHot (MTXRifft *x, t_symbol *s,
                               int argc, t_atom *argv)
{
  int rows = atom_getint (argv++);
  int columns_re = atom_getint (argv++);
  int columns = x->columns;
  int size = x->rows * x->columns;
  int in_size = argc-2;
  int list_size = x->list_size;
  t_complex *c_in = x->c_in;

  t_float renorm_fac = 1. / (t_float)columns;
  (void)s; /* unused */

  /* ifftsize check */
  if ((rows != x->rows) || (columns != ((columns_re-1)<<1))) {
    pd_error(x, "[mtx_rifft]: matrix dimensions do not match: expected %dx%d, got %dx%d[%d]", x->rows, x->columns, rows, ((columns_re-1)<<1), columns_re);
    return;
  } else if (in_size<list_size) {
    pd_error(x, "[mtx_rifft]: sparse matrix not yet supported: use [mtx_check]");
    return;
  } else if (!x->list_size) {
    pd_error(x, "[mtx_rifft]: invalid right side matrix");
    return;
  }

  /* main part */
  for (int r=0; r<rows; r++) {
    list2real(columns_re, argv, c_in);
    iemfft_execute(x->plan[r]);

    c_in+=columns_re;
    argv += columns_re;
  }

  multiplyVector(size, x->f_out, renorm_fac);
  iemmatrix_floats2list(x->list_re + 2, x->f_out, size);

  SETFLOAT(x->list_re, rows);
  SETFLOAT(x->list_re+1, x->columns);


  outlet_anything(x->list_re_out, gensym("matrix"), size+2, x->list_re);
  x->size = size;
}

static void mTXRifftBang (MTXRifft *x)
{
  if (x->list_re)
    outlet_anything(x->list_re_out, gensym("matrix"),
                    x->size+2, x->list_re);
}


static void deleteMTXRifft (MTXRifft *x)
{
  if (x->plan) {
    int n;
    for (n=0; n<x->rows; n++) {
      iemfft_destroy_plan(x->plan[n]);
    }
    free(x->plan);
  }

  free(x->c_in);
  free(x->f_out);

  free(x->list_re);
}

void mtx_rifft_setup (void)
{
  mtx_rifft_class = class_new
                    (gensym("mtx_rifft"),
                     (t_newmethod) newMTXRifft,
                     (t_method) deleteMTXRifft,
                     sizeof (MTXRifft),
                     CLASS_DEFAULT, 0);
  class_addbang (mtx_rifft_class, (t_method) mTXRifftBang);
  class_addmethod (mtx_rifft_class, (t_method) mTXRifftMatrixHot,
                   gensym("matrix"), A_GIMME,0);
  class_addmethod (mtx_rifft_class, (t_method) mTXRifftMatrixCold,
                   gensym(""), A_GIMME,0);

  iemfft_init(mtx_rifft_class);
}

void iemtx_rifft_setup(void)
{
  mtx_rifft_setup();
}
