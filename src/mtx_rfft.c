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

#include <stdlib.h>

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

#ifndef HAVE_FFTW
typedef struct _fftw_plan_private_tag *fftw_plan;
typedef double fftw_complex[2];
#define FFTW_ESTIMATE (1U << 6)
#endif
typedef void(*t_fftw_destroy_plan)(fftw_plan);
typedef void(*t_fftw_execute)(const fftw_plan);
typedef fftw_plan(*t_fftw_plan_dft_r2c_1d)(int, double*, fftw_complex*, unsigned);

static t_fftw_destroy_plan my_destroy_plan = 0;
static t_fftw_execute my_execute = 0;
static t_fftw_plan_dft_r2c_1d my_plan_dft_r2c_1d = 0;

static int have_fftw = 0;

static t_class *mtx_rfft_class;

enum ComplexPart { REALPART=0,  IMAGPART=1};

typedef struct _MTXRfft_ MTXRfft;
struct _MTXRfft_ {
  t_object x_obj;
  int size;
  int size2;

  /* fftw */
  int fftn;
  int rows;
  fftw_plan *fftplan;
  fftw_complex *f_out;
  double *f_in;

  /* mayer */
  t_float *f_re;
  t_float *f_im;


  t_outlet *list_re_out;
  t_outlet *list_im_out;

  t_atom *list_re;
  t_atom *list_im;
};

static void deleteMTXRfft (MTXRfft *x)
{
  if(have_fftw) {
    int n;
    if (x->fftplan) {
      for (n=0; n<x->rows; n++) {
        my_destroy_plan(x->fftplan[n]);
      }
      free(x->fftplan);
    }
    if (x->f_out) {
      free(x->f_out);
    }
    if (x->f_in) {
      free(x->f_in);
    }
  } else {
    if (x->f_re) {
      free (x->f_re);
    }
    if (x->f_im) {
      free (x->f_im);
    }
  }
  if (x->list_re) {
    free (x->list_re);
  }
  if (x->list_im) {
    free (x->list_im);
  }
}

static void *newMTXRfft (t_symbol *s, int argc, t_atom *argv)
{
  MTXRfft *x = (MTXRfft *) pd_new (mtx_rfft_class);
  (void)argc; (void)argv; /* unused */

  x->list_re_out = outlet_new (&x->x_obj, gensym("matrix"));
  x->list_im_out = outlet_new (&x->x_obj, gensym("matrix"));
  if (!have_fftw) {
    static int warn_fftw = 1;
    if(warn_fftw)
#ifdef HAVE_FFTW
      pd_error(x, "[%s] couldn't find (recent enough) FFTW", s->s_name);
#else
      pd_error(x, "[%s] compiled without FFTW", s->s_name);
#endif
    warn_fftw = 0;
  }
  return ((void *) x);
}

static void mTXRfftBang (MTXRfft *x)
{
  if (x->list_im) {
    outlet_anything(x->list_im_out, gensym("matrix"), x->size2, x->list_im);
    outlet_anything(x->list_re_out, gensym("matrix"), x->size2, x->list_re);
  }
}

static void fftRestoreImag (int n, t_float *re, t_float *im)
{
  t_float *im2;
  n >>= 1;
  *im=0;
  re += n;
  im += n;
  im2 = im;
  *im=0;
  while (--n) {
    *--im = -*++re;
    *++im2 = 0;
    *re = 0;
  }
}

static void writeFFTWComplexPartIntoList (int n, t_atom *l,
    fftw_complex *c, enum ComplexPart p)
{
  t_float f;
  while (n--) {
    f=(t_float)c[n][p];
    SETFLOAT (l+n, f);
  }
}
static void readDoubleFromList (int n, t_atom *l, double *f)
{
  while (n--) {
    *f++ = (double)atom_getfloat (l++);
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
  int size2 = columns_re * rows +
              2; /* +2 since the list also contains matrix row+col */
  int fft_count;
  t_atom *list_re = x->list_re;
  t_atom *list_im = x->list_im;

  fftw_complex *f_out = x->f_out;
  double *f_in = x->f_in;

  t_float *f_re = x->f_re;
  t_float *f_im = x->f_im;

  const int use_fftw = have_fftw;
  (void)s; /* unused */



  /* fftsize check */
  if (!size) {
    pd_error(x, "[mtx_rfft]: invalid dimensions");
  } else if (in_size<size) {
    pd_error(x, "[mtx_rfft]: sparse matrix not yet supported: use \"mtx_check\"");
  } else if (columns < 4) {
    pd_error(x, "[mtx_rfft]: matrix must have at least 4 columns");
  } else if (columns == (1 << ilog2(columns))) {
    /* ok, do the FFT! */

    /* memory things */
    if(use_fftw) {
      if ((x->rows!=rows)||(columns!=x->fftn)) {
        f_out=(fftw_complex*)realloc(f_out, sizeof(fftw_complex)*(size2-2));
        f_in=(double*)realloc(f_in, sizeof(double)*size);
        x->f_in = f_in;
        x->f_out = f_out;
        for (fft_count=0; fft_count<x->rows; fft_count++) {
          my_destroy_plan(x->fftplan[fft_count]);
        }
        x->fftplan = (fftw_plan*)realloc(x->fftplan, sizeof(fftw_plan)*rows);
        for (fft_count=0; fft_count<rows;
             fft_count++, f_in+=columns, f_out+=columns_re) {
          x->fftplan[fft_count] = my_plan_dft_r2c_1d (columns,f_in,f_out, FFTW_ESTIMATE);
        }
        x->fftn=columns;
        x->rows=rows;
        f_in=x->f_in;
        f_out=x->f_out;
      }
    } else {
      f_re=(t_float*)realloc(f_re, sizeof(t_float)*size);
      f_im=(t_float*)realloc(f_im, sizeof(t_float)*size);
      x->f_re = f_re;
      x->f_im = f_im;
    }

    list_re=(t_atom*)realloc(list_re, sizeof(t_atom)*size2);
    list_im=(t_atom*)realloc(list_im, sizeof(t_atom)*size2);

    x->size = size;
    x->size2 = size2;
    x->list_im = list_im;
    x->list_re = list_re;

    /* main part */
    if(use_fftw) {
      readDoubleFromList (size, argv, f_in);
    } else {
      iemmatrix_list2floats(f_re, argv, size);
    }

    list_re += 2;
    list_im += 2;
    for (fft_count=0; fft_count<rows; fft_count++) {
      if(use_fftw) {
        my_execute(x->fftplan[fft_count]);
        writeFFTWComplexPartIntoList(columns_re,list_re,f_out,REALPART);
        writeFFTWComplexPartIntoList(columns_re,list_im,f_out,IMAGPART);
        f_out+=columns_re;
      } else {
        mayer_realfft (columns, f_re);
        fftRestoreImag (columns, f_re, f_im);
        iemmatrix_floats2list(list_re, f_re, columns_re);
        iemmatrix_floats2list(list_im, f_im, columns_re);
        f_im += columns;
        f_re += columns;
      }
      list_re += columns_re;
      list_im += columns_re;
    }

    list_re = x->list_re;
    list_im = x->list_im;

    SETSYMBOL(list_re, gensym("matrix"));
    SETSYMBOL(list_im, gensym("matrix"));
    SETFLOAT(list_re, rows);
    SETFLOAT(list_im, rows);
    SETFLOAT(list_re+1, columns_re);
    SETFLOAT(list_im+1, columns_re);
    outlet_anything(x->list_im_out, gensym("matrix"),
                    x->size2, list_im);
    outlet_anything(x->list_re_out, gensym("matrix"),
                    x->size2, list_re);
  } else {
    pd_error(x, "[mtx_rfft]: rowvector size no power of 2!");
  }

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

#ifdef HAVE_FFTW
  my_destroy_plan = iemmatrix_get_stub("fftw_destroy_plan", mtx_rfft_class);
  my_execute = iemmatrix_get_stub("fftw_execute", mtx_rfft_class);
  my_plan_dft_r2c_1d = iemmatrix_get_stub("fftw_plan_dft_r2c_1d", mtx_rfft_class);
#endif
  have_fftw = (1
               && my_destroy_plan
               && my_execute
               && my_plan_dft_r2c_1d
               );}

void iemtx_rfft_setup(void)
{
  mtx_rfft_setup();
}
