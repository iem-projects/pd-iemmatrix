/*
 *  iemmatrix
 *
 *  objects for manipulating simple matrices
 *  mostly referring to matlab/octave matrix functions
 *  this functions depends on the GNU scientific library
 *
 * Copyright (c) 2009, Franz Zotter
 * IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */

#include "iemmatrix.h"
#include "iemmatrix_stub.h"
#include <stdlib.h>

#ifdef HAVE_LIBGSL
#include <gsl/gsl_linalg.h>
#else
# include "stub/gsl.h"
#endif

IEMMATRIX_DECLARE_ALLOCFREE_STUB(my_vector);
IEMMATRIX_DECLARE_ALLOCFREE2_STUB(my_matrix);
typedef int(*t_linalg_QR_decomp)(gsl_matrix*, gsl_vector*);
typedef int(*t_linalg_QR_unpack)(gsl_matrix*, gsl_vector*, gsl_matrix*, gsl_matrix*);

static t_linalg_QR_decomp my_linalg_QR_decomp = 0;
static t_linalg_QR_unpack my_linalg_QR_unpack = 0;

static int have_gsl = 0;


static t_class *mtx_qr_class;

typedef struct _MTXQr_ MTXQr;
struct _MTXQr_ {
  t_object x_obj;

  gsl_matrix *a;
  gsl_vector *tau;
  gsl_matrix *q;
  gsl_matrix *r;

  t_outlet *list_q_out;
  t_outlet *list_r_out;
  t_atom *list_q;
  t_atom *list_r;
  int rows;
  int columns;
};

static void allocMTXqr (MTXQr *x)
{
  x->a=(gsl_matrix*)my_matrix_alloc(x->rows,x->columns);
  x->tau=(gsl_vector*)my_vector_alloc(
           ((x->columns<x->rows)?x->columns:x->rows));
  x->q=(gsl_matrix*)my_matrix_alloc(x->rows,x->rows);
  x->r=(gsl_matrix*)my_matrix_alloc(x->rows,x->columns);
  x->list_q=(t_atom*)calloc(x->rows*x->rows+2, sizeof(t_atom));
  x->list_r=(t_atom*)calloc(x->rows*x->columns+2, sizeof(t_atom));
}

static void deleteMTXqr (MTXQr *x)
{
  if (x->list_q!=0) {
    free(x->list_q);
  }
  if (x->list_r!=0) {
    free(x->list_r);
  }

  x->list_q = x->list_r = 0;

  if(have_gsl) {
    if (x->a!=0) {
      my_matrix_free(x->a);
    }
    if (x->tau!=0) {
      my_vector_free(x->tau);
    }

    if (x->q!=0) {
      my_matrix_free(x->q);
    }
    if (x->r!=0) {
      my_matrix_free(x->r);
    }
  }

  x->list_q = 0;
  x->list_r = 0;

  x->a = 0;
  x->tau = 0;
  x->q = 0;
  x->r = 0;
}

static void deleteMTXQr (MTXQr *x)
{
  deleteMTXqr(x);
}

static void *newMTXQr (t_symbol *s, int argc, t_atom *argv)
{
  MTXQr *x = (MTXQr *) pd_new (mtx_qr_class);
  x->list_q_out = outlet_new (&x->x_obj, gensym("matrix"));
  x->list_r_out = outlet_new (&x->x_obj, gensym("matrix"));

  if (!have_gsl) {
    static int warn_gsl = 1;
    if(warn_gsl)
#ifdef HAVE_LIBGSL
      pd_error(x, "[%s] couldn't find (recent enough) GSL", s->s_name);
#else
      pd_error(x, "[%s] compiled without GSL", s->s_name);
#endif
    warn_gsl = 0;
  }
  return ((void *) x);
}

static void mTXQrBang (MTXQr *x)
{
  if (x->list_q) {
     outlet_anything(x->list_r_out, gensym("matrix"), x->rows*x->columns+2,
                    x->list_r);
     outlet_anything(x->list_q_out, gensym("matrix"), x->rows*x->rows+2, x->list_q);
  }
}

static void mTXQrMatrix (MTXQr *x, t_symbol *s,
                         int argc, t_atom *argv)
{
  int rows, columns, size;
  int m,n;
  if(!have_gsl)
    return;

  if(iemmatrix_check(x, argc, argv, 0))return;
  rows = atom_getint (argv++);
  columns = atom_getint (argv++);
  size = rows * columns;

  x->rows=rows;
  x->columns=columns;

  deleteMTXqr(x);
  allocMTXqr(x);

  for (n=0; n<size; n++) {
    x->a->data[n]=(double) atom_getfloat(argv++);
  }

  my_linalg_QR_decomp(x->a,x->tau);
  my_linalg_QR_unpack(x->a,x->tau, x->q, x->r);

  SETFLOAT((x->list_q),(float) x->rows);
  SETFLOAT((x->list_q+1),(float) x->rows);
  for (n=0,size=0; n<x->rows; n++) {
    for (m=0; m<x->rows; m++) {
      SETFLOAT((x->list_q+2+size), x->q->data[size]);
      size++;
    }
  }

  SETFLOAT((x->list_r),(float) x->rows);
  SETFLOAT((x->list_r+1),(float) x->columns);
  for (n=0,size=0; n<x->rows; n++) {
    for (m=0; m<x->columns; m++) {
      SETFLOAT((x->list_r+2+size), x->r->data[size]);
      size++;
    }
  }

  mTXQrBang(x);
}

void mtx_qr_setup (void)
{
  mtx_qr_class = class_new
                 (gensym("mtx_qr"),
                  (t_newmethod) newMTXQr,
                  (t_method) deleteMTXQr,
                  sizeof (MTXQr),
                  CLASS_DEFAULT, A_GIMME, 0);
  class_addbang (mtx_qr_class, (t_method) mTXQrBang);
  class_addmethod (mtx_qr_class, (t_method) mTXQrMatrix, gensym("matrix"),
                   A_GIMME,0);

#ifdef HAVE_LIBGSL
  my_matrix_alloc = iemmatrix_get_stub("gsl_matrix_alloc", mtx_qr_class);
  my_matrix_free = iemmatrix_get_stub("gsl_matrix_free", mtx_qr_class);
  my_vector_alloc = iemmatrix_get_stub("gsl_vector_complex_alloc", mtx_qr_class);
  my_vector_free = iemmatrix_get_stub("gsl_vector_complex_free", mtx_qr_class);
  my_linalg_QR_decomp = iemmatrix_get_stub("gsl_linalg_QR_decomp", mtx_qr_class);
  my_linalg_QR_unpack = iemmatrix_get_stub("gsl_linalg_QR_unpack", mtx_qr_class);
#endif
  have_gsl = (
              my_matrix_alloc && my_matrix_free
              && my_vector_alloc && my_vector_free
              && my_linalg_QR_decomp
              && my_linalg_QR_unpack
              );
}

void iemtx_qr_setup(void)
{
  mtx_qr_setup();
}
