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

//#undef HAVE_GSL_EIGEN_NONSYMM

#ifdef HAVE_GSL_EIGEN_NONSYMM
#include <gsl/gsl_eigen.h>
#else
#include "stub/gsl.h"
#endif

typedef int(*t_eigen_nonsymm)(gsl_matrix*, gsl_vector_complex*, gsl_eigen_nonsymm_workspace*);
typedef int(*t_eigen_nonsymmv)(gsl_matrix*, gsl_vector_complex*, gsl_matrix_complex*, gsl_eigen_nonsymmv_workspace*);

IEMMATRIX_DECLARE_ALLOCFREE2_STUB(my_matrix);
IEMMATRIX_DECLARE_ALLOCFREE2_STUB(my_matrix_complex);
IEMMATRIX_DECLARE_ALLOCFREE_STUB(my_vector_complex);
IEMMATRIX_DECLARE_ALLOCFREE_STUB(my_eigen_nonsymmv);
IEMMATRIX_DECLARE_ALLOCFREE_STUB(my_eigen_nonsymm);
static t_eigen_nonsymm my_eigen_nonsymm = 0;
static t_eigen_nonsymmv my_eigen_nonsymmv = 0;

static int have_gsl = 0;

static t_class *mtx_eig_class;
enum WithEigenVectors {WITHEVS=1, WITHOUTEVS=0};
typedef struct _MTXEig_ MTXEig;
struct _MTXEig_ {
  t_object x_obj;
  gsl_matrix *a;
  gsl_matrix_complex *q;
  gsl_vector_complex *l;
  gsl_eigen_nonsymm_workspace *w;
  gsl_eigen_nonsymmv_workspace *wv;
  t_outlet *list_q_out_re;
  t_outlet *list_q_out_im;
  t_outlet *list_l_out_re;
  t_outlet *list_l_out_im;
  t_atom *list_q_re;
  t_atom *list_q_im;
  t_atom *list_l_re;
  t_atom *list_l_im;
  int size;
  enum WithEigenVectors withevs;
};

static void allocMTXqlw (MTXEig *x)
{
  if (!have_gsl)return;
  x->a=(gsl_matrix*)my_matrix_alloc(x->size,x->size);
  x->l=(gsl_vector_complex*)my_vector_complex_alloc(x->size);

  switch (x->withevs) {
  case WITHEVS:
    x->wv=(gsl_eigen_nonsymmv_workspace*)my_eigen_nonsymmv_alloc(x->size);
    x->q=(gsl_matrix_complex*)my_matrix_complex_alloc(x->size,x->size);
    break;
  case WITHOUTEVS:
    x->w=(gsl_eigen_nonsymm_workspace*)my_eigen_nonsymm_alloc(x->size);
  }

  x->list_q_re=(t_atom*)calloc(sizeof(t_atom),x->size*x->size+2);
  x->list_q_im=(t_atom*)calloc(sizeof(t_atom),x->size*x->size+2);
  x->list_l_re=(t_atom*)calloc(sizeof(t_atom),x->size);
  x->list_l_im=(t_atom*)calloc(sizeof(t_atom),x->size);
}

static void deleteMTXqlw (MTXEig *x)
{
  if (x->list_q_re!=0) {
    free(x->list_q_re);
  }
  if (x->list_q_im!=0) {
    free(x->list_q_im);
  }
  if (x->list_l_re!=0) {
    free(x->list_l_re);
  }
  if (x->list_l_im!=0) {
    free(x->list_l_im);
  }

  x->list_q_re = 0;
  x->list_q_im = 0;
  x->list_l_re = 0;
  x->list_l_im = 0;

  if (have_gsl) {
    if (x->a!=0) {
      my_matrix_free(x->a);
    }
    if (x->q!=0) {
      my_matrix_complex_free(x->q);
    }
    if (x->l!=0) {
      my_vector_complex_free(x->l);
    }
    if (x->w!=0) {
      my_eigen_nonsymm_free(x->w);
    }
    if (x->wv!=0) {
      my_eigen_nonsymmv_free(x->wv);
    }
  }

  x->a = 0;
  x->q = 0;
  x->l = 0;
  x->w = 0;
  x->wv = 0;
}

static void deleteMTXEig (MTXEig *x)
{
  deleteMTXqlw(x);
}

static void *newMTXEig (t_symbol *s, int argc, t_atom *argv)
{
  MTXEig *x = (MTXEig *) pd_new (mtx_eig_class);

  x->list_l_out_re = outlet_new (&x->x_obj, gensym("list"));
  x->list_l_out_im = outlet_new (&x->x_obj, gensym("list"));
  if (argc && atom_getsymbol(argv)==gensym("v")) {
    x->withevs=1;
    x->list_q_out_re = outlet_new (&x->x_obj, gensym("matrix"));
    x->list_q_out_im = outlet_new (&x->x_obj, gensym("matrix"));
  }
  if (!have_gsl) {
    static int warn_gsl = 1;
    if(warn_gsl)
#ifdef HAVE_GSL_EIGEN_NONSYMM
      pd_error(x, "[%s] couldn't find (recent enough) GSL", s->s_name);
#else
      pd_error(x, "[%s] compiled without GSL", s->s_name);
#endif
    warn_gsl = 0;
  }

  return ((void *) x);
}

static void mTXEigBang (MTXEig *x)
{
  if (!have_gsl)return;
  if (x->list_l_re) {
    switch (x->withevs) {
    case WITHEVS:
      outlet_anything(x->list_q_out_im, gensym("matrix"), x->size*x->size+2,
                      x->list_q_im);
      outlet_anything(x->list_q_out_re, gensym("matrix"), x->size*x->size+2,
                      x->list_q_re);
    case WITHOUTEVS:
      outlet_anything(x->list_l_out_im, gensym("list"), x->size, x->list_l_im);
      outlet_anything(x->list_l_out_re, gensym("list"), x->size, x->list_l_re);
    }
  }
}

static void mTXEigMatrix (MTXEig *x, t_symbol *s,
                          int argc, t_atom *argv)
{
  int rows, columns, size;
  int n,m;
  float f;

  gsl_complex c;
  if (!have_gsl) {
    return;
  }
  /* size check */
  if(iemmatrix_check(x, argc, argv, 0))return;
  rows = atom_getint (argv++);
  columns = atom_getint (argv++);
  size=rows*columns;
  if (rows!=columns) {
    pd_error(x, "[mtx_eig]: Eigendecomposition works for square matrices only!");
    return;
  }
  size=rows;
  x->size=size;

  deleteMTXqlw(x);
  allocMTXqlw(x);

  for (n=0; n<size; n++) {
    x->a->data[n]=(double) atom_getfloat(argv++);
  }

  switch (x->withevs) {
  case WITHOUTEVS:
    my_eigen_nonsymm(x->a,x->l,x->w);
    break;
  case WITHEVS:
    my_eigen_nonsymmv(x->a,x->l,x->q,x->wv);
    SETFLOAT((x->list_q_re),(float) x->size);
    SETFLOAT((x->list_q_im),(float) x->size);
    SETFLOAT((x->list_q_re+1),(float) x->size);
    SETFLOAT((x->list_q_im+1),(float) x->size);
    for (n=0; n<size; n++) {
      SETFLOAT((x->list_q_im+2+n), (float) x->q->data[2*n+1]);
      SETFLOAT((x->list_q_re+2+n), (float) x->q->data[2*n]);
    }
    break;
  }

  for (n=0; n<x->size; n++) {
    f=(float) GSL_VECTOR_IMAG(x->l, n);
    SETFLOAT((x->list_l_im+n), f);
    f=(float) GSL_VECTOR_REAL(x->l, n);
    SETFLOAT((x->list_l_re+n), f);
  }

  mTXEigBang(x);
}

void mtx_eig_setup (void)
{
  mtx_eig_class = class_new
                  (gensym("mtx_eig"),
                   (t_newmethod) newMTXEig,
                   (t_method) deleteMTXEig,
                   sizeof (MTXEig),
                   CLASS_DEFAULT, A_GIMME, 0);
  class_addbang (mtx_eig_class, (t_method) mTXEigBang);
  class_addmethod (mtx_eig_class, (t_method) mTXEigMatrix, gensym("matrix"),
                   A_GIMME,0);

#ifdef HAVE_GSL_EIGEN_NONSYMM
  my_matrix_alloc = iemmatrix_get_stub("gsl_matrix_alloc", mtx_eig_class);
  my_matrix_free = iemmatrix_get_stub("gsl_matrix_free", mtx_eig_class);
  my_matrix_complex_alloc = iemmatrix_get_stub("gsl_matrix_complex_alloc", mtx_eig_class);
  my_matrix_complex_free = iemmatrix_get_stub("gsl_matrix_complex_free", mtx_eig_class);
  my_vector_complex_alloc = iemmatrix_get_stub("gsl_vector_complex_alloc", mtx_eig_class);
  my_vector_complex_free = iemmatrix_get_stub("gsl_vector_complex_free", mtx_eig_class);
  my_eigen_nonsymmv_alloc = iemmatrix_get_stub("gsl_eigen_nonsymmv_alloc", mtx_eig_class);
  my_eigen_nonsymmv_free = iemmatrix_get_stub("gsl_eigen_nonsymmv_free", mtx_eig_class);
  my_eigen_nonsymm_alloc = iemmatrix_get_stub("gsl_eigen_nonsymm_alloc", mtx_eig_class);
  my_eigen_nonsymm_free = iemmatrix_get_stub("gsl_eigen_nonsymm_free", mtx_eig_class);
  my_eigen_nonsymm = iemmatrix_get_stub("gsl_eigen_nonsymm", mtx_eig_class);
  my_eigen_nonsymmv = iemmatrix_get_stub("gsl_eigen_nonsymmv", mtx_eig_class);
#endif
  have_gsl = (
              my_matrix_alloc && my_matrix_free
              && my_vector_complex_alloc && my_vector_complex_free
              && my_eigen_nonsymmv_alloc && my_eigen_nonsymmv_free
              && my_matrix_complex_alloc && my_matrix_complex_free
              && my_eigen_nonsymm_alloc && my_eigen_nonsymm_free
              && my_eigen_nonsymm_free
              && my_eigen_nonsymmv_free
              );
}

void iemtx_eig_setup(void)
{
  mtx_eig_setup();
}
