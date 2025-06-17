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
# include <gsl/gsl_linalg.h>
#else
# include "stub/gsl.h"
#endif

IEMMATRIX_DECLARE_ALLOCFREE_STUB(my_vector);
IEMMATRIX_DECLARE_ALLOCFREE2_STUB(my_matrix);
typedef int(*t_linalg_SV_decomp)(gsl_matrix*, gsl_matrix*, gsl_vector*, gsl_vector*);
static t_linalg_SV_decomp my_linalg_SV_decomp = 0;

static int have_gsl = 0;


static t_class *mtx_svd_class;

typedef struct _MTXSvd_ MTXSvd;
struct _MTXSvd_ {
  t_object x_obj;

  gsl_matrix *u;
  gsl_vector *s;
  gsl_matrix *v;
  gsl_vector *w;

  t_outlet *list_u_out;
  t_outlet *list_s_out;
  t_outlet *list_v_out;
  t_atom *list_u;
  t_atom *list_s;
  t_atom *list_v;
  int rows;
  int columns;
};

static void allocMTXusvw (MTXSvd *x)
{
  if(have_gsl) {
    x->u=(gsl_matrix*)my_matrix_alloc(x->rows,x->columns);
    x->s=(gsl_vector*)my_vector_alloc(x->columns);
    x->v=(gsl_matrix*)my_matrix_alloc(x->columns,x->columns);
    x->w=(gsl_vector*)my_vector_alloc(x->columns);

    x->list_u=(t_atom*)calloc(x->rows*x->columns+2, sizeof(t_atom));
    x->list_s=(t_atom*)calloc(x->columns, sizeof(t_atom));
    x->list_v=(t_atom*)calloc(x->columns*x->columns+2, sizeof(t_atom));
  }
}

static void deleteMTXusvw (MTXSvd *x)
{
  if (x->list_u!=0) {
    free(x->list_u);
  }
  if (x->list_s!=0) {
    free(x->list_s);
  }
  if (x->list_v!=0) {
    free(x->list_v);
  }

  x->list_u = x->list_s = x->list_v = 0;

  if(have_gsl) {
    if (x->u!=0) {
      my_matrix_free(x->u);
    }
    if (x->s!=0) {
      my_vector_free(x->s);
    }
    if (x->v!=0) {
      my_matrix_free(x->v);
    }
    if (x->w!=0) {
      my_vector_free(x->w);
    }
  }

  x->u = 0;
  x->s = 0;
  x->v = 0;
  x->w = 0;
}

static void deleteMTXSvd (MTXSvd *x)
{
  deleteMTXusvw(x);
}

static void *newMTXSvd (t_symbol *s, int argc, t_atom *argv)
{
  MTXSvd *x = (MTXSvd *) pd_new (mtx_svd_class);
  (void)argc; (void)argv; /* unused */
  x->list_u_out = outlet_new (&x->x_obj, gensym("matrix"));
  x->list_s_out = outlet_new (&x->x_obj, gensym("list"));
  x->list_v_out = outlet_new (&x->x_obj, gensym("matrix"));
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

static void mTXSvdBang (MTXSvd *x)
{
  if (x->list_u) {
    outlet_anything(x->list_v_out, gensym("matrix"), x->columns*x->columns+2,
                    x->list_v);
    outlet_anything(x->list_s_out, gensym("list"), x->columns, x->list_s);
    outlet_anything(x->list_u_out, gensym("matrix"), x->rows*x->columns+2,
                    x->list_u);
  }
}

static void mTXSvdMatrix (MTXSvd *x, t_symbol *s,
                          int argc, t_atom *argv)
{
  int rows, columns, size;
  int n;
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  rows = atom_getint (argv++);
  columns = atom_getint (argv++);
  size=rows*columns;

  if(!have_gsl)
    return;
  /* size check */
  if (rows<columns) {
    pd_error(x, "[mtx_svd]: gsl_linalg_SV_decomp does not support M<N");
    return;
  }
  x->rows=rows;
  x->columns=columns;

  deleteMTXusvw(x);
  allocMTXusvw(x);

  for (n=0; n<size; n++) {
    x->u->data[n]=(double) atom_getfloat(argv++);
  }

  my_linalg_SV_decomp(x->u,x->v,x->s,x->w);

  SETFLOAT((x->list_u),(float) x->rows);
  SETFLOAT((x->list_u+1),(float) x->columns);
  for (n=0; n<size; n++) {
    SETFLOAT((x->list_u+2+n), (float) x->u->data[n]);
  }

  for (n=0; n<x->columns; n++) {
    SETFLOAT((x->list_s+n),(float) x->s->data[n]);
  }

  SETFLOAT((x->list_v),(float) x->columns);
  SETFLOAT((x->list_v+1),(float) x->columns);
  size=x->columns*x->columns;
  for (n=0; n<size; n++) {
    SETFLOAT((x->list_v+n+2), (float) x->v->data[n]);
  }

  mTXSvdBang(x);
}

void mtx_svd_setup (void)
{
  mtx_svd_class = class_new
                  (gensym("mtx_svd"),
                   (t_newmethod) newMTXSvd,
                   (t_method) deleteMTXSvd,
                   sizeof (MTXSvd),
                   CLASS_DEFAULT, A_GIMME, 0);
  class_addbang (mtx_svd_class, (t_method) mTXSvdBang);
  class_addmethod (mtx_svd_class, (t_method) mTXSvdMatrix, gensym("matrix"),
                   A_GIMME,0);

#ifdef HAVE_LIBGSL
  my_matrix_alloc = iemmatrix_get_stub("gsl_matrix_alloc", mtx_svd_class);
  my_matrix_free = iemmatrix_get_stub("gsl_matrix_free", mtx_svd_class);
  my_vector_alloc = iemmatrix_get_stub("gsl_vector_complex_alloc", mtx_svd_class);
  my_vector_free = iemmatrix_get_stub("gsl_vector_complex_free", mtx_svd_class);
  my_linalg_SV_decomp = iemmatrix_get_stub("gsl_linalg_SV_decomp", mtx_svd_class);
#endif
  have_gsl = (
              my_matrix_alloc && my_matrix_free
              && my_vector_alloc && my_vector_free
              && my_linalg_SV_decomp
              );
}

void iemtx_svd_setup(void)
{
  mtx_svd_setup();
}
