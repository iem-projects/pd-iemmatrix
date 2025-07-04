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
#include <math.h>
#include <stdlib.h>
#include "iemmatrix_stub.h"

static t_class *mtx_bessel_class;

static t_besselfn my_jn = 0;
static t_besselfn my_yn = 0;

typedef struct _MTXBessel_ MTXBessel;
struct _MTXBessel_ {
  t_object x_obj;
  t_outlet *list_h_re_out;
  t_outlet *list_h_im_out;
  t_atom *list_h_re;
  t_atom *list_h_im;
  double *kr;
  double *h_re;
  double *h_im;
  size_t nmax;
  size_t l;
};

static void allocMTXBesseldata (MTXBessel *x)
{
  x->kr=(double*)calloc(x->l,sizeof(double));
  if (x->list_h_re_out!=0) {
    x->list_h_re=(t_atom*)calloc(x->l*(x->nmax+1)+2,sizeof(t_atom));
    x->h_re=(double*)calloc(x->l*(x->nmax+1),sizeof(double));
  }
  if (x->list_h_im_out!=0) {
    x->list_h_im=(t_atom*)calloc(x->l*(x->nmax+1)+2,sizeof(t_atom));
    x->h_im=(double*)calloc(x->l*(x->nmax+1),sizeof(double));
  }
}

static void deleteMTXBesseldata (MTXBessel *x)
{
  if (x->kr!=0) {
    free(x->kr);
  }
  if (x->h_re!=0) {
    free(x->h_re);
  }
  if (x->h_im!=0) {
    free(x->h_im);
  }
  if (x->list_h_re!=0) {
    free(x->list_h_re);
  }
  if (x->list_h_im!=0) {
    free(x->list_h_im);
  }
  x->list_h_re=0;
  x->list_h_im=0;
  x->h_re=0;
  x->h_im=0;
  x->kr=0;
}

static void *newMTXBessel (t_symbol *s, int argc, t_atom *argv)
{
  int nmax;
  char whichfunction = 'j';
  t_symbol *fsym = 0;
  MTXBessel *x = (MTXBessel *) pd_new (mtx_bessel_class);
  (void)s; /* unused */
  x->list_h_re = 0;
  x->list_h_im = 0;
  x->list_h_im_out = 0;
  x->list_h_re_out = 0;
  x->kr = 0;
  x->h_re = 0;
  x->h_im = 0;
  x->l=0;
  if(argc)
    fsym=atom_getsymbol(argv);
  if (fsym && fsym->s_name!=0) {
    whichfunction=fsym->s_name[0];
  }
  switch (whichfunction) {
  default:
  case 'J':
  case 'j':
    x->list_h_re_out = outlet_new (&x->x_obj, gensym("matrix"));
    break;
  case 'H':
  case 'h':
    x->list_h_re_out = outlet_new (&x->x_obj, gensym("matrix"));
  /* coverity[unterminated_case]: H has both real&imag outlet */ /* fall through */
  case 'Y':
  case 'y':
    x->list_h_im_out = outlet_new (&x->x_obj, gensym("matrix"));
  }
  nmax=(int) atom_getfloat(argv+1);
  if (nmax<0) {
    nmax=0;
  }
  x->nmax=nmax;

  return ((void *) x);
}

static void mTXBesselBang (MTXBessel *x)
{
  if (x->list_h_im!=0) {
    outlet_anything(x->list_h_im_out, gensym("matrix"), x->l*(x->nmax+1)+2,
                    x->list_h_im);
  }
  if (x->list_h_re!=0) {
    outlet_anything(x->list_h_re_out, gensym("matrix"), x->l*(x->nmax+1)+2,
                    x->list_h_re);
  }
}

static void mTXBesselMatrix (MTXBessel *x, t_symbol *s,
                             int argc, t_atom *argv)
{
  size_t columns;
  size_t n, m;
  (void)s; /* unused */

  /* size check */
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  argv++; /* rows */
  columns = atom_getint (argv++);

  if(my_jn && my_yn) {
    if (x->l != columns) {
      deleteMTXBesseldata(x);
      x->l=columns;
      allocMTXBesseldata(x);
    }
    for (n=0; n < x->l; n++) {
      x->kr[n]=(double) atom_getfloat(argv+n);
    }

    if (x->h_re != 0)
      for (m=0; m < x->l; m++)
        for (n=0; n < x->nmax+1; n++) {
          x->h_re[n+m*(x->nmax+1)]=my_jn(n,x->kr[m]);
        }

    if (x->h_im != 0)
      for (m=0; m < x->l; m++)
        for (n=0; n < x->nmax+1; n++) {
          x->h_im[n+m*(x->nmax+1)]=my_yn(n,x->kr[m]);
        }

    if (x->h_re != 0) {
      SETFLOAT(x->list_h_re+1,(t_float)(x->nmax+1));
      SETFLOAT(x->list_h_re,(t_float)x->l);
      for (n=0; n < x->l*(x->nmax+1); n++) {
        SETFLOAT(x->list_h_re+n+2,(t_float)x->h_re[n]);
      }
    }

    if (x->h_im != 0) {
      SETFLOAT(x->list_h_im+1,(t_float)(x->nmax+1));
      SETFLOAT(x->list_h_im,(t_float)x->l);
      for (n=0; n < x->l*(x->nmax+1); n++) {
        SETFLOAT(x->list_h_im+n+2,(t_float)x->h_im[n]);
      }
    }

    mTXBesselBang(x);
  } else {
    pd_error(x, "[mtx_bessel]: implementation requires math.h that implements jn and yn Bessel functions or gsl_sf_bessel.h");
  }
}

void mtx_bessel_setup (void)
{
  mtx_bessel_class = class_new
                     (gensym("mtx_bessel"),
                      (t_newmethod) newMTXBessel,
                      (t_method) deleteMTXBesseldata,
                      sizeof (MTXBessel),
                      CLASS_DEFAULT, A_GIMME, 0);
  class_addbang (mtx_bessel_class, (t_method) mTXBesselBang);
  class_addmethod (mtx_bessel_class, (t_method) mTXBesselMatrix,
                   gensym("matrix"), A_GIMME,0);

  my_jn = iemmatrix_get_stub("jn", mtx_bessel_class);
  my_yn = iemmatrix_get_stub("yn", mtx_bessel_class);
}

void iemtx_bessel_setup(void)
{
  mtx_bessel_setup();
}
