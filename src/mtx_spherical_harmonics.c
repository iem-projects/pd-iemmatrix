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
#include <stdlib.h>
#include "mtx_spherical_harmonics/sharmonics.h"
#include "mtx_spherical_harmonics/legendre_a.h"
#include "mtx_spherical_harmonics/chebyshev12.h"
#include "mtx_spherical_harmonics/sharmonics_normalization.h"

static t_class *mtx_spherical_harmonics_class;

typedef struct _MTXSh_ MTXSh;
struct _MTXSh_ {
  t_object x_obj;
  t_outlet *list_sh_out;
  t_atom *list_sh;
  SHNormType ntype;
  int legacy_azi_reverse;
  double *phi;
  double *theta;
  SHWorkSpace *ws;
  size_t nmax;
  size_t l;
};

static t_class *mtx_circular_harmonics_class;

typedef struct _MTXCh_ MTXCh;
struct _MTXCh_ {
  t_object x_obj;
  t_outlet *list_ch_out;
  t_atom *list_ch;
  CHNormType ntype;
  double *phi;
  Cheby12WorkSpace *wc;
  size_t nmax;
  size_t l;
};


static void allocMTXShdata (MTXSh *x)
{
  x->phi=(double*)calloc(x->l,sizeof(double));
  x->theta=(double*)calloc(x->l,sizeof(double));
  x->ws=sharmonics_alloc(x->nmax,x->l,x->ntype);
  x->list_sh=(t_atom*)calloc(x->l*(x->nmax+1)*(x->nmax+1)+2,sizeof(t_atom));
}

static void deleteMTXShdata (MTXSh *x)
{
  if (x->phi!=0) {
    free(x->phi);
  }
  if (x->theta!=0) {
    free(x->theta);
  }
  if (x->list_sh!=0) {
    free(x->list_sh);
  }
  sharmonics_free(x->ws);
  x->ws=0;
  x->list_sh=0;
  x->theta=0;
  x->phi=0;
}

static void *newMTXSh (t_symbol *s, int argc, t_atom *argv)
{
  int nmax;
  MTXSh *x = (MTXSh *) pd_new (mtx_spherical_harmonics_class);
  x->list_sh_out = outlet_new (&x->x_obj, gensym("matrix"));
  x->list_sh = 0;
  x->phi = 0;
  x->theta = 0;
  x->ws = 0;
  x->l=0;
  x->ntype=N3D;
  x->legacy_azi_reverse=1;
  t_symbol *nt;
  if (argc<1) {
     nmax=1;
  } else {
     if (argc>1) {
         nt=atom_getsymbol(argv+1);
	 x->legacy_azi_reverse=0;
	 if (nt==gensym("N3D")) {
	    x->ntype=N3D;
	 } else if (nt==gensym("N3D4PI")) {
	    x->ntype=N3D4PI;
	 } else if (nt==gensym("SN3D")) {
	    x->ntype=SN3D;
	 } else {
	    x->legacy_azi_reverse=1;
	    x->ntype=N3D;
	 }
     }
     nmax=(int) atom_getfloat(argv);
     if (nmax<0) {
         nmax=0;
     }
  }
  x->nmax=nmax;
  return ((void *) x);
}

static void mTXShBang (MTXSh *x)
{
  if (x->list_sh!=0) {
    outlet_anything(x->list_sh_out, gensym("matrix"),
                    x->l*(x->nmax+1)*(x->nmax+1)+2, x->list_sh);
  }
}

static void mTXShMatrix (MTXSh *x, t_symbol *s,
                         int argc, t_atom *argv)
{
  int rows, columns, size;

  /* size check */
  if(iemmatrix_check(x, argc, argv, 0))return;

  rows = atom_getint (argv++);
  columns = atom_getint (argv++);
  size = rows * columns;

  if ((rows!=2)||(columns<1)) {
    pd_error(x, "[mtx_spherical_harmonics]: 2 X L matrix expected with phi and theta vector, but got more rows/no entries");
    return;
  }
  if (x->l!=columns) {
    deleteMTXShdata(x);
    x->l=columns;
    allocMTXShdata(x);
  }
  if (1) {
    unsigned int n;

    switch (x->legacy_azi_reverse) {
    case 0:
      for (n=0; n<x->l; n++) {
        x->phi[n]=(double) atom_getfloat(argv+n);
        x->theta[n]=(double) atom_getfloat(argv+columns+n);
      }
      break;
    default:
      for (n=0; n<x->l; n++) {
        x->phi[n]=-(double) atom_getfloat(argv+n);
        x->theta[n]=(double) atom_getfloat(argv+columns+n);
      }
      break;
    }
  }

  if (x->ws!=0) {
    int n;
    sharmonics(x->phi, x->theta, x->ws);
    size=x->l*(x->nmax+1)*(x->nmax+1);
    SETFLOAT(x->list_sh,(t_float)x->l);
    SETFLOAT(x->list_sh+1,(t_float)(x->nmax+1)*(x->nmax+1));
    for (n=0; n<size; n++) {
      SETFLOAT(x->list_sh+n+2,(t_float)x->ws->y[n]);
    }
    mTXShBang(x);
  } else {
    pd_error(x, "[mtx_spherical_harmonics]: memory error, no operation");
  }
}

static void allocMTXChdata (MTXCh *x)
{
  x->phi=(double*)calloc(x->l,sizeof(double));
  x->wc=chebyshev12_alloc(x->nmax,x->l,x->ntype);
  x->list_ch=(t_atom*)calloc(x->l*(2*x->nmax+1)+2,sizeof(t_atom));
}

static void deleteMTXChdata (MTXCh *x)
{
  if (x->phi!=0) {
    free(x->phi);
  }
  if (x->list_ch!=0) {
    free(x->list_ch);
  }
  chebyshev12_free(x->wc);
  x->wc=0;
  x->list_ch=0;
  x->phi=0;
}

static void *newMTXCh (t_symbol *s, int argc, t_atom *argv)
{
  int nmax;
  MTXCh *x = (MTXCh *) pd_new (mtx_circular_harmonics_class);
  x->list_ch_out = outlet_new (&x->x_obj, gensym("matrix"));
  x->list_ch = 0;
  x->phi = 0;
  x->wc = 0;
  x->l=0;
  x->ntype=N2D;
  t_symbol *nt;
  if (argc<1) {
     nmax=1;
  } else {
     if (argc>1) {
         nt=atom_getsymbol(argv+1);
	 if (nt==gensym("N2D")) {
            x->ntype=N2D;
	 } else if (nt==gensym("N2D2PI")) {
	    x->ntype=N2D2PI;
	 } else if (nt==gensym("SN2D")) {
	    x->ntype=SN2D;
	 } else {
	    x->ntype=N2D;
	 }
     }
     nmax=(int) atom_getfloat(argv);
     if (nmax<0) {
         nmax=0;
     }
  }
  x->nmax=nmax;

  return ((void *) x);
}

static void mTXChBang (MTXCh *x)
{
  if (x->list_ch!=0) {
    outlet_anything(x->list_ch_out, gensym("matrix"), x->l*(2*x->nmax+1)+2,
                    x->list_ch);
  }
}

static void mTXChMatrix (MTXCh *x, t_symbol *s,
                         int argc, t_atom *argv)
{
  int rows, columns, size;

  /* size check */
  if(iemmatrix_check(x, argc, argv, 0))return;
  rows = atom_getint (argv++);
  columns = atom_getint (argv++);
  size = rows * columns;

  if ((rows!=1)||(columns<1)) {
    pd_error(x, "[mtx_circular_harmonics]: 1*L matrix expected with phi vector, but got more rows/no entries");
  } else {
    if (x->l!=columns) {
      deleteMTXChdata(x);
      x->l=columns;
      allocMTXChdata(x);
    }
    if (1) {
      unsigned int n;
      for (n=0; n<x->l; n++) {
        x->phi[n]=(double) atom_getfloat(argv+n);
      }
    }

    if (x->wc!=0) {
      int n;
      chebyshev12(x->phi, x->wc);
      size=x->l*(2*x->nmax+1);
      SETFLOAT(x->list_ch,(t_float)x->l);
      SETFLOAT(x->list_ch+1,(t_float)(2*x->nmax+1));
      for (n=0; n<size; n++) {
        SETFLOAT(x->list_ch+n+2,(t_float)x->wc->t[n]);
      }
      mTXChBang(x);
    } else {
      pd_error(x, "[mtx_circular_harmonics]: memory error, no operation");
    }
  }


}

void mtx_spherical_harmonics_setup (void)
{
  mtx_spherical_harmonics_class = class_new
                                  (gensym("mtx_spherical_harmonics"),
                                   (t_newmethod) newMTXSh,
                                   (t_method) deleteMTXShdata,
                                   sizeof (MTXSh),
                                   CLASS_DEFAULT, A_GIMME, 0);
  class_addbang (mtx_spherical_harmonics_class, (t_method) mTXShBang);
  class_addmethod (mtx_spherical_harmonics_class, (t_method) mTXShMatrix,
                   gensym("matrix"), A_GIMME,0);
}


void mtx_circular_harmonics_setup (void)
{
  mtx_circular_harmonics_class = class_new
                                 (gensym("mtx_circular_harmonics"),
                                  (t_newmethod) newMTXCh,
                                  (t_method) deleteMTXChdata,
                                  sizeof (MTXCh),
                                  CLASS_DEFAULT, A_GIMME, 0);
  class_addbang (mtx_circular_harmonics_class, (t_method) mTXChBang);
  class_addmethod (mtx_circular_harmonics_class, (t_method) mTXChMatrix,
                   gensym("matrix"), A_GIMME,0);
}

void iemtx_spherical_harmonics_setup(void)
{
  mtx_circular_harmonics_setup();
  mtx_spherical_harmonics_setup();
}
