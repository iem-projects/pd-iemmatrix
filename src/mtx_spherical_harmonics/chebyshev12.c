/*
 Evaluates all fully normalized circular harmonics
 at the angles phi up to the order nmax.
 using the recurrence for the Chebyshev
 polynomials of the first and second kind
 T has the dimensions length(phi) x 2nmax+1

 Implementation by Franz Zotter, Institute of Electronic Music and Acoustics
 (IEM), University of Music and Dramatic Arts (KUG), Graz, Austria
 http://iem.at/Members/zotter, 2008.

 This code is published under the Gnu General Public License, see
 "LICENSE.txt"
*/

#include "mtx_spherical_harmonics/chebyshev12.h"

Cheby12WorkSpace *chebyshev12_alloc(const size_t nmax, const size_t l, CHNormType type)
{
  Cheby12WorkSpace *wc;
  // memory allocation
  if ((wc=(Cheby12WorkSpace*)calloc(1,sizeof(Cheby12WorkSpace)))!=0) {
    wc->l=l;
    wc->nmax=nmax;
    if ((wc->t=(double*)calloc(l*(2*nmax+1),sizeof(double)))==0) {
      free(wc);
      return 0;
    }
    switch(type) {
       case N2D2PI:
	   wc->n0=1.0;
	   wc->nm=sqrt(2.0);
	   break;
       case SN2D:
	   wc->n0=1.0;
	   wc->nm=1.0;
	   break;
       case N2D:
       default:
           wc->n0=1.0/sqrt(2.0*M_PI);
	   wc->nm=1.0/sqrt(M_PI);
    }
    return wc;
  }
  return 0;
}

void chebyshev12_free(Cheby12WorkSpace *wc)
{
  if (wc!=0) {
    free(wc->t);
    free(wc);
  }
}

void chebyshev12(double *phi, Cheby12WorkSpace *wc)
{
  unsigned int l,l0,n;
  const int incr=wc?(2*wc->nmax+1):0;
  double *cosphi;
  double *sinphi;
  // memory allocation
  if ((wc!=0)&&(phi!=0)) {
    // constants and initialization
    for (l=0, l0=wc->nmax; l<wc->l; l++, l0+=incr) {
      // initial value T_0=1
      wc->t[l0]=wc->n0;
    }
    if (wc->nmax==0) {
      return;
    } 
    if ((cosphi=(double*)calloc(wc->l,sizeof(double)))==0) {
      return;
    }
    if ((sinphi=(double*)calloc(wc->l,sizeof(double)))==0) {
      free(cosphi);
      return;
    }
    for (l=0, l0=wc->nmax; l<wc->l; l++, l0+=incr) {
      // initial value T_0=1
      cosphi[l]=cos(phi[l]);
      sinphi[l]=sin(phi[l]);
      wc->t[l0+1]=cosphi[l]*wc->nm;
      wc->t[l0-1]=sinphi[l]*wc->nm;
    }
    // recurrence for n>1
    for (n=2; n<=wc->nmax; n++) {
      for (l=0, l0=wc->nmax; l<wc->l; l++, l0+=incr) {
        wc->t[l0+n]=cosphi[l]* wc->t[l0+n-1] - sinphi[l]* wc->t[l0-n+1];
        wc->t[l0-n]=sinphi[l]* wc->t[l0+n-1] + cosphi[l]* wc->t[l0-n+1];
      }
    }
    free(cosphi);
    free(sinphi);
  }
}
