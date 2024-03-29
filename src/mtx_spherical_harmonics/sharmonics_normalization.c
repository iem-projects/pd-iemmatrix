/*
 * Recursive computation of (arbitrary degree) normalization constants
 * for spherical harmonics, according to Gumerov and Duraiswami,
 * "The Fast Multipole Methods for the Helmholtz Equation in Three Dimensions",
 * Elsevier, 2005.
 *
 * Implementation by Franz Zotter, Institute of Electronic Music and Acoustics
 * (IEM), University of Music and Dramatic Arts (KUG), Graz, Austria
 * http://iem.at/Members/zotter, 2007.
 *
 * This code is published under the Gnu General Public License, see
 * "LICENSE.txt"
 */

#include "mtx_spherical_harmonics/sharmonics_normalization.h"

SHNorml *sharmonics_normalization_new (const size_t nmax, const SHNormType type)
{
  SHNorml *wn;
  unsigned int n,n0,m;
  const double oneoversqrt2 = 1.0/sqrt(2);

  // memory allocation
  if ((wn=(SHNorml*)calloc(1,sizeof(SHNorml)))!=0) {
    wn->nmax=nmax;
    if ((wn->n=(double*)calloc((nmax+1)*(nmax+2)/2,sizeof(double)))==0) {
      free(wn);
      wn=0;
    } else {
      // computing N_n^m for m=0,
      switch(type) {
       	 case N3D4PI:
	 case SN3D:
	     wn->n[0]=1.0;
	     break;
	 case N3D:
	 default:
	     wn->n[0]=oneoversqrt2;
      }
      switch(type) {
	  case N3D:
	  case N3D4PI:
             for (n=1,n0=1; n<=nmax; n++) {
	        wn->n[n0]=wn->n[0] * sqrt(2*n+1);
		n0+=n+1;
	     }
	     break;
	  case SN3D:
	  default:
             for (n=1,n0=1; n<=nmax; n++) {
	        wn->n[n0]=wn->n[0];
		n0+=n+1;
	     }
      }
      // computing N_n^m for 0<m<=n (incl. Condon-Shortley sign)
      for (n=1,n0=1; n<=nmax; n++) {
        for (m=1; m<=n; m++) {
          wn->n[n0+m]= - wn->n[n0+m-1] / sqrt((n+m)*(n-m+1));
        }
        n0+=n+1;
      }
    }
  }
  return wn;
}

void sharmonics_normalization_free(SHNorml *wn)
{
  if (wn!=0) {
    free(wn->n);
    free(wn);
  }
}
