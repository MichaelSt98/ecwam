#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "frcutindex_c.h"

__device__ void frcutindex_c(int kijs, int kijl, const double * fm, const double * fmws, 
  const double * ufric, const double * cicover, int * mij, double * rhowgdfth, 
  double cithrsh_tail, double epsmin, double flogsprdm1, const double * fr, double fric, 
  double g, int nfre, const double * rhowg_dfim, double tailfactor, 
  double tailfactor_pm, const double * zpifr, int ichnk, int nchnk, int ij) {
  
  
  
  const int nfre_loki_param = 36;
  int m;
  
  double fpmh;
  double fppm;
  double fm2;
  double fpm;
  double fpm4;
  fpmh = tailfactor / fr[1 - 1];
  fppm = tailfactor_pm*g / (fric*zpifr[1 - 1]);
  

  if (cicover[ij - 1 + kijl*(ichnk - 1)] <= cithrsh_tail) {
    fm2 = max((double) (fmws[ij - 1]), (double) (fm[ij - 1]))*fpmh;
    fpm = fppm / max((double) (ufric[ij - 1 + kijl*(ichnk - 1)]), (double) (epsmin));
    fpm4 = max((double) (fm2), (double) (fpm));
    // mij[ij - 1 + kijl*(ichnk - 1)] = nint(log10(fpm4)*flogsprdm1) + 1;
    mij[ij - 1 + kijl*(ichnk - 1)] = rint(log10(fpm4)*flogsprdm1) + 1;
    mij[ij - 1 + kijl*(ichnk - 1)] = min((double) (max((double) (1), (double) (mij[ij - 1
       + kijl*(ichnk - 1)]))), (double) (nfre));
  } else {
    mij[ij - 1 + kijl*(ichnk - 1)] = nfre;
  }
  for (m = 1; m <= mij[ij - 1 + kijl*(ichnk - 1)]; m += 1) {
    rhowgdfth[ij - 1 + kijl*(m - 1)] = rhowg_dfim[m - 1];
  }
  if (mij[ij - 1 + kijl*(ichnk - 1)] != nfre) {
    rhowgdfth[ij - 1 + kijl*(mij[ij - 1 + kijl*(ichnk - 1)] - 1)] = 
      (double) 0.5*rhowgdfth[ij - 1 + kijl*(mij[ij - 1 + kijl*(ichnk - 1)] - 1)];
  }
  for (m = mij[ij - 1 + kijl*(ichnk - 1)] + 1; m <= nfre; m += 1) {
    rhowgdfth[ij - 1 + kijl*(m - 1)] = (double) 0.0;
  }

  
  
}
