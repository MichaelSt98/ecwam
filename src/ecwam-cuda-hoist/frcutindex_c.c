#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "frcutindex_c.h"


__device__ void frcutindex_c(int kijs, int kijl, const double * __restrict__ fm, 
  const double * __restrict__ fmws, const double * __restrict__ ufric, 
  const double * __restrict__ cicover, int * __restrict__ mij, 
  double * __restrict__ rhowgdfth, double cithrsh_tail, double epsmin, 
  double flogsprdm1, const double * __restrict__ fr, double fric, double g, int nfre, 
  const double * __restrict__ rhowg_dfim, double tailfactor, double tailfactor_pm, 
  const double * __restrict__ zpifr, int ichnk, int nchnk, int ij) {
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  
  int m;
  
  double fpmh;
  double fppm;
  double fm2;
  double fpm;
  double fpm4;
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  
  //*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
  //*    FREQUENCIES LE MAX(TAILFACTOR*MAX(FMNWS,FM),TAILFACTOR_PM*FPM),
  //*    WHERE FPM IS THE PIERSON-MOSKOWITZ FREQUENCY BASED ON FRICTION
  //*    VELOCITY. (FPM=G/(FRIC*ZPI*USTAR))
  //     ------------------------------------------------------------
  
  fpmh = tailfactor / fr[1 - 1];
  fppm = tailfactor_pm*g / (fric*zpifr[1 - 1]);
  
  
  if (cicover[ij - 1 + kijl*(ichnk - 1)] <= cithrsh_tail) {
    fm2 = fmax((double) (fmws[ij - 1]), (double) (fm[ij - 1]))*fpmh;
    fpm = fppm / fmax((double) (ufric[ij - 1 + kijl*(ichnk - 1)]), (double) (epsmin));
    fpm4 = fmax((double) (fm2), (double) (fpm));
    mij[ij - 1 + kijl*(ichnk - 1)] = rint((double) (log10(fpm4)*flogsprdm1)) + 1;
    mij[ij - 1 + kijl*(ichnk - 1)] = fmin((double) (fmax((double) (1), (double) (mij[ij -
       1 + kijl*(ichnk - 1)]))), (double) (nfre));
  } else {
    mij[ij - 1 + kijl*(ichnk - 1)] = nfre;
  }
  
  //     SET RHOWGDFTH
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
