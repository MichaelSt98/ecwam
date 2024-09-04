#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "halphap_c.h"
#include "meansqs_lf_c.h"
#include "femean_c.h"


__device__ void halphap_c(int kijs, int kijl, const double * __restrict__ wavnum, 
  const double * __restrict__ coswdif, const double * __restrict__ fl1, 
  double * __restrict__ halp, double alphapmax, double delth, 
  const double * __restrict__ dfim, const double * __restrict__ dfimofr, double epsmin, 
  const double * __restrict__ fr, const double * __restrict__ fr5, double frtail, 
  int nang, int nfre, double wetail, double zpi4gm2, int ichnk, int nchnk, int ij, 
  double * __restrict__ alphap, double * __restrict__ xmss, double * __restrict__ em, 
  double * __restrict__ fm, double * __restrict__ f1d, double * __restrict__ wd, 
  double * __restrict__ flwd, double * __restrict__ femean_temp2, 
  double * __restrict__ meansqs_lf_fd, double * __restrict__ meansqs_lf_temp1, 
  double * __restrict__ meansqs_lf_temp2) {
  
  
  
  // ----------------------------------------------------------------------
  


  
  int k;
  int m;
  
  double zlnfrnfre;
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  
  zlnfrnfre = log(fr[nfre - 1]);
  
  // Find spectrum in wind direction
  
  for (k = 1; k <= nang; k += 1) {
    wd[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] = (double) 0.5 + (double) 
      0.5*copysign((double) ((double) 1.0), (double) (coswdif[ij - 1 + kijl*(k - 1)]));
  }
  
  for (m = 1; m <= nfre; m += 1) {
    for (k = 1; k <= nang; k += 1) {
      flwd[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = fl1[ij - 1 + 
        kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))]*wd[ij - 1 + kijl*(k - 1 + 
        nang*(ichnk - 1))];
    }
  }
  
  
  meansqs_lf_c(nfre, kijs, kijl,  (&flwd[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), 
    wavnum,  (&xmss[ + kijl*(ichnk - 1)]), dfim, nang, nfre, ichnk, nchnk, ij, 
    meansqs_lf_fd, meansqs_lf_temp1, meansqs_lf_temp2);
  
  femean_c(kijs, kijl,  (&flwd[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), 
     (&em[ + kijl*(ichnk - 1)]),  (&fm[ + kijl*(ichnk - 1)]), delth, dfim, dfimofr, 
    epsmin, fr, frtail, nang, nfre, wetail, ichnk, nchnk, ij, femean_temp2);
  
  
  if (em[ij - 1 + kijl*(ichnk - 1)] > (double) 0.0 && fm[ij - 1 + kijl*(ichnk - 1)] < 
    fr[nfre - 2 - 1]) {
    alphap[ij - 1 + kijl*(ichnk - 1)] = 
      xmss[ij - 1 + kijl*(ichnk - 1)] / (zlnfrnfre - log(fm[ij - 1 + kijl*(ichnk - 1)]));
    if (alphap[ij - 1 + kijl*(ichnk - 1)] > alphapmax) {
      // some odd cases, revert to tail value
      f1d[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
      for (k = 1; k <= nang; k += 1) {
        f1d[ij - 1 + kijl*(ichnk - 1)] = f1d[ij - 1 + kijl*(ichnk - 1)] + flwd[ij - 1 + 
          kijl*(k - 1 + nang*(nfre - 1 + nfre*(ichnk - 1)))]*delth;
      }
      alphap[ij - 1 + kijl*(ichnk - 1)] = 
        zpi4gm2*fr5[nfre - 1]*f1d[ij - 1 + kijl*(ichnk - 1)];
    }
  } else {
    f1d[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
    for (k = 1; k <= nang; k += 1) {
      f1d[ij - 1 + kijl*(ichnk - 1)] = f1d[ij - 1 + kijl*(ichnk - 1)] + flwd[ij - 1 + 
        kijl*(k - 1 + nang*(nfre - 1 + nfre*(ichnk - 1)))]*delth;
    }
    alphap[ij - 1 + kijl*(ichnk - 1)] = 
      zpi4gm2*fr5[nfre - 1]*f1d[ij - 1 + kijl*(ichnk - 1)];
  }
  
  //     1/2 ALPHAP:
  halp[ij - 1] = 
    (double) 0.5*fmin((double) (alphap[ij - 1 + kijl*(ichnk - 1)]), (double) (alphapmax))
    ;
  
  
  
}
