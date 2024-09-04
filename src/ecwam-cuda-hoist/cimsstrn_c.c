#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "cimsstrn_c.h"
#include "aki_ice_c.h"


__device__ void cimsstrn_c(int kijs, int kijl, const double * __restrict__ fl1, 
  const double * __restrict__ wavnum, const double * __restrict__ depth, 
  const double * __restrict__ cithick, double * __restrict__ strn, double delth, 
  const double * __restrict__ dfim, double flmin, double g, int nang, int nfre, 
  double rowater, int ichnk, int nchnk, int ij, double * __restrict__ xki, 
  double * __restrict__ e, double * __restrict__ sume) {
  
  
  
  // ----------------------------------------------------------------------
  

  
  
  int m;
  int k;
  double f1lim;
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  
  //*    1. INITIALISE
  //        ----------
  
  f1lim = flmin / delth;
  
  
  strn[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  
  // ----------------------------------------------------------------------
  
  //*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
  //        ------------------------------------------
  
  for (m = 1; m <= nfre; m += 1) {
    xki[ij - 1 + kijl*(ichnk - 1)] = aki_ice_c(g, wavnum[ij - 1 + kijl*(m - 1 + 
      nfre*(ichnk - 1))], depth[ij - 1 + kijl*(ichnk - 1)], rowater, cithick[ij - 1 + 
      kijl*(ichnk - 1)]);
    e[ij - 1 + kijl*(ichnk - 1)] = (double) 0.5*cithick[ij - 1 + kijl*(ichnk - 1)
      ]*(pow(xki[ij - 1 + kijl*(ichnk - 1)], 3)) / wavnum[ij - 1 + kijl*(m - 1 + 
      nfre*(ichnk - 1))];
    
    sume[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
    for (k = 1; k <= nang; k += 1) {
      sume[ij - 1 + kijl*(ichnk - 1)] = sume[ij - 1 + kijl*(ichnk - 1)] + fl1[ij - 1 + 
        kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
    }
    
    if (sume[ij - 1 + kijl*(ichnk - 1)] > f1lim) {
      strn[ij - 1 + kijl*(ichnk - 1)] = strn[ij - 1 + kijl*(ichnk - 1)] + (pow(e[ij - 1 +
         kijl*(ichnk - 1)], 2))*sume[ij - 1 + kijl*(ichnk - 1)]*dfim[m - 1];
    }
    
  }
  
  
  
}
