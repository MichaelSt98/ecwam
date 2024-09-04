#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "semean_c.h"


__device__ void semean_c(const double * __restrict__ fl1, int kijs, int kijl, 
  double * __restrict__ em, int llepsmin, double delth, 
  const double * __restrict__ dfim, double epsmin, const double * __restrict__ fr, 
  int nang, int nfre, double wetail, int ichnk, int nchnk, int ij, 
  double * __restrict__ temp) {
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  
  int m;
  int k;
  
  double delt25;
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  
  //*    1. INITIALISE ENERGY ARRAY.
  //        ------------------------
  
  
  if (llepsmin) {
    em[ij - 1] = epsmin;
  } else {
    em[ij - 1] = (double) 0.0;
  }
  
  // ----------------------------------------------------------------------
  
  //*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
  //        -----------------------------------------
  
  for (m = 1; m <= nfre; m += 1) {
    k = 1;
    temp[ij - 1 + kijl*(ichnk - 1)] = 
      fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
    for (k = 2; k <= nang; k += 1) {
      temp[ij - 1 + kijl*(ichnk - 1)] = temp[ij - 1 + kijl*(ichnk - 1)] + fl1[ij - 1 + 
        kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
    }
    em[ij - 1] = em[ij - 1] + dfim[m - 1]*temp[ij - 1 + kijl*(ichnk - 1)];
  }
  
  // ----------------------------------------------------------------------
  
  //*    3. ADD TAIL ENERGY.
  //        ----------------
  
  delt25 = wetail*fr[nfre - 1]*delth;
  em[ij - 1] = em[ij - 1] + delt25*temp[ij - 1 + kijl*(ichnk - 1)];
  
  
  
}
