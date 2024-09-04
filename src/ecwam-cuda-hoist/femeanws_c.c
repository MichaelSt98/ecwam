#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "femeanws_c.h"


__device__ void femeanws_c(int kijs, int kijl, const double * __restrict__ fl1, 
  const double * __restrict__ xllws, double * __restrict__ fm, double * __restrict__ em, 
  double delth, const double * __restrict__ dfim, const double * __restrict__ dfimofr, 
  double epsmin, const double * __restrict__ fr, double frtail, int nang, int nfre, 
  double wetail, int ichnk, int nchnk, int ij, double * __restrict__ temp2, 
  double * __restrict__ em_loc) {
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  
  int m;
  int k;
  
  double delt25;
  double delt2;
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  
  //*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
  //        ------------------------------------------------
  
  
  em_loc[ij - 1 + kijl*(ichnk - 1)] = epsmin;
  fm[ij - 1] = epsmin;
  
  delt25 = wetail*fr[nfre - 1]*delth;
  delt2 = frtail*delth;
  
  
  //*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
  //        ------------------------------------------
  
  for (m = 1; m <= nfre; m += 1) {
    temp2[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
    for (k = 1; k <= nang; k += 1) {
      temp2[ij - 1 + kijl*(ichnk - 1)] = temp2[ij - 1 + kijl*(ichnk - 1)] + xllws[ij - 1 
        + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))]*fl1[ij - 1 + kijl*(k - 1 + 
        nang*(m - 1 + nfre*(ichnk - 1)))];
    }
    em_loc[ij - 1 + kijl*(ichnk - 1)] = 
      em_loc[ij - 1 + kijl*(ichnk - 1)] + dfim[m - 1]*temp2[ij - 1 + kijl*(ichnk - 1)];
    fm[ij - 1] = fm[ij - 1] + dfimofr[m - 1]*temp2[ij - 1 + kijl*(ichnk - 1)];
  }
  
  //*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
  //*       NORMALIZE WITH TOTAL ENERGY.
  //        ------------------------------------------
  
  em_loc[ij - 1 + kijl*(ichnk - 1)] = 
    em_loc[ij - 1 + kijl*(ichnk - 1)] + delt25*temp2[ij - 1 + kijl*(ichnk - 1)];
  fm[ij - 1] = fm[ij - 1] + delt2*temp2[ij - 1 + kijl*(ichnk - 1)];
  fm[ij - 1] = em_loc[ij - 1 + kijl*(ichnk - 1)] / fm[ij - 1];
  
  // IF(PRESENT(EM))THEN
  em[ij - 1] = em_loc[ij - 1 + kijl*(ichnk - 1)];
  
  // ENDIF
  
  
}
