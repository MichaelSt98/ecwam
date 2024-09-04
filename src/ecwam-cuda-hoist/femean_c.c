#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "femean_c.h"


__device__ void femean_c(int kijs, int kijl, const double * __restrict__ f, 
  double * __restrict__ em, double * __restrict__ fm, double delth, 
  const double * __restrict__ dfim, const double * __restrict__ dfimofr, double epsmin, 
  const double * __restrict__ fr, double frtail, int nang, int nfre, double wetail, 
  int ichnk, int nchnk, int ij, double * __restrict__ temp2) {
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  
  
  int m;
  int k;
  double delt25;
  double delt2;
  double del2;
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  
  //*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
  //        ------------------------------------------------
  
  
  em[ij - 1] = (double) 0.0;
  fm[ij - 1] = (double) 0.0;
  
  delt25 = wetail*fr[nfre - 1]*delth;
  delt2 = frtail*delth;
  
  //*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
  //        ------------------------------------------
  
  for (m = 1; m <= nfre; m += 1) {
    k = 1;
    temp2[ij - 1 + kijl*(ichnk - 1)] = 
      fmax((double) (f[ij - 1 + kijl*(k - 1 + nang*(m - 1))]), (double) (epsmin));
    for (k = 2; k <= nang; k += 1) {
      temp2[ij - 1 + kijl*(ichnk - 1)] = temp2[ij - 1 + kijl*(ichnk - 1)] + fmax((double)
         (f[ij - 1 + kijl*(k - 1 + nang*(m - 1))]), (double) (epsmin));
    }
    em[ij - 1] = em[ij - 1] + temp2[ij - 1 + kijl*(ichnk - 1)]*dfim[m - 1];
    fm[ij - 1] = fm[ij - 1] + dfimofr[m - 1]*temp2[ij - 1 + kijl*(ichnk - 1)];
  }
  
  
  //*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
  //*       NORMALIZE WITH TOTAL ENERGY.
  //        ------------------------------------------
  
  em[ij - 1] = em[ij - 1] + delt25*temp2[ij - 1 + kijl*(ichnk - 1)];
  fm[ij - 1] = fm[ij - 1] + delt2*temp2[ij - 1 + kijl*(ichnk - 1)];
  fm[ij - 1] = em[ij - 1] / fm[ij - 1];
  fm[ij - 1] = fmax((double) (fm[ij - 1]), (double) (fr[1 - 1]));
  
  
  
}
