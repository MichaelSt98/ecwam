#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sbottom_c.h"


__device__ void sbottom_c(int kijs, int kijl, const double * __restrict__ fl1, 
  double * __restrict__ fld, double * __restrict__ sl, 
  const double * __restrict__ wavnum, const double * __restrict__ depth, 
  double bathymax, double gm1, int nang, int nfre, int nfre_red, int ichnk, int nchnk, 
  int ij, double * __restrict__ sbo) {
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  
  int k;
  int m;
  double const_var;
  double arg;
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  
  const_var = -(double) 2.0*(double) 0.038*gm1;
  
  for (m = 1; m <= nfre_red; m += 1) {
    if (depth[ij - 1 + kijl*(ichnk - 1)] < bathymax) {
      arg = (double) 2.0*depth[ij - 1 + kijl*(ichnk - 1)]*wavnum[ij - 1 + kijl*(m - 1 + 
        nfre*(ichnk - 1))];
      arg = fmin((double) (arg), (double) ((double) 50.0));
      sbo[ij - 1 + kijl*(ichnk - 1)] = 
        const_var*wavnum[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))] / sinh(arg);
    } else {
      sbo[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
    }
    
    for (k = 1; k <= nang; k += 1) {
      sl[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = sl[ij - 1 + kijl*(k - 1 + nang*(m - 1))]
         + sbo[ij - 1 + kijl*(ichnk - 1)]*fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + 
        nfre*(ichnk - 1)))];
      fld[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = 
        fld[ij - 1 + kijl*(k - 1 + nang*(m - 1))] + sbo[ij - 1 + kijl*(ichnk - 1)];
    }
  }
  
  
  
}
