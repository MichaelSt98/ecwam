#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "imphftail_c.h"


__device__ void imphftail_c(int kijs, int kijl, const int * __restrict__ mij, 
  const double * __restrict__ flm, const double * __restrict__ wavnum, 
  const double * __restrict__ xk2cg, double * __restrict__ fl1, int nang, int nfre, 
  int ichnk, int nchnk, int ij, double * __restrict__ temp1, double * __restrict__ temp2)
   {
  
  
  // ----------------------------------------------------------------------
  
  
  
  
  int k;
  int m;
  
  double akm1;
  double tfac;
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  
  //*    DIAGNOSTIC TAIL.
  //     ----------------
  
  
  temp1[ij - 1 + kijl*(ichnk - 1)] = (double) 1.0 / xk2cg[ij - 1 + kijl*(mij[ij - 1 + 
    kijl*(ichnk - 1)] - 1 + nfre*(ichnk - 1))] / wavnum[ij - 1 + kijl*(mij[ij - 1 + 
    kijl*(ichnk - 1)] - 1 + nfre*(ichnk - 1))];
  
  for (m = mij[ij - 1 + kijl*(ichnk - 1)] + 1; m <= nfre; m += 1) {
    temp2[ij - 1 + kijl*(ichnk - 1)] = (double) 1.0 / xk2cg[ij - 1 + kijl*(m - 1 + 
      nfre*(ichnk - 1))] / wavnum[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))];
    temp2[ij - 1 + kijl*(ichnk - 1)] = 
      temp2[ij - 1 + kijl*(ichnk - 1)] / temp1[ij - 1 + kijl*(ichnk - 1)];
    
    //*    MERGE TAIL INTO SPECTRA.
    //     ------------------------
    for (k = 1; k <= nang; k += 1) {
      tfac = fl1[ij - 1 + kijl*(k - 1 + nang*(mij[ij - 1 + kijl*(ichnk - 1)] - 1 + 
        nfre*(ichnk - 1)))];
      fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = fmax((double) 
        (temp2[ij - 1 + kijl*(ichnk - 1)]*tfac), (double) (flm[ij - 1 + kijl*(k - 1)]));
    }
  }
  
  
  // ----------------------------------------------------------------------
  
  
}
