#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "setice_c.h"


__device__ void setice_c(int kijs, int kijl, double * __restrict__ fl1, 
  const double * __restrict__ cicover, const double * __restrict__ coswdif, 
  double cithrsh, double epsmin, double flmin, int nang, int nfre, int ichnk, int nchnk, 
  int ij, double * __restrict__ cireduc, double * __restrict__ temp, 
  double * __restrict__ icefree) {
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  int m;
  int k;
  
  double zhook_handle;
  // ----------------------------------------------------------------------
  
  
  //*    1. SET SPECTRA TO NOISE LEVEL OVER ICE POINTS.
  //     ----------------------------------------------
  
  
  if (cicover[ij - 1 + kijl*(ichnk - 1)] > cithrsh) {
    cireduc[ij - 1 + kijl*(ichnk - 1)] = fmax((double) (epsmin), (double) (((double) 1.0 
      - cicover[ij - 1 + kijl*(ichnk - 1)])));
    icefree[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  } else {
    cireduc[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
    icefree[ij - 1 + kijl*(ichnk - 1)] = (double) 1.0;
  }
  
  temp[ij - 1 + kijl*(ichnk - 1)] = cireduc[ij - 1 + kijl*(ichnk - 1)]*flmin;
  for (m = 1; m <= nfre; m += 1) {
    for (k = 1; k <= nang; k += 1) {
      fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = fl1[ij - 1 + kijl*(k
         - 1 + nang*(m - 1 + nfre*(ichnk - 1)))]*icefree[ij - 1 + kijl*(ichnk - 1)] + 
        temp[ij - 1 + kijl*(ichnk - 1)]*(pow(fmax((double) ((double) 0.0), (double) 
        (coswdif[ij - 1 + kijl*(k - 1)])), 2));
    }
  }
  
  
  
}
