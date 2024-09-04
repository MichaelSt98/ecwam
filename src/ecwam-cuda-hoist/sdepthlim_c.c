#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sdepthlim_c.h"
#include "semean_c.h"


__device__ void sdepthlim_c(int kijs, int kijl, const double * __restrict__ emaxdpt, 
  double * __restrict__ fl1, double delth, const double * __restrict__ dfim, 
  double epsmin, const double * __restrict__ fr, int nang, int nfre, double wetail, 
  int ichnk, int nchnk, int ij, double * __restrict__ em, 
  double * __restrict__ semean_temp) {
  
  
  
  // ----------------------------------------------------------------------

  
  int k;
  int m;
  double zhook_handle;
  int llepsmin;
  
  // ----------------------------------------------------------------------
  
  
  llepsmin = true;
  semean_c(fl1, kijs, kijl,  (&em[ + kijl*(ichnk - 1)]), llepsmin, delth, dfim, epsmin, 
    fr, nang, nfre, wetail, ichnk, nchnk, ij, semean_temp);
  
  
  em[ij - 1 + kijl*(ichnk - 1)] = fmin((double) (emaxdpt[ij - 1 + kijl*(ichnk - 1)] / 
    em[ij - 1 + kijl*(ichnk - 1)]), (double) ((double) 1.0));
  
  for (m = 1; m <= nfre; m += 1) {
    for (k = 1; k <= nang; k += 1) {
      fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = fmax((double) 
        (fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))]*em[ij - 1 + 
        kijl*(ichnk - 1)]), (double) (epsmin));
    }
  }
  
  
  
}
