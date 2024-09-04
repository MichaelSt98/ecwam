#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "ciwabr_c.h"


__device__ void ciwabr_c(int kijs, int kijl, const double * __restrict__ cicover, 
  const double * __restrict__ fl1, const double * __restrict__ wavnum, 
  const double * __restrict__ cgroup, double * __restrict__ ciwab, double cdicwa, 
  const double * __restrict__ dfim, double epsmin, int idelt, int licerun, int lmaskice, 
  int nang, int nfre, int ichnk, int nchnk, int ij) {
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  int k;
  int m;
  double ewh;
  double x;
  double alp;
  double xk2;
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  
  
  if (!licerun || lmaskice) {
    
    for (m = 1; m <= nfre; m += 1) {
      for (k = 1; k <= nang; k += 1) {
        ciwab[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = (double) 1.0;
      }
    }
    
  } else {
    
    for (m = 1; m <= nfre; m += 1) {
      for (k = 1; k <= nang; k += 1) {
        ewh = (double) 4.0*sqrt((double) (fmax((double) (epsmin), (double) (fl1[ij - 1 + 
          kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))]*dfim[m - 1]))));
        xk2 = pow(wavnum[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))], 2);
        alp = cdicwa*xk2*ewh;
        x = alp*cgroup[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))]*idelt;
        ciwab[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = (double) 1.0 - cicover[ij - 1 + 
          kijl*(ichnk - 1)]*((double) 1.0 - exp((double) (-fmin((double) (x), (double) 
          ((double) 50.0)))));
      }
    }
    
  }
  
  
  
}
