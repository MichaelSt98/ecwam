#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "meansqs_lf_c.h"


__device__ void meansqs_lf_c(int nfre_eff, int kijs, int kijl, 
  const double * __restrict__ f, const double * __restrict__ wavnum, 
  double * __restrict__ xmss, const double * __restrict__ dfim, int nang, int nfre, 
  int ichnk, int nchnk, int ij, double * __restrict__ fd, double * __restrict__ temp1, 
  double * __restrict__ temp2) {
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  int m;
  int k;
  int kfre;
  
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  //*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
  //        ------------------------------------------
  
  kfre = fmin((double) (nfre_eff), (double) (nfre));
  
  
  xmss[ij - 1] = (double) 0.0;
  
  //*    2.2 SHALLOW WATER INTEGRATION.
  //         --------------------------
  
  for (m = 1; m <= kfre; m += 1) {
    temp1[ij - 1 + kijl*(ichnk - 1)] = 
      dfim[m - 1]*(pow(wavnum[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))], 2));
    temp2[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
    for (k = 1; k <= nang; k += 1) {
      temp2[ij - 1 + kijl*(ichnk - 1)] = 
        temp2[ij - 1 + kijl*(ichnk - 1)] + f[ij - 1 + kijl*(k - 1 + nang*(m - 1))];
    }
    xmss[ij - 1] = 
      xmss[ij - 1] + temp1[ij - 1 + kijl*(ichnk - 1)]*temp2[ij - 1 + kijl*(ichnk - 1)];
  }
  
  
  
}
