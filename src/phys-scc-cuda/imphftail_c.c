#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "imphftail_c.h"

__device__ void imphftail_c(int kijs, int kijl, const int * mij, const double * flm, 
  const double * wavnum, const double * xk2cg, double * fl1, int nang, int nfre, 
  int ichnk, int nchnk, int ij) {
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  int k;
  int m;
  
  double akm1;
  double tfac;
  double temp1;
  double temp2;

  temp1 = (double) 1.0 / xk2cg[ij - 1 + kijl*(mij[ij - 1 + kijl*(ichnk - 1)] - 1 + 
    nfre_loki_param*(ichnk - 1))] / wavnum[ij - 1 + kijl*(mij[ij - 1 + kijl*(ichnk - 1)] 
    - 1 + nfre_loki_param*(ichnk - 1))];
  
  for (m = mij[ij - 1 + kijl*(ichnk - 1)] + 1; m <= nfre; m += 1) {
    temp2 = (double) 1.0 / xk2cg[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))] / 
      wavnum[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))];
    temp2 = temp2 / temp1;
    for (k = 1; k <= nang; k += 1) {
      tfac = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(mij[ij - 1 + kijl*(ichnk - 1)] -
         1 + nfre_loki_param*(ichnk - 1)))];
      fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))]
         = max((double) (temp2*tfac), (double) (flm[ij - 1 + kijl*(k - 1)]));
    }
  }

  
}
