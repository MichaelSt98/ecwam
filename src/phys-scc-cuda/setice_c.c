#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "setice_c.h"

__device__ void setice_c(int kijs, int kijl, double * fl1, const double * cicover, 
  const double * coswdif, double cithrsh, double epsmin, double flmin, int nang, 
  int nfre, int ichnk, int nchnk, int ij) {
  
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  int m;
  int k;
  
  double cireduc;
  double temp;
  double icefree;

  if (cicover[ij - 1 + kijl*(ichnk - 1)] > cithrsh) {
    cireduc = max((double) (epsmin), (double) (((double) 1.0 - cicover[ij - 1 + 
      kijl*(ichnk - 1)])));
    icefree = (double) 0.0;
  } else {
    cireduc = (double) 0.0;
    icefree = (double) 1.0;
  }
  
  temp = cireduc*flmin;
  for (m = 1; m <= nfre; m += 1) {
    for (k = 1; k <= nang; k += 1) {
      fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))]
         = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1
        )))]*icefree + temp*(pow(max((double) ((double) 0.0), (double) (coswdif[ij - 1 + 
        kijl*(k - 1)])), 2));
    }
  }

  
  
}
