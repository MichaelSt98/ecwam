#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "cimsstrn_c.h"
#include "aki_ice_c.h"

__device__ void cimsstrn_c(int kijs, int kijl, const double * fl1, 
  const double * wavnum, const double * depth, const double * cithick, double * strn, 
  double delth, const double * dfim, double flmin, double g, int nang, int nfre, 
  double rowater, int ichnk, int nchnk, int ij) {
  
  

  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  int m;
  int k;
  double f1lim;
  double xki;
  double e;
  double sume;
  f1lim = flmin / delth;
  

  strn[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  for (m = 1; m <= nfre; m += 1) {
    xki = aki_ice_c(g, wavnum[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))], 
      depth[ij - 1 + kijl*(ichnk - 1)], rowater, cithick[ij - 1 + kijl*(ichnk - 1)]);
    e = (double) 0.5*cithick[ij - 1 + kijl*(ichnk - 1)]*(pow(xki, 3)) / wavnum[ij - 1 + 
      kijl*(m - 1 + nfre_loki_param*(ichnk - 1))];
    
    sume = (double) 0.0;
    for (k = 1; k <= nang; k += 1) {
      sume = sume + fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
        nfre_loki_param*(ichnk - 1)))];
    }
    
    if (sume > f1lim) {
      strn[ij - 1 + kijl*(ichnk - 1)] = 
        strn[ij - 1 + kijl*(ichnk - 1)] + (pow(e, 2))*sume*dfim[m - 1];
    }
    
  }

  
  
}
