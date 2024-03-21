#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sdissip_jan_c.h"

__device__ void sdissip_jan_c(int kijs, int kijl, const double * fl1, double * fld, 
  double * sl, const double * wavnum, const double * emean, const double * f1mean, 
  const double * xkmean, double cdis, double cdisvis, double delta_sdis, int nang, 
  int nfre, double rnu, double zpi, int ichnk, int nchnk, int ij) {
  
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  
  int k;
  int m;
  
  double scdfm;
  double consd;
  double conss;
  double delta_sdism1;
  double cvis;
  double temp1;
  double sds;
  double x;
  double xk2;
  delta_sdism1 = (double) 1.0 - delta_sdis;
  
  conss = cdis*zpi;

  sds = conss*f1mean[ij - 1]*(pow(emean[ij - 1], 2))*(pow(xkmean[ij - 1], 4));
  
  for (m = 1; m <= nfre; m += 1) {
    x = wavnum[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))] / xkmean[ij - 1];
    xk2 = pow(wavnum[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))], 2);
    
    cvis = rnu*cdisvis;
    temp1 = sds*x*(delta_sdism1 + delta_sdis*x) + cvis*xk2;
    
    for (k = 1; k <= nang; k += 1) {
      fld[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] = 
        fld[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] + temp1;
      sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] = sl[ij - 1 + kijl*(k - 1 + 
        nang_loki_param*(m - 1))] + temp1*fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m -
         1 + nfre_loki_param*(ichnk - 1)))];
    }
    
  }

  
  
}
