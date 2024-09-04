#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sdissip_jan_c.h"


__device__ void sdissip_jan_c(int kijs, int kijl, const double * __restrict__ fl1, 
  double * __restrict__ fld, double * __restrict__ sl, 
  const double * __restrict__ wavnum, const double * __restrict__ emean, 
  const double * __restrict__ f1mean, const double * __restrict__ xkmean, double cdis, 
  double cdisvis, double delta_sdis, int nang, int nfre, double rnu, double zpi, 
  int ichnk, int nchnk, int ij, double * __restrict__ temp1, double * __restrict__ sds, 
  double * __restrict__ x, double * __restrict__ xk2) {
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  
  
  int k;
  int m;
  
  double scdfm;
  double consd;
  double conss;
  double delta_sdism1;
  double cvis;
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  
  //*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
  //*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
  //        --------------------------------------------------------------
  
  delta_sdism1 = (double) 1.0 - delta_sdis;
  
  conss = cdis*zpi;
  
  sds[ij - 1 + kijl*(ichnk - 1)] = 
    conss*f1mean[ij - 1]*(pow(emean[ij - 1], 2))*(pow(xkmean[ij - 1], 4));
  
  for (m = 1; m <= nfre; m += 1) {
    x[ij - 1 + kijl*(ichnk - 1)] = 
      wavnum[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))] / xkmean[ij - 1];
    xk2[ij - 1 + kijl*(ichnk - 1)] = 
      pow(wavnum[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))], 2);
    
    cvis = rnu*cdisvis;
    temp1[ij - 1 + kijl*(ichnk - 1)] = sds[ij - 1 + kijl*(ichnk - 1)]*x[ij - 1 + 
      kijl*(ichnk - 1)]*(delta_sdism1 + delta_sdis*x[ij - 1 + kijl*(ichnk - 1)]) + 
      cvis*xk2[ij - 1 + kijl*(ichnk - 1)];
    
    for (k = 1; k <= nang; k += 1) {
      fld[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = 
        fld[ij - 1 + kijl*(k - 1 + nang*(m - 1))] + temp1[ij - 1 + kijl*(ichnk - 1)];
      sl[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = sl[ij - 1 + kijl*(k - 1 + nang*(m - 1))]
         + temp1[ij - 1 + kijl*(ichnk - 1)]*fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + 
        nfre*(ichnk - 1)))];
    }
    
  }
  
  
  
}
