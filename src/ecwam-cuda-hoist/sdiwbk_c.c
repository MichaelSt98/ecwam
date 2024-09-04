#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sdiwbk_c.h"


__device__ void sdiwbk_c(int kijs, int kijl, const double * __restrict__ fl1, 
  double * __restrict__ fld, double * __restrict__ sl, 
  const double * __restrict__ depth, const double * __restrict__ emaxdpt, 
  const double * __restrict__ emean, const double * __restrict__ f1mean, int lbiwbk, 
  int nang, int nfre, int nfre_red, int ichnk, int nchnk, int ij, 
  double * __restrict__ sds) {
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  
  int k;
  int m;
  int ic;
  double zhook_handle;
  double alph;
  double arg;
  double q;
  double q_old;
  double rel_err;
  double expq;
  
  double alph_b_j = (double) 1.0;
  double coef_b_j = 2*alph_b_j;
  double depthtrs = (double) 50.0;
  
  // ----------------------------------------------------------------------
  
  
  //*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
  //*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
  //        --------------------------------------------------------------
  
  
  if (lbiwbk) {
    //       (FOLLOWING BATTJES-JANSSEN AND BEJI)
    if (depth[ij - 1 + kijl*(ichnk - 1)] < depthtrs) {
      alph = (double) 2.0*emaxdpt[ij - 1 + kijl*(ichnk - 1)] / emean[ij - 1];
      arg = fmin((double) (alph), (double) ((double) 50.0));
      q_old = exp((double) (-arg));
      //            USE NEWTON-RAPHSON METHOD
      for (ic = 1; ic <= 15; ic += 1) {
        expq = exp((double) (-arg*((double) 1.0 - q_old)));
        q = q_old - (expq - q_old) / (arg*expq - (double) 1.0);
        rel_err = fabs((double) (q - q_old)) / q_old;
        // IF(REL_ERR.LT.0.00001_JWRB) EXIT
        q_old = q;
      }
      q = fmin((double) (q), (double) ((double) 1.0));
      sds[ij - 1 + kijl*(ichnk - 1)] = coef_b_j*alph*q*f1mean[ij - 1];
    }
    
    for (m = 1; m <= nfre_red; m += 1) {
      for (k = 1; k <= nang; k += 1) {
        if (depth[ij - 1 + kijl*(ichnk - 1)] < depthtrs) {
          sl[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = sl[ij - 1 + kijl*(k - 1 + nang*(m - 
            1))] - sds[ij - 1 + kijl*(ichnk - 1)]*fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 
            + nfre*(ichnk - 1)))];
          fld[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = 
            fld[ij - 1 + kijl*(k - 1 + nang*(m - 1))] - sds[ij - 1 + kijl*(ichnk - 1)];
        }
      }
    }
    
  }
  
  
  
}
