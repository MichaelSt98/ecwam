#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sdiwbk_c.h"

__device__ void sdiwbk_c(int kijs, int kijl, const double * fl1, double * fld, 
  double * sl, const double * depth, const double * emaxdpt, const double * emean, 
  const double * f1mean, int lbiwbk, int nang, int nfre_red, int ichnk, int nchnk, int ij
  ) {
  
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  int k;
  int m;
  int ic;
  double alph;
  double arg;
  double q;
  double q_old;
  double rel_err;
  double expq;
  double sds;
  
  double alph_b_j = (double) 1.0;
  double coef_b_j = 2*alph_b_j;
  double depthtrs = (double) 50.0;

  if (lbiwbk) {
    //       (FOLLOWING BATTJES-JANSSEN AND BEJI)
    if (depth[ij - 1 + kijl*(ichnk - 1)] < depthtrs) {
      alph = (double) 2.0*emaxdpt[ij - 1 + kijl*(ichnk - 1)] / emean[ij - 1];
      arg = min((double) (alph), (double) ((double) 50.0));
      q_old = exp((double) (-arg));
      //            USE NEWTON-RAPHSON METHOD
      for (ic = 1; ic <= 15; ic += 1) {
        expq = exp((double) (-arg*((double) 1.0 - q_old)));
        q = q_old - (expq - q_old) / (arg*expq - (double) 1.0);
        rel_err = abs((double) (q - q_old)) / q_old;
        if (rel_err < (double) 0.00001) {
          // EXIT
        }
        q_old = q;
      }
      q = min((double) (q), (double) ((double) 1.0));
      sds = coef_b_j*alph*q*f1mean[ij - 1];
    }
    
    for (m = 1; m <= nfre_red; m += 1) {
      for (k = 1; k <= nang; k += 1) {
        if (depth[ij - 1 + kijl*(ichnk - 1)] < depthtrs) {
          sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] = sl[ij - 1 + kijl*(k - 1 +
             nang_loki_param*(m - 1))] - sds*fl1[ij - 1 + kijl*(k - 1 + 
            nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))];
          fld[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] = 
            fld[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] - sds;
        }
      }
    }
    
  }

  
  
}
