#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sdepthlim_c.h"

__device__ void sdepthlim_c(int kijs, int kijl, const double * emaxdpt, double * fl1, 
  double delth, const double * dfim, double epsmin, const double * fr, int nang, 
  int nfre, double wetail, int ichnk, int nchnk, int ij) {
  
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  
  int k;
  int m;
  double delt25;
  double em;
  double temp;
  int llepsmin;
  

  em = epsmin;
  for (m = 1; m <= nfre; m += 1) {
    k = 1;
    temp = 
      fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))];
    for (k = 2; k <= nang; k += 1) {
      temp = temp + fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
        nfre_loki_param*(ichnk - 1)))];
    }
    em = em + dfim[m - 1]*temp;
  }
  delt25 = wetail*fr[nfre - 1]*delth;
  em = em + delt25*temp;
  
  em = min((double) (emaxdpt[ij - 1 + kijl*(ichnk - 1)] / em), (double) ((double) 1.0));
  
  for (m = 1; m <= nfre; m += 1) {
    for (k = 1; k <= nang; k += 1) {
      fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))]
         = max((double) (fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
        nfre_loki_param*(ichnk - 1)))]*em), (double) (epsmin));
    }
  }

  
  
}
