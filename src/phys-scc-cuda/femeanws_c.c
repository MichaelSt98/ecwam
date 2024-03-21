#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "femeanws_c.h"

__device__ void femeanws_c(int kijs, int kijl, const double * fl1, const double * xllws, 
  double * fm, double * em, double delth, const double * dfim, const double * dfimofr, 
  double epsmin, const double * fr, double frtail, int nang, int nfre, double wetail, 
  int ichnk, int nchnk, int ij) {
  
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  int m;
  int k;
  
  double delt25;
  double delt2;
  double temp2;
  double em_loc;

  em_loc = epsmin;
  fm[ij - 1] = epsmin;
  
  delt25 = wetail*fr[nfre - 1]*delth;
  delt2 = frtail*delth;
  for (m = 1; m <= nfre; m += 1) {
    temp2 = (double) 0.0;
    for (k = 1; k <= nang; k += 1) {
      temp2 = temp2 + xllws[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
        nfre_loki_param*(ichnk - 1)))]*fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 
        + nfre_loki_param*(ichnk - 1)))];
    }
    em_loc = em_loc + dfim[m - 1]*temp2;
    fm[ij - 1] = fm[ij - 1] + dfimofr[m - 1]*temp2;
  }
  em_loc = em_loc + delt25*temp2;
  fm[ij - 1] = fm[ij - 1] + delt2*temp2;
  fm[ij - 1] = em_loc / fm[ij - 1];
  
  // if (present( (*em))) {
    em[ij - 1] = em_loc;
  // }

  
  
}
