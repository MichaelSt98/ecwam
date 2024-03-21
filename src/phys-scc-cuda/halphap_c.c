#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "halphap_c.h"

__device__ void halphap_c(int kijs, int kijl, const double * wavnum, 
  const double * coswdif, const double * fl1, double * halp, double alphapmax, 
  double delth, const double * dfim, const double * dfimofr, double epsmin, 
  const double * fr, const double * fr5, double frtail, int nang, int nfre, 
  double wetail, double zpi4gm2, int ichnk, int nchnk, int ij) {
  
  // Loki: parameters from YOWPARAM inlined
  
  
  const int nfre_loki_param = 36;
  const int nang_loki_param = 24;
  
  int k;
  int m;
  
  double zlnfrnfre;
  double delt25;
  double delt2;
  double del2;
  double temp1;
  double temp2;
  double alphap;
  double xmss;
  double em;
  double fm;
  double f1d;
  double flwd[36];
  
  zlnfrnfre = log(fr[nfre - 1]);
  
  delt25 = wetail*fr[nfre - 1]*delth;
  delt2 = frtail*delth;

  for (m = 1; m <= nfre; m += 1) {
    for (k = 1; k <= nang; k += 1) {
      flwd[k - 1] = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
        nfre_loki_param*(ichnk - 1)))]*(double) 0.5 + (double) 0.5*copysign((double) ((double
        ) 1.0), (double) (coswdif[ij - 1 + kijl*(k - 1)]));
    }
    
    xmss = (double) 0.;
    temp1 = 
      dfim[m - 1]*(pow(wavnum[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))], 2));
    temp2 = (double) 0.0;
    for (k = 1; k <= nang; k += 1) {
      temp2 = temp2 + flwd[k - 1];
    }
    xmss = xmss + temp1*temp2;
    
    k = 1;
    em = (double) 0.;
    fm = (double) 0.;
    temp2 = max((double) (flwd[k - 1]), (double) (epsmin));
    for (k = 2; k <= nang; k += 1) {
      temp2 = temp2 + max((double) (flwd[k - 1]), (double) (epsmin));
    }
    em = em + temp2*dfim[m - 1];
    fm = fm + dfimofr[m - 1]*temp2;
  }
  
  for (k = 1; k <= nang; k += 1) {
    flwd[k - 1] = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(nfre - 1 + 
      nfre_loki_param*(ichnk - 1)))]*(double) 0.5 + (double) 0.5*copysign((double) ((double) 
      1.0), (double) (coswdif[ij - 1 + kijl*(k - 1)]));
  }
  
  em = em + delt25*temp2;
  fm = fm + delt2*temp2;
  fm = em / fm;
  fm = max((double) (fm), (double) (fr[1 - 1]));
  
  if (em > (double) 0.0 && fm < fr[-2 + nfre - 1]) {
    alphap = xmss / (zlnfrnfre - log(fm));
    if (alphap > alphapmax) {
      // some odd cases, revert to tail value
      f1d = (double) 0.0;
      for (k = 1; k <= nang; k += 1) {
        f1d = f1d + flwd[k - 1]*delth;
      }
      alphap = zpi4gm2*fr5[nfre - 1]*f1d;
    }
  } else {
    f1d = (double) 0.0;
    for (k = 1; k <= nang; k += 1) {
      f1d = f1d + flwd[k - 1]*delth;
    }
    alphap = zpi4gm2*fr5[nfre - 1]*f1d;
  }
  halp[ij - 1] = (double) 0.5*min((double) (alphap), (double) (alphapmax));

  
  
}
