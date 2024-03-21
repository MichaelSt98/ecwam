#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "fkmean_c.h"

__device__ void fkmean_c(int kijs, int kijl, const double * fl1, const double * wavnum, 
  double * em, double * fm1, double * f1, double * ak, double * xk, double delth, 
  const double * dfim, const double * dfimfr, const double * dfimofr, double epsmin, 
  const double * fr, double frtail, double g, int nang, int nfre, double wetail, 
  double wp1tail, double zpi, int ichnk, int nchnk, int ij) {
  
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  
  
  int m;
  int k;
  double delt25;
  double coefm1;
  double coef1;
  double coefa;
  double coefx;
  double sqrtk;
  double tempa;
  double tempx;
  double temp2;

  em[ij - 1] = epsmin;
  fm1[ij - 1] = epsmin;
  f1[ij - 1] = epsmin;
  ak[ij - 1] = epsmin;
  xk[ij - 1] = epsmin;
  
  delt25 = wetail*fr[nfre - 1]*delth;
  coefm1 = frtail*delth;
  coef1 = wp1tail*delth*(pow(fr[nfre - 1], 2));
  coefa = coefm1*sqrt((double) (g)) / zpi;
  coefx = coef1*(zpi / sqrt((double) (g)));
  for (m = 1; m <= nfre; m += 1) {
    sqrtk = sqrt((double) (wavnum[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))]));
    tempa = dfim[m - 1] / sqrtk;
    tempx = sqrtk*dfim[m - 1];
    k = 1;
    temp2 = 
      fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))];
    for (k = 2; k <= nang; k += 1) {
      temp2 = temp2 + fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
        nfre_loki_param*(ichnk - 1)))];
    }
    em[ij - 1] = em[ij - 1] + dfim[m - 1]*temp2;
    fm1[ij - 1] = fm1[ij - 1] + dfimofr[m - 1]*temp2;
    f1[ij - 1] = f1[ij - 1] + dfimfr[m - 1]*temp2;
    ak[ij - 1] = ak[ij - 1] + tempa*temp2;
    xk[ij - 1] = xk[ij - 1] + tempx*temp2;
  }
  em[ij - 1] = em[ij - 1] + delt25*temp2;
  fm1[ij - 1] = fm1[ij - 1] + coefm1*temp2;
  fm1[ij - 1] = em[ij - 1] / fm1[ij - 1];
  f1[ij - 1] = f1[ij - 1] + coef1*temp2;
  f1[ij - 1] = f1[ij - 1] / em[ij - 1];
  ak[ij - 1] = ak[ij - 1] + coefa*temp2;
  ak[ij - 1] = pow((em[ij - 1] / ak[ij - 1]), 2);
  xk[ij - 1] = xk[ij - 1] + coefx*temp2;
  xk[ij - 1] = pow((xk[ij - 1] / em[ij - 1]), 2);

  
  
}
