#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "fkmean_c.h"


__device__ void fkmean_c(int kijs, int kijl, const double * __restrict__ fl1, 
  const double * __restrict__ wavnum, double * __restrict__ em, 
  double * __restrict__ fm1, double * __restrict__ f1, double * __restrict__ ak, 
  double * __restrict__ xk, double delth, const double * __restrict__ dfim, 
  const double * __restrict__ dfimfr, const double * __restrict__ dfimofr, 
  double epsmin, const double * __restrict__ fr, double frtail, double g, int nang, 
  int nfre, double wetail, double wp1tail, double zpi, int ichnk, int nchnk, int ij, 
  double * __restrict__ tempa, double * __restrict__ tempx, double * __restrict__ temp2)
   {
  
  
  // ----------------------------------------------------------------------
  
  
  
  
  
  
  int m;
  int k;
  double zhook_handle;
  double delt25;
  double coefm1;
  double coef1;
  double coefa;
  double coefx;
  double sqrtk;
  
  // ----------------------------------------------------------------------
  
  
  
  //*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
  //        ------------------------------------------------
  
  
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
  
  //*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
  //        ------------------------------------------
  
  //*    2.2 SHALLOW WATER INTEGRATION.
  //         --------------------------
  
  for (m = 1; m <= nfre; m += 1) {
    sqrtk = sqrt((double) (wavnum[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))]));
    tempa[ij - 1 + kijl*(ichnk - 1)] = dfim[m - 1] / sqrtk;
    tempx[ij - 1 + kijl*(ichnk - 1)] = sqrtk*dfim[m - 1];
    k = 1;
    temp2[ij - 1 + kijl*(ichnk - 1)] = 
      fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
    for (k = 2; k <= nang; k += 1) {
      temp2[ij - 1 + kijl*(ichnk - 1)] = temp2[ij - 1 + kijl*(ichnk - 1)] + fl1[ij - 1 + 
        kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
    }
    em[ij - 1] = em[ij - 1] + dfim[m - 1]*temp2[ij - 1 + kijl*(ichnk - 1)];
    fm1[ij - 1] = fm1[ij - 1] + dfimofr[m - 1]*temp2[ij - 1 + kijl*(ichnk - 1)];
    f1[ij - 1] = f1[ij - 1] + dfimfr[m - 1]*temp2[ij - 1 + kijl*(ichnk - 1)];
    ak[ij - 1] = 
      ak[ij - 1] + tempa[ij - 1 + kijl*(ichnk - 1)]*temp2[ij - 1 + kijl*(ichnk - 1)];
    xk[ij - 1] = 
      xk[ij - 1] + tempx[ij - 1 + kijl*(ichnk - 1)]*temp2[ij - 1 + kijl*(ichnk - 1)];
  }
  
  //*      ADD TAIL CORRECTION TO MEAN FREQUENCY AND
  //*      NORMALIZE WITH TOTAL ENERGY.
  em[ij - 1] = em[ij - 1] + delt25*temp2[ij - 1 + kijl*(ichnk - 1)];
  fm1[ij - 1] = fm1[ij - 1] + coefm1*temp2[ij - 1 + kijl*(ichnk - 1)];
  fm1[ij - 1] = em[ij - 1] / fm1[ij - 1];
  f1[ij - 1] = f1[ij - 1] + coef1*temp2[ij - 1 + kijl*(ichnk - 1)];
  f1[ij - 1] = f1[ij - 1] / em[ij - 1];
  ak[ij - 1] = ak[ij - 1] + coefa*temp2[ij - 1 + kijl*(ichnk - 1)];
  ak[ij - 1] = pow((em[ij - 1] / ak[ij - 1]), 2);
  xk[ij - 1] = xk[ij - 1] + coefx*temp2[ij - 1 + kijl*(ichnk - 1)];
  xk[ij - 1] = pow((xk[ij - 1] / em[ij - 1]), 2);
  
  
  
}
