#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "peak_ang_c.h"


__device__ void peak_ang_c(int kijs, int kijl, const double * __restrict__ fl1, 
  double * __restrict__ xnu, double * __restrict__ sig_th, 
  const double * __restrict__ costh, double delth, const double * __restrict__ dfim, 
  const double * __restrict__ dfimfr, const double * __restrict__ dfimfr2, 
  const double * __restrict__ fr, double fratio, int nang, int nfre, 
  const double * __restrict__ sinth, const double * __restrict__ th, double wetail, 
  double wp1tail, double wp2tail, int ichnk, int nchnk, int ij, int * __restrict__ mmax, 
  double * __restrict__ sum0, double * __restrict__ sum1, double * __restrict__ sum2, 
  double * __restrict__ xmax, double * __restrict__ temp, double * __restrict__ thmean, 
  double * __restrict__ sum_s, double * __restrict__ sum_c) {
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  int nsh;
  int m;
  int k;
  int mmstart;
  int mmstop;
  double const_sig = (double) 1.0;
  double r1;
  double delt25;
  double coef_fr;
  double coef_fr2;
  double zhook_handle;
  double zepsilon;
  
  // ----------------------------------------------------------------------
  
  //***  1. DETERMINE L-H SPECTRAL WIDTH OF THE 2-D SPECTRUM.
  //     ---------------------------------------------------
  
  zepsilon = (double) 10.*DBL_EPSILON;
  
  nsh = 1 + (int) (log((double) 1.5) / log(fratio));
  
  
  sum0[ij - 1 + kijl*(ichnk - 1)] = zepsilon;
  sum1[ij - 1 + kijl*(ichnk - 1)] = (double) 0.;
  sum2[ij - 1 + kijl*(ichnk - 1)] = (double) 0.;
  
  for (m = 1; m <= nfre; m += 1) {
    k = 1;
    temp[ij - 1 + kijl*(ichnk - 1)] = 
      fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
    for (k = 2; k <= nang; k += 1) {
      temp[ij - 1 + kijl*(ichnk - 1)] = temp[ij - 1 + kijl*(ichnk - 1)] + fl1[ij - 1 + 
        kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
    }
    sum0[ij - 1 + kijl*(ichnk - 1)] = 
      sum0[ij - 1 + kijl*(ichnk - 1)] + temp[ij - 1 + kijl*(ichnk - 1)]*dfim[m - 1];
    sum1[ij - 1 + kijl*(ichnk - 1)] = 
      sum1[ij - 1 + kijl*(ichnk - 1)] + temp[ij - 1 + kijl*(ichnk - 1)]*dfimfr[m - 1];
    sum2[ij - 1 + kijl*(ichnk - 1)] = 
      sum2[ij - 1 + kijl*(ichnk - 1)] + temp[ij - 1 + kijl*(ichnk - 1)]*dfimfr2[m - 1];
  }
  
  //     ADD TAIL CORRECTIONS
  delt25 = wetail*fr[nfre - 1]*delth;
  coef_fr = wp1tail*delth*(pow(fr[nfre - 1], 2));
  coef_fr2 = wp2tail*delth*(pow(fr[nfre - 1], 3));
  sum0[ij - 1 + kijl*(ichnk - 1)] = 
    sum0[ij - 1 + kijl*(ichnk - 1)] + delt25*temp[ij - 1 + kijl*(ichnk - 1)];
  sum1[ij - 1 + kijl*(ichnk - 1)] = 
    sum1[ij - 1 + kijl*(ichnk - 1)] + coef_fr*temp[ij - 1 + kijl*(ichnk - 1)];
  sum2[ij - 1 + kijl*(ichnk - 1)] = 
    sum2[ij - 1 + kijl*(ichnk - 1)] + coef_fr2*temp[ij - 1 + kijl*(ichnk - 1)];
  
  if (sum0[ij - 1 + kijl*(ichnk - 1)] > zepsilon) {
    xnu[ij - 1] = sqrt((double) (fmax((double) (zepsilon), (double) (sum2[ij - 1 + 
      kijl*(ichnk - 1)]*sum0[ij - 1 + kijl*(ichnk - 1)] / (pow(sum1[ij - 1 + kijl*(ichnk 
      - 1)], 2)) - (double) 1.))));
  } else {
    xnu[ij - 1] = zepsilon;
  }
  
  //***  2. DETERMINE ANGULAR WIDTH OF THE 2-D SPECTRUM.
  //     ----------------------------------------------
  
  //     MAX OF 2-D SPECTRUM
  xmax[ij - 1 + kijl*(ichnk - 1)] = (double) 0.;
  mmax[ij - 1 + kijl*(ichnk - 1)] = 2;
  
  for (m = 2; m <= nfre - 1; m += 1) {
    for (k = 1; k <= nang; k += 1) {
      if (fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] > xmax[ij - 1 + 
        kijl*(ichnk - 1)]) {
        mmax[ij - 1 + kijl*(ichnk - 1)] = m;
        xmax[ij - 1 + kijl*(ichnk - 1)] = 
          fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
      }
    }
  }
  
  sum1[ij - 1 + kijl*(ichnk - 1)] = zepsilon;
  sum2[ij - 1 + kijl*(ichnk - 1)] = (double) 0.;
  
  mmstart = fmax((double) (1), (double) (mmax[ij - 1 + kijl*(ichnk - 1)] - nsh));
  mmstop = fmin((double) (nfre), (double) (mmax[ij - 1 + kijl*(ichnk - 1)] + nsh));
  
  sum_s[ij - 1 + kijl*(ichnk - 1)] = (double) 0.;
  sum_c[ij - 1 + kijl*(ichnk - 1)] = zepsilon;
  for (m = mmstart; m <= mmstop; m += 1) {
    for (k = 1; k <= nang; k += 1) {
      sum_s[ij - 1 + kijl*(ichnk - 1)] = sum_s[ij - 1 + kijl*(ichnk - 1)] + sinth[k - 
        1]*fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
      sum_c[ij - 1 + kijl*(ichnk - 1)] = sum_c[ij - 1 + kijl*(ichnk - 1)] + costh[k - 
        1]*fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
    }
    thmean[ij - 1 + kijl*(ichnk - 1)] = 
      atan2(sum_s[ij - 1 + kijl*(ichnk - 1)], sum_c[ij - 1 + kijl*(ichnk - 1)]);
    for (k = 1; k <= nang; k += 1) {
      sum1[ij - 1 + kijl*(ichnk - 1)] = sum1[ij - 1 + kijl*(ichnk - 1)] + fl1[ij - 1 + 
        kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))]*dfim[m - 1];
      sum2[ij - 1 + kijl*(ichnk - 1)] = sum2[ij - 1 + kijl*(ichnk - 1)] + cos(th[k - 1] -
         thmean[ij - 1 + kijl*(ichnk - 1)])*fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + 
        nfre*(ichnk - 1)))]*dfim[m - 1];
    }
  }
  
  if (sum1[ij - 1 + kijl*(ichnk - 1)] > zepsilon) {
    r1 = sum2[ij - 1 + kijl*(ichnk - 1)] / sum1[ij - 1 + kijl*(ichnk - 1)];
    sig_th[ij - 1] = const_sig*sqrt((double) ((double) 2.*((double) 1. - r1)));
  } else {
    sig_th[ij - 1] = (double) 0.;
  }
  
  
  
}
