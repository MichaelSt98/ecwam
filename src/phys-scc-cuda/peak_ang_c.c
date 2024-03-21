#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "peak_ang_c.h"

__device__ void peak_ang_c(int kijs, int kijl, const double * fl1, double * xnu, 
  double * sig_th, const double * costh, double delth, const double * dfim, 
  const double * dfimfr, const double * dfimfr2, const double * fr, double fratio, 
  int nang, int nfre, const double * sinth, const double * th, double wetail, 
  double wp1tail, double wp2tail, int ichnk, int nchnk, int ij) {
  
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  int nsh;
  int m;
  int k;
  int mmax;
  int mmstart;
  int mmstop;
  double const_sig = (double) 1.0;
  double r1;
  double delt25;
  double coef_fr;
  double coef_fr2;
  double zepsilon;
  double sum0;
  double sum1;
  double sum2;
  double xmax;
  double temp;
  double thmean;
  double sum_s;
  double sum_c;
  zepsilon = (double) 10.* DBL_EPSILON; // epsilon(zepsilon);
  nsh = 1 + (int) (log((double) 1.5) / log(fratio));
  

  sum0 = zepsilon;
  sum1 = (double) 0.;
  sum2 = (double) 0.;
  
  for (m = 1; m <= nfre; m += 1) {
    k = 1;
    temp = 
      fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))];
    for (k = 2; k <= nang; k += 1) {
      temp = temp + fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
        nfre_loki_param*(ichnk - 1)))];
    }
    sum0 = sum0 + temp*dfim[m - 1];
    sum1 = sum1 + temp*dfimfr[m - 1];
    sum2 = sum2 + temp*dfimfr2[m - 1];
  }
  delt25 = wetail*fr[nfre - 1]*delth;
  coef_fr = wp1tail*delth*(pow(fr[nfre - 1], 2));
  coef_fr2 = wp2tail*delth*(pow(fr[nfre - 1], 3));
  sum0 = sum0 + delt25*temp;
  sum1 = sum1 + coef_fr*temp;
  sum2 = sum2 + coef_fr2*temp;
  
  if (sum0 > zepsilon) {
    xnu[ij - 1] = sqrt((double) (max((double) (zepsilon), (double) (sum2*sum0 / 
      (pow(sum1, 2)) - (double) 1.))));
  } else {
    xnu[ij - 1] = zepsilon;
  }
  xmax = (double) 0.;
  mmax = 2;
  
  for (m = 2; m <= nfre - 1; m += 1) {
    for (k = 1; k <= nang; k += 1) {
      if (fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)
        ))] > xmax) {
        mmax = m;
        xmax = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk
           - 1)))];
      }
    }
  }
  
  sum1 = zepsilon;
  sum2 = (double) 0.;
  
  mmstart = max((double) (1), (double) (mmax - nsh));
  mmstop = min((double) (nfre), (double) (mmax + nsh));
  for (m = mmstart; m <= mmstop; m += 1) {
    sum_s = (double) 0.;
    sum_c = zepsilon;
    for (k = 1; k <= nang; k += 1) {
      sum_s = sum_s + sinth[k - 1]*fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
        nfre_loki_param*(ichnk - 1)))];
      sum_c = sum_c + costh[k - 1]*fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
        nfre_loki_param*(ichnk - 1)))];
    }
    thmean = atan2(sum_s, sum_c);
    for (k = 1; k <= nang; k += 1) {
      sum1 = sum1 + fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
        nfre_loki_param*(ichnk - 1)))]*dfim[m - 1];
      sum2 = sum2 + cos(th[k - 1] - thmean)*fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m
         - 1 + nfre_loki_param*(ichnk - 1)))]*dfim[m - 1];
    }
  }
  
  if (sum1 > zepsilon) {
    r1 = sum2 / sum1;
    sig_th[ij - 1] = const_sig*sqrt((double) ((double) 2.*((double) 1. - r1)));
  } else {
    sig_th[ij - 1] = (double) 0.;
  }

  
  
}
