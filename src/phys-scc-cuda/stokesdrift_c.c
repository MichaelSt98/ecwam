#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "stokesdrift_c.h"

__device__ void stokesdrift_c(int kijs, int kijl, const double * fl1, 
  const double * stokfac, const double * wswave, const double * wdwave, 
  const double * cicover, double * ustokes, double * vstokes, double cithrsh, 
  const double * costh, double delth, const double * dfim_sim, const double * fr, 
  double g, int licerun, int lwamrsetci, int nang, int nfre_odd, const double * sinth, 
  double zpi, int ichnk, int nchnk, int ij) {
  
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  
  int m;
  int k;
  
  double stmax = (double) 1.5;  // maximum magnitude (this is for safety when coupled)
  double const_var;
  double fac;
  double fac1;
  double fac2;
  double fac3;
  double stfac;
  const_var = (double) 2.0*delth*(pow(zpi, 3)) / g*(pow(fr[nfre_odd - 1], 4));

  ustokes[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  vstokes[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  
  for (m = 1; m <= nfre_odd; m += 1) {
    stfac = stokfac[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))]*dfim_sim[m - 1];
    for (k = 1; k <= nang; k += 1) {
      fac3 = stfac*fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
        nfre_loki_param*(ichnk - 1)))];
      ustokes[ij - 1 + kijl*(ichnk - 1)] = 
        ustokes[ij - 1 + kijl*(ichnk - 1)] + fac3*sinth[k - 1];
      vstokes[ij - 1 + kijl*(ichnk - 1)] = 
        vstokes[ij - 1 + kijl*(ichnk - 1)] + fac3*costh[k - 1];
    }
  }
  for (k = 1; k <= nang; k += 1) {
    fac1 = const_var*sinth[k - 1];
    fac2 = const_var*costh[k - 1];
    ustokes[ij - 1 + kijl*(ichnk - 1)] = ustokes[ij - 1 + kijl*(ichnk - 1)] + fac1*fl1[ij
       - 1 + kijl*(k - 1 + nang_loki_param*(nfre_odd - 1 + nfre_loki_param*(ichnk - 1)))]
      ;
    vstokes[ij - 1 + kijl*(ichnk - 1)] = vstokes[ij - 1 + kijl*(ichnk - 1)] + fac2*fl1[ij
       - 1 + kijl*(k - 1 + nang_loki_param*(nfre_odd - 1 + nfre_loki_param*(ichnk - 1)))]
      ;
  }
  if (licerun && lwamrsetci) {
    if (cicover[ij - 1 + kijl*(ichnk - 1)] > cithrsh) {
      ustokes[ij - 1 + kijl*(ichnk - 1)] = (double) 0.016*wswave[ij - 1 + kijl*(ichnk - 1
        )]*sin(wdwave[ij - 1 + kijl*(ichnk - 1)])*((double) 1.0 - cicover[ij - 1 + 
        kijl*(ichnk - 1)]);
      vstokes[ij - 1 + kijl*(ichnk - 1)] = (double) 0.016*wswave[ij - 1 + kijl*(ichnk - 1
        )]*cos(wdwave[ij - 1 + kijl*(ichnk - 1)])*((double) 1.0 - cicover[ij - 1 + 
        kijl*(ichnk - 1)]);
    }
  }
  ustokes[ij - 1 + kijl*(ichnk - 1)] = min((double) (max((double) (ustokes[ij - 1 + 
    kijl*(ichnk - 1)]), (double) (-stmax))), (double) (stmax));
  vstokes[ij - 1 + kijl*(ichnk - 1)] = min((double) (max((double) (vstokes[ij - 1 + 
    kijl*(ichnk - 1)]), (double) (-stmax))), (double) (stmax));

  
  
}
