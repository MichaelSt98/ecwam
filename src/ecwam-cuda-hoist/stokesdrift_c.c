#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "stokesdrift_c.h"


__device__ void stokesdrift_c(int kijs, int kijl, const double * __restrict__ fl1, 
  const double * __restrict__ stokfac, const double * __restrict__ wswave, 
  const double * __restrict__ wdwave, const double * __restrict__ cicover, 
  double * __restrict__ ustokes, double * __restrict__ vstokes, double cithrsh, 
  const double * __restrict__ costh, double delth, const double * __restrict__ dfim_sim, 
  const double * __restrict__ fr, double g, int licerun, int lwamrsetci, int nang, 
  int nfre, int nfre_odd, const double * __restrict__ sinth, double zpi, int ichnk, 
  int nchnk, int ij, double * __restrict__ stfac) {
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  
  int m;
  int k;
  
  double stmax = (double) 1.5;  // maximum magnitude (this is for safety when coupled)
  double const_var;
  double fac;
  double fac1;
  double fac2;
  double fac3;
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  
  
  //***  1. DETERMINE STOKE DRIFT VECTOR.
  //     --------------------------------
  
  const_var = (double) 2.0*delth*(pow(zpi, 3)) / g*(pow(fr[nfre_odd - 1], 4));
  
  //***  1.1 PERFORM INTEGRATION.
  //     ------------------------
  
  
  ustokes[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  vstokes[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  
  for (m = 1; m <= nfre_odd; m += 1) {
    stfac[ij - 1 + kijl*(ichnk - 1)] = 
      stokfac[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))]*dfim_sim[m - 1];
    for (k = 1; k <= nang; k += 1) {
      fac3 = stfac[ij - 1 + kijl*(ichnk - 1)]*fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + 
        nfre*(ichnk - 1)))];
      ustokes[ij - 1 + kijl*(ichnk - 1)] = 
        ustokes[ij - 1 + kijl*(ichnk - 1)] + fac3*sinth[k - 1];
      vstokes[ij - 1 + kijl*(ichnk - 1)] = 
        vstokes[ij - 1 + kijl*(ichnk - 1)] + fac3*costh[k - 1];
    }
  }
  
  //***  1.2 ADD CONTRIBUTION OF UNRESOLVED WAVES.
  //     -----------------------------------------
  
  for (k = 1; k <= nang; k += 1) {
    fac1 = const_var*sinth[k - 1];
    fac2 = const_var*costh[k - 1];
    ustokes[ij - 1 + kijl*(ichnk - 1)] = ustokes[ij - 1 + kijl*(ichnk - 1)] + fac1*fl1[ij
       - 1 + kijl*(k - 1 + nang*(nfre_odd - 1 + nfre*(ichnk - 1)))];
    vstokes[ij - 1 + kijl*(ichnk - 1)] = vstokes[ij - 1 + kijl*(ichnk - 1)] + fac2*fl1[ij
       - 1 + kijl*(k - 1 + nang*(nfre_odd - 1 + nfre*(ichnk - 1)))];
  }
  
  
  //***  1.3 Sea Ice exception
  //     ---------------------
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
  
  //***  1.4 Protection
  //     --------------
  
  ustokes[ij - 1 + kijl*(ichnk - 1)] = fmin((double) (fmax((double) (ustokes[ij - 1 + 
    kijl*(ichnk - 1)]), (double) (-stmax))), (double) (stmax));
  vstokes[ij - 1 + kijl*(ichnk - 1)] = fmin((double) (fmax((double) (vstokes[ij - 1 + 
    kijl*(ichnk - 1)]), (double) (-stmax))), (double) (stmax));
  
  
  
}
