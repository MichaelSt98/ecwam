#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "ciwabr_c.h"

__device__ void ciwabr_c(int kijs, int kijl, const double * cicover, const double * fl1, 
  const double * wavnum, const double * cgroup, double * ciwab, double cdicwa, 
  const double * dfim, double epsmin, int idelt, int licerun, int lmaskice, int nang, 
  int nfre, int ichnk, int nchnk, int ij) {
  
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  int k;
  int m;
  double ewh;
  double x;
  double alp;
  double xk2;
  

  if (!licerun || lmaskice) {
    
    for (m = 1; m <= nfre; m += 1) {
      for (k = 1; k <= nang; k += 1) {
        ciwab[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] = (double) 1.0;
      }
    }
    
  } else {
    
    for (m = 1; m <= nfre; m += 1) {
      for (k = 1; k <= nang; k += 1) {
        ewh = (double) 4.0*sqrt((double) (max((double) (epsmin), (double) (fl1[ij - 1 + 
          kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))]*dfim[m - 
          1]))));
        xk2 = pow(wavnum[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))], 2);
        alp = cdicwa*xk2*ewh;
        x = alp*cgroup[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))]*idelt;
        ciwab[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] = (double) 1.0 - 
          cicover[ij - 1 + kijl*(ichnk - 1)]*((double) 1.0 - exp((double) (-min((double) 
          (x), (double) ((double) 50.0)))));
      }
    }
    
  }

  
  
}
