#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sdissip_ard_c.h"

__device__ void sdissip_ard_c(int kijs, int kijl, const double * fl1, double * fld, 
  double * sl, const int * indep, const double * wavnum, const double * xk2cg, 
  const double * ufric, const double * coswdif, const double * raorw, 
  const double * cumulw, double g, const int * indicessat, int ipsat, double miche, 
  int nang, int ndepth, int ndikcumul, int nfre, int nsdsnth, const double * satweights, 
  double sdsbr, double ssdsc2, double ssdsc3, double ssdsc4, double ssdsc5, 
  double ssdsc6, double zpi, const double * zpifr, int ichnk, int nchnk, int ij) {
  
  // Loki: parameters from YOWPARAM inlined
  // Loki: parameters from YOWPHYS inlined
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  
  int k;
  int m;
  int i;
  int j;
  int m2;
  int k2;
  int kk;
  
  double tpiinv;
  double tpiinvh;
  double tmp01;
  double tmp03;
  double epsr;
  double ssdsc6m1;
  double zcoef;
  double zcoefm1;
  double ssdsc2_sig;
  double facturb;
  double bth;
  double bth0;
  double scumul[36];
  double d[36];
  
  double renewalfreq;
  int foo;
  foo = ndepth;    // necessary for Loki ...
  epsr = sqrt((double) (sdsbr));
  
  tpiinv = (double) 1.0 / zpi;
  tpiinvh = (double) 0.5*tpiinv;
  tmp03 = (double) 1.0 / (sdsbr*miche);
  ssdsc6m1 = (double) 1. - ssdsc6;
  

  for (m = 1; m <= nfre; m += 1) {
    ssdsc2_sig = ssdsc2*zpifr[m - 1];
    zcoef = ssdsc2_sig*ssdsc6;
    zcoefm1 = ssdsc2_sig*ssdsc6m1;
    bth0 = (double) 0.0;
    
    for (k = 1; k <= nang; k += 1) {
      bth = (double) 0.0;
      // integrates in directional sector
      for (k2 = 1; k2 <= nsdsnth*2 + 1; k2 += 1) {
        kk = indicessat[k - 1 + nang_loki_param*(k2 - 1)];
        bth = bth + satweights[k - 1 + nang_loki_param*(k2 - 1)]*fl1[ij - 1 + kijl*(kk - 
          1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))];
      }
      bth = bth*wavnum[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))
        ]*tpiinv*xk2cg[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))];
      bth0 = max((double) (bth0), (double) (bth));
      
      d[k - 1] = 
        zcoefm1*(pow(max((double) ((double) 0.), (double) (bth*tmp03 - ssdsc4)), ipsat));
      
      scumul[k - 1] = pow(max((double) (sqrt((double) (abs((double) (bth)))) - epsr), 
        (double) ((double) 0.)), 2);
    }
    
    for (k = 1; k <= nang; k += 1) {
      // cumulative term
      d[k - 1] = d[k - 1] + zcoef*(pow(max((double) ((double) 0.), (double) (bth0*tmp03 -
         ssdsc4)), ipsat));
      if (bth0 <= sdsbr) {
        scumul[k - 1] = (double) 0.;
      }
      
    }
    
    if (m > ndikcumul) {
      // CUMULATIVE TERM
      if (ssdsc3 != (double) 0.0) {
        
        for (k = 1; k <= nang; k += 1) {
          renewalfreq = (double) 0.0;
          
          for (m2 = 1; m2 <= m - ndikcumul; m2 += 1) {
            for (k2 = 1; k2 <= nang; k2 += 1) {
              kk = abs((double) (k2 - k));
              if (kk > nang / 2) {
                kk = kk - nang / 2;
              }
              renewalfreq = renewalfreq + cumulw[indep[ij - 1 + kijl*(ichnk - 1)] - 1 + 
                ndepth*(1 + kk - 1 + (1 + nang / 2)*(m2 - 1 + nfre_loki_param*(m - 1)))
                ]*scumul[k2 - 1];
            }
          }
          
          d[k - 1] = d[k - 1] + renewalfreq;
        }
      }
    }
    if (ssdsc5 != (double) 0.0) {
      tmp01 = (double) 2.*ssdsc5 / g;
      facturb = tmp01*raorw[ij - 1]*ufric[ij - 1 + kijl*(ichnk - 1)]*ufric[ij - 1 + 
        kijl*(ichnk - 1)];
      for (k = 1; k <= nang; k += 1) {
        d[k - 1] = d[k - 1] - zpifr[m - 1]*wavnum[ij - 1 + kijl*(m - 1 + 
          nfre_loki_param*(ichnk - 1))]*facturb*coswdif[ij - 1 + kijl*(k - 1)];
      }
    }
    for (k = 1; k <= nang; k += 1) {
      sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] = sl[ij - 1 + kijl*(k - 1 + 
        nang_loki_param*(m - 1))] + d[k - 1]*fl1[ij - 1 + kijl*(k - 1 + 
        nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))];
      fld[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] = 
        fld[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] + d[k - 1];
    }
  }

  
  
}
