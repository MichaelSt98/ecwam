#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sdissip_ard_c.h"


__device__ void sdissip_ard_c(int kijs, int kijl, const double * __restrict__ fl1, 
  double * __restrict__ fld, double * __restrict__ sl, 
  const double * __restrict__ wavnum, const double * __restrict__ cgroup, 
  const double * __restrict__ xk2cg, const double * __restrict__ ufric, 
  const double * __restrict__ coswdif, const double * __restrict__ raorw, 
  double brkpbcoef, double delth, double fratio, double g, 
  const int * __restrict__ indicessat, int ipsat, double miche, int nang, int nfre, 
  int nsdsnth, const double * __restrict__ satweights, double sdsbr, double ssdsbrf1, 
  double ssdsc2, double ssdsc3, double ssdsc4, double ssdsc5, double ssdsc6, double zpi, 
  const double * __restrict__ zpifr, int ichnk, int nchnk, int ij, 
  double * __restrict__ facturb, double * __restrict__ facsat, 
  double * __restrict__ facwtrb, double * __restrict__ temp1, 
  double * __restrict__ bth0, double * __restrict__ c_, double * __restrict__ c_c, 
  double * __restrict__ dsip, double * __restrict__ trpz_dsip, 
  double * __restrict__ bth, double * __restrict__ temp2, double * __restrict__ d, 
  double * __restrict__ scumul, double * __restrict__ renewalfreq, 
  double * __restrict__ wcumul) {
  
  // Loki: parameters from YOWPHYS inlined
  
  
  // ----------------------------------------------------------------------
  
  
  
  
  
  int k;
  int m;
  int i;
  int j;
  int m2;
  int k2;
  int kk;
  int nangd;
  //     NDIKCUMUL is the  integer difference in frequency bands
  int ndikcumul;
  
  double tpiinv;
  double tpiinvh;
  double tmp01;
  double tmp02;
  double tmp03;
  double epsr;
  double ssdsc6m1;
  double zcoef;
  double zcoefm1;
  double xlogdfrth;
  double brlambda;
  
  double zhook_handle;
  
  
  double ssdsc2_sig;
  
  // ----------------------------------------------------------------------
  
  
  // INITIALISATION
  
  epsr = sqrt((double) (sdsbr));
  
  tpiinv = (double) 1.0 / zpi;
  tpiinvh = (double) 0.5*tpiinv;
  
  
  tmp03 = (double) 1.0 / (sdsbr*miche);
  
  ssdsc6m1 = (double) 1. - ssdsc6;
  
  
  for (m = 1; m <= nfre; m += 1) {
    facsat[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))] = wavnum[ij - 1 + kijl*(m - 1 + 
      nfre*(ichnk - 1))]*tpiinv*xk2cg[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))];
  }
  
  // COMPUTE SATURATION SPECTRUM
  for (m = 1; m <= nfre; m += 1) {
    bth0[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))] = (double) 0.0;
    for (k = 1; k <= nang; k += 1) {
      bth[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = (double) 0.0;
    }
  }
  
  for (m = 1; m <= nfre; m += 1) {
    for (k = 1; k <= nang; k += 1) {
      // integrates in directional sector
      for (k2 = 1; k2 <= nsdsnth*2 + 1; k2 += 1) {
        kk = indicessat[k - 1 + nang*(k2 - 1)];
        bth[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = bth[ij - 1 + 
          kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] + satweights[k - 1 + nang*(k2 -
           1)]*fl1[ij - 1 + kijl*(kk - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
      }
      bth[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = bth[ij - 1 + kijl*(k
         - 1 + nang*(m - 1 + nfre*(ichnk - 1)))]*facsat[ij - 1 + kijl*(m - 1 + 
        nfre*(ichnk - 1))];
      bth0[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))] = fmax((double) (bth0[ij - 1 + 
        kijl*(m - 1 + nfre*(ichnk - 1))]), (double) (bth[ij - 1 + kijl*(k - 1 + nang*(m -
         1 + nfre*(ichnk - 1)))]));
    }
  }
  
  
  // SATURATION TERM
  
  for (m = 1; m <= nfre; m += 1) {
    ssdsc2_sig = ssdsc2*zpifr[m - 1];
    zcoef = ssdsc2_sig*ssdsc6;
    zcoefm1 = ssdsc2_sig*ssdsc6m1;
    for (k = 1; k <= nang; k += 1) {
      d[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = 
        zcoef*(pow(fmax((double) ((double) 0.), (double) (bth0[ij - 1 + kijl*(m - 1 + 
        nfre*(ichnk - 1))]*tmp03 - ssdsc4)), ipsat)) + zcoefm1*(pow(fmax((double) 
        ((double) 0.), (double) (bth[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1
        )))]*tmp03 - ssdsc4)), ipsat));
    }
  }
  
  
  // CUMULATIVE TERM
  if (ssdsc3 != (double) 0.0) {
    
    nangd = nang / 2;
    xlogdfrth = log(fratio)*delth;
    //       l(k,th)=1/(2*pi²)= the breaking crest density
    brlambda = brkpbcoef / ((double) 2.0*(pow(zpi, 2)));
    tmp02 = ssdsc3*brlambda;
    
    //       NDIKCUMUL is the  integer difference in frequency bands
    //       between the "large breakers" and short "wiped-out waves"
    //!! wrong !!???        NDIKCUMUL = NINT(SSDSBRF1/(FRATIO-1.))
    ndikcumul = rint((double) (-log(ssdsbrf1) / log(fratio)));
    
    for (m = 1; m <= nfre; m += 1) {
      c_[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))] = 
        zpifr[m - 1] / wavnum[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))];
      c_c[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))] = 
        pow(c_[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))], 2);
      dsip[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))] = 
        tmp02*zpifr[m - 1]*xlogdfrth / cgroup[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))];        //  coef*dtheta*dk = coef*dtheta*dsigma/cg
    }
    
    for (m2 = 1; m2 <= nfre - ndikcumul; m2 += 1) {
      if (bth0[ij - 1 + kijl*(m2 - 1 + nfre*(ichnk - 1))] > sdsbr) {
        temp1[ij - 1 + kijl*(m2 - 1 + nfre*(ichnk - 1))] = (double) 1.0;
      } else {
        temp1[ij - 1 + kijl*(m2 - 1 + nfre*(ichnk - 1))] = (double) 0.0;
      }
    }
    for (m2 = 1; m2 <= nfre - ndikcumul; m2 += 1) {
      for (k2 = 1; k2 <= nang; k2 += 1) {
        scumul[ij - 1 + kijl*(k2 - 1 + nang*(m2 - 1 + nfre*(ichnk - 1)))] = temp1[ij - 1 
          + kijl*(m2 - 1 + nfre*(ichnk - 1))]*(pow(fmax((double) (sqrt((double) (bth[ij -
           1 + kijl*(k2 - 1 + nang*(m2 - 1 + nfre*(ichnk - 1)))])) - epsr), (double) 
          ((double) 0.0)), 2));
      }
    }
    
    for (m = ndikcumul + 1; m <= nfre; m += 1) {
      for (k = 1; k <= nang; k += 1) {
        renewalfreq[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = 
          (double) 0.0;
      }
    }
    
    
    for (m = ndikcumul + 1; m <= nfre; m += 1) {
      
      if (m - ndikcumul >= 3) {
        trpz_dsip[ij - 1 + kijl*(1 - 1 + nfre*(ichnk - 1))] = 
          (double) 0.5*dsip[ij - 1 + kijl*(1 - 1 + nfre*(ichnk - 1))];
        for (m2 = 2; m2 <= m - ndikcumul - 1; m2 += 1) {
          trpz_dsip[ij - 1 + kijl*(m2 - 1 + nfre*(ichnk - 1))] = 
            dsip[ij - 1 + kijl*(m2 - 1 + nfre*(ichnk - 1))];
        }
        trpz_dsip[ij - 1 + kijl*(m - ndikcumul - 1 + nfre*(ichnk - 1))] = 
          (double) 0.5*dsip[ij - 1 + kijl*(m - ndikcumul - 1 + nfre*(ichnk - 1))];
      } else {
        for (m2 = 1; m2 <= m - ndikcumul; m2 += 1) {
          trpz_dsip[ij - 1 + kijl*(m2 - 1 + nfre*(ichnk - 1))] = 
            dsip[ij - 1 + kijl*(m2 - 1 + nfre*(ichnk - 1))];
        }
      }
      
      for (m2 = 1; m2 <= m - ndikcumul; m2 += 1) {
        for (kk = 0; kk <= nangd; kk += 1) {
          wcumul[ij - 1 + kijl*(1 + kk - 1 + (1 + nang / 2)*(m2 - 1 + nfre*(ichnk - 1)))]
             = sqrt((double) (fabs((double) (c_c[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))
            ] + c_c[ij - 1 + kijl*(m2 - 1 + nfre*(ichnk - 1))] - (double) 2.0*c_[ij - 1 +
             kijl*(m - 1 + nfre*(ichnk - 1))]*c_[ij - 1 + kijl*(m2 - 1 + nfre*(ichnk - 1)
            )]*cos(kk*delth)))))*trpz_dsip[ij - 1 + kijl*(m2 - 1 + nfre*(ichnk - 1))];
        }
      }
      
      for (k = 1; k <= nang; k += 1) {
        // Correction of saturation level for shallow-water kinematics
        // Cumulative effect based on lambda   (breaking probability is
        // the expected rate of sweeping by larger breaking waves)
        
        for (k2 = 1; k2 <= nang; k2 += 1) {
        }
        
        for (m2 = 1; m2 <= m - ndikcumul; m2 += 1) {
          for (k2 = 1; k2 <= nang; k2 += 1) {
            kk = fabs((double) (k2 - k));
            if (kk > nangd) {
              kk = kk - nangd;
            }
            // Integrates over frequencies M2 and directions K2 to
            // Integration is performed from M2=1 to a frequency lower than M: IK-NDIKCUMUL
            renewalfreq[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = 
              renewalfreq[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] + 
              wcumul[ij - 1 + kijl*(1 + kk - 1 + (1 + nang / 2)*(m2 - 1 + nfre*(ichnk - 1
              )))]*scumul[ij - 1 + kijl*(k2 - 1 + nang*(m2 - 1 + nfre*(ichnk - 1)))];
          }
        }
      }
    }
    
    for (m = ndikcumul + 1; m <= nfre; m += 1) {
      for (k = 1; k <= nang; k += 1) {
        d[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = d[ij - 1 + kijl*(k -
           1 + nang*(m - 1 + nfre*(ichnk - 1)))] + renewalfreq[ij - 1 + kijl*(k - 1 + 
          nang*(m - 1 + nfre*(ichnk - 1)))];
      }
    }
  }
  
  
  //     WAVE-TURBULENCE INTERACTION TERM
  if (ssdsc5 != (double) 0.0) {
    tmp01 = (double) 2.*ssdsc5 / g;
    facturb[ij - 1 + kijl*(ichnk - 1)] = tmp01*raorw[ij - 1]*ufric[ij - 1 + kijl*(ichnk -
       1)]*ufric[ij - 1 + kijl*(ichnk - 1)];
    for (m = 1; m <= nfre; m += 1) {
      facwtrb[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))] = zpifr[m - 1]*wavnum[ij - 1 + 
        kijl*(m - 1 + nfre*(ichnk - 1))]*facturb[ij - 1 + kijl*(ichnk - 1)];
      for (k = 1; k <= nang; k += 1) {
        d[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = d[ij - 1 + kijl*(k -
           1 + nang*(m - 1 + nfre*(ichnk - 1)))] - facwtrb[ij - 1 + kijl*(m - 1 + 
          nfre*(ichnk - 1))]*coswdif[ij - 1 + kijl*(k - 1)];
      }
    }
  }
  
  
  // ADD ALL CONTRIBUTIONS TO SOURCE TERM
  for (m = 1; m <= nfre; m += 1) {
    for (k = 1; k <= nang; k += 1) {
      sl[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = sl[ij - 1 + kijl*(k - 1 + nang*(m - 1))]
         + d[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))]*fl1[ij - 1 + 
        kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
      fld[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = fld[ij - 1 + kijl*(k - 1 + nang*(m - 1)
        )] + d[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
    }
  }
  
  
  
}
