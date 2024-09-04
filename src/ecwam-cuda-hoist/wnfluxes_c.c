#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "wnfluxes_c.h"


__device__ void wnfluxes_c(int kijs, int kijl, const int * __restrict__ mij, 
  const double * __restrict__ rhowgdfth, const double * __restrict__ cinv, 
  const double * __restrict__ ssurf, const double * __restrict__ cicover, 
  const double * __restrict__ phiwa, const double * __restrict__ em, 
  const double * __restrict__ f1, const double * __restrict__ wswave, 
  const double * __restrict__ wdwave, const double * __restrict__ ustra, 
  const double * __restrict__ vstra, const double * __restrict__ ufric, 
  const double * __restrict__ aird, double * __restrict__ nphieps, 
  double * __restrict__ ntauoc, double * __restrict__ nswh, double * __restrict__ nmwp, 
  double * __restrict__ nemotaux, double * __restrict__ nemotauy, 
  double * __restrict__ nemowswave, double * __restrict__ nemophif, 
  double * __restrict__ tauxd, double * __restrict__ tauyd, 
  double * __restrict__ tauocxd, double * __restrict__ tauocyd, 
  double * __restrict__ tauoc, double * __restrict__ phiocd, 
  double * __restrict__ phieps, double * __restrict__ phiaw, int lnupd, double afcrv, 
  double bfcrv, double ciblock, double cithrsh, const double * __restrict__ costh, 
  double egrcrv, double epsu10, double epsus, const double * __restrict__ fr, double g, 
  int licerun, int lwamrsetci, int lwcouast, int lwnemocou, int lwnemotauoc, int nang, 
  int nfre, double phiepsmax, double phiepsmin, const double * __restrict__ sinth, 
  double tauocmax, double tauocmin, int ichnk, int nchnk, int ij, 
  double * __restrict__ xstress, double * __restrict__ ystress, 
  double * __restrict__ ustar, double * __restrict__ philf, double * __restrict__ ooval, 
  double * __restrict__ em_oc, double * __restrict__ f1_oc, 
  double * __restrict__ cmrhowgdfth, double * __restrict__ sumt, 
  double * __restrict__ sumx, double * __restrict__ sumy) {
  
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  
  
  int k;
  int m;
  
  //     FICTITIOUS VALUE OF THE NORMALISED WAVE ENERGY FLUX UNDER THE SEA ICE
  //     (negative because it is defined as leaving the waves)
  double phioc_ice = -(double) 3.75;
  double phiaw_ice = (double) 3.75;
  
  //     USE HERSBACH 2011 FOR CD(U10) (SEE ALSO EDSON et al. 2013)
  double c1 = (double) 1.03E-3;
  double c2 = (double) 0.04E-3;
  double p1 = (double) 1.48;
  double p2 = -(double) 0.21;
  double cdmax = (double) 0.003;
  
  double efd_min = (double) 0.0625;  // corresponds to min Hs=1m under sea ice
  double efd_max = (double) 6.25;  // corresponds to max Hs=10m under sea ice
  
  double tau;
  double xn;
  double tauo;
  double u10p;
  double cd_bulk;
  double cd_wave;
  double cd_ice;
  double cnst;
  double epsus3;
  double cithrsh_inv;
  double efd;
  double ffd;
  double efd_fac;
  double ffd_fac;
  double zhook_handle;
  
  
  // ----------------------------------------------------------------------
  
  
  epsus3 = epsus*sqrt((double) (epsus));
  
  cithrsh_inv = (double) 1. / fmax((double) (cithrsh), (double) ((double) 0.01));
  
  efd_fac = (double) 4.0*egrcrv / (pow(g, 2));
  ffd_fac = (pow((egrcrv / afcrv), ((double) 1.0 / bfcrv)))*g;
  
  //*    DETERMINE NORMALIZED FLUXES FROM AIR TO WAVE AND FROM WAVE TO OCEAN.
  //     -------------------------------------------------------------------
  
  //     ENERGY FLUX from SSURF
  //     MOMENTUM FLUX FROM SSURF
  
  philf[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  xstress[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  ystress[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  
  //     THE INTEGRATION ONLY UP TO FR=MIJ
  for (m = 1; m <= nfre; m += 1) {
    k = 1;
    sumt[ij - 1 + kijl*(ichnk - 1)] = ssurf[ij - 1 + kijl*(k - 1 + nang*(m - 1))];
    sumx[ij - 1 + kijl*(ichnk - 1)] = 
      sinth[k - 1]*ssurf[ij - 1 + kijl*(k - 1 + nang*(m - 1))];
    sumy[ij - 1 + kijl*(ichnk - 1)] = 
      costh[k - 1]*ssurf[ij - 1 + kijl*(k - 1 + nang*(m - 1))];
    for (k = 2; k <= nang; k += 1) {
      sumt[ij - 1 + kijl*(ichnk - 1)] = 
        sumt[ij - 1 + kijl*(ichnk - 1)] + ssurf[ij - 1 + kijl*(k - 1 + nang*(m - 1))];
      sumx[ij - 1 + kijl*(ichnk - 1)] = sumx[ij - 1 + kijl*(ichnk - 1)] + sinth[k - 
        1]*ssurf[ij - 1 + kijl*(k - 1 + nang*(m - 1))];
      sumy[ij - 1 + kijl*(ichnk - 1)] = sumy[ij - 1 + kijl*(ichnk - 1)] + costh[k - 
        1]*ssurf[ij - 1 + kijl*(k - 1 + nang*(m - 1))];
    }
    philf[ij - 1 + kijl*(ichnk - 1)] = philf[ij - 1 + kijl*(ichnk - 1)] + sumt[ij - 1 + 
      kijl*(ichnk - 1)]*rhowgdfth[ij - 1 + kijl*(m - 1)];
    cmrhowgdfth[ij - 1 + kijl*(ichnk - 1)] = 
      cinv[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))]*rhowgdfth[ij - 1 + kijl*(m - 1)];
    xstress[ij - 1 + kijl*(ichnk - 1)] = xstress[ij - 1 + kijl*(ichnk - 1)] + sumx[ij - 1
       + kijl*(ichnk - 1)]*cmrhowgdfth[ij - 1 + kijl*(ichnk - 1)];
    ystress[ij - 1 + kijl*(ichnk - 1)] = ystress[ij - 1 + kijl*(ichnk - 1)] + sumy[ij - 1
       + kijl*(ichnk - 1)]*cmrhowgdfth[ij - 1 + kijl*(ichnk - 1)];
  }
  
  if (licerun && lwamrsetci) {
    if (cicover[ij - 1 + kijl*(ichnk - 1)] > ciblock) {
      ooval[ij - 1 + kijl*(ichnk - 1)] = exp((double) (-fmin((double) (pow((cicover[ij - 
        1 + kijl*(ichnk - 1)]*cithrsh_inv), 4)), (double) ((double) 10.))));
      //           ADJUST USTAR FOR THE PRESENCE OF SEA ICE
      u10p = fmax((double) (wswave[ij - 1 + kijl*(ichnk - 1)]), (double) (epsu10));
      cd_bulk = 
        fmin((double) ((c1 + c2*(pow(u10p, p1)))*(pow(u10p, p2))), (double) (cdmax));
      cd_wave = pow((ufric[ij - 1 + kijl*(ichnk - 1)] / u10p), 2);
      cd_ice = ooval[ij - 1 + kijl*(ichnk - 1)]*cd_wave + ((double) 1.0 - ooval[ij - 1 + 
        kijl*(ichnk - 1)])*cd_bulk;
      ustar[ij - 1 + kijl*(ichnk - 1)] = 
        fmax((double) (sqrt((double) (cd_ice))*u10p), (double) (epsus));
      
      // EM_OC and F1_OC with fully developed model ENERGY
      // The significant wave height derived from EM_OC will be used
      // by NEMO as a scaling factor as if it was open ocean
      efd = fmin((double) (efd_fac*(pow(ustar[ij - 1 + kijl*(ichnk - 1)], 4))), (double) 
        (efd_max));
      em_oc[ij - 1 + kijl*(ichnk - 1)] = fmax((double) (ooval[ij - 1 + kijl*(ichnk - 1)
        ]*em[ij - 1] + ((double) 1.0 - ooval[ij - 1 + kijl*(ichnk - 1)])*efd), (double) 
        (efd_min));
      ffd = ffd_fac / ustar[ij - 1 + kijl*(ichnk - 1)];
      f1_oc[ij - 1 + kijl*(ichnk - 1)] = ooval[ij - 1 + kijl*(ichnk - 1)]*f1[ij - 1] + 
        ((double) 1.0 - ooval[ij - 1 + kijl*(ichnk - 1)])*ffd;
      f1_oc[ij - 1 + kijl*(ichnk - 1)] = fmin((double) (fmax((double) (f1_oc[ij - 1 + 
        kijl*(ichnk - 1)]), (double) (fr[2 - 1]))), (double) (fr[nfre - 1]));
    } else {
      ooval[ij - 1 + kijl*(ichnk - 1)] = (double) 1.0;
      ustar[ij - 1 + kijl*(ichnk - 1)] = ufric[ij - 1 + kijl*(ichnk - 1)];
      em_oc[ij - 1 + kijl*(ichnk - 1)] = em[ij - 1];
      f1_oc[ij - 1 + kijl*(ichnk - 1)] = f1[ij - 1];
    }
  } else {
    ooval[ij - 1 + kijl*(ichnk - 1)] = (double) 1.0;
    ustar[ij - 1 + kijl*(ichnk - 1)] = ufric[ij - 1 + kijl*(ichnk - 1)];
    em_oc[ij - 1 + kijl*(ichnk - 1)] = em[ij - 1];
    f1_oc[ij - 1 + kijl*(ichnk - 1)] = f1[ij - 1];
  }
  
  tau = aird[ij - 1 + kijl*(ichnk - 1)]*fmax((double) (pow(ustar[ij - 1 + kijl*(ichnk - 1
    )], 2)), (double) (epsus));
  tauxd[ij - 1 + kijl*(ichnk - 1)] = tau*sin(wdwave[ij - 1 + kijl*(ichnk - 1)]);
  tauyd[ij - 1 + kijl*(ichnk - 1)] = tau*cos(wdwave[ij - 1 + kijl*(ichnk - 1)]);
  
  tauocxd[ij - 1 + kijl*(ichnk - 1)] = tauxd[ij - 1 + kijl*(ichnk - 1)] - ooval[ij - 1 + 
    kijl*(ichnk - 1)]*xstress[ij - 1 + kijl*(ichnk - 1)];
  tauocyd[ij - 1 + kijl*(ichnk - 1)] = tauyd[ij - 1 + kijl*(ichnk - 1)] - ooval[ij - 1 + 
    kijl*(ichnk - 1)]*ystress[ij - 1 + kijl*(ichnk - 1)];
  tauo = sqrt((double) ((pow(tauocxd[ij - 1 + kijl*(ichnk - 1)], 2)) + (pow(tauocyd[ij - 
    1 + kijl*(ichnk - 1)], 2))));
  tauoc[ij - 1 + kijl*(ichnk - 1)] = fmin((double) (fmax((double) (tauo / tau), (double) 
    (tauocmin))), (double) (tauocmax));
  
  if (lwcouast) {
    //       USE THE SURFACE STRESSES AS PROVIDED BY THE ATMOSPHERE AND
    //       ONLY RESCALE THEM WITH TAUOC
    if (ustra[ij - 1 + kijl*(ichnk - 1)] != (double) 0.0 || vstra[ij - 1 + kijl*(ichnk - 
      1)] != (double) 0.0) {
      tauxd[ij - 1 + kijl*(ichnk - 1)] = ustra[ij - 1 + kijl*(ichnk - 1)];
      tauocxd[ij - 1 + kijl*(ichnk - 1)] = 
        ustra[ij - 1 + kijl*(ichnk - 1)]*tauoc[ij - 1 + kijl*(ichnk - 1)];
      tauyd[ij - 1 + kijl*(ichnk - 1)] = vstra[ij - 1 + kijl*(ichnk - 1)];
      tauocyd[ij - 1 + kijl*(ichnk - 1)] = 
        vstra[ij - 1 + kijl*(ichnk - 1)]*tauoc[ij - 1 + kijl*(ichnk - 1)];
    }
  }
  
  
  xn = aird[ij - 1 + kijl*(ichnk - 1)]*fmax((double) (pow(ustar[ij - 1 + kijl*(ichnk - 1)
    ], 3)), (double) (epsus3));
  phiocd[ij - 1 + kijl*(ichnk - 1)] = ooval[ij - 1 + kijl*(ichnk - 1)]*(philf[ij - 1 + 
    kijl*(ichnk - 1)] - phiwa[ij - 1]) + ((double) 1.0 - ooval[ij - 1 + kijl*(ichnk - 1)]
    )*phioc_ice*xn;
  
  phieps[ij - 1 + kijl*(ichnk - 1)] = phiocd[ij - 1 + kijl*(ichnk - 1)] / xn;
  phieps[ij - 1 + kijl*(ichnk - 1)] = fmin((double) (fmax((double) (phieps[ij - 1 + 
    kijl*(ichnk - 1)]), (double) (phiepsmin))), (double) (phiepsmax));
  
  phiocd[ij - 1 + kijl*(ichnk - 1)] = phieps[ij - 1 + kijl*(ichnk - 1)]*xn;
  
  phiaw[ij - 1 + kijl*(ichnk - 1)] = phiwa[ij - 1] / xn;
  phiaw[ij - 1 + kijl*(ichnk - 1)] = ooval[ij - 1 + kijl*(ichnk - 1)]*phiwa[ij - 1] / xn 
    + ((double) 1.0 - ooval[ij - 1 + kijl*(ichnk - 1)])*phiaw_ice;
  
  if (lwnemocou && lnupd) {
    nphieps[ij - 1 + kijl*(ichnk - 1)] = phieps[ij - 1 + kijl*(ichnk - 1)];
    ntauoc[ij - 1 + kijl*(ichnk - 1)] = tauoc[ij - 1 + kijl*(ichnk - 1)];
    if (em_oc[ij - 1 + kijl*(ichnk - 1)] != (double) 0.0) {
      nswh[ij - 1 + kijl*(ichnk - 1)] = 
        (double) 4.0*sqrt((double) (em_oc[ij - 1 + kijl*(ichnk - 1)]));
    } else {
      nswh[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
    }
    if (f1_oc[ij - 1 + kijl*(ichnk - 1)] != (double) 0.0) {
      nmwp[ij - 1 + kijl*(ichnk - 1)] = (double) 1.0 / f1_oc[ij - 1 + kijl*(ichnk - 1)];
    } else {
      nmwp[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
    }
    
    if (lwnemotauoc) {
      nemotaux[ij - 1 + kijl*(ichnk - 1)] = 
        nemotaux[ij - 1 + kijl*(ichnk - 1)] + tauocxd[ij - 1 + kijl*(ichnk - 1)];
      nemotauy[ij - 1 + kijl*(ichnk - 1)] = 
        nemotauy[ij - 1 + kijl*(ichnk - 1)] + tauocyd[ij - 1 + kijl*(ichnk - 1)];
    } else {
      nemotaux[ij - 1 + kijl*(ichnk - 1)] = 
        nemotaux[ij - 1 + kijl*(ichnk - 1)] + tauxd[ij - 1 + kijl*(ichnk - 1)];
      nemotauy[ij - 1 + kijl*(ichnk - 1)] = 
        nemotauy[ij - 1 + kijl*(ichnk - 1)] + tauyd[ij - 1 + kijl*(ichnk - 1)];
    }
    nemowswave[ij - 1 + kijl*(ichnk - 1)] = 
      nemowswave[ij - 1 + kijl*(ichnk - 1)] + wswave[ij - 1 + kijl*(ichnk - 1)];
    nemophif[ij - 1 + kijl*(ichnk - 1)] = 
      nemophif[ij - 1 + kijl*(ichnk - 1)] + phiocd[ij - 1 + kijl*(ichnk - 1)];
  }
  
  
  // ----------------------------------------------------------------------
  
}
