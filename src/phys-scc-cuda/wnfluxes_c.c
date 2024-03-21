#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "wnfluxes_c.h"

__device__ void wnfluxes_c(int kijs, int kijl, const int * mij, 
  const double * rhowgdfth, const double * cinv, const double * ssurf, 
  const double * cicover, const double * phiwa, const double * em, const double * f1, 
  const double * wswave, const double * wdwave, const double * ufric, 
  const double * aird, double * nphieps, double * ntauoc, double * nswh, double * nmwp, 
  double * nemotaux, double * nemotauy, double * nemowswave, double * nemophif, 
  double * tauxd, double * tauyd, double * tauocxd, double * tauocyd, double * tauoc, 
  double * phiocd, double * phieps, double * phiaw, int lnupd, double afcrv, 
  double bfcrv, double ciblock, double cithrsh, const double * costh, double egrcrv, 
  double epsu10, double epsus, const double * fr, double g, int licerun, int lwamrsetci, 
  int lwnemocou, int lwnemotauoc, int nang, int nfre, double phiepsmax, 
  double phiepsmin, const double * sinth, double tauocmax, double tauocmin, int ichnk, 
  int nchnk, int ij) {
  
  
  
  
  const int nfre_loki_param = 36;
  const int nang_loki_param = 24;
  
  int k;
  int m;
  double phioc_ice = -(double) 3.75;
  double phiaw_ice = (double) 3.75;
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
  
  double xstress;
  double ystress;
  double ustar;
  double philf;
  double ooval;
  double em_oc;
  double f1_oc;
  double cmrhowgdfth;
  double sumt;
  double sumx;
  double sumy;
  
  epsus3 = epsus*sqrt((double) (epsus));
  
  cithrsh_inv = (double) 1. / max((double) (cithrsh), (double) ((double) 0.01));
  
  efd_fac = (double) 4.0*egrcrv / (pow(g, 2));
  ffd_fac = (pow((egrcrv / afcrv), ((double) 1.0 / bfcrv)))*g;

  philf = (double) 0.0;
  xstress = (double) 0.0;
  ystress = (double) 0.0;
  for (m = 1; m <= nfre; m += 1) {
    k = 1;
    sumt = ssurf[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))];
    sumx = sinth[k - 1]*ssurf[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))];
    sumy = costh[k - 1]*ssurf[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))];
    for (k = 2; k <= nang; k += 1) {
      sumt = sumt + ssurf[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))];
      sumx = sumx + sinth[k - 1]*ssurf[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))];
      sumy = sumy + costh[k - 1]*ssurf[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))];
    }
    philf = philf + sumt*rhowgdfth[ij - 1 + kijl*(m - 1)];
    cmrhowgdfth = cinv[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))]*rhowgdfth[ij 
      - 1 + kijl*(m - 1)];
    xstress = xstress + sumx*cmrhowgdfth;
    ystress = ystress + sumy*cmrhowgdfth;
  }
  
  if (licerun && lwamrsetci) {
    if (cicover[ij - 1 + kijl*(ichnk - 1)] > ciblock) {
      ooval = exp((double) (-min((double) (pow((cicover[ij - 1 + kijl*(ichnk - 1)
        ]*cithrsh_inv), 4)), (double) ((double) 10.))));
      //           ADJUST USTAR FOR THE PRESENCE OF SEA ICE
      u10p = max((double) (wswave[ij - 1 + kijl*(ichnk - 1)]), (double) (epsu10));
      cd_bulk = 
        min((double) ((c1 + c2*(pow(u10p, p1)))*(pow(u10p, p2))), (double) (cdmax));
      cd_wave = pow((ufric[ij - 1 + kijl*(ichnk - 1)] / u10p), 2);
      cd_ice = ooval*cd_wave + ((double) 1.0 - ooval)*cd_bulk;
      ustar = max((double) (sqrt((double) (cd_ice))*u10p), (double) (epsus));
      efd = min((double) (efd_fac*(pow(ustar, 4))), (double) (efd_max));
      em_oc = 
        max((double) (ooval*em[ij - 1] + ((double) 1.0 - ooval)*efd), (double) (efd_min))
        ;
      ffd = ffd_fac / ustar;
      f1_oc = ooval*f1[ij - 1] + ((double) 1.0 - ooval)*ffd;
      f1_oc = min((double) (max((double) (f1_oc), (double) (fr[2 - 1]))), (double) 
        (fr[nfre - 1]));
    } else {
      ooval = (double) 1.0;
      ustar = ufric[ij - 1 + kijl*(ichnk - 1)];
      em_oc = em[ij - 1];
      f1_oc = f1[ij - 1];
    }
  } else {
    ooval = (double) 1.0;
    ustar = ufric[ij - 1 + kijl*(ichnk - 1)];
    em_oc = em[ij - 1];
    f1_oc = f1[ij - 1];
  }
  
  
  tau = aird[ij - 1 + kijl*(ichnk - 1)]*max((double) (pow(ustar, 2)), (double) (epsus));
  tauxd[ij - 1 + kijl*(ichnk - 1)] = tau*sin(wdwave[ij - 1 + kijl*(ichnk - 1)]);
  tauyd[ij - 1 + kijl*(ichnk - 1)] = tau*cos(wdwave[ij - 1 + kijl*(ichnk - 1)]);
  
  tauocxd[ij - 1 + kijl*(ichnk - 1)] = tauxd[ij - 1 + kijl*(ichnk - 1)] - ooval*xstress;
  tauocyd[ij - 1 + kijl*(ichnk - 1)] = tauyd[ij - 1 + kijl*(ichnk - 1)] - ooval*ystress;
  tauo = sqrt((double) ((pow(tauocxd[ij - 1 + kijl*(ichnk - 1)], 2)) + (pow(tauocyd[ij - 
    1 + kijl*(ichnk - 1)], 2))));
  tauoc[ij - 1 + kijl*(ichnk - 1)] = 
    min((double) (max((double) (tauo / tau), (double) (tauocmin))), (double) (tauocmax));
  
  xn = aird[ij - 1 + kijl*(ichnk - 1)]*max((double) (pow(ustar, 3)), (double) (epsus3));
  phiocd[ij - 1 + kijl*(ichnk - 1)] = 
    ooval*(philf - phiwa[ij - 1]) + ((double) 1.0 - ooval)*phioc_ice*xn;
  
  phieps[ij - 1 + kijl*(ichnk - 1)] = phiocd[ij - 1 + kijl*(ichnk - 1)] / xn;
  phieps[ij - 1 + kijl*(ichnk - 1)] = min((double) (max((double) (phieps[ij - 1 + 
    kijl*(ichnk - 1)]), (double) (phiepsmin))), (double) (phiepsmax));
  
  phiocd[ij - 1 + kijl*(ichnk - 1)] = phieps[ij - 1 + kijl*(ichnk - 1)]*xn;
  
  phiaw[ij - 1 + kijl*(ichnk - 1)] = phiwa[ij - 1] / xn;
  phiaw[ij - 1 + kijl*(ichnk - 1)] = 
    ooval*phiwa[ij - 1] / xn + ((double) 1.0 - ooval)*phiaw_ice;
  
  if (lwnemocou && lnupd) {
    nphieps[ij - 1 + kijl*(ichnk - 1)] = phieps[ij - 1 + kijl*(ichnk - 1)];
    ntauoc[ij - 1 + kijl*(ichnk - 1)] = tauoc[ij - 1 + kijl*(ichnk - 1)];
    if (em_oc != (double) 0.0) {
      nswh[ij - 1 + kijl*(ichnk - 1)] = (double) 4.0*sqrt((double) (em_oc));
    } else {
      nswh[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
    }
    if (f1_oc != (double) 0.0) {
      nmwp[ij - 1 + kijl*(ichnk - 1)] = (double) 1.0 / f1_oc;
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

  
}
