#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "stresso_c.h"
#include "tau_phi_hf_c.h"

__device__ void stresso_c(int kijs, int kijl, const int * mij, const double * rhowgdfth, 
  const double * fl1, const double * sl, const double * spos, const double * cinv, 
  const double * wdwave, const double * ufric, const double * z0m, const double * aird, 
  const double * rnfac, const double * coswdif, const double * sinwdif2, double * tauw, 
  double * tauwdir, double * phiwa, int llphiwa, const double * costh, double delth, 
  double eps1, const double * fr5, double g, double gamnconst, double gm1, int iphys, 
  int jtot_tauhf, int llgcbz0, int llnormagam, int nang, int nfre, int nwav_gc, 
  const double * omega_gc, const double * rhowg_dfim, const double * sinth, 
  double sqrtgosurft, double tauwshelter, const double * wtauhf, double x0tauhf, 
  double xkappa, const double * xkm_gc, const double * xk_gc, double xlogkratiom1_gc, 
  double zalp, double zpi4gm1, double zpi4gm2, const double * zpifr, int ichnk, 
  int nchnk, int ij, double * xstress, double * ystress, double * tauhf, double * phihf, 
  double * usdirp, double * ust) {
  
  
  
  const int nfre_loki_param = 36;
  const int nang_loki_param = 24;
  int m;
  int k;
  int i;
  int j;
  int ii;
  
  double tautous2;
  double cosw;
  double fcosw2;
  
  double cmrhowgdfth;
  double taux;
  double tauy;
  double taupx;
  double taupy;
  double sumt;
  double sumx;
  double sumy;
  
  int ltauwshelter;
  

  phiwa[ij - 1] = (double) 0.0;
  xstress[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  ystress[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  if (llphiwa) {
    for (m = 1; m <= nfre; m += 1) {
      for (k = 1; k <= nang; k += 1) {
        phiwa[ij - 1] = phiwa[ij - 1] + (sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1
          ))] - spos[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))])*rhowg_dfim[m - 1];
      }
    }
  }
  for (m = 1; m <= nfre; m += 1) {
    //     THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
    k = 1;
    sumx = spos[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))]*sinth[k - 1];
    sumy = spos[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))]*costh[k - 1];
    for (k = 2; k <= nang; k += 1) {
      sumx = sumx + spos[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))]*sinth[k - 1];
      sumy = sumy + spos[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))]*costh[k - 1];
    }
    cmrhowgdfth = rhowgdfth[ij - 1 + kijl*(m - 1)]*cinv[ij - 1 + kijl*(m - 1 + 
      nfre_loki_param*(ichnk - 1))];
    xstress[ij - 1 + kijl*(ichnk - 1)] = 
      xstress[ij - 1 + kijl*(ichnk - 1)] + cmrhowgdfth*sumx;
    ystress[ij - 1 + kijl*(ichnk - 1)] = 
      ystress[ij - 1 + kijl*(ichnk - 1)] + cmrhowgdfth*sumy;
  }
  xstress[ij - 1 + kijl*(ichnk - 1)] = xstress[ij - 1 + kijl*(ichnk - 1)] / max((double) 
    (aird[ij - 1 + kijl*(ichnk - 1)]), (double) ((double) 1.0));
  ystress[ij - 1 + kijl*(ichnk - 1)] = ystress[ij - 1 + kijl*(ichnk - 1)] / max((double) 
    (aird[ij - 1 + kijl*(ichnk - 1)]), (double) ((double) 1.0));
  
  if (llphiwa) {
    for (m = 1; m <= nfre; m += 1) {
      //       THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
      k = 1;
      sumt = spos[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))];
      for (k = 2; k <= nang; k += 1) {
        sumt = sumt + spos[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))];
      }
      phiwa[ij - 1] = phiwa[ij - 1] + rhowgdfth[ij - 1 + kijl*(m - 1)]*sumt;
    }
  }
  if (iphys == 0 || tauwshelter == (double) 0.0) {
    ltauwshelter = false;
    usdirp[ij - 1 + kijl*(ichnk - 1)] = wdwave[ij - 1 + kijl*(ichnk - 1)];
    ust[ij - 1 + kijl*(ichnk - 1)] = ufric[ij - 1 + kijl*(ichnk - 1)];
  } else {
    ltauwshelter = true;
    taux = 
      (pow(ufric[ij - 1 + kijl*(ichnk - 1)], 2))*sin(wdwave[ij - 1 + kijl*(ichnk - 1)]);
    tauy = 
      (pow(ufric[ij - 1 + kijl*(ichnk - 1)], 2))*cos(wdwave[ij - 1 + kijl*(ichnk - 1)]);
    taupx = taux - tauwshelter*xstress[ij - 1 + kijl*(ichnk - 1)];
    taupy = tauy - tauwshelter*ystress[ij - 1 + kijl*(ichnk - 1)];
    usdirp[ij - 1 + kijl*(ichnk - 1)] = atan2(taupx, taupy);
    ust[ij - 1 + kijl*(ichnk - 1)] = 
      pow(((pow(taupx, 2)) + (pow(taupy, 2))), (double) 0.25);
  }

  
  tau_phi_hf_c(kijs, kijl, mij, ltauwshelter, ufric, z0m, fl1, aird, rnfac, coswdif, 
    sinwdif2,  (&ust[ + kijl*(ichnk - 1)]),  (&tauhf[ + kijl*(ichnk - 1)]), 
     (&phihf[ + kijl*(ichnk - 1)]), llphiwa, delth, fr5, g, gamnconst, gm1, jtot_tauhf, 
    llgcbz0, llnormagam, nang, nwav_gc, omega_gc, sqrtgosurft, tauwshelter, wtauhf, 
    x0tauhf, xkappa, xkm_gc, xk_gc, xlogkratiom1_gc, zalp, zpi4gm1, zpi4gm2, zpifr, 
    ichnk, nchnk, ij);
  

  xstress[ij - 1 + kijl*(ichnk - 1)] = xstress[ij - 1 + kijl*(ichnk - 1)] + tauhf[ij - 1 
    + kijl*(ichnk - 1)]*sin(usdirp[ij - 1 + kijl*(ichnk - 1)]);
  ystress[ij - 1 + kijl*(ichnk - 1)] = ystress[ij - 1 + kijl*(ichnk - 1)] + tauhf[ij - 1 
    + kijl*(ichnk - 1)]*cos(usdirp[ij - 1 + kijl*(ichnk - 1)]);
  tauw[ij - 1 + kijl*(ichnk - 1)] = sqrt((double) ((pow(xstress[ij - 1 + kijl*(ichnk - 1)
    ], 2)) + (pow(ystress[ij - 1 + kijl*(ichnk - 1)], 2))));
  tauw[ij - 1 + kijl*(ichnk - 1)] = 
    max((double) (tauw[ij - 1 + kijl*(ichnk - 1)]), (double) ((double) 0.0));
  tauwdir[ij - 1 + kijl*(ichnk - 1)] = 
    atan2(xstress[ij - 1 + kijl*(ichnk - 1)], ystress[ij - 1 + kijl*(ichnk - 1)]);
  
  if (!llgcbz0) {
    tautous2 = (double) 1.0 / ((double) 1.0 + eps1);
    tauw[ij - 1 + kijl*(ichnk - 1)] = min((double) (tauw[ij - 1 + kijl*(ichnk - 1)]), 
      (double) ((pow(ufric[ij - 1 + kijl*(ichnk - 1)], 2))*tautous2));
  }
  
  if (llphiwa) {
    phiwa[ij - 1] = phiwa[ij - 1] + phihf[ij - 1 + kijl*(ichnk - 1)];
  }

  
  
}
