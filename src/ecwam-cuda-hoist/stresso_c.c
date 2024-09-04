#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "stresso_c.h"
#include "tau_phi_hf_c.h"


__device__ void stresso_c(int kijs, int kijl, const int * __restrict__ mij, 
  const double * __restrict__ rhowgdfth, const double * __restrict__ fl1, 
  const double * __restrict__ sl, const double * __restrict__ spos, 
  const double * __restrict__ cinv, const double * __restrict__ wdwave, 
  const double * __restrict__ ufric, const double * __restrict__ z0m, 
  const double * __restrict__ aird, const double * __restrict__ rnfac, 
  const double * __restrict__ coswdif, const double * __restrict__ sinwdif2, 
  double * __restrict__ tauw, double * __restrict__ tauwdir, 
  double * __restrict__ phiwa, int llphiwa, const double * __restrict__ costh, 
  double delth, double eps1, const double * __restrict__ fr5, double g, 
  double gamnconst, double gm1, int iphys, int jtot_tauhf, int llgcbz0, int llnormagam, 
  int nang, int nfre, int nwav_gc, const double * __restrict__ omega_gc, 
  const double * __restrict__ rhowg_dfim, const double * __restrict__ sinth, 
  double sqrtgosurft, double tauwshelter, const double * __restrict__ wtauhf, 
  double x0tauhf, double xkappa, const double * __restrict__ xkm_gc, 
  const double * __restrict__ xk_gc, double xlogkratiom1_gc, double zalp, 
  double zpi4gm1, double zpi4gm2, const double * __restrict__ zpifr, int ichnk, 
  int nchnk, int ij, double * __restrict__ xstress, double * __restrict__ ystress, 
  double * __restrict__ tauhf, double * __restrict__ phihf, 
  double * __restrict__ cmrhowgdfth, double * __restrict__ us2, 
  double * __restrict__ taux, double * __restrict__ tauy, double * __restrict__ taupx, 
  double * __restrict__ taupy, double * __restrict__ usdirp, double * __restrict__ ust, 
  double * __restrict__ sumt, double * __restrict__ sumx, double * __restrict__ sumy, 
  int * __restrict__ tau_phi_hf_ns, double * __restrict__ tau_phi_hf_xks, 
  double * __restrict__ tau_phi_hf_oms, double * __restrict__ tau_phi_hf_sqrtz0og, 
  double * __restrict__ tau_phi_hf_zsup, double * __restrict__ tau_phi_hf_zinf, 
  double * __restrict__ tau_phi_hf_delz, double * __restrict__ tau_phi_hf_taul, 
  double * __restrict__ tau_phi_hf_xloggz0, double * __restrict__ tau_phi_hf_sqrtgz0, 
  double * __restrict__ tau_phi_hf_ustph, double * __restrict__ tau_phi_hf_const1, 
  double * __restrict__ tau_phi_hf_const2, double * __restrict__ tau_phi_hf_consttau, 
  double * __restrict__ tau_phi_hf_constphi, double * __restrict__ tau_phi_hf_f1dcos2, 
  double * __restrict__ tau_phi_hf_f1dcos3, double * __restrict__ tau_phi_hf_f1d, 
  double * __restrict__ tau_phi_hf_f1dsin2) {
  
  
  
  // ----------------------------------------------------------------------
  

  
  
  int m;
  int k;
  int i;
  int j;
  int ii;
  
  double tautous2;
  double cosw;
  double fcosw2;
  double zhook_handle;
  
  int ltauwshelter;
  
  // ----------------------------------------------------------------------
  
  
  
  phiwa[ij - 1] = (double) 0.0;
  xstress[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  ystress[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  
  //*    CONTRIBUTION TO THE WAVE STRESS FROM THE NEGATIVE PART OF THE WIND INPUT
  //     ------------------------------------------------------------------------
  
  if (llphiwa) {
    //     full energy flux due to negative Sinput (SL-SPOS)
    //     we assume that above NFRE, the contibutions can be neglected
    for (m = 1; m <= nfre; m += 1) {
      for (k = 1; k <= nang; k += 1) {
        phiwa[ij - 1] = phiwa[ij - 1] + (sl[ij - 1 + kijl*(k - 1 + nang*(m - 1))] - 
          spos[ij - 1 + kijl*(k - 1 + nang*(m - 1))])*rhowg_dfim[m - 1];
      }
    }
  }
  
  //*    CALCULATE LOW-FREQUENCY CONTRIBUTION TO STRESS AND ENERGY FLUX (positive sinput).
  //     ---------------------------------------------------------------------------------
  for (m = 1; m <= nfre; m += 1) {
    //     THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
    k = 1;
    sumx[ij - 1 + kijl*(ichnk - 1)] = 
      spos[ij - 1 + kijl*(k - 1 + nang*(m - 1))]*sinth[k - 1];
    sumy[ij - 1 + kijl*(ichnk - 1)] = 
      spos[ij - 1 + kijl*(k - 1 + nang*(m - 1))]*costh[k - 1];
    for (k = 2; k <= nang; k += 1) {
      sumx[ij - 1 + kijl*(ichnk - 1)] = sumx[ij - 1 + kijl*(ichnk - 1)] + spos[ij - 1 + 
        kijl*(k - 1 + nang*(m - 1))]*sinth[k - 1];
      sumy[ij - 1 + kijl*(ichnk - 1)] = sumy[ij - 1 + kijl*(ichnk - 1)] + spos[ij - 1 + 
        kijl*(k - 1 + nang*(m - 1))]*costh[k - 1];
    }
    cmrhowgdfth[ij - 1 + kijl*(ichnk - 1)] = 
      rhowgdfth[ij - 1 + kijl*(m - 1)]*cinv[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))];
    xstress[ij - 1 + kijl*(ichnk - 1)] = xstress[ij - 1 + kijl*(ichnk - 1)] + 
      cmrhowgdfth[ij - 1 + kijl*(ichnk - 1)]*sumx[ij - 1 + kijl*(ichnk - 1)];
    ystress[ij - 1 + kijl*(ichnk - 1)] = ystress[ij - 1 + kijl*(ichnk - 1)] + 
      cmrhowgdfth[ij - 1 + kijl*(ichnk - 1)]*sumy[ij - 1 + kijl*(ichnk - 1)];
  }
  
  //     TAUW is the kinematic wave stress !
  xstress[ij - 1 + kijl*(ichnk - 1)] = xstress[ij - 1 + kijl*(ichnk - 1)] / fmax((double)
     (aird[ij - 1 + kijl*(ichnk - 1)]), (double) ((double) 1.0));
  ystress[ij - 1 + kijl*(ichnk - 1)] = ystress[ij - 1 + kijl*(ichnk - 1)] / fmax((double)
     (aird[ij - 1 + kijl*(ichnk - 1)]), (double) ((double) 1.0));
  
  if (llphiwa) {
    for (m = 1; m <= nfre; m += 1) {
      //       THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
      k = 1;
      sumt[ij - 1 + kijl*(ichnk - 1)] = spos[ij - 1 + kijl*(k - 1 + nang*(m - 1))];
      for (k = 2; k <= nang; k += 1) {
        sumt[ij - 1 + kijl*(ichnk - 1)] = 
          sumt[ij - 1 + kijl*(ichnk - 1)] + spos[ij - 1 + kijl*(k - 1 + nang*(m - 1))];
      }
      phiwa[ij - 1] = 
        phiwa[ij - 1] + rhowgdfth[ij - 1 + kijl*(m - 1)]*sumt[ij - 1 + kijl*(ichnk - 1)];
    }
  }
  
  //*    CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS and energy flux (positive sinput).
  //     ----------------------------------------------------------------------------------
  
  if (iphys == 0 || tauwshelter == (double) 0.0) {
    ltauwshelter = false;
    usdirp[ij - 1 + kijl*(ichnk - 1)] = wdwave[ij - 1 + kijl*(ichnk - 1)];
    ust[ij - 1 + kijl*(ichnk - 1)] = ufric[ij - 1 + kijl*(ichnk - 1)];
  } else {
    ltauwshelter = true;
    taux[ij - 1 + kijl*(ichnk - 1)] = 
      (pow(ufric[ij - 1 + kijl*(ichnk - 1)], 2))*sin(wdwave[ij - 1 + kijl*(ichnk - 1)]);
    tauy[ij - 1 + kijl*(ichnk - 1)] = 
      (pow(ufric[ij - 1 + kijl*(ichnk - 1)], 2))*cos(wdwave[ij - 1 + kijl*(ichnk - 1)]);
    taupx[ij - 1 + kijl*(ichnk - 1)] = 
      taux[ij - 1 + kijl*(ichnk - 1)] - tauwshelter*xstress[ij - 1 + kijl*(ichnk - 1)];
    taupy[ij - 1 + kijl*(ichnk - 1)] = 
      tauy[ij - 1 + kijl*(ichnk - 1)] - tauwshelter*ystress[ij - 1 + kijl*(ichnk - 1)];
    usdirp[ij - 1 + kijl*(ichnk - 1)] = 
      atan2(taupx[ij - 1 + kijl*(ichnk - 1)], taupy[ij - 1 + kijl*(ichnk - 1)]);
    ust[ij - 1 + kijl*(ichnk - 1)] = pow(((pow(taupx[ij - 1 + kijl*(ichnk - 1)], 2)) + 
      (pow(taupy[ij - 1 + kijl*(ichnk - 1)], 2))), (double) 0.25);
  }
  
  
  tau_phi_hf_c(kijs, kijl, mij, ltauwshelter, ufric, z0m, fl1, aird, rnfac, coswdif, 
    sinwdif2,  (&ust[ + kijl*(ichnk - 1)]),  (&tauhf[ + kijl*(ichnk - 1)]), 
     (&phihf[ + kijl*(ichnk - 1)]), llphiwa, delth, fr5, g, gamnconst, gm1, jtot_tauhf, 
    llgcbz0, llnormagam, nang, nfre, nwav_gc, omega_gc, sqrtgosurft, tauwshelter, 
    wtauhf, x0tauhf, xkappa, xkm_gc, xk_gc, xlogkratiom1_gc, zalp, zpi4gm1, zpi4gm2, 
    zpifr, ichnk, nchnk, ij, tau_phi_hf_ns, tau_phi_hf_xks, tau_phi_hf_oms, 
    tau_phi_hf_sqrtz0og, tau_phi_hf_zsup, tau_phi_hf_zinf, tau_phi_hf_delz, 
    tau_phi_hf_taul, tau_phi_hf_xloggz0, tau_phi_hf_sqrtgz0, tau_phi_hf_ustph, 
    tau_phi_hf_const1, tau_phi_hf_const2, tau_phi_hf_consttau, tau_phi_hf_constphi, 
    tau_phi_hf_f1dcos2, tau_phi_hf_f1dcos3, tau_phi_hf_f1d, tau_phi_hf_f1dsin2);
  
  
  xstress[ij - 1 + kijl*(ichnk - 1)] = xstress[ij - 1 + kijl*(ichnk - 1)] + tauhf[ij - 1 
    + kijl*(ichnk - 1)]*sin(usdirp[ij - 1 + kijl*(ichnk - 1)]);
  ystress[ij - 1 + kijl*(ichnk - 1)] = ystress[ij - 1 + kijl*(ichnk - 1)] + tauhf[ij - 1 
    + kijl*(ichnk - 1)]*cos(usdirp[ij - 1 + kijl*(ichnk - 1)]);
  tauw[ij - 1 + kijl*(ichnk - 1)] = sqrt((double) ((pow(xstress[ij - 1 + kijl*(ichnk - 1)
    ], 2)) + (pow(ystress[ij - 1 + kijl*(ichnk - 1)], 2))));
  tauw[ij - 1 + kijl*(ichnk - 1)] = 
    fmax((double) (tauw[ij - 1 + kijl*(ichnk - 1)]), (double) ((double) 0.0));
  tauwdir[ij - 1 + kijl*(ichnk - 1)] = 
    atan2(xstress[ij - 1 + kijl*(ichnk - 1)], ystress[ij - 1 + kijl*(ichnk - 1)]);
  
  if (!llgcbz0) {
    tautous2 = (double) 1.0 / ((double) 1.0 + eps1);
    tauw[ij - 1 + kijl*(ichnk - 1)] = fmin((double) (tauw[ij - 1 + kijl*(ichnk - 1)]), 
      (double) ((pow(ufric[ij - 1 + kijl*(ichnk - 1)], 2))*tautous2));
  }
  
  if (llphiwa) {
    phiwa[ij - 1] = phiwa[ij - 1] + phihf[ij - 1 + kijl*(ichnk - 1)];
  }
  
  
  
}
