#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "tau_phi_hf_c.h"
#include "omegagc_c.h"


__device__ void tau_phi_hf_c(int kijs, int kijl, const int * __restrict__ mij, 
  int ltauwshelter, const double * __restrict__ ufric, const double * __restrict__ z0m, 
  const double * __restrict__ fl1, const double * __restrict__ aird, 
  const double * __restrict__ rnfac, const double * __restrict__ coswdif, 
  const double * __restrict__ sinwdif2, double * __restrict__ ust, 
  double * __restrict__ tauhf, double * __restrict__ phihf, int llphihf, double delth, 
  const double * __restrict__ fr5, double g, double gamnconst, double gm1, 
  int jtot_tauhf, int llgcbz0, int llnormagam, int nang, int nfre, int nwav_gc, 
  const double * __restrict__ omega_gc, double sqrtgosurft, double tauwshelter, 
  const double * __restrict__ wtauhf, double x0tauhf, double xkappa, 
  const double * __restrict__ xkm_gc, const double * __restrict__ xk_gc, 
  double xlogkratiom1_gc, double zalp, double zpi4gm1, double zpi4gm2, 
  const double * __restrict__ zpifr, int ichnk, int nchnk, int ij, 
  int * __restrict__ ns, double * __restrict__ xks, double * __restrict__ oms, 
  double * __restrict__ sqrtz0og, double * __restrict__ zsup, 
  double * __restrict__ zinf, double * __restrict__ delz, double * __restrict__ taul, 
  double * __restrict__ xloggz0, double * __restrict__ sqrtgz0, 
  double * __restrict__ ustph, double * __restrict__ const1, 
  double * __restrict__ const2, double * __restrict__ consttau, 
  double * __restrict__ constphi, double * __restrict__ f1dcos2, 
  double * __restrict__ f1dcos3, double * __restrict__ f1d, double * __restrict__ f1dsin2
  ) {
  
  
  
  // ----------------------------------------------------------------------
  

  
  
  int j;
  int k;
  
  double zsupmax = (double) 0.0;  //  LOG(1.)
  double omega;
  double omegacc;
  double x0g;
  double yc;
  double y;
  double cm1;
  double zx;
  double zarg;
  double zlog;
  double zbeta;
  double fnc;
  double fnc2;
  double gamnorma;  // RENORMALISATION FACTOR OF THE GROWTH RATE
  double znz;
  double confg;
  double cosw;
  double fcosw2;
  double zhook_handle;
  
  
  // ----------------------------------------------------------------------
  
  
  if (llgcbz0) {
    omegagc_c(kijs, kijl, ufric,  (&ns[ + kijl*(ichnk - 1)]), 
       (&xks[ + kijl*(ichnk - 1)]),  (&oms[ + kijl*(ichnk - 1)]), nwav_gc, omega_gc, 
      sqrtgosurft, xkm_gc, xk_gc, xlogkratiom1_gc, ichnk, nchnk, ij);
  }
  
  //     See INIT_X0TAUHF
  x0g = x0tauhf*g;
  
  
  if (llphihf) {
    ustph[ij - 1 + kijl*(ichnk - 1)] = ust[ij - 1];
  }
  
  //*    COMPUTE THE INTEGRALS
  //     ---------------------
  
  xloggz0[ij - 1 + kijl*(ichnk - 1)] = log(g*z0m[ij - 1 + kijl*(ichnk - 1)]);
  omegacc = fmax((double) (zpifr[mij[ij - 1 + kijl*(ichnk - 1)] - 1]), (double) (x0g / 
    ust[ij - 1]));
  sqrtz0og[ij - 1 + kijl*(ichnk - 1)] = 
    sqrt((double) (z0m[ij - 1 + kijl*(ichnk - 1)]*gm1));
  sqrtgz0[ij - 1 + kijl*(ichnk - 1)] = (double) 1.0 / sqrtz0og[ij - 1 + kijl*(ichnk - 1)]
    ;
  yc = omegacc*sqrtz0og[ij - 1 + kijl*(ichnk - 1)];
  zinf[ij - 1 + kijl*(ichnk - 1)] = log(yc);
  
  consttau[ij - 1 + kijl*(ichnk - 1)] = zpi4gm2*fr5[mij[ij - 1 + kijl*(ichnk - 1)] - 1];
  
  k = 1;
  cosw = fmax((double) (coswdif[ij - 1 + kijl*(k - 1)]), (double) ((double) 0.0));
  fcosw2 = fl1[ij - 1 + kijl*(k - 1 + nang*(mij[ij - 1 + kijl*(ichnk - 1)] - 1 + 
    nfre*(ichnk - 1)))]*(pow(cosw, 2));
  f1dcos3[ij - 1 + kijl*(ichnk - 1)] = fcosw2*cosw;
  f1dcos2[ij - 1 + kijl*(ichnk - 1)] = fcosw2;
  f1dsin2[ij - 1 + kijl*(ichnk - 1)] = fl1[ij - 1 + kijl*(k - 1 + nang*(mij[ij - 1 + 
    kijl*(ichnk - 1)] - 1 + nfre*(ichnk - 1)))]*sinwdif2[ij - 1 + kijl*(k - 1)];
  f1d[ij - 1 + kijl*(ichnk - 1)] = fl1[ij - 1 + kijl*(k - 1 + nang*(mij[ij - 1 + 
    kijl*(ichnk - 1)] - 1 + nfre*(ichnk - 1)))];
  for (k = 2; k <= nang; k += 1) {
    cosw = fmax((double) (coswdif[ij - 1 + kijl*(k - 1)]), (double) ((double) 0.0));
    fcosw2 = fl1[ij - 1 + kijl*(k - 1 + nang*(mij[ij - 1 + kijl*(ichnk - 1)] - 1 + 
      nfre*(ichnk - 1)))]*(pow(cosw, 2));
    f1dcos3[ij - 1 + kijl*(ichnk - 1)] = f1dcos3[ij - 1 + kijl*(ichnk - 1)] + fcosw2*cosw
      ;
    f1dcos2[ij - 1 + kijl*(ichnk - 1)] = f1dcos2[ij - 1 + kijl*(ichnk - 1)] + fcosw2;
    f1dsin2[ij - 1 + kijl*(ichnk - 1)] = f1dsin2[ij - 1 + kijl*(ichnk - 1)] + fl1[ij - 1 
      + kijl*(k - 1 + nang*(mij[ij - 1 + kijl*(ichnk - 1)] - 1 + nfre*(ichnk - 1)))
      ]*sinwdif2[ij - 1 + kijl*(k - 1)];
    f1d[ij - 1 + kijl*(ichnk - 1)] = f1d[ij - 1 + kijl*(ichnk - 1)] + fl1[ij - 1 + 
      kijl*(k - 1 + nang*(mij[ij - 1 + kijl*(ichnk - 1)] - 1 + nfre*(ichnk - 1)))];
  }
  f1dcos3[ij - 1 + kijl*(ichnk - 1)] = delth*f1dcos3[ij - 1 + kijl*(ichnk - 1)];
  f1dcos2[ij - 1 + kijl*(ichnk - 1)] = delth*f1dcos2[ij - 1 + kijl*(ichnk - 1)];
  f1dsin2[ij - 1 + kijl*(ichnk - 1)] = delth*f1dsin2[ij - 1 + kijl*(ichnk - 1)];
  f1d[ij - 1 + kijl*(ichnk - 1)] = delth*f1d[ij - 1 + kijl*(ichnk - 1)];
  
  if (llnormagam) {
    confg = gamnconst*fr5[mij[ij - 1 + kijl*(ichnk - 1)] - 1]*rnfac[ij - 1]*sqrtgz0[ij - 
      1 + kijl*(ichnk - 1)];
    const1[ij - 1 + kijl*(ichnk - 1)] = confg*f1dsin2[ij - 1 + kijl*(ichnk - 1)];
    const2[ij - 1 + kijl*(ichnk - 1)] = confg*f1d[ij - 1 + kijl*(ichnk - 1)];
  } else {
    const1[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
    const2[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
  }
  
  
  //     TAUHF :
  if (llgcbz0) {
    zsup[ij - 1 + kijl*(ichnk - 1)] = fmin((double) (log(oms[ij - 1 + kijl*(ichnk - 1)
      ]*sqrtz0og[ij - 1 + kijl*(ichnk - 1)])), (double) (zsupmax));
  } else {
    zsup[ij - 1 + kijl*(ichnk - 1)] = zsupmax;
  }
  
  taul[ij - 1 + kijl*(ichnk - 1)] = pow(ust[ij - 1], 2);
  delz[ij - 1 + kijl*(ichnk - 1)] = fmax((double) ((zsup[ij - 1 + kijl*(ichnk - 1)] - 
    zinf[ij - 1 + kijl*(ichnk - 1)]) / (double) (jtot_tauhf - 1)), (double) ((double) 0.0
    ));
  tauhf[ij - 1] = (double) 0.0;
  
  // Intergrals are integrated following a change of variable : Z=LOG(Y)
  if (ltauwshelter) {
    for (j = 1; j <= jtot_tauhf; j += 1) {
      y = exp((double) (zinf[ij - 1 + kijl*(ichnk - 1)] + (double) (j - 1)*delz[ij - 1 + 
        kijl*(ichnk - 1)]));
      omega = y*sqrtgz0[ij - 1 + kijl*(ichnk - 1)];
      cm1 = omega*gm1;
      zx = ust[ij - 1]*cm1 + zalp;
      zarg = xkappa / zx;
      zlog = xloggz0[ij - 1 + kijl*(ichnk - 1)] + (double) 2.0*log(cm1) + zarg;
      zlog = fmin((double) (zlog), (double) ((double) 0.0));
      zbeta = (pow(zlog, 4))*exp((double) (zlog));
      znz = zbeta*ust[ij - 1]*y;
      gamnorma = ((double) 1.0 + const1[ij - 1 + kijl*(ichnk - 1)]*znz) / ((double) 1.0 +
         const2[ij - 1 + kijl*(ichnk - 1)]*znz);
      fnc2 = f1dcos3[ij - 1 + kijl*(ichnk - 1)]*consttau[ij - 1 + kijl*(ichnk - 1)
        ]*zbeta*taul[ij - 1 + kijl*(ichnk - 1)]*wtauhf[j - 1]*delz[ij - 1 + kijl*(ichnk -
         1)]*gamnorma;
      taul[ij - 1 + kijl*(ichnk - 1)] = fmax((double) (taul[ij - 1 + kijl*(ichnk - 1)] - 
        tauwshelter*fnc2), (double) ((double) 0.0));
      
      ust[ij - 1] = sqrt((double) (taul[ij - 1 + kijl*(ichnk - 1)]));
      tauhf[ij - 1] = tauhf[ij - 1] + fnc2;
    }
  } else {
    for (j = 1; j <= jtot_tauhf; j += 1) {
      y = exp((double) (zinf[ij - 1 + kijl*(ichnk - 1)] + (double) (j - 1)*delz[ij - 1 + 
        kijl*(ichnk - 1)]));
      omega = y*sqrtgz0[ij - 1 + kijl*(ichnk - 1)];
      cm1 = omega*gm1;
      zx = ust[ij - 1]*cm1 + zalp;
      zarg = xkappa / zx;
      zlog = xloggz0[ij - 1 + kijl*(ichnk - 1)] + (double) 2.0*log(cm1) + zarg;
      zlog = fmin((double) (zlog), (double) ((double) 0.0));
      zbeta = (pow(zlog, 4))*exp((double) (zlog));
      fnc2 = zbeta*wtauhf[j - 1];
      znz = zbeta*ust[ij - 1]*y;
      gamnorma = ((double) 1.0 + const1[ij - 1 + kijl*(ichnk - 1)]*znz) / ((double) 1.0 +
         const2[ij - 1 + kijl*(ichnk - 1)]*znz);
      tauhf[ij - 1] = tauhf[ij - 1] + fnc2*gamnorma;
    }
    tauhf[ij - 1] = f1dcos3[ij - 1 + kijl*(ichnk - 1)]*consttau[ij - 1 + kijl*(ichnk - 1)
      ]*taul[ij - 1 + kijl*(ichnk - 1)]*tauhf[ij - 1]*delz[ij - 1 + kijl*(ichnk - 1)];
  }
  
  
  phihf[ij - 1] = (double) 0.0;
  if (llphihf) {
    //       PHIHF:
    //       We are neglecting the gravity-capillary contribution
    //       Recompute DELZ over the full interval
    taul[ij - 1 + kijl*(ichnk - 1)] = pow(ustph[ij - 1 + kijl*(ichnk - 1)], 2);
    zsup[ij - 1 + kijl*(ichnk - 1)] = zsupmax;
    delz[ij - 1 + kijl*(ichnk - 1)] = fmax((double) ((zsup[ij - 1 + kijl*(ichnk - 1)] - 
      zinf[ij - 1 + kijl*(ichnk - 1)]) / (double) (jtot_tauhf - 1)), (double) ((double) 
      0.0));
    
    constphi[ij - 1 + kijl*(ichnk - 1)] = 
      aird[ij - 1 + kijl*(ichnk - 1)]*zpi4gm1*fr5[mij[ij - 1 + kijl*(ichnk - 1)] - 1];
    
    // Intergrals are integrated following a change of variable : Z=LOG(Y)
    if (ltauwshelter) {
      for (j = 1; j <= jtot_tauhf; j += 1) {
        y = exp((double) (zinf[ij - 1 + kijl*(ichnk - 1)] + (double) (j - 1)*delz[ij - 1 
          + kijl*(ichnk - 1)]));
        omega = y*sqrtgz0[ij - 1 + kijl*(ichnk - 1)];
        cm1 = omega*gm1;
        zx = ustph[ij - 1 + kijl*(ichnk - 1)]*cm1 + zalp;
        zarg = xkappa / zx;
        zlog = xloggz0[ij - 1 + kijl*(ichnk - 1)] + (double) 2.0*log(cm1) + zarg;
        zlog = fmin((double) (zlog), (double) ((double) 0.0));
        zbeta = (pow(zlog, 4))*exp((double) (zlog));
        znz = zbeta*ust[ij - 1]*y;
        gamnorma = ((double) 1.0 + const1[ij - 1 + kijl*(ichnk - 1)]*znz) / ((double) 1.0
           + const2[ij - 1 + kijl*(ichnk - 1)]*znz);
        fnc2 = zbeta*taul[ij - 1 + kijl*(ichnk - 1)]*wtauhf[j - 1]*delz[ij - 1 + 
          kijl*(ichnk - 1)]*gamnorma;
        taul[ij - 1 + kijl*(ichnk - 1)] = fmax((double) (taul[ij - 1 + kijl*(ichnk - 1)] 
          - tauwshelter*f1dcos3[ij - 1 + kijl*(ichnk - 1)]*consttau[ij - 1 + kijl*(ichnk 
          - 1)]*fnc2), (double) ((double) 0.0));
        ustph[ij - 1 + kijl*(ichnk - 1)] = 
          sqrt((double) (taul[ij - 1 + kijl*(ichnk - 1)]));
        phihf[ij - 1] = phihf[ij - 1] + fnc2 / y;
      }
      phihf[ij - 1] = f1dcos2[ij - 1 + kijl*(ichnk - 1)]*constphi[ij - 1 + kijl*(ichnk - 
        1)]*sqrtz0og[ij - 1 + kijl*(ichnk - 1)]*phihf[ij - 1];
    } else {
      for (j = 1; j <= jtot_tauhf; j += 1) {
        y = exp((double) (zinf[ij - 1 + kijl*(ichnk - 1)] + (double) (j - 1)*delz[ij - 1 
          + kijl*(ichnk - 1)]));
        omega = y*sqrtgz0[ij - 1 + kijl*(ichnk - 1)];
        cm1 = omega*gm1;
        zx = ustph[ij - 1 + kijl*(ichnk - 1)]*cm1 + zalp;
        zarg = xkappa / zx;
        zlog = xloggz0[ij - 1 + kijl*(ichnk - 1)] + (double) 2.0*log(cm1) + zarg;
        zlog = fmin((double) (zlog), (double) ((double) 0.0));
        zbeta = (pow(zlog, 4))*exp((double) (zlog));
        znz = zbeta*ust[ij - 1]*y;
        gamnorma = ((double) 1.0 + const1[ij - 1 + kijl*(ichnk - 1)]*znz) / ((double) 1.0
           + const2[ij - 1 + kijl*(ichnk - 1)]*znz);
        fnc2 = zbeta*wtauhf[j - 1]*gamnorma;
        phihf[ij - 1] = phihf[ij - 1] + fnc2 / y;
      }
      phihf[ij - 1] = f1dcos2[ij - 1 + kijl*(ichnk - 1)]*constphi[ij - 1 + kijl*(ichnk - 
        1)]*sqrtz0og[ij - 1 + kijl*(ichnk - 1)]*taul[ij - 1 + kijl*(ichnk - 1)]*phihf[ij 
        - 1]*delz[ij - 1 + kijl*(ichnk - 1)];
    }
  }
  
  
  
}
