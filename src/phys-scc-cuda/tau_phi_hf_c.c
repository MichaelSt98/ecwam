#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "tau_phi_hf_c.h"
#include "omegagc_c.c"

__device__ void tau_phi_hf_c(int kijs, int kijl, const int * mij, int ltauwshelter, 
  const double * ufric, const double * z0m, const double * fl1, const double * aird, 
  const double * rnfac, const double * coswdif, const double * sinwdif2, double * ust, 
  double * tauhf, double * phihf, int llphihf, double delth, const double * fr5, 
  double g, double gamnconst, double gm1, int jtot_tauhf, int llgcbz0, int llnormagam, 
  int nang, int nwav_gc, const double * omega_gc, double sqrtgosurft, 
  double tauwshelter, const double * wtauhf, double x0tauhf, double xkappa, 
  const double * xkm_gc, const double * xk_gc, double xlogkratiom1_gc, double zalp, 
  double zpi4gm1, double zpi4gm2, const double * zpifr, int ichnk, int nchnk, int ij) {
  
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  int j;
  int k;
  int ns;
  
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
  
  double xks;
  double oms;
  double sqrtz0og;
  double zsup;
  double zinf;
  double delz;
  double taul;
  double xloggz0;
  double sqrtgz0;
  double ustph;
  double const1;
  double const2;
  double consttau;
  double constphi;
  double f1dcos2;
  double f1dcos3;
  double f1d;
  double f1dsin2;
  

  if (llgcbz0) {
    omegagc_c(ufric[ij - 1 + kijl*(ichnk - 1)],  (&ns),  (&xks),  (&oms), nwav_gc, 
      omega_gc, sqrtgosurft, xkm_gc, xk_gc, xlogkratiom1_gc);
  }
  x0g = x0tauhf*g;
  
  if (llphihf) {
    ustph = ust[ij - 1];
  }
  xloggz0 = log(g*z0m[ij - 1 + kijl*(ichnk - 1)]);
  omegacc = max((double) (zpifr[mij[ij - 1 + kijl*(ichnk - 1)] - 1]), (double) (x0g / 
    ust[ij - 1]));
  sqrtz0og = sqrt((double) (z0m[ij - 1 + kijl*(ichnk - 1)]*gm1));
  sqrtgz0 = (double) 1.0 / sqrtz0og;
  yc = omegacc*sqrtz0og;
  zinf = log(yc);
  
  consttau = zpi4gm2*fr5[mij[ij - 1 + kijl*(ichnk - 1)] - 1];
  
  k = 1;
  cosw = max((double) (coswdif[ij - 1 + kijl*(k - 1)]), (double) ((double) 0.0));
  fcosw2 = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(mij[ij - 1 + kijl*(ichnk - 1)] - 1
     + nfre_loki_param*(ichnk - 1)))]*(pow(cosw, 2));
  f1dcos3 = fcosw2*cosw;
  f1dcos2 = fcosw2;
  f1dsin2 = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(mij[ij - 1 + kijl*(ichnk - 1)] - 
    1 + nfre_loki_param*(ichnk - 1)))]*sinwdif2[ij - 1 + kijl*(k - 1)];
  f1d = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(mij[ij - 1 + kijl*(ichnk - 1)] - 1 + 
    nfre_loki_param*(ichnk - 1)))];
  for (k = 2; k <= nang; k += 1) {
    cosw = max((double) (coswdif[ij - 1 + kijl*(k - 1)]), (double) ((double) 0.0));
    fcosw2 = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(mij[ij - 1 + kijl*(ichnk - 1)] -
       1 + nfre_loki_param*(ichnk - 1)))]*(pow(cosw, 2));
    f1dcos3 = f1dcos3 + fcosw2*cosw;
    f1dcos2 = f1dcos2 + fcosw2;
    f1dsin2 = f1dsin2 + fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(mij[ij - 1 + 
      kijl*(ichnk - 1)] - 1 + nfre_loki_param*(ichnk - 1)))]*sinwdif2[ij - 1 + kijl*(k - 
      1)];
    f1d = f1d + fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(mij[ij - 1 + kijl*(ichnk - 1)
      ] - 1 + nfre_loki_param*(ichnk - 1)))];
  }
  f1dcos3 = delth*f1dcos3;
  f1dcos2 = delth*f1dcos2;
  f1dsin2 = delth*f1dsin2;
  f1d = delth*f1d;
  
  if (llnormagam) {
    confg = gamnconst*fr5[mij[ij - 1 + kijl*(ichnk - 1)] - 1]*rnfac[ij - 1]*sqrtgz0;
    const1 = confg*f1dsin2;
    const2 = confg*f1d;
  } else {
    const1 = (double) 0.0;
    const2 = (double) 0.0;
  }
  if (llgcbz0) {
    zsup = min((double) (log(oms*sqrtz0og)), (double) (zsupmax));
  } else {
    zsup = zsupmax;
  }
  
  taul = pow(ust[ij - 1], 2);
  delz = 
    max((double) ((zsup - zinf) / (double) (jtot_tauhf - 1)), (double) ((double) 0.0));
  tauhf[ij - 1] = (double) 0.0;
  if (ltauwshelter) {
    for (j = 1; j <= jtot_tauhf; j += 1) {
      y = exp((double) (zinf + (double) (j - 1)*delz));
      omega = y*sqrtgz0;
      cm1 = omega*gm1;
      zx = ust[ij - 1]*cm1 + zalp;
      zarg = xkappa / zx;
      zlog = xloggz0 + (double) 2.0*log(cm1) + zarg;
      zlog = min((double) (zlog), (double) ((double) 0.0));
      zbeta = (pow(zlog, 4))*exp((double) (zlog));
      znz = zbeta*ust[ij - 1]*y;
      gamnorma = ((double) 1.0 + const1*znz) / ((double) 1.0 + const2*znz);
      fnc2 = f1dcos3*consttau*zbeta*taul*wtauhf[j - 1]*delz*gamnorma;
      taul = max((double) (taul - tauwshelter*fnc2), (double) ((double) 0.0));
      
      ust[ij - 1] = sqrt((double) (taul));
      tauhf[ij - 1] = tauhf[ij - 1] + fnc2;
    }
  } else {
    for (j = 1; j <= jtot_tauhf; j += 1) {
      y = exp((double) (zinf + (double) (j - 1)*delz));
      omega = y*sqrtgz0;
      cm1 = omega*gm1;
      zx = ust[ij - 1]*cm1 + zalp;
      zarg = xkappa / zx;
      zlog = xloggz0 + (double) 2.0*log(cm1) + zarg;
      zlog = min((double) (zlog), (double) ((double) 0.0));
      zbeta = (pow(zlog, 4))*exp((double) (zlog));
      fnc2 = zbeta*wtauhf[j - 1];
      znz = zbeta*ust[ij - 1]*y;
      gamnorma = ((double) 1.0 + const1*znz) / ((double) 1.0 + const2*znz);
      tauhf[ij - 1] = tauhf[ij - 1] + fnc2*gamnorma;
    }
    tauhf[ij - 1] = f1dcos3*consttau*taul*tauhf[ij - 1]*delz;
  }
  phihf[ij - 1] = (double) 0.0;
  if (llphihf) {
    taul = pow(ustph, 2);
    zsup = zsupmax;
    delz = 
      max((double) ((zsup - zinf) / (double) (jtot_tauhf - 1)), (double) ((double) 0.0));
    
    constphi = 
      aird[ij - 1 + kijl*(ichnk - 1)]*zpi4gm1*fr5[mij[ij - 1 + kijl*(ichnk - 1)] - 1];
    if (ltauwshelter) {
      for (j = 1; j <= jtot_tauhf; j += 1) {
        y = exp((double) (zinf + (double) (j - 1)*delz));
        omega = y*sqrtgz0;
        cm1 = omega*gm1;
        zx = ustph*cm1 + zalp;
        zarg = xkappa / zx;
        zlog = xloggz0 + (double) 2.0*log(cm1) + zarg;
        zlog = min((double) (zlog), (double) ((double) 0.0));
        zbeta = (pow(zlog, 4))*exp((double) (zlog));
        znz = zbeta*ust[ij - 1]*y;
        gamnorma = ((double) 1.0 + const1*znz) / ((double) 1.0 + const2*znz);
        fnc2 = zbeta*taul*wtauhf[j - 1]*delz*gamnorma;
        taul = max((double) (taul - tauwshelter*f1dcos3*consttau*fnc2), (double) ((double
          ) 0.0));
        ustph = sqrt((double) (taul));
        phihf[ij - 1] = phihf[ij - 1] + fnc2 / y;
      }
      phihf[ij - 1] = f1dcos2*constphi*sqrtz0og*phihf[ij - 1];
    } else {
      for (j = 1; j <= jtot_tauhf; j += 1) {
        y = exp((double) (zinf + (double) (j - 1)*delz));
        omega = y*sqrtgz0;
        cm1 = omega*gm1;
        zx = ustph*cm1 + zalp;
        zarg = xkappa / zx;
        zlog = xloggz0 + (double) 2.0*log(cm1) + zarg;
        zlog = min((double) (zlog), (double) ((double) 0.0));
        zbeta = (pow(zlog, 4))*exp((double) (zlog));
        znz = zbeta*ust[ij - 1]*y;
        gamnorma = ((double) 1.0 + const1*znz) / ((double) 1.0 + const2*znz);
        fnc2 = zbeta*wtauhf[j - 1]*gamnorma;
        phihf[ij - 1] = phihf[ij - 1] + fnc2 / y;
      }
      phihf[ij - 1] = f1dcos2*constphi*sqrtz0og*taul*phihf[ij - 1]*delz;
    }
  }

  
  
}
