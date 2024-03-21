#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sinput_jan_c.h"
#include "wsigstar_c.c"

__device__ void sinput_jan_c(int ngst, int llsneg, int kijs, int kijl, 
  const double * fl1, const double * wavnum, const double * cinv, const double * xk2cg, 
  const double * wswave, const double * ufric, const double * z0m, 
  const double * coswdif, const double * sinwdif2, const double * raorw, 
  const double * wstar, const double * rnfac, double * fld, double * sl, double * spos, 
  double * xllws, double acdlin, double alphamax, double alphamin, double bcdlin, 
  double betamaxoxkappa2, double delth, double epsus, double g, int idamping, 
  int llgcbz0, int llnormagam, int nang, int nfre, double rnum, double wspmin, 
  double xkappa, double zalp, double zpi, const double * zpifr, int ichnk, int nchnk, 
  int ij) {
  
  // Loki: parameters from YOWPARAM inlined
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  int ig;
  int k;
  int m;
  int igst;
  
  double const1;
  double const3;
  double xkappad;
  double constn;
  double znz;
  double x;
  double zlog;
  double zlog2x;
  double zbeta;
  double tempd;
  
  double wsin[2];
  double ztanhkd;
  double sig_n;
  double cnsn;
  double sumf;
  double sumfsin2;
  double cstrnfac;
  double ufac1;
  double ufac2;
  double gamnorma[2];  // ! RENORMALISATION FACTOR OF THE GROWTH RATE
  double sigdev[2];
  double us[2];
  double z0[2];
  double ucn[2];
  double zcn[2];
  double ustpm1[2];
  double xvd[2];
  double ucnd[2];
  double const3_ucn2[2];
  
  int lz;
  double gam0[36*2];
  
  const1 = betamaxoxkappa2;
  const3 = (double) 2.0*xkappa / const1;    // SEE IDAMPING
  xkappad = (double) 1.E0 / xkappa;
  
  const3 = idamping*const3;
  
  constn = delth / (xkappa*zpi);

  if (ngst > 1) {
    wsigstar_c(wswave[ij - 1 + kijl*(ichnk - 1)], ufric[ij - 1 + kijl*(ichnk - 1)], 
      z0m[ij - 1 + kijl*(ichnk - 1)], wstar[ij - 1 + kijl*(ichnk - 1)],  (&sig_n), 
      acdlin, alphamax, alphamin, bcdlin, epsus, g, llgcbz0, rnum, wspmin, xkappa);
  }
  if (ngst == 1) {
    wsin[1 - 1] = (double) 1.0;
    sigdev[1 - 1] = (double) 1.0;
  } else {
    wsin[1 - 1] = (double) 0.5;
    wsin[2 - 1] = (double) 0.5;
    sigdev[1 - 1] = (double) 1.0 - sig_n;
    sigdev[2 - 1] = (double) 1.0 + sig_n;
  }
  if (ngst == 1) {
    us[1 - 1] = ufric[ij - 1 + kijl*(ichnk - 1)];
    z0[1 - 1] = z0m[ij - 1 + kijl*(ichnk - 1)];
  } else {
    for (igst = 1; igst <= ngst; igst += 1) {
      us[igst - 1] = ufric[ij - 1 + kijl*(ichnk - 1)]*sigdev[igst - 1];
      z0[igst - 1] = z0m[ij - 1 + kijl*(ichnk - 1)];
    }
  }
  
  for (igst = 1; igst <= ngst; igst += 1) {
    ustpm1[igst - 1] = (double) 1.0 / max((double) (us[igst - 1]), (double) (epsus));
  }
  for (m = 1; m <= nfre; m += 1) {
    ztanhkd = (pow(zpifr[m - 1], 2)) / (g*wavnum[ij - 1 + kijl*(m - 1 + 
      nfre_loki_param*(ichnk - 1))]);
    cnsn = const1*zpifr[m - 1]*ztanhkd*raorw[ij - 1];
    
    for (igst = 1; igst <= ngst; igst += 1) {
      ucn[igst - 1] = 
        us[igst - 1]*cinv[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))] + zalp;
      const3_ucn2[igst - 1] = const3*(pow(ucn[igst - 1], 2));
      ucnd[igst - 1] = (double) 1.0 / ucn[igst - 1];
      zcn[igst - 1] = 
        log(wavnum[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))]*z0[igst - 1]);
      xvd[igst - 1] = (double) 1.0 / (-us[igst - 1]*xkappad*zcn[igst - 1]*cinv[ij - 1 + 
        kijl*(m - 1 + nfre_loki_param*(ichnk - 1))]);
    }
    for (k = 1; k <= nang; k += 1) {
      xllws[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))
        ] = (double) 0.0;
      
      for (igst = 1; igst <= ngst; igst += 1) {
        
        if (coswdif[ij - 1 + kijl*(k - 1)] > (double) 0.01) {
          lz = true;
          tempd = xkappa / coswdif[ij - 1 + kijl*(k - 1)];
        } else {
          lz = false;
          tempd = xkappa;
        }
        
        gam0[igst - 1 + 2*(k - 1)] = (double) 0.0;
        if (lz) {
          zlog = zcn[igst - 1] + tempd*ucnd[igst - 1];
          if (zlog < (double) 0.0) {
            x = coswdif[ij - 1 + kijl*(k - 1)]*ucn[igst - 1];
            zlog2x = zlog*zlog*x;
            gam0[igst - 1 + 2*(k - 1)] = zlog2x*zlog2x*exp((double) (zlog))*cnsn;
            xllws[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk 
              - 1)))] = (double) 1.0;
          }
        }
      }
      
    }
    if (llnormagam) {
      
      sumf = (double) 0.0;
      sumfsin2 = (double) 0.0;
      for (k = 1; k <= nang; k += 1) {
        for (igst = 1; igst <= ngst; igst += 1) {
          sumf = sumf + gam0[igst - 1 + 2*(k - 1)]*fl1[ij - 1 + kijl*(k - 1 + 
            nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))];
          sumfsin2 = sumfsin2 + gam0[igst - 1 + 2*(k - 1)]*fl1[ij - 1 + kijl*(k - 1 + 
            nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))]*sinwdif2[ij - 1 + 
            kijl*(k - 1)];
        }
        
        cstrnfac = constn*rnfac[ij - 1] / raorw[ij - 1];
        znz = cstrnfac*xk2cg[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))
          ]*ustpm1[igst - 1];
        gamnorma[igst - 1] = ((double) 1.0 + znz*sumfsin2) / ((double) 1.0 + znz*sumf);
        
      }
    } else {
      for (igst = 1; igst <= ngst; igst += 1) {
        gamnorma[igst - 1] = (double) 1.0;
      }
    }
    
    for (k = 1; k <= nang; k += 1) {
      ufac1 = wsin[1 - 1]*gam0[1 - 1 + 2*(k - 1)]*gamnorma[1 - 1];
      for (igst = 2; igst <= ngst; igst += 1) {
        ufac1 = ufac1 + wsin[igst - 1]*gam0[igst - 1 + 2*(k - 1)]*gamnorma[igst - 1];
      }
      
      ufac2 = (double) 0.0;
      if (llsneg) {
        //         SWELL DAMPING:
        zbeta = const3_ucn2[1 - 1]*(coswdif[ij - 1 + kijl*(k - 1)] - xvd[1 - 1]);
        ufac2 = wsin[1 - 1]*zbeta;
        for (igst = 2; igst <= ngst; igst += 1) {
          zbeta = const3_ucn2[igst - 1]*(coswdif[ij - 1 + kijl*(k - 1)] - xvd[igst - 1]);
          ufac2 = ufac2 + wsin[igst - 1]*zbeta;
        }
      }
      
      fld[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] = ufac1 + ufac2*cnsn;
      spos[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] = ufac1*fl1[ij - 1 + kijl*(k 
        - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))];
      sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] = fld[ij - 1 + kijl*(k - 1 + 
        nang_loki_param*(m - 1))]*fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
        nfre_loki_param*(ichnk - 1)))];
    }
  }

  
  
}
