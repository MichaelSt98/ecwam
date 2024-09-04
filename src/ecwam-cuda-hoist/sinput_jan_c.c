#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sinput_jan_c.h"
#include "wsigstar_c.h"


__device__ void sinput_jan_c(int ngst, int llsneg, int kijs, int kijl, 
  const double * __restrict__ fl1, const double * __restrict__ wavnum, 
  const double * __restrict__ cinv, const double * __restrict__ xk2cg, 
  const double * __restrict__ wswave, const double * __restrict__ ufric, 
  const double * __restrict__ z0m, const double * __restrict__ coswdif, 
  const double * __restrict__ sinwdif2, const double * __restrict__ raorw, 
  const double * __restrict__ wstar, const double * __restrict__ rnfac, 
  double * __restrict__ fld, double * __restrict__ sl, double * __restrict__ spos, 
  double * __restrict__ xllws, double acdlin, double alphamax, double alphamin, 
  double bcdlin, double betamaxoxkappa2, double delth, double epsus, double g, 
  int idamping, int llgcbz0, int llnormagam, int nang, int nfre, double rnum, 
  double wspmin, double xkappa, double zalp, double zpi, 
  const double * __restrict__ zpifr, int ichnk, int nchnk, int ij, 
  double * __restrict__ ztanhkd, double * __restrict__ sig_n, 
  double * __restrict__ cnsn, double * __restrict__ sumf, 
  double * __restrict__ sumfsin2, double * __restrict__ cstrnfac, 
  double * __restrict__ xngamconst, double * __restrict__ gamnorma, 
  double * __restrict__ sigdev, double * __restrict__ us, double * __restrict__ z0, 
  double * __restrict__ ucn, double * __restrict__ zcn, double * __restrict__ ustpm1, 
  double * __restrict__ xvd, double * __restrict__ ucnd, 
  double * __restrict__ const3_ucn2, double * __restrict__ ufac1, 
  double * __restrict__ ufac2, double * __restrict__ tempd, double * __restrict__ gam0, 
  int * __restrict__ lz) {
  
  
  
  // ----------------------------------------------------------------------
  


  
  
  int ig;
  int k;
  int m;
  int igst;
  
  double const1;
  double const3;
  double xkappad;
  double constn;
  double const_var;
  double znz;
  double x;
  double zlog;
  double zlog2x;
  double zbeta;
  double zhook_handle;
  double wsin[2];
  
  int i_gamnorma_1;
  int i_ufac2_1;
  
  // ----------------------------------------------------------------------
  
  
  const1 = betamaxoxkappa2;
  const3 = (double) 2.0*xkappa / const1;    // SEE IDAMPING
  xkappad = (double) 1.E0 / xkappa;
  
  const3 = idamping*const3;
  
  constn = delth / (xkappa*zpi);
  
  //     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
  if (ngst > 1) {
    wsigstar_c(kijs, kijl, wswave, ufric, z0m, wstar,  (&sig_n[ + kijl*(ichnk - 1)]), 
      acdlin, alphamax, alphamin, bcdlin, epsus, g, llgcbz0, rnum, wspmin, xkappa, 
      ichnk, nchnk, ij);
  }
  
  //*    1. PRECALCULATED ANGULAR DEPENDENCE.
  //        ---------------------------------
  
  
  for (k = 1; k <= nang; k += 1) {
    if (coswdif[ij - 1 + kijl*(k - 1)] > (double) 0.01) {
      lz[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] = true;
      tempd[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] = 
        xkappa / coswdif[ij - 1 + kijl*(k - 1)];
    } else {
      lz[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] = false;
      tempd[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] = xkappa;
    }
  }
  
  if (llnormagam) {
    cstrnfac[ij - 1 + kijl*(ichnk - 1)] = constn*rnfac[ij - 1] / raorw[ij - 1];
  }
  
  
  //     DEFINE WHERE SINPUT WILL BE EVALUATED IN RELATIVE TERM WRT USTAR
  //     DEFINE ALSO THE RELATIVE WEIGHT OF EACH.
  
  if (ngst == 1) {
    wsin[1 - 1] = (double) 1.0;
    sigdev[ij - 1 + kijl*(1 - 1 + 2*(ichnk - 1))] = (double) 1.0;
  } else if (ngst == 2) {
    wsin[1 - 1] = (double) 0.5;
    wsin[2 - 1] = (double) 0.5;
    sigdev[ij - 1 + kijl*(1 - 1 + 2*(ichnk - 1))] = 
      (double) 1.0 - sig_n[ij - 1 + kijl*(ichnk - 1)];
    sigdev[ij - 1 + kijl*(2 - 1 + 2*(ichnk - 1))] = 
      (double) 1.0 + sig_n[ij - 1 + kijl*(ichnk - 1)];
  }
  
  
  if (ngst == 1) {
    us[ij - 1 + kijl*(1 - 1 + 2*(ichnk - 1))] = ufric[ij - 1 + kijl*(ichnk - 1)];
    z0[ij - 1 + kijl*(1 - 1 + 2*(ichnk - 1))] = z0m[ij - 1 + kijl*(ichnk - 1)];
  } else {
    for (igst = 1; igst <= ngst; igst += 1) {
      us[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = 
        ufric[ij - 1 + kijl*(ichnk - 1)]*sigdev[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))]
        ;
      z0[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = z0m[ij - 1 + kijl*(ichnk - 1)];
    }
  }
  
  for (igst = 1; igst <= ngst; igst += 1) {
    ustpm1[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = (double) 1.0 / fmax((double) 
      (us[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))]), (double) (epsus));
  }
  
  // ----------------------------------------------------------------------
  
  //*    2. LOOP OVER FREQUENCIES.
  //        ----------------------
  
  if (!llnormagam) {
    for (i_gamnorma_1 = 1; i_gamnorma_1 <= 2; i_gamnorma_1 += 1) {
      gamnorma[ij - 1 + kijl*(i_gamnorma_1 - 1 + 2*(ichnk - 1))] = (double) 1.0;
    }
  }
  
  if (!llsneg) {
    for (i_ufac2_1 = 1; i_ufac2_1 <= nang; i_ufac2_1 += 1) {
      ufac2[ij - 1 + kijl*(i_ufac2_1 - 1 + nang*(ichnk - 1))] = (double) 0.0;
    }
  }
  
  for (m = 1; m <= nfre; m += 1) {
    
    const_var = zpifr[m - 1]*const1;
    
    //*      PRECALCULATE FREQUENCY DEPENDENCE.
    //       ----------------------------------
    
    ztanhkd[ij - 1 + kijl*(ichnk - 1)] = 
      (pow(zpifr[m - 1], 2)) / (g*wavnum[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))]);
    
    cnsn[ij - 1 + kijl*(ichnk - 1)] = 
      const_var*ztanhkd[ij - 1 + kijl*(ichnk - 1)]*raorw[ij - 1];
    
    for (igst = 1; igst <= ngst; igst += 1) {
      ucn[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = us[ij - 1 + kijl*(igst - 1 + 
        2*(ichnk - 1))]*cinv[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))] + zalp;
      const3_ucn2[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = 
        const3*(pow(ucn[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))], 2));
      ucnd[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = 
        (double) 1.0 / ucn[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))];
      zcn[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = log(wavnum[ij - 1 + kijl*(m - 1 + 
        nfre*(ichnk - 1))]*z0[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))]);
      xvd[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = (double) 1.0 / (-us[ij - 1 + 
        kijl*(igst - 1 + 2*(ichnk - 1))]*xkappad*zcn[ij - 1 + kijl*(igst - 1 + 2*(ichnk -
         1))]*cinv[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))]);
    }
    
    //*    2.1 LOOP OVER DIRECTIONS.
    //         ---------------------
    
    //       WIND INPUT:
    for (k = 1; k <= nang; k += 1) {
      xllws[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = (double) 0.0;
      
      for (igst = 1; igst <= ngst; igst += 1) {
        if (lz[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))]) {
          zlog = zcn[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] + tempd[ij - 1 + kijl*(k -
             1 + nang*(ichnk - 1))]*ucnd[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))];
          if (zlog < (double) 0.0) {
            x = coswdif[ij - 1 + kijl*(k - 1)]*ucn[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1
              ))];
            zlog2x = zlog*zlog*x;
            gam0[ij - 1 + kijl*(k - 1 + nang*(igst - 1 + 2*(ichnk - 1)))] = 
              zlog2x*zlog2x*exp((double) (zlog))*cnsn[ij - 1 + kijl*(ichnk - 1)];
            xllws[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = (double) 1.0
              ;
          } else {
            gam0[ij - 1 + kijl*(k - 1 + nang*(igst - 1 + 2*(ichnk - 1)))] = (double) 0.0;
          }
        } else {
          gam0[ij - 1 + kijl*(k - 1 + nang*(igst - 1 + 2*(ichnk - 1)))] = (double) 0.0;
        }
      }
      
    }
    
    
    if (llnormagam) {
      
      xngamconst[ij - 1 + kijl*(ichnk - 1)] = cstrnfac[ij - 1 + kijl*(ichnk - 1)
        ]*xk2cg[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))];
      
      for (igst = 1; igst <= ngst; igst += 1) {
        
        sumf[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
        sumfsin2[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
        for (k = 1; k <= nang; k += 1) {
          sumf[ij - 1 + kijl*(ichnk - 1)] = sumf[ij - 1 + kijl*(ichnk - 1)] + gam0[ij - 1
             + kijl*(k - 1 + nang*(igst - 1 + 2*(ichnk - 1)))]*fl1[ij - 1 + kijl*(k - 1 +
             nang*(m - 1 + nfre*(ichnk - 1)))];
          sumfsin2[ij - 1 + kijl*(ichnk - 1)] = sumfsin2[ij - 1 + kijl*(ichnk - 1)] + 
            gam0[ij - 1 + kijl*(k - 1 + nang*(igst - 1 + 2*(ichnk - 1)))]*fl1[ij - 1 + 
            kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))]*sinwdif2[ij - 1 + kijl*(k - 1
            )];
        }
        
        znz = xngamconst[ij - 1 + kijl*(ichnk - 1)]*ustpm1[ij - 1 + kijl*(igst - 1 + 
          2*(ichnk - 1))];
        gamnorma[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = ((double) 1.0 + 
          znz*sumfsin2[ij - 1 + kijl*(ichnk - 1)]) / ((double) 1.0 + znz*sumf[ij - 1 + 
          kijl*(ichnk - 1)]);
        
      }
      
    }
    
    
    for (k = 1; k <= nang; k += 1) {
      ufac1[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] = wsin[1 - 1]*gam0[ij - 1 + kijl*(k
         - 1 + nang*(1 - 1 + 2*(ichnk - 1)))]*gamnorma[ij - 1 + kijl*(1 - 1 + 2*(ichnk - 
        1))];
      for (igst = 2; igst <= ngst; igst += 1) {
        ufac1[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] = ufac1[ij - 1 + kijl*(k - 1 + 
          nang*(ichnk - 1))] + wsin[igst - 1]*gam0[ij - 1 + kijl*(k - 1 + nang*(igst - 1 
          + 2*(ichnk - 1)))]*gamnorma[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))];
      }
    }
    
    if (llsneg) {
      //         SWELL DAMPING:
      for (k = 1; k <= nang; k += 1) {
        for (igst = 1; igst <= 1; igst += 1) {
          zbeta = const3_ucn2[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))]*(coswdif[ij - 1 +
             kijl*(k - 1)] - xvd[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))]);
          ufac2[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] = wsin[igst - 1]*zbeta;
        }
        for (igst = 2; igst <= ngst; igst += 1) {
          zbeta = const3_ucn2[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))]*(coswdif[ij - 1 +
             kijl*(k - 1)] - xvd[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))]);
          ufac2[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] = 
            ufac2[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] + wsin[igst - 1]*zbeta;
        }
      }
    }
    
    //*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
    //         ------------------------------------------------
    
    for (k = 1; k <= nang; k += 1) {
      fld[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = ufac1[ij - 1 + kijl*(k - 1 + 
        nang*(ichnk - 1))] + ufac2[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))]*cnsn[ij - 1 
        + kijl*(ichnk - 1)];
      spos[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = ufac1[ij - 1 + kijl*(k - 1 + 
        nang*(ichnk - 1))]*fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
      sl[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = fld[ij - 1 + kijl*(k - 1 + nang*(m - 1))
        ]*fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
    }
    
  }
  
  
  
}
