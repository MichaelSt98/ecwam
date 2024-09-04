#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sinput_ard_c.h"
#include "wsigstar_c.h"


__device__ void sinput_ard_c(int ngst, int llsneg, int kijs, int kijl, 
  const double * __restrict__ fl1, const double * __restrict__ wavnum, 
  const double * __restrict__ cinv, const double * __restrict__ xk2cg, 
  const double * __restrict__ wdwave, const double * __restrict__ wswave, 
  const double * __restrict__ ufric, const double * __restrict__ z0m, 
  const double * __restrict__ coswdif, const double * __restrict__ sinwdif2, 
  const double * __restrict__ raorw, const double * __restrict__ wstar, 
  const double * __restrict__ rnfac, double * __restrict__ fld, 
  double * __restrict__ sl, double * __restrict__ spos, double * __restrict__ xllws, 
  double abmax, double abmin, double acdlin, double alphamax, double alphamin, 
  double bcdlin, double betamaxoxkappa2, const double * __restrict__ costh, 
  double delth, const double * __restrict__ dfim, double epsmin, double epsus, double g, 
  int iab, int llgcbz0, int llnormagam, int nang, int nfre, double rnu, double rnum, 
  const double * __restrict__ sinth, double swellf, double swellf2, double swellf3, 
  double swellf4, double swellf5, double swellf6, double swellf7, double swellf7m1, 
  const double * __restrict__ swellft, double tauwshelter, 
  const double * __restrict__ th, double wspmin, double xkappa, double z0rat, 
  double z0tubmax, double zalp, double zpi, const double * __restrict__ zpifr, 
  int ichnk, int nchnk, int ij, double * __restrict__ constf, 
  double * __restrict__ const11, double * __restrict__ const22, 
  double * __restrict__ z0vis, double * __restrict__ z0noz, double * __restrict__ fww, 
  double * __restrict__ pvisc, double * __restrict__ pturb, double * __restrict__ zcn, 
  double * __restrict__ sig_n, double * __restrict__ uorbt, double * __restrict__ aorb, 
  double * __restrict__ temp, double * __restrict__ re, double * __restrict__ re_c, 
  double * __restrict__ zorb, double * __restrict__ cnsn, double * __restrict__ sumf, 
  double * __restrict__ sumfsin2, double * __restrict__ cstrnfac, 
  double * __restrict__ flp_avg, double * __restrict__ slp_avg, 
  double * __restrict__ rogoroair, double * __restrict__ aird_pvisc, 
  double * __restrict__ xstress, double * __restrict__ ystress, 
  double * __restrict__ flp, double * __restrict__ slp, double * __restrict__ usg2, 
  double * __restrict__ taux, double * __restrict__ tauy, double * __restrict__ ustp, 
  double * __restrict__ ustpm1, double * __restrict__ usdirp, double * __restrict__ ucn, 
  double * __restrict__ ucnzalpd, double * __restrict__ xngamconst, 
  double * __restrict__ gamnorma, double * __restrict__ dstab1, 
  double * __restrict__ temp1, double * __restrict__ temp2, double * __restrict__ gam0, 
  double * __restrict__ dstab, double * __restrict__ coslp) {
  
  
  
  // ----------------------------------------------------------------------
  

  
  
  int k;
  int m;
  int ind;
  int igst;
  
  double constn;
  double avg_gst;
  double abs_tauwshelter;
  double const1;
  double znz;
  double x1;
  double x2;
  double zlog;
  double zlog1;
  double zlog2;
  double zlog2x;
  double xv1;
  double xv2;
  double zbeta1;
  double zbeta2;
  double xi;
  double x;
  double deli1;
  double deli2;
  double fu;
  double fud;
  double nu_air;
  double smooth;
  double hftswellf6;
  double z0tub;
  double fac_nu_air;
  double facm1_nu_air;
  double arg;
  double delabm1;
  double taupx;
  double taupy;
  double dstab2;
  double zhook_handle;
  double const_var;
  double sig;
  double sig2;
  double coef;
  double coef5;
  double dfim_sig2;
  
  int ltauwshelter;
  int i_dstab_2;
  int i_dstab_1;
  int i_gamnorma_1;
  // ----------------------------------------------------------------------
  
  
  avg_gst = (double) 1.0 / ngst;
  const1 = betamaxoxkappa2;
  constn = delth / (xkappa*zpi);
  
  abs_tauwshelter = fabs((double) (tauwshelter));
  if (abs_tauwshelter == (double) 0.0) {
    ltauwshelter = false;
  } else {
    ltauwshelter = true;
  }
  
  
  //     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
  if (ngst > 1) {
    wsigstar_c(kijs, kijl, wswave, ufric, z0m, wstar,  (&sig_n[ + kijl*(ichnk - 1)]), 
      acdlin, alphamax, alphamin, bcdlin, epsus, g, llgcbz0, rnum, wspmin, xkappa, 
      ichnk, nchnk, ij);
  }
  
  
  
  if (llnormagam) {
    cstrnfac[ij - 1 + kijl*(ichnk - 1)] = constn*rnfac[ij - 1] / raorw[ij - 1];
  }
  // ----------------------------------------------------------------------
  
  if (llsneg) {
    //!!!  only for the negative sinput
    nu_air = rnu;
    facm1_nu_air = (double) 4.0 / nu_air;
    
    fac_nu_air = rnum;
    
    fu = fabs((double) (swellf3));
    fud = swellf2;
    delabm1 = (double) (iab) / (abmax - abmin);
    
    
    //       computation of Uorb and Aorb
    uorbt[ij - 1 + kijl*(ichnk - 1)] = epsmin;
    aorb[ij - 1 + kijl*(ichnk - 1)] = epsmin;
    
    for (m = 1; m <= nfre; m += 1) {
      sig = zpifr[m - 1];
      sig2 = pow(sig, 2);
      dfim_sig2 = dfim[m - 1]*sig2;
      
      k = 1;
      temp[ij - 1 + kijl*(ichnk - 1)] = 
        fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
      for (k = 2; k <= nang; k += 1) {
        temp[ij - 1 + kijl*(ichnk - 1)] = temp[ij - 1 + kijl*(ichnk - 1)] + fl1[ij - 1 + 
          kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
      }
      
      uorbt[ij - 1 + kijl*(ichnk - 1)] = 
        uorbt[ij - 1 + kijl*(ichnk - 1)] + dfim_sig2*temp[ij - 1 + kijl*(ichnk - 1)];
      aorb[ij - 1 + kijl*(ichnk - 1)] = 
        aorb[ij - 1 + kijl*(ichnk - 1)] + dfim[m - 1]*temp[ij - 1 + kijl*(ichnk - 1)];
    }
    
    uorbt[ij - 1 + kijl*(ichnk - 1)] = 
      (double) 2.0*sqrt((double) (uorbt[ij - 1 + kijl*(ichnk - 1)]));      // this is the significant orbital amplitude
    aorb[ij - 1 + kijl*(ichnk - 1)] = 
      (double) 2.0*sqrt((double) (aorb[ij - 1 + kijl*(ichnk - 1)]));      // this 1/2 Hs
    re[ij - 1 + kijl*(ichnk - 1)] = 
      facm1_nu_air*uorbt[ij - 1 + kijl*(ichnk - 1)]*aorb[ij - 1 + kijl*(ichnk - 1)];      // this is the Reynolds number
    z0vis[ij - 1 + kijl*(ichnk - 1)] = fac_nu_air / fmax((double) (ufric[ij - 1 + 
      kijl*(ichnk - 1)]), (double) ((double) 0.0001));
    z0tub = z0rat*fmin((double) (z0tubmax), (double) (z0m[ij - 1 + kijl*(ichnk - 1)]));
    z0noz[ij - 1 + kijl*(ichnk - 1)] = 
      fmax((double) (z0vis[ij - 1 + kijl*(ichnk - 1)]), (double) (z0tub));
    zorb[ij - 1 + kijl*(ichnk - 1)] = 
      aorb[ij - 1 + kijl*(ichnk - 1)] / z0noz[ij - 1 + kijl*(ichnk - 1)];
    
    //         compute fww
    xi = (log10(fmax((double) (zorb[ij - 1 + kijl*(ichnk - 1)]), (double) ((double) 3.0))
      ) - abmin)*delabm1;
    ind = fmin((double) (iab - 1), (double) ((int) (xi)));
    deli1 = fmin((double) ((double) 1.0), (double) (xi - (double) (ind)));
    deli2 = (double) 1.0 - deli1;
    fww[ij - 1 + kijl*(ichnk - 1)] = swellft[ind - 1]*deli2 + swellft[ind + 1 - 1]*deli1;
    temp2[ij - 1 + kijl*(ichnk - 1)] = 
      fww[ij - 1 + kijl*(ichnk - 1)]*uorbt[ij - 1 + kijl*(ichnk - 1)];
    
    //       Define the critical Reynolds number
    if (swellf6 == (double) 1.0) {
      re_c[ij - 1 + kijl*(ichnk - 1)] = swellf4;
    } else {
      hftswellf6 = (double) 1.0 - swellf6;
      re_c[ij - 1 + kijl*(ichnk - 1)] = 
        swellf4*(pow(((double) 2.0 / aorb[ij - 1 + kijl*(ichnk - 1)]), hftswellf6));
    }
    
    //       Swell damping weight between viscous and turbulent boundary layer
    if (swellf7 > (double) 0.0) {
      smooth = (double) 0.5*tanh((re[ij - 1 + kijl*(ichnk - 1)] - re_c[ij - 1 + 
        kijl*(ichnk - 1)])*swellf7m1);
      pturb[ij - 1 + kijl*(ichnk - 1)] = (double) 0.5 + smooth;
      pvisc[ij - 1 + kijl*(ichnk - 1)] = (double) 0.5 - smooth;
    } else {
      if (re[ij - 1 + kijl*(ichnk - 1)] <= re_c[ij - 1 + kijl*(ichnk - 1)]) {
        pturb[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
        pvisc[ij - 1 + kijl*(ichnk - 1)] = (double) 0.5;
      } else {
        pturb[ij - 1 + kijl*(ichnk - 1)] = (double) 0.5;
        pvisc[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
      }
    }
    
    aird_pvisc[ij - 1 + kijl*(ichnk - 1)] = 
      pvisc[ij - 1 + kijl*(ichnk - 1)]*raorw[ij - 1];
    
  }
  
  
  
  // Initialisation
  
  if (ngst == 1) {
    ustp[ij - 1 + kijl*(1 - 1 + 2*(ichnk - 1))] = ufric[ij - 1 + kijl*(ichnk - 1)];
  } else if (ngst == 2) {
    ustp[ij - 1 + kijl*(1 - 1 + 2*(ichnk - 1))] = 
      ufric[ij - 1 + kijl*(ichnk - 1)]*((double) 1.0 + sig_n[ij - 1 + kijl*(ichnk - 1)]);
    ustp[ij - 1 + kijl*(2 - 1 + 2*(ichnk - 1))] = 
      ufric[ij - 1 + kijl*(ichnk - 1)]*((double) 1.0 - sig_n[ij - 1 + kijl*(ichnk - 1)]);
  }
  
  for (igst = 1; igst <= ngst; igst += 1) {
    ustpm1[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = (double) 1.0 / fmax((double) 
      (ustp[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))]), (double) (epsus));
  }
  
  if (ltauwshelter) {
    for (igst = 1; igst <= ngst; igst += 1) {
      xstress[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = (double) 0.0;
      ystress[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = (double) 0.0;
      usg2[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = 
        pow(ustp[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))], 2);
      taux[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = usg2[ij - 1 + kijl*(igst - 1 + 
        2*(ichnk - 1))]*sin(wdwave[ij - 1 + kijl*(ichnk - 1)]);
      tauy[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = usg2[ij - 1 + kijl*(igst - 1 + 
        2*(ichnk - 1))]*cos(wdwave[ij - 1 + kijl*(ichnk - 1)]);
    }
    
    rogoroair[ij - 1 + kijl*(ichnk - 1)] = g / raorw[ij - 1];
    
  } else {
    for (k = 1; k <= nang; k += 1) {
      coslp[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] = coswdif[ij - 1 + kijl*(k - 1)];
    }
  }
  
  
  //*    2. MAIN LOOP OVER FREQUENCIES.
  //        ---------------------------
  
  if (!llnormagam) {
    for (i_gamnorma_1 = 1; i_gamnorma_1 <= 2; i_gamnorma_1 += 1) {
      gamnorma[ij - 1 + kijl*(i_gamnorma_1 - 1 + 2*(ichnk - 1))] = (double) 1.0;
    }
  }
  
  if (!llsneg) {
    for (i_dstab_2 = 1; i_dstab_2 <= 2; i_dstab_2 += 1) {
      for (i_dstab_1 = 1; i_dstab_1 <= nang; i_dstab_1 += 1) {
        dstab[ij - 1 + kijl*(i_dstab_1 - 1 + nang*(i_dstab_2 - 1 + 2*(ichnk - 1)))] = 
          (double) 0.0;
      }
    }
  }
  
  for (m = 1; m <= nfre; m += 1) {
    sig = zpifr[m - 1];
    sig2 = pow(sig, 2);
    const_var = sig*const1;
    
    if (llsneg) {
      coef = -swellf*(double) 16.*sig2 / g;
      coef5 = -swellf5*(double) 2.*sqrt((double) ((double) 2.*nu_air*sig));
    }
    
    if (ltauwshelter) {
      for (igst = 1; igst <= ngst; igst += 1) {
        taupx = taux[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] - 
          abs_tauwshelter*xstress[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))];
        taupy = tauy[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] - 
          abs_tauwshelter*ystress[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))];
        usdirp[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = atan2(taupx, taupy);
        ustp[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = 
          pow(((pow(taupx, 2)) + (pow(taupy, 2))), (double) 0.25);
        ustpm1[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = (double) 1.0 / fmax((double) 
          (ustp[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))]), (double) (epsus));
      }
      
      constf[ij - 1 + kijl*(ichnk - 1)] = rogoroair[ij - 1 + kijl*(ichnk - 1)]*cinv[ij - 
        1 + kijl*(m - 1 + nfre*(ichnk - 1))]*dfim[m - 1];
    }
    
    
    //*      PRECALCULATE FREQUENCY DEPENDENCE.
    //       ----------------------------------
    
    for (igst = 1; igst <= ngst; igst += 1) {
      ucn[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = ustp[ij - 1 + kijl*(igst - 1 + 
        2*(ichnk - 1))]*cinv[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))];
      ucnzalpd[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = 
        xkappa / (ucn[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] + zalp);
    }
    zcn[ij - 1 + kijl*(ichnk - 1)] = log(wavnum[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))
      ]*z0m[ij - 1 + kijl*(ichnk - 1)]);
    cnsn[ij - 1 + kijl*(ichnk - 1)] = const_var*raorw[ij - 1];
    
    //*    2.1 LOOP OVER DIRECTIONS.
    //         ---------------------
    
    for (k = 1; k <= nang; k += 1) {
      xllws[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = (double) 0.0;
    }
    
    if (llnormagam) {
      xngamconst[ij - 1 + kijl*(ichnk - 1)] = cstrnfac[ij - 1 + kijl*(ichnk - 1)
        ]*xk2cg[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))];
    }
    
    if (llsneg) {
      //       SWELL DAMPING:
      dstab1[ij - 1 + kijl*(ichnk - 1)] = coef5*aird_pvisc[ij - 1 + kijl*(ichnk - 1)
        ]*wavnum[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))];
      temp1[ij - 1 + kijl*(ichnk - 1)] = coef*raorw[ij - 1];
    }
    
    for (igst = 1; igst <= ngst; igst += 1) {
      for (k = 1; k <= nang; k += 1) {
        if (ltauwshelter) {
          coslp[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] = 
            cos(th[k - 1] - usdirp[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))]);
        }
        gam0[ij - 1 + kijl*(k - 1 + nang*(igst - 1 + 2*(ichnk - 1)))] = (double) 0.0;
        if (coslp[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] > (double) 0.01) {
          x = coslp[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))]*ucn[ij - 1 + kijl*(igst - 1
             + 2*(ichnk - 1))];
          zlog = zcn[ij - 1 + kijl*(ichnk - 1)] + ucnzalpd[ij - 1 + kijl*(igst - 1 + 
            2*(ichnk - 1))] / coslp[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))];
          if (zlog < (double) 0.0) {
            zlog2x = zlog*zlog*x;
            gam0[ij - 1 + kijl*(k - 1 + nang*(igst - 1 + 2*(ichnk - 1)))] = 
              exp((double) (zlog))*zlog2x*zlog2x*cnsn[ij - 1 + kijl*(ichnk - 1)];
            xllws[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = (double) 1.0
              ;
          }
        }
      }
      
      if (llnormagam) {
        
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
      
      if (llsneg) {
        for (k = 1; k <= nang; k += 1) {
          dstab2 = temp1[ij - 1 + kijl*(ichnk - 1)]*(temp2[ij - 1 + kijl*(ichnk - 1)] + 
            (fu + fud*coslp[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))])*ustp[ij - 1 + 
            kijl*(igst - 1 + 2*(ichnk - 1))]);
          dstab[ij - 1 + kijl*(k - 1 + nang*(igst - 1 + 2*(ichnk - 1)))] = 
            dstab1[ij - 1 + kijl*(ichnk - 1)] + pturb[ij - 1 + kijl*(ichnk - 1)]*dstab2;
        }
      }
    }
    
    
    //*    2.2 UPDATE THE SHELTERING STRESS (in any),
    //         AND THEN ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
    //         ---------------------------------------------------------
    
    for (k = 1; k <= nang; k += 1) {
      
      for (igst = 1; igst <= ngst; igst += 1) {
        // SLP: only the positive contributions
        slp[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = gam0[ij - 1 + kijl*(k - 1 + 
          nang*(igst - 1 + 2*(ichnk - 1)))]*gamnorma[ij - 1 + kijl*(igst - 1 + 2*(ichnk -
           1))];
        flp[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = slp[ij - 1 + kijl*(igst - 1 + 
          2*(ichnk - 1))] + dstab[ij - 1 + kijl*(k - 1 + nang*(igst - 1 + 2*(ichnk - 1)))
          ];
      }
      
      for (igst = 1; igst <= ngst; igst += 1) {
        slp[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = slp[ij - 1 + kijl*(igst - 1 + 
          2*(ichnk - 1))]*fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
      }
      
      if (ltauwshelter) {
        const11[ij - 1 + kijl*(ichnk - 1)] = 
          constf[ij - 1 + kijl*(ichnk - 1)]*sinth[k - 1];
        const22[ij - 1 + kijl*(ichnk - 1)] = 
          constf[ij - 1 + kijl*(ichnk - 1)]*costh[k - 1];
        for (igst = 1; igst <= ngst; igst += 1) {
          xstress[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = xstress[ij - 1 + kijl*(igst
             - 1 + 2*(ichnk - 1))] + slp[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))
            ]*const11[ij - 1 + kijl*(ichnk - 1)];
          ystress[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))] = ystress[ij - 1 + kijl*(igst
             - 1 + 2*(ichnk - 1))] + slp[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))
            ]*const22[ij - 1 + kijl*(ichnk - 1)];
        }
      }
      
      igst = 1;
      slp_avg[ij - 1 + kijl*(ichnk - 1)] = slp[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))];
      flp_avg[ij - 1 + kijl*(ichnk - 1)] = flp[ij - 1 + kijl*(igst - 1 + 2*(ichnk - 1))];
      for (igst = 2; igst <= ngst; igst += 1) {
        slp_avg[ij - 1 + kijl*(ichnk - 1)] = slp_avg[ij - 1 + kijl*(ichnk - 1)] + slp[ij 
          - 1 + kijl*(igst - 1 + 2*(ichnk - 1))];
        flp_avg[ij - 1 + kijl*(ichnk - 1)] = flp_avg[ij - 1 + kijl*(ichnk - 1)] + flp[ij 
          - 1 + kijl*(igst - 1 + 2*(ichnk - 1))];
      }
      
      spos[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = 
        avg_gst*slp_avg[ij - 1 + kijl*(ichnk - 1)];
      
      fld[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = 
        avg_gst*flp_avg[ij - 1 + kijl*(ichnk - 1)];
      sl[ij - 1 + kijl*(k - 1 + nang*(m - 1))] = fld[ij - 1 + kijl*(k - 1 + nang*(m - 1))
        ]*fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
      
    }
    
  }
  
  // END LOOP OVER FREQUENCIES
  
  
}
