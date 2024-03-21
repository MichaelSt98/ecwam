#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sinput_ard_c.h"
#include "wsigstar_c.h"

__device__ void sinput_ard_c(int ngst, int llsneg, int kijs, int kijl, 
  const double * fl1, const double * wavnum, const double * cinv, const double * xk2cg, 
  const double * wdwave, const double * wswave, const double * ufric, 
  const double * z0m, const double * coswdif, const double * sinwdif2, 
  const double * raorw, const double * wstar, const double * rnfac, double * fld, 
  double * sl, double * spos, double * xllws, double abmax, double abmin, double acdlin, 
  double alphamax, double alphamin, double bcdlin, double betamaxoxkappa2, 
  const double * costh, double delth, const double * dfim, double epsmin, double epsus, 
  double g, int iab, int llgcbz0, int llnormagam, int nang, int nfre, double rnu, 
  double rnum, const double * sinth, double swellf, double swellf2, double swellf3, 
  double swellf4, double swellf5, double swellf6, double swellf7, double swellf7m1, 
  const double * swellft, double tauwshelter, const double * th, double wspmin, 
  double xkappa, double z0rat, double z0tubmax, double zalp, double zpi, 
  const double * zpifr, int ichnk, int nchnk, int ij) {
  
  // Loki: parameters from YOWPARAM inlined
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
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
  
  double sig2;
  double coef;
  double coef5;
  double dfim_sig2;
  double coslp;
  
  double xngamconst;
  double constf;
  double const11;
  double const22;
  double z0vis;
  double z0noz;
  double fww;
  double pvisc;
  double pturb;
  double zcn;
  double sig_n;
  double uorbt;
  double aorb;
  double temp;
  double re;
  double re_c;
  double zorb;
  double cnsn;
  double sumf;
  double sumfsin2;
  double cstrnfac;
  double flp_avg;
  double slp_avg;
  double rogoroair;
  double aird_pvisc;
  double dstab1;
  double temp1;
  double temp2;
  
  double xstress[2];
  double ystress[2];
  double flp[2];
  double slp[2];
  double usg2[2];
  double taux[2];
  double tauy[2];
  double ustp[2];
  double ustpm1[2];
  double usdirp[2];
  double ucn[2];
  double ucnzalpd[2];
  double gamnorma[2];  // ! RENORMALISATION FACTOR OF THE GROWTH RATE
  
  int ltauwshelter;
  double gam0[36*2];
  double dstab[36*2];
  
  avg_gst = (double) 1.0 / ngst;
  const1 = betamaxoxkappa2;
  constn = delth / (xkappa*zpi);
  
  abs_tauwshelter = abs((double) (tauwshelter));
  if (abs_tauwshelter == (double) 0.0) {
    ltauwshelter = false;
  } else {
    ltauwshelter = true;
  }
  

  if (ngst > 1) {
    wsigstar_c(wswave[ij - 1 + kijl*(ichnk - 1)], ufric[ij - 1 + kijl*(ichnk - 1)], 
      z0m[ij - 1 + kijl*(ichnk - 1)], wstar[ij - 1 + kijl*(ichnk - 1)],  (&sig_n), 
      acdlin, alphamax, alphamin, bcdlin, epsus, g, llgcbz0, rnum, wspmin, xkappa);
  }
  if (llnormagam) {
    cstrnfac = constn*rnfac[ij - 1] / raorw[ij - 1];
  }
  if (llsneg) {
    //!!!  only for the negative sinput
    nu_air = rnu;
    facm1_nu_air = (double) 4.0 / nu_air;
    
    fac_nu_air = rnum;
    
    fu = abs((double) (swellf3));
    fud = swellf2;
    delabm1 = (double) (iab) / (abmax - abmin);
    uorbt = epsmin;
    aorb = epsmin;
    
    for (m = 1; m <= nfre; m += 1) {
      sig2 = pow(zpifr[m - 1], 2);
      dfim_sig2 = dfim[m - 1]*sig2;
      
      k = 1;
      temp = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk -
         1)))];
      for (k = 2; k <= nang; k += 1) {
        temp = temp + fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
          nfre_loki_param*(ichnk - 1)))];
      }
      
      uorbt = uorbt + dfim_sig2*temp;
      aorb = aorb + dfim[m - 1]*temp;
    }
    
    uorbt = (double) 2.0*sqrt((double) (uorbt));      // this is the significant orbital amplitude
    aorb = (double) 2.0*sqrt((double) (aorb));      // this 1/2 Hs
    re = facm1_nu_air*uorbt*aorb;      // this is the Reynolds number
    z0vis = fac_nu_air / max((double) (ufric[ij - 1 + kijl*(ichnk - 1)]), (double) 
      ((double) 0.0001));
    z0tub = z0rat*min((double) (z0tubmax), (double) (z0m[ij - 1 + kijl*(ichnk - 1)]));
    z0noz = max((double) (z0vis), (double) (z0tub));
    zorb = aorb / z0noz;
    xi = (log10(max((double) (zorb), (double) ((double) 3.0))) - abmin)*delabm1;
    ind = min((double) (iab - 1), (double) ((int) (xi)));
    deli1 = min((double) ((double) 1.0), (double) (xi - (double) (ind)));
    deli2 = (double) 1.0 - deli1;
    fww = swellft[ind - 1]*deli2 + swellft[ind + 1 - 1]*deli1;
    temp2 = fww*uorbt;
    if (swellf6 == (double) 1.0) {
      re_c = swellf4;
    } else {
      hftswellf6 = (double) 1.0 - swellf6;
      re_c = swellf4*(pow(((double) 2.0 / aorb), hftswellf6));
    }
    if (swellf7 > (double) 0.0) {
      smooth = (double) 0.5*tanh((re - re_c)*swellf7m1);
      pturb = (double) 0.5 + smooth;
      pvisc = (double) 0.5 - smooth;
    } else {
      if (re <= re_c) {
        pturb = (double) 0.0;
        pvisc = (double) 0.5;
      } else {
        pturb = (double) 0.5;
        pvisc = (double) 0.0;
      }
    }
    
    aird_pvisc = pvisc*raorw[ij - 1];
    
  }
  if (ngst == 1) {
    ustp[1 - 1] = ufric[ij - 1 + kijl*(ichnk - 1)];
  } else {
    ustp[1 - 1] = ufric[ij - 1 + kijl*(ichnk - 1)]*((double) 1.0 + sig_n);
    ustp[2 - 1] = ufric[ij - 1 + kijl*(ichnk - 1)]*((double) 1.0 - sig_n);
  }
  
  for (igst = 1; igst <= ngst; igst += 1) {
    ustpm1[igst - 1] = (double) 1.0 / max((double) (ustp[igst - 1]), (double) (epsus));
  }
  
  if (ltauwshelter) {
    for (igst = 1; igst <= ngst; igst += 1) {
      xstress[igst - 1] = (double) 0.0;
      ystress[igst - 1] = (double) 0.0;
      usg2[igst - 1] = pow(ustp[igst - 1], 2);
      taux[igst - 1] = usg2[igst - 1]*sin(wdwave[ij - 1 + kijl*(ichnk - 1)]);
      tauy[igst - 1] = usg2[igst - 1]*cos(wdwave[ij - 1 + kijl*(ichnk - 1)]);
    }
    
    rogoroair = g / raorw[ij - 1];
  }
  if (!llnormagam) {
    for (igst = 1; igst <= ngst; igst += 1) {
      gamnorma[igst - 1] = (double) 1.0;
    }
  }
  
  if (!llsneg) {
    for (k = 1; k <= nang; k += 1) {
      for (igst = 1; igst <= ngst; igst += 1) {
        dstab[igst - 1 + 2*(k - 1)] = (double) 0.0;
      }
    }
  }
  
  for (m = 1; m <= nfre; m += 1) {
    
    if (ltauwshelter) {
      for (igst = 1; igst <= ngst; igst += 1) {
        taupx = taux[igst - 1] - abs_tauwshelter*xstress[igst - 1];
        taupy = tauy[igst - 1] - abs_tauwshelter*ystress[igst - 1];
        usdirp[igst - 1] = atan2(taupx, taupy);
        ustp[igst - 1] = pow(((pow(taupx, 2)) + (pow(taupy, 2))), (double) 0.25);
        ustpm1[igst - 1] = 
          (double) 1.0 / max((double) (ustp[igst - 1]), (double) (epsus));
      }
      
      constf = 
        rogoroair*cinv[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))]*dfim[m - 1];
    }
    for (igst = 1; igst <= ngst; igst += 1) {
      ucn[igst - 1] = 
        ustp[igst - 1]*cinv[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))];
      ucnzalpd[igst - 1] = xkappa / (ucn[igst - 1] + zalp);
    }
    zcn = log(wavnum[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))]*z0m[ij - 1 + 
      kijl*(ichnk - 1)]);
    cnsn = zpifr[m - 1]*const1*raorw[ij - 1];
    for (k = 1; k <= nang; k += 1) {
      xllws[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))
        ] = (double) 0.0;
    }
    
    if (llsneg) {
      sig2 = pow(zpifr[m - 1], 2);
      dfim_sig2 = dfim[m - 1]*sig2;
      
      coef = -swellf*(double) 16.*sig2 / g;
      coef5 = -swellf5*(double) 2.*sqrt((double) ((double) 2.*nu_air*zpifr[m - 1]));
      
      dstab1 = 
        coef5*aird_pvisc*wavnum[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))];
      temp1 = coef*raorw[ij - 1];
    }
    
    for (k = 1; k <= nang; k += 1) {
      for (igst = 1; igst <= ngst; igst += 1) {
        
        sumf = (double) 0.0;
        sumfsin2 = (double) 0.0;
        
        if (ltauwshelter) {
          coslp = cos(th[k - 1] - usdirp[igst - 1]);
        } else {
          coslp = coswdif[ij - 1 + kijl*(k - 1)];
        }
        
        gam0[igst - 1 + 2*(k - 1)] = (double) 0.;
        if (coslp > (double) 0.01) {
          x = coslp*ucn[igst - 1];
          zlog = zcn + ucnzalpd[igst - 1] / coslp;
          if (zlog < (double) 0.0) {
            zlog2x = zlog*zlog*x;
            gam0[igst - 1 + 2*(k - 1)] = exp((double) (zlog))*zlog2x*zlog2x*cnsn;
            xllws[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk 
              - 1)))] = (double) 1.0;
          }
        }
        
        if (llsneg) {
          dstab2 = temp1*(temp2 + (fu + fud*coslp)*ustp[igst - 1]);
          dstab[igst - 1 + 2*(k - 1)] = dstab1 + pturb*dstab2;
        }
        
        sumf = sumf + gam0[igst - 1 + 2*(k - 1)]*fl1[ij - 1 + kijl*(k - 1 + 
          nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))];
        sumfsin2 = sumfsin2 + gam0[igst - 1 + 2*(k - 1)]*fl1[ij - 1 + kijl*(k - 1 + 
          nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))]*sinwdif2[ij - 1 + 
          kijl*(k - 1)];
      }
    }
    
    if (llnormagam) {
      
      xngamconst = cstrnfac*xk2cg[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))];
      for (igst = 1; igst <= ngst; igst += 1) {
        znz = xngamconst*ustpm1[igst - 1];
        gamnorma[igst - 1] = ((double) 1.0 + znz*sumfsin2) / ((double) 1.0 + znz*sumf);
      }
      
    }
    for (k = 1; k <= nang; k += 1) {
      
      for (igst = 1; igst <= ngst; igst += 1) {
        // SLP: only the positive contributions
        slp[igst - 1] = gam0[igst - 1 + 2*(k - 1)]*gamnorma[igst - 1];
        flp[igst - 1] = slp[igst - 1] + dstab[igst - 1 + 2*(k - 1)];
      }
      
      for (igst = 1; igst <= ngst; igst += 1) {
        slp[igst - 1] = slp[igst - 1]*fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 +
           nfre_loki_param*(ichnk - 1)))];
      }
      
      if (ltauwshelter) {
        const11 = constf*sinth[k - 1];
        const22 = constf*costh[k - 1];
        for (igst = 1; igst <= ngst; igst += 1) {
          xstress[igst - 1] = xstress[igst - 1] + slp[igst - 1]*const11;
          ystress[igst - 1] = ystress[igst - 1] + slp[igst - 1]*const22;
        }
      }
      
      igst = 1;
      slp_avg = slp[igst - 1];
      flp_avg = flp[igst - 1];
      for (igst = 2; igst <= ngst; igst += 1) {
        slp_avg = slp_avg + slp[igst - 1];
        flp_avg = flp_avg + flp[igst - 1];
      }
      
      spos[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] = avg_gst*slp_avg;
      fld[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] = avg_gst*flp_avg;
      sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1))] = fld[ij - 1 + kijl*(k - 1 + 
        nang_loki_param*(m - 1))]*fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
        nfre_loki_param*(ichnk - 1)))];
      
    }
    
  }

  
}
