#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sinflx_c.h"
#include "stresso_c.h"
#include "frcutindex_c.h"
#include "femeanws_c.h"
#include "sinput_c.h"
#include "airsea_c.h"
#include "halphap_c.h"

__device__ void sinflx_c(int icall, int kijs, int kijl, int lupdtus, double * fl1, 
  const double * wavnum, const double * cinv, const double * xk2cg, double * wswave, 
  const double * wdwave, const double * aird, const double * raorw, 
  const double * wstar, const double * cicover, const double * coswdif, 
  const double * sinwdif2, const double * fmean, double * halp, double * fmeanws, 
  const double * flm, double * ufric, double * tauw, double * tauwdir, double * z0m, 
  double * z0b, double * chrnck, double * phiwa, double * fld, double * sl, 
  double * spos, int * mij, double * rhowgdfth, double * xllws, double abmax, 
  double abmin, double acd, double acdlin, double alpha, double alphamax, 
  double alphamin, double alphapmax, double ang_gc_a, double ang_gc_b, double ang_gc_c, 
  double bcd, double bcdlin, double betamaxoxkappa2, double bmaxokap, 
  const double * c2osqrtvg_gc, double cdmax, double chnkmin_u, double cithrsh_tail, 
  const double * cm_gc, const double * costh, const double * delkcc_gc_ns, 
  const double * delkcc_omxkm3_gc, double delth, const double * dfim, 
  const double * dfimofr, double dthrn_a, double dthrn_u, double eps1, double epsmin, 
  double epsus, double flogsprdm1, const double * fr, const double * fr5, double fric, 
  double frtail, double g, double gamnconst, double gm1, int iab, int icode, 
  int icode_cpl, int idamping, int iphys, int jtot_tauhf, int llcapchnk, int llgcbz0, 
  int llnormagam, int lwcou, int nang, int nfre, int nwav_gc, const double * om3gmkm_gc, 
  const double * omega_gc, const double * omxkm3_gc, const double * rhowg_dfim, 
  double rn1_rn, double rnu, double rnum, const double * sinth, double sqrtgosurft, 
  double swellf, double swellf2, double swellf3, double swellf4, double swellf5, 
  double swellf6, double swellf7, double swellf7m1, const double * swellft, 
  double tailfactor, double tailfactor_pm, double tauwshelter, const double * th, 
  double wetail, double wspmin, const double * wtauhf, double x0tauhf, double xkappa, 
  const double * xkmsqrtvgoc2_gc, const double * xkm_gc, const double * xk_gc, 
  double xlogkratiom1_gc, double xnlev, double z0rat, double z0tubmax, double zalp, 
  double zpi, double zpi4gm1, double zpi4gm2, const double * zpifr, int ichnk, 
  int nchnk, int ij, double * rnfac, double * tmp_em, double * stresso_xstress, 
  double * stresso_ystress, double * stresso_tauhf, double * stresso_phihf, 
  double * stresso_usdirp, double * stresso_ust) {
  
  






  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  
  
  
  
  
  
  int k;
  int iusfg;
  int icode_wnd;
  int ngst;
  
  
  int llphiwa;
  int llsneg;
  if (icall == 1) {
    iusfg = 0;
    if (lwcou) {
      icode_wnd = icode_cpl;
    } else {
      icode_wnd = icode;
    }
    
    llphiwa = false;
    llsneg = false;
  } else {
    iusfg = 1;
    icode_wnd = 3;
    
    llphiwa = true;
    llsneg = true;
  }
  

  if (llnormagam && llcapchnk) {
    rnfac[ij - 1 + kijl*(ichnk - 1)] = (double) 1.0 + dthrn_a*((double) 1.0 + 
      tanh(wswave[ij - 1 + kijl*(ichnk - 1)] - dthrn_u));
  } else {
    rnfac[ij - 1 + kijl*(ichnk - 1)] = (double) 1.0;
  }

  if (lupdtus) {
    // increase noise level in the tail
    if (icall == 1) {

      for (k = 1; k <= nang; k += 1) {
        fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(nfre - 1 + nfre_loki_param*(ichnk - 1
          )))] = max((double) (fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(nfre - 1 + 
          nfre_loki_param*(ichnk - 1)))]), (double) (flm[ij - 1 + kijl*(k - 1)]));
      }

      
      if (llgcbz0) {
        halphap_c(kijs, kijl, wavnum, coswdif, fl1, halp, alphapmax, delth, dfim, 
          dfimofr, epsmin, fr, fr5, frtail, nang, nfre, wetail, zpi4gm2, ichnk, nchnk, ij
          );
      } else {

        halp[ij - 1] = (double) 0.0;

      }
      
    }
    
    airsea_c(kijs, kijl, halp, wswave, wdwave, tauw, tauwdir, 
       (&rnfac[ + kijl*(ichnk - 1)]), ufric, z0m, z0b, chrnck, icode_wnd, iusfg, acd, 
      alpha, alphamax, alphamin, ang_gc_a, ang_gc_b, ang_gc_c, bcd, betamaxoxkappa2, 
      bmaxokap, c2osqrtvg_gc, cdmax, chnkmin_u, cm_gc, delkcc_gc_ns, delkcc_omxkm3_gc, 
      eps1, epsmin, epsus, g, gm1, llcapchnk, llgcbz0, llnormagam, nwav_gc, om3gmkm_gc, 
      omxkm3_gc, rn1_rn, rnu, rnum, sqrtgosurft, wspmin, xkappa, xkmsqrtvgoc2_gc, 
      xkm_gc, xk_gc, xlogkratiom1_gc, xnlev, zalp, ichnk, nchnk, ij);
    
  }
  sinput_c(icall, llsneg, kijs, kijl, fl1, wavnum, cinv, xk2cg, wdwave, wswave, ufric, 
    z0m, coswdif, sinwdif2, raorw, wstar,  (&rnfac[ + kijl*(ichnk - 1)]), fld, sl, spos, 
    xllws, abmax, abmin, acdlin, alphamax, alphamin, bcdlin, betamaxoxkappa2, costh, 
    delth, dfim, epsmin, epsus, g, iab, idamping, iphys, llgcbz0, llnormagam, nang, 
    nfre, rnu, rnum, sinth, swellf, swellf2, swellf3, swellf4, swellf5, swellf6, 
    swellf7, swellf7m1, swellft, tauwshelter, th, wspmin, xkappa, z0rat, z0tubmax, zalp, 
    zpi, zpifr, ichnk, nchnk, ij);
  femeanws_c(kijs, kijl, fl1, xllws, fmeanws,  (&tmp_em[ + kijl*(ichnk - 1)]), delth, 
    dfim, dfimofr, epsmin, fr, frtail, nang, nfre, wetail, ichnk, nchnk, ij);
  frcutindex_c(kijs, kijl, fmean, fmeanws, ufric, cicover, mij, rhowgdfth, cithrsh_tail, 
    epsmin, flogsprdm1, fr, fric, g, nfre, rhowg_dfim, tailfactor, tailfactor_pm, zpifr, 
    ichnk, nchnk, ij);
  stresso_c(kijs, kijl, mij, rhowgdfth, fl1, sl, spos, cinv, wdwave, ufric, z0m, aird, 
     (&rnfac[ + kijl*(ichnk - 1)]), coswdif, sinwdif2, tauw, tauwdir, phiwa, llphiwa, 
    costh, delth, eps1, fr5, g, gamnconst, gm1, iphys, jtot_tauhf, llgcbz0, llnormagam, 
    nang, nfre, nwav_gc, omega_gc, rhowg_dfim, sinth, sqrtgosurft, tauwshelter, wtauhf, 
    x0tauhf, xkappa, xkm_gc, xk_gc, xlogkratiom1_gc, zalp, zpi4gm1, zpi4gm2, zpifr, 
    ichnk, nchnk, ij, stresso_xstress, stresso_ystress, stresso_tauhf, stresso_phihf, 
    stresso_usdirp, stresso_ust);
  
}
