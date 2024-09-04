#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sinflx_c.h"
#include "stresso_c.h"
#include "sinput_c.h"
#include "halphap_c.h"
#include "frcutindex_c.h"
#include "femeanws_c.h"
#include "airsea_c.h"


__device__ void sinflx_c(int icall, int ncall, int kijs, int kijl, int lupdtus, 
  double * __restrict__ fl1, const double * __restrict__ wavnum, 
  const double * __restrict__ cinv, const double * __restrict__ xk2cg, 
  double * __restrict__ wswave, const double * __restrict__ wdwave, 
  const double * __restrict__ aird, const double * __restrict__ raorw, 
  const double * __restrict__ wstar, const double * __restrict__ cicover, 
  const double * __restrict__ coswdif, const double * __restrict__ sinwdif2, 
  const double * __restrict__ fmean, double * __restrict__ halp, 
  double * __restrict__ fmeanws, const double * __restrict__ flm, 
  double * __restrict__ ufric, double * __restrict__ tauw, 
  double * __restrict__ tauwdir, double * __restrict__ z0m, double * __restrict__ z0b, 
  double * __restrict__ chrnck, double * __restrict__ phiwa, double * __restrict__ fld, 
  double * __restrict__ sl, double * __restrict__ spos, int * __restrict__ mij, 
  double * __restrict__ rhowgdfth, double * __restrict__ xllws, double abmax, 
  double abmin, double acd, double acdlin, double alpha, double alphamax, 
  double alphamin, double alphapmax, double ang_gc_a, double ang_gc_b, double ang_gc_c, 
  double bcd, double bcdlin, double betamaxoxkappa2, double bmaxokap, 
  const double * __restrict__ c2osqrtvg_gc, double cdmax, double chnkmin_u, 
  double cithrsh_tail, const double * __restrict__ cm_gc, 
  const double * __restrict__ costh, const double * __restrict__ delkcc_gc_ns, 
  const double * __restrict__ delkcc_omxkm3_gc, double delth, 
  const double * __restrict__ dfim, const double * __restrict__ dfimofr, double dthrn_a, 
  double dthrn_u, double eps1, double epsmin, double epsus, double flogsprdm1, 
  const double * __restrict__ fr, const double * __restrict__ fr5, double fric, 
  double frtail, double g, double gamnconst, double gm1, int iab, int icode, 
  int icode_cpl, int idamping, int iphys, int jtot_tauhf, int llcapchnk, int llgcbz0, 
  int llnormagam, int lwcou, int nang, int nfre, int nwav_gc, 
  const double * __restrict__ om3gmkm_gc, const double * __restrict__ omega_gc, 
  const double * __restrict__ omxkm3_gc, const double * __restrict__ rhowg_dfim, 
  double rn1_rn, double rnu, double rnum, const double * __restrict__ sinth, 
  double sqrtgosurft, double swellf, double swellf2, double swellf3, double swellf4, 
  double swellf5, double swellf6, double swellf7, double swellf7m1, 
  const double * __restrict__ swellft, double tailfactor, double tailfactor_pm, 
  double tauwshelter, const double * __restrict__ th, double wetail, double wspmin, 
  const double * __restrict__ wtauhf, double x0tauhf, double xkappa, 
  const double * __restrict__ xkmsqrtvgoc2_gc, const double * __restrict__ xkm_gc, 
  const double * __restrict__ xk_gc, double xlogkratiom1_gc, double xnlev, double z0rat, 
  double z0tubmax, double zalp, double zpi, double zpi4gm1, double zpi4gm2, 
  const double * __restrict__ zpifr, int ichnk, int nchnk, int ij, 
  double * __restrict__ rnfac, double * __restrict__ tmp_em, 
  double * __restrict__ taut_z0_alphaog, double * __restrict__ taut_z0_xmin, 
  double * __restrict__ taut_z0_w1, double * __restrict__ taut_z0_tauwact, 
  double * __restrict__ taut_z0_tauweff, double * __restrict__ taut_z0_ang_gc, 
  double * __restrict__ taut_z0_tauunr, int * __restrict__ taut_z0_llcosdiff, 
  double * __restrict__ stress_gc_gam_w, double * __restrict__ z0wave_alphaog, 
  double * __restrict__ femeanws_temp2, double * __restrict__ femeanws_em_loc, 
  double * __restrict__ halphap_alphap, double * __restrict__ halphap_xmss, 
  double * __restrict__ halphap_em, double * __restrict__ halphap_fm, 
  double * __restrict__ halphap_f1d, double * __restrict__ halphap_wd, 
  double * __restrict__ halphap_flwd, double * __restrict__ femean_temp2, 
  double * __restrict__ meansqs_lf_fd, double * __restrict__ meansqs_lf_temp1, 
  double * __restrict__ meansqs_lf_temp2, double * __restrict__ sinput_ard_constf, 
  double * __restrict__ sinput_ard_const11, double * __restrict__ sinput_ard_const22, 
  double * __restrict__ sinput_ard_z0vis, double * __restrict__ sinput_ard_z0noz, 
  double * __restrict__ sinput_ard_fww, double * __restrict__ sinput_ard_pvisc, 
  double * __restrict__ sinput_ard_pturb, double * __restrict__ sinput_ard_zcn, 
  double * __restrict__ sinput_ard_sig_n, double * __restrict__ sinput_ard_uorbt, 
  double * __restrict__ sinput_ard_aorb, double * __restrict__ sinput_ard_temp, 
  double * __restrict__ sinput_ard_re, double * __restrict__ sinput_ard_re_c, 
  double * __restrict__ sinput_ard_zorb, double * __restrict__ sinput_ard_cnsn, 
  double * __restrict__ sinput_ard_sumf, double * __restrict__ sinput_ard_sumfsin2, 
  double * __restrict__ sinput_ard_cstrnfac, double * __restrict__ sinput_ard_flp_avg, 
  double * __restrict__ sinput_ard_slp_avg, double * __restrict__ sinput_ard_rogoroair, 
  double * __restrict__ sinput_ard_aird_pvisc, double * __restrict__ sinput_ard_xstress, 
  double * __restrict__ sinput_ard_ystress, double * __restrict__ sinput_ard_flp, 
  double * __restrict__ sinput_ard_slp, double * __restrict__ sinput_ard_usg2, 
  double * __restrict__ sinput_ard_taux, double * __restrict__ sinput_ard_tauy, 
  double * __restrict__ sinput_ard_ustp, double * __restrict__ sinput_ard_ustpm1, 
  double * __restrict__ sinput_ard_usdirp, double * __restrict__ sinput_ard_ucn, 
  double * __restrict__ sinput_ard_ucnzalpd, 
  double * __restrict__ sinput_ard_xngamconst, 
  double * __restrict__ sinput_ard_gamnorma, double * __restrict__ sinput_ard_dstab1, 
  double * __restrict__ sinput_ard_temp1, double * __restrict__ sinput_ard_temp2, 
  double * __restrict__ sinput_ard_gam0, double * __restrict__ sinput_ard_dstab, 
  double * __restrict__ sinput_ard_coslp, double * __restrict__ sinput_jan_ztanhkd, 
  double * __restrict__ sinput_jan_sig_n, double * __restrict__ sinput_jan_cnsn, 
  double * __restrict__ sinput_jan_sumf, double * __restrict__ sinput_jan_sumfsin2, 
  double * __restrict__ sinput_jan_cstrnfac, 
  double * __restrict__ sinput_jan_xngamconst, 
  double * __restrict__ sinput_jan_gamnorma, double * __restrict__ sinput_jan_sigdev, 
  double * __restrict__ sinput_jan_us, double * __restrict__ sinput_jan_z0, 
  double * __restrict__ sinput_jan_ucn, double * __restrict__ sinput_jan_zcn, 
  double * __restrict__ sinput_jan_ustpm1, double * __restrict__ sinput_jan_xvd, 
  double * __restrict__ sinput_jan_ucnd, double * __restrict__ sinput_jan_const3_ucn2, 
  double * __restrict__ sinput_jan_ufac1, double * __restrict__ sinput_jan_ufac2, 
  double * __restrict__ sinput_jan_tempd, double * __restrict__ sinput_jan_gam0, 
  int * __restrict__ sinput_jan_lz, double * __restrict__ stresso_xstress, 
  double * __restrict__ stresso_ystress, double * __restrict__ stresso_tauhf, 
  double * __restrict__ stresso_phihf, double * __restrict__ stresso_cmrhowgdfth, 
  double * __restrict__ stresso_us2, double * __restrict__ stresso_taux, 
  double * __restrict__ stresso_tauy, double * __restrict__ stresso_taupx, 
  double * __restrict__ stresso_taupy, double * __restrict__ stresso_usdirp, 
  double * __restrict__ stresso_ust, double * __restrict__ stresso_sumt, 
  double * __restrict__ stresso_sumx, double * __restrict__ stresso_sumy, 
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
  
  // needed for Loki
  
  
  // ----------------------------------------------------------------------
  






  
  
  
  
  
  
  
  
  int k;
  int iusfg;
  int icode_wnd;
  int ngst;
  
  double zhook_handle;
  
  int llphiwa;
  int llsneg;
  
  // ----------------------------------------------------------------------
  
  
  // UPDATE UFRIC AND Z0M
  if (icall == 1) {
    iusfg = 0;
    if (lwcou) {
      icode_wnd = icode_cpl;
    } else {
      icode_wnd = icode;
    }
  } else {
    iusfg = 1;
    icode_wnd = 3;
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
        fl1[ij - 1 + kijl*(k - 1 + nang*(nfre - 1 + nfre*(ichnk - 1)))] = fmax((double) 
          (fl1[ij - 1 + kijl*(k - 1 + nang*(nfre - 1 + nfre*(ichnk - 1)))]), (double) 
          (flm[ij - 1 + kijl*(k - 1)]));
      }
      
      
      if (llgcbz0) {
        halphap_c(kijs, kijl, wavnum, coswdif, fl1, halp, alphapmax, delth, dfim, 
          dfimofr, epsmin, fr, fr5, frtail, nang, nfre, wetail, zpi4gm2, ichnk, nchnk, 
          ij, halphap_alphap, halphap_xmss, halphap_em, halphap_fm, halphap_f1d, 
          halphap_wd, halphap_flwd, femean_temp2, meansqs_lf_fd, meansqs_lf_temp1, 
          meansqs_lf_temp2);
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
      xkm_gc, xk_gc, xlogkratiom1_gc, xnlev, zalp, ichnk, nchnk, ij, taut_z0_alphaog, 
      taut_z0_xmin, taut_z0_w1, taut_z0_tauwact, taut_z0_tauweff, taut_z0_ang_gc, 
      taut_z0_tauunr, taut_z0_llcosdiff, stress_gc_gam_w, z0wave_alphaog);
    
  }
  
  // COMPUTE WIND INPUT
  //!  FLD AND SL ARE INITIALISED IN SINPUT
  if (icall < ncall) {
    ngst = 1;
    llphiwa = false;
    llsneg = false;
  } else {
    ngst = 2;
    llphiwa = true;
    llsneg = true;
  }
  
  sinput_c(ngst, llsneg, kijs, kijl, fl1, wavnum, cinv, xk2cg, wdwave, wswave, ufric, 
    z0m, coswdif, sinwdif2, raorw, wstar,  (&rnfac[ + kijl*(ichnk - 1)]), fld, sl, spos, 
    xllws, abmax, abmin, acdlin, alphamax, alphamin, bcdlin, betamaxoxkappa2, costh, 
    delth, dfim, epsmin, epsus, g, iab, idamping, iphys, llgcbz0, llnormagam, nang, 
    nfre, rnu, rnum, sinth, swellf, swellf2, swellf3, swellf4, swellf5, swellf6, 
    swellf7, swellf7m1, swellft, tauwshelter, th, wspmin, xkappa, z0rat, z0tubmax, zalp, 
    zpi, zpifr, ichnk, nchnk, ij, sinput_ard_constf, sinput_ard_const11, 
    sinput_ard_const22, sinput_ard_z0vis, sinput_ard_z0noz, sinput_ard_fww, 
    sinput_ard_pvisc, sinput_ard_pturb, sinput_ard_zcn, sinput_ard_sig_n, 
    sinput_ard_uorbt, sinput_ard_aorb, sinput_ard_temp, sinput_ard_re, sinput_ard_re_c, 
    sinput_ard_zorb, sinput_ard_cnsn, sinput_ard_sumf, sinput_ard_sumfsin2, 
    sinput_ard_cstrnfac, sinput_ard_flp_avg, sinput_ard_slp_avg, sinput_ard_rogoroair, 
    sinput_ard_aird_pvisc, sinput_ard_xstress, sinput_ard_ystress, sinput_ard_flp, 
    sinput_ard_slp, sinput_ard_usg2, sinput_ard_taux, sinput_ard_tauy, sinput_ard_ustp, 
    sinput_ard_ustpm1, sinput_ard_usdirp, sinput_ard_ucn, sinput_ard_ucnzalpd, 
    sinput_ard_xngamconst, sinput_ard_gamnorma, sinput_ard_dstab1, sinput_ard_temp1, 
    sinput_ard_temp2, sinput_ard_gam0, sinput_ard_dstab, sinput_ard_coslp, 
    sinput_jan_ztanhkd, sinput_jan_sig_n, sinput_jan_cnsn, sinput_jan_sumf, 
    sinput_jan_sumfsin2, sinput_jan_cstrnfac, sinput_jan_xngamconst, 
    sinput_jan_gamnorma, sinput_jan_sigdev, sinput_jan_us, sinput_jan_z0, 
    sinput_jan_ucn, sinput_jan_zcn, sinput_jan_ustpm1, sinput_jan_xvd, sinput_jan_ucnd, 
    sinput_jan_const3_ucn2, sinput_jan_ufac1, sinput_jan_ufac2, sinput_jan_tempd, 
    sinput_jan_gam0, sinput_jan_lz);
  
  
  // MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
  femeanws_c(kijs, kijl, fl1, xllws, fmeanws,  (&tmp_em[ + kijl*(ichnk - 1)]), delth, 
    dfim, dfimofr, epsmin, fr, frtail, nang, nfre, wetail, ichnk, nchnk, ij, 
    femeanws_temp2, femeanws_em_loc);
  
  // COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
  frcutindex_c(kijs, kijl, fmean, fmeanws, ufric, cicover, mij, rhowgdfth, cithrsh_tail, 
    epsmin, flogsprdm1, fr, fric, g, nfre, rhowg_dfim, tailfactor, tailfactor_pm, zpifr, 
    ichnk, nchnk, ij);
  
  // UPDATE TAUW
  stresso_c(kijs, kijl, mij, rhowgdfth, fl1, sl, spos, cinv, wdwave, ufric, z0m, aird, 
     (&rnfac[ + kijl*(ichnk - 1)]), coswdif, sinwdif2, tauw, tauwdir, phiwa, llphiwa, 
    costh, delth, eps1, fr5, g, gamnconst, gm1, iphys, jtot_tauhf, llgcbz0, llnormagam, 
    nang, nfre, nwav_gc, omega_gc, rhowg_dfim, sinth, sqrtgosurft, tauwshelter, wtauhf, 
    x0tauhf, xkappa, xkm_gc, xk_gc, xlogkratiom1_gc, zalp, zpi4gm1, zpi4gm2, zpifr, 
    ichnk, nchnk, ij, stresso_xstress, stresso_ystress, stresso_tauhf, stresso_phihf, 
    stresso_cmrhowgdfth, stresso_us2, stresso_taux, stresso_tauy, stresso_taupx, 
    stresso_taupy, stresso_usdirp, stresso_ust, stresso_sumt, stresso_sumx, 
    stresso_sumy, tau_phi_hf_ns, tau_phi_hf_xks, tau_phi_hf_oms, tau_phi_hf_sqrtz0og, 
    tau_phi_hf_zsup, tau_phi_hf_zinf, tau_phi_hf_delz, tau_phi_hf_taul, 
    tau_phi_hf_xloggz0, tau_phi_hf_sqrtgz0, tau_phi_hf_ustph, tau_phi_hf_const1, 
    tau_phi_hf_const2, tau_phi_hf_consttau, tau_phi_hf_constphi, tau_phi_hf_f1dcos2, 
    tau_phi_hf_f1dcos3, tau_phi_hf_f1d, tau_phi_hf_f1dsin2);
  
  // ----------------------------------------------------------------------
  
  
}
