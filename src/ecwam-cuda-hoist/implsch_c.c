#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "wnfluxes_c.h"
#include "stokestrn_c.h"
#include "snonlin_c.h"
#include "sinflx_c.h"
#include "setice_c.h"
#include "sdiwbk_c.h"
#include "sdissip_c.h"
#include "sdepthlim_c.h"
#include "imphftail_c.h"
#include "sbottom_c.h"
#include "fkmean_c.h"
#include "femeanws_c.h"
#include "ciwabr_c.h"

#include "implsch_c.h"
#include "implsch_c_launch.h"


__global__ void implsch_c(int kijs, int kijl, double * __restrict__ fl1, 
  const double * __restrict__ wavnum, const double * __restrict__ cgroup, 
  const double * __restrict__ ciwa, const double * __restrict__ cinv, 
  const double * __restrict__ xk2cg, const double * __restrict__ stokfac, 
  const double * __restrict__ emaxdpt, const double * __restrict__ depth, 
  const int * __restrict__ iobnd, const int * __restrict__ iodp, 
  double * __restrict__ aird, double * __restrict__ wdwave, 
  double * __restrict__ cicover, double * __restrict__ wswave, 
  double * __restrict__ wstar, double * __restrict__ ustra, double * __restrict__ vstra, 
  double * __restrict__ ufric, double * __restrict__ tauw, 
  double * __restrict__ tauwdir, double * __restrict__ z0m, double * __restrict__ z0b, 
  double * __restrict__ chrnck, double * __restrict__ cithick, 
  double * __restrict__ nemoustokes, double * __restrict__ nemovstokes, 
  double * __restrict__ nemostrn, double * __restrict__ nphieps, 
  double * __restrict__ ntauoc, double * __restrict__ nswh, double * __restrict__ nmwp, 
  double * __restrict__ nemotaux, double * __restrict__ nemotauy, 
  double * __restrict__ nemowswave, double * __restrict__ nemophif, 
  double * __restrict__ wsemean, double * __restrict__ wsfmean, 
  double * __restrict__ ustokes, double * __restrict__ vstokes, 
  double * __restrict__ strnms, double * __restrict__ tauxd, 
  double * __restrict__ tauyd, double * __restrict__ tauocxd, 
  double * __restrict__ tauocyd, double * __restrict__ tauoc, 
  double * __restrict__ phiocd, double * __restrict__ phieps, 
  double * __restrict__ phiaw, int * __restrict__ mij, double * __restrict__ xllws, 
  double abmax, double abmin, double acd, double acdlin, 
  const double * __restrict__ af11, double afcrv, double alpha, double alphamax, 
  double alphamin, double alphapmax, double ang_gc_a, double ang_gc_b, double ang_gc_c, 
  double bathymax, double bcd, double bcdlin, double betamaxoxkappa2, double bfcrv, 
  double bmaxokap, double brkpbcoef, const double * __restrict__ c2osqrtvg_gc, 
  double cdicwa, double cdis, double cdisvis, double cdmax, double chnkmin_u, 
  double ciblock, double cithrsh, double cithrsh_tail, 
  const double * __restrict__ cm_gc, const double * __restrict__ cofrm4, 
  const double * __restrict__ costh, double dal1, double dal2, 
  const double * __restrict__ delkcc_gc_ns, 
  const double * __restrict__ delkcc_omxkm3_gc, double delta_sdis, double delth, 
  const double * __restrict__ dfim, const double * __restrict__ dfimfr, 
  const double * __restrict__ dfimfr2, const double * __restrict__ dfimofr, 
  const double * __restrict__ dfim_sim, double dkmax, double dthrn_a, double dthrn_u, 
  double egrcrv, double eps1, double epsmin, double epsu10, double epsus, 
  const double * __restrict__ fklam, const double * __restrict__ fklam1, 
  const double * __restrict__ fklap, const double * __restrict__ fklap1, 
  const double * __restrict__ flmax, double flmin, double flogsprdm1, 
  const double * __restrict__ fr, const double * __restrict__ fr5, double fratio, 
  double fric, double frtail, double g, double gamnconst, double gm1, int iab, 
  int icode, int icode_cpl, int idamping, int idelt, const int * __restrict__ ikm, 
  const int * __restrict__ ikm1, const int * __restrict__ ikp, 
  const int * __restrict__ ikp1, const int * __restrict__ indicessat, 
  const int * __restrict__ inlcoef, int iphys, int ipsat, int isnonlin, int jtot_tauhf, 
  const int * __restrict__ k11w, const int * __restrict__ k1w, 
  const int * __restrict__ k21w, const int * __restrict__ k2w, int kfrh, int lbiwbk, 
  int lciwabr, int licerun, int llcapchnk, int llgcbz0, int llnormagam, int llunstr, 
  int lmaskice, int lwamrsetci, int lwcou, int lwcouast, int lwflux, int lwfluxout, 
  int lwnemocou, int lwnemocousend, int lwnemocoustk, int lwnemocoustrn, 
  int lwnemotauoc, int lwvflx_snl, int mfrstlw, double miche, int mlsthg, int nang, 
  int nfre, int nfre_odd, int nfre_red, int nsdsnth, int nwav_gc, 
  const double * __restrict__ om3gmkm_gc, const double * __restrict__ omega_gc, 
  const double * __restrict__ omxkm3_gc, double phiepsmax, double phiepsmin, 
  const double * __restrict__ rhowg_dfim, double rn1_rn, 
  const double * __restrict__ rnlcoef, double rnu, double rnum, double rowater, 
  double rowaterm1, const double * __restrict__ satweights, double sdsbr, 
  const double * __restrict__ sinth, double sqrtgosurft, double ssdsbrf1, double ssdsc2, 
  double ssdsc3, double ssdsc4, double ssdsc5, double ssdsc6, double swellf, 
  double swellf2, double swellf3, double swellf4, double swellf5, double swellf6, 
  double swellf7, double swellf7m1, const double * __restrict__ swellft, 
  double tailfactor, double tailfactor_pm, double tauocmax, double tauocmin, 
  double tauwshelter, const double * __restrict__ th, double wetail, double wp1tail, 
  double wp2tail, double wsemean_min, double wspmin, const double * __restrict__ wtauhf, 
  double x0tauhf, double xkappa, double xkdmin, 
  const double * __restrict__ xkmsqrtvgoc2_gc, const double * __restrict__ xkm_gc, 
  const double * __restrict__ xk_gc, double xlogkratiom1_gc, double xnlev, double z0rat, 
  double z0tubmax, double zalp, double zpi, double zpi4gm1, double zpi4gm2, 
  const double * __restrict__ zpifr, int nchnk, int nproma_wam, 
  double * __restrict__ delfl, double * __restrict__ raorw, double * __restrict__ emean, 
  double * __restrict__ fmean, double * __restrict__ halp, 
  double * __restrict__ emeanws, double * __restrict__ fmeanws, 
  double * __restrict__ usfm, double * __restrict__ f1mean, 
  double * __restrict__ akmean, double * __restrict__ xkmean, 
  double * __restrict__ phiwa, double * __restrict__ flm, double * __restrict__ coswdif, 
  double * __restrict__ sinwdif2, double * __restrict__ temp, 
  double * __restrict__ rhowgdfth, double * __restrict__ fld, double * __restrict__ sl, 
  double * __restrict__ spos, double * __restrict__ cireduc, 
  double * __restrict__ ssource, double * __restrict__ femeanws_temp2, 
  double * __restrict__ femeanws_em_loc, double * __restrict__ fkmean_tempa, 
  double * __restrict__ fkmean_tempx, double * __restrict__ fkmean_temp2, 
  double * __restrict__ sbottom_sbo, double * __restrict__ imphftail_temp1, 
  double * __restrict__ imphftail_temp2, double * __restrict__ sdepthlim_em, 
  double * __restrict__ semean_temp, double * __restrict__ sdissip_ard_facturb, 
  double * __restrict__ sdissip_ard_facsat, double * __restrict__ sdissip_ard_facwtrb, 
  double * __restrict__ sdissip_ard_temp1, double * __restrict__ sdissip_ard_bth0, 
  double * __restrict__ sdissip_ard_c_, double * __restrict__ sdissip_ard_c_c, 
  double * __restrict__ sdissip_ard_dsip, double * __restrict__ sdissip_ard_trpz_dsip, 
  double * __restrict__ sdissip_ard_bth, double * __restrict__ sdissip_ard_temp2, 
  double * __restrict__ sdissip_ard_d, double * __restrict__ sdissip_ard_scumul, 
  double * __restrict__ sdissip_ard_renewalfreq, 
  double * __restrict__ sdissip_ard_wcumul, double * __restrict__ sdissip_jan_temp1, 
  double * __restrict__ sdissip_jan_sds, double * __restrict__ sdissip_jan_x, 
  double * __restrict__ sdissip_jan_xk2, double * __restrict__ sdiwbk_sds, 
  double * __restrict__ setice_cireduc, double * __restrict__ setice_temp, 
  double * __restrict__ setice_icefree, double * __restrict__ sinflx_rnfac, 
  double * __restrict__ sinflx_tmp_em, double * __restrict__ taut_z0_alphaog, 
  double * __restrict__ taut_z0_xmin, double * __restrict__ taut_z0_w1, 
  double * __restrict__ taut_z0_tauwact, double * __restrict__ taut_z0_tauweff, 
  double * __restrict__ taut_z0_ang_gc, double * __restrict__ taut_z0_tauunr, 
  int * __restrict__ taut_z0_llcosdiff, double * __restrict__ stress_gc_gam_w, 
  double * __restrict__ z0wave_alphaog, double * __restrict__ halphap_alphap, 
  double * __restrict__ halphap_xmss, double * __restrict__ halphap_em, 
  double * __restrict__ halphap_fm, double * __restrict__ halphap_f1d, 
  double * __restrict__ halphap_wd, double * __restrict__ halphap_flwd, 
  double * __restrict__ femean_temp2, double * __restrict__ meansqs_lf_fd, 
  double * __restrict__ meansqs_lf_temp1, double * __restrict__ meansqs_lf_temp2, 
  double * __restrict__ sinput_ard_constf, double * __restrict__ sinput_ard_const11, 
  double * __restrict__ sinput_ard_const22, double * __restrict__ sinput_ard_z0vis, 
  double * __restrict__ sinput_ard_z0noz, double * __restrict__ sinput_ard_fww, 
  double * __restrict__ sinput_ard_pvisc, double * __restrict__ sinput_ard_pturb, 
  double * __restrict__ sinput_ard_zcn, double * __restrict__ sinput_ard_sig_n, 
  double * __restrict__ sinput_ard_uorbt, double * __restrict__ sinput_ard_aorb, 
  double * __restrict__ sinput_ard_temp, double * __restrict__ sinput_ard_re, 
  double * __restrict__ sinput_ard_re_c, double * __restrict__ sinput_ard_zorb, 
  double * __restrict__ sinput_ard_cnsn, double * __restrict__ sinput_ard_sumf, 
  double * __restrict__ sinput_ard_sumfsin2, double * __restrict__ sinput_ard_cstrnfac, 
  double * __restrict__ sinput_ard_flp_avg, double * __restrict__ sinput_ard_slp_avg, 
  double * __restrict__ sinput_ard_rogoroair, 
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
  double * __restrict__ tau_phi_hf_f1dsin2, double * __restrict__ snonlin_enh, 
  double * __restrict__ snonlin_xnu, double * __restrict__ snonlin_sig_th, 
  double * __restrict__ snonlin_ftemp, double * __restrict__ snonlin_ad, 
  double * __restrict__ snonlin_delad, double * __restrict__ snonlin_delap, 
  double * __restrict__ snonlin_delam, double * __restrict__ snonlin_enhfr, 
  int * __restrict__ peak_ang_mmax, double * __restrict__ peak_ang_sum0, 
  double * __restrict__ peak_ang_sum1, double * __restrict__ peak_ang_sum2, 
  double * __restrict__ peak_ang_xmax, double * __restrict__ peak_ang_temp, 
  double * __restrict__ peak_ang_thmean, double * __restrict__ peak_ang_sum_s, 
  double * __restrict__ peak_ang_sum_c, double * __restrict__ cimsstrn_xki, 
  double * __restrict__ cimsstrn_e, double * __restrict__ cimsstrn_sume, 
  double * __restrict__ stokesdrift_stfac, double * __restrict__ wnfluxes_xstress, 
  double * __restrict__ wnfluxes_ystress, double * __restrict__ wnfluxes_ustar, 
  double * __restrict__ wnfluxes_philf, double * __restrict__ wnfluxes_ooval, 
  double * __restrict__ wnfluxes_em_oc, double * __restrict__ wnfluxes_f1_oc, 
  double * __restrict__ wnfluxes_cmrhowgdfth, double * __restrict__ wnfluxes_sumt, 
  double * __restrict__ wnfluxes_sumx, double * __restrict__ wnfluxes_sumy) {
  
  // needed for Loki
  // needed for Loki
  
  
  // ----------------------------------------------------------------------
  













  
  
  
  
  int ij;
  int k;
  int m;
  int icall;
  int ncall;
  
  double delt;
  double deltm;
  double ximp;
  double delt5;
  double gtemp1;
  double gtemp2;
  double flhab;
  double zhook_handle;
  
  //     *FLD* DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
  //     *SL*  TOTAL SOURCE FUNCTION ARRAY.
  //     *SPOS* : POSITIVE SINPUT ONLY
  
  int lcflx;
  int lupdtus;
  int ichnk;
  struct dim3 blockdim;
  struct dim3 griddim;
  ij = threadIdx.x;
  ichnk = blockIdx.x;
  
  if (ichnk <= nchnk && ij <= kijl) {
    ij = ij + 1;
    ichnk = ichnk + 1;
    
    // ----------------------------------------------------------------------
    
    
    
    //*    1. INITIALISATION.
    //        ---------------
    
    delt = idelt;
    deltm = (double) 1.0 / delt;
    ximp = (double) 1.0;
    delt5 = ximp*delt;
    
    lcflx = lwflux || lwfluxout || lwnemocou;
    
    
    
    raorw[ij - 1 + kijl*(ichnk - 1)] = 
      fmax((double) (aird[ij - 1 + kijl*(ichnk - 1)]), (double) ((double) 1.0))*rowaterm1
      ;
    
    for (k = 1; k <= nang; k += 1) {
      coswdif[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] = 
        cos(th[k - 1] - wdwave[ij - 1 + kijl*(ichnk - 1)]);
      sinwdif2[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] = 
        pow(sin(th[k - 1] - wdwave[ij - 1 + kijl*(ichnk - 1)]), 2);
    }
    
    
    // ----------------------------------------------------------------------
    
    //*    2. COMPUTATION OF IMPLICIT INTEGRATION.
    //        ------------------------------------
    
    //         INTEGRATION IS DONE FROM CDATE UNTIL CDTPRO FOR A BLOCK
    //         OF LATITUDES BETWEEN PROPAGATION CALLS.
    
    
    //     REDUCE WAVE ENERGY IF LARGER THAN DEPTH LIMITED WAVE HEIGHT
    if (lbiwbk) {
      sdepthlim_c(kijs, kijl, emaxdpt, fl1, delth, dfim, epsmin, fr, nang, nfre, wetail, 
        ichnk, nchnk, ij, sdepthlim_em, semean_temp);
    }
    
    //*    2.2 COMPUTE MEAN PARAMETERS.
    //        ------------------------
    
    fkmean_c(kijs, kijl, fl1, wavnum,  (&emean[ + kijl*(ichnk - 1)]), 
       (&fmean[ + kijl*(ichnk - 1)]),  (&f1mean[ + kijl*(ichnk - 1)]), 
       (&akmean[ + kijl*(ichnk - 1)]),  (&xkmean[ + kijl*(ichnk - 1)]), delth, dfim, 
      dfimfr, dfimofr, epsmin, fr, frtail, g, nang, nfre, wetail, wp1tail, zpi, ichnk, 
      nchnk, ij, fkmean_tempa, fkmean_tempx, fkmean_temp2);
    
    
    for (k = 1; k <= nang; k += 1) {
      flm[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))] = flmin*(pow(fmax((double) ((double) 
        0.0), (double) (coswdif[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))])), 2));
    }
    
    
    //     COMPUTE DAMPING COEFFICIENT DUE TO FRICTION ON BOTTOM OF THE SEA ICE.
    //!! testing sea ice attenuation (might need to restrict usage when needed)
    if (lciwabr) {
      ciwabr_c(kijs, kijl, cicover, fl1, wavnum, cgroup, 
         (&cireduc[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), cdicwa, dfim, epsmin, 
        idelt, licerun, lmaskice, nang, nfre, ichnk, nchnk, ij);
      
      for (m = 1; m <= nfre; m += 1) {
        for (k = 1; k <= nang; k += 1) {
          cireduc[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = ciwa[ij - 1 
            + kijl*(m - 1 + nfre*(ichnk - 1))]*cireduc[ij - 1 + kijl*(k - 1 + nang*(m - 1
             + nfre*(ichnk - 1)))];
        }
      }
      
    } else {
      
      for (m = 1; m <= nfre; m += 1) {
        for (k = 1; k <= nang; k += 1) {
          cireduc[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = 
            ciwa[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))];
        }
      }
      
    }
    
    // ----------------------------------------------------------------------
    
    //*    2.3 COMPUTATION OF SOURCE FUNCTIONS.
    //         --------------------------------
    
    //*    2.3.1 ITERATIVELY UPDATE STRESS AND COMPUTE WIND INPUT TERMS.
    //           -------------------------------------------------------
    
    lupdtus = true;
    ncall = 2;
    for (icall = 1; icall <= ncall; icall += 1) {
      sinflx_c(icall, ncall, kijs, kijl, lupdtus, fl1, wavnum, cinv, xk2cg, wswave, 
        wdwave, aird,  (&raorw[ + kijl*(ichnk - 1)]), wstar, cicover, 
         (&coswdif[ + kijl*( + nang*(ichnk - 1))]), 
         (&sinwdif2[ + kijl*( + nang*(ichnk - 1))]),  (&fmean[ + kijl*(ichnk - 1)]), 
         (&halp[ + kijl*(ichnk - 1)]),  (&fmeanws[ + kijl*(ichnk - 1)]), 
         (&flm[ + kijl*( + nang*(ichnk - 1))]), ufric, tauw, tauwdir, z0m, z0b, chrnck, 
         (&phiwa[ + kijl*(ichnk - 1)]), 
         (&fld[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), 
         (&sl[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), 
         (&spos[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), mij, 
         (&rhowgdfth[ + kijl*( + nfre*(ichnk - 1))]), xllws, abmax, abmin, acd, acdlin, 
        alpha, alphamax, alphamin, alphapmax, ang_gc_a, ang_gc_b, ang_gc_c, bcd, bcdlin, 
        betamaxoxkappa2, bmaxokap, c2osqrtvg_gc, cdmax, chnkmin_u, cithrsh_tail, cm_gc, 
        costh, delkcc_gc_ns, delkcc_omxkm3_gc, delth, dfim, dfimofr, dthrn_a, dthrn_u, 
        eps1, epsmin, epsus, flogsprdm1, fr, fr5, fric, frtail, g, gamnconst, gm1, iab, 
        icode, icode_cpl, idamping, iphys, jtot_tauhf, llcapchnk, llgcbz0, llnormagam, 
        lwcou, nang, nfre, nwav_gc, om3gmkm_gc, omega_gc, omxkm3_gc, rhowg_dfim, rn1_rn, 
        rnu, rnum, sinth, sqrtgosurft, swellf, swellf2, swellf3, swellf4, swellf5, 
        swellf6, swellf7, swellf7m1, swellft, tailfactor, tailfactor_pm, tauwshelter, 
        th, wetail, wspmin, wtauhf, x0tauhf, xkappa, xkmsqrtvgoc2_gc, xkm_gc, xk_gc, 
        xlogkratiom1_gc, xnlev, z0rat, z0tubmax, zalp, zpi, zpi4gm1, zpi4gm2, zpifr, 
        ichnk, nchnk, ij, sinflx_rnfac, sinflx_tmp_em, taut_z0_alphaog, taut_z0_xmin, 
        taut_z0_w1, taut_z0_tauwact, taut_z0_tauweff, taut_z0_ang_gc, taut_z0_tauunr, 
        taut_z0_llcosdiff, stress_gc_gam_w, z0wave_alphaog, femeanws_temp2, 
        femeanws_em_loc, halphap_alphap, halphap_xmss, halphap_em, halphap_fm, 
        halphap_f1d, halphap_wd, halphap_flwd, femean_temp2, meansqs_lf_fd, 
        meansqs_lf_temp1, meansqs_lf_temp2, sinput_ard_constf, sinput_ard_const11, 
        sinput_ard_const22, sinput_ard_z0vis, sinput_ard_z0noz, sinput_ard_fww, 
        sinput_ard_pvisc, sinput_ard_pturb, sinput_ard_zcn, sinput_ard_sig_n, 
        sinput_ard_uorbt, sinput_ard_aorb, sinput_ard_temp, sinput_ard_re, 
        sinput_ard_re_c, sinput_ard_zorb, sinput_ard_cnsn, sinput_ard_sumf, 
        sinput_ard_sumfsin2, sinput_ard_cstrnfac, sinput_ard_flp_avg, 
        sinput_ard_slp_avg, sinput_ard_rogoroair, sinput_ard_aird_pvisc, 
        sinput_ard_xstress, sinput_ard_ystress, sinput_ard_flp, sinput_ard_slp, 
        sinput_ard_usg2, sinput_ard_taux, sinput_ard_tauy, sinput_ard_ustp, 
        sinput_ard_ustpm1, sinput_ard_usdirp, sinput_ard_ucn, sinput_ard_ucnzalpd, 
        sinput_ard_xngamconst, sinput_ard_gamnorma, sinput_ard_dstab1, sinput_ard_temp1, 
        sinput_ard_temp2, sinput_ard_gam0, sinput_ard_dstab, sinput_ard_coslp, 
        sinput_jan_ztanhkd, sinput_jan_sig_n, sinput_jan_cnsn, sinput_jan_sumf, 
        sinput_jan_sumfsin2, sinput_jan_cstrnfac, sinput_jan_xngamconst, 
        sinput_jan_gamnorma, sinput_jan_sigdev, sinput_jan_us, sinput_jan_z0, 
        sinput_jan_ucn, sinput_jan_zcn, sinput_jan_ustpm1, sinput_jan_xvd, 
        sinput_jan_ucnd, sinput_jan_const3_ucn2, sinput_jan_ufac1, sinput_jan_ufac2, 
        sinput_jan_tempd, sinput_jan_gam0, sinput_jan_lz, stresso_xstress, 
        stresso_ystress, stresso_tauhf, stresso_phihf, stresso_cmrhowgdfth, stresso_us2, 
        stresso_taux, stresso_tauy, stresso_taupx, stresso_taupy, stresso_usdirp, 
        stresso_ust, stresso_sumt, stresso_sumx, stresso_sumy, tau_phi_hf_ns, 
        tau_phi_hf_xks, tau_phi_hf_oms, tau_phi_hf_sqrtz0og, tau_phi_hf_zsup, 
        tau_phi_hf_zinf, tau_phi_hf_delz, tau_phi_hf_taul, tau_phi_hf_xloggz0, 
        tau_phi_hf_sqrtgz0, tau_phi_hf_ustph, tau_phi_hf_const1, tau_phi_hf_const2, 
        tau_phi_hf_consttau, tau_phi_hf_constphi, tau_phi_hf_f1dcos2, 
        tau_phi_hf_f1dcos3, tau_phi_hf_f1d, tau_phi_hf_f1dsin2);
      
    }
    
    //     2.3.3 ADD THE OTHER SOURCE TERMS.
    //           ---------------------------
    
    sdissip_c(kijs, kijl, fl1,  (&fld[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), 
       (&sl[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), wavnum, cgroup, xk2cg, 
       (&emean[ + kijl*(ichnk - 1)]),  (&f1mean[ + kijl*(ichnk - 1)]), 
       (&xkmean[ + kijl*(ichnk - 1)]), ufric, 
       (&coswdif[ + kijl*( + nang*(ichnk - 1))]),  (&raorw[ + kijl*(ichnk - 1)]), 
      brkpbcoef, cdis, cdisvis, delta_sdis, delth, fratio, g, indicessat, iphys, ipsat, 
      miche, nang, nfre, nsdsnth, rnu, satweights, sdsbr, ssdsbrf1, ssdsc2, ssdsc3, 
      ssdsc4, ssdsc5, ssdsc6, zpi, zpifr, ichnk, nchnk, ij, sdissip_ard_facturb, 
      sdissip_ard_facsat, sdissip_ard_facwtrb, sdissip_ard_temp1, sdissip_ard_bth0, 
      sdissip_ard_c_, sdissip_ard_c_c, sdissip_ard_dsip, sdissip_ard_trpz_dsip, 
      sdissip_ard_bth, sdissip_ard_temp2, sdissip_ard_d, sdissip_ard_scumul, 
      sdissip_ard_renewalfreq, sdissip_ard_wcumul, sdissip_jan_temp1, sdissip_jan_sds, 
      sdissip_jan_x, sdissip_jan_xk2);
    
    //     Save source term contributions relevant for the calculation of ocean fluxes
    
    if (lcflx && !lwvflx_snl) {
      for (m = 1; m <= nfre; m += 1) {
        for (k = 1; k <= nang; k += 1) {
          ssource[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = 
            sl[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))];
        }
      }
    }
    
    
    snonlin_c(kijs, kijl, fl1,  (&fld[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), 
       (&sl[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), wavnum, depth, 
       (&akmean[ + kijl*(ichnk - 1)]), af11, bathymax, costh, dal1, dal2, delth, dfim, 
      dfimfr, dfimfr2, dkmax, fklam, fklam1, fklap, fklap1, fr, fratio, g, gm1, ikm, 
      ikm1, ikp, ikp1, inlcoef, isnonlin, k11w, k1w, k21w, k2w, kfrh, mfrstlw, mlsthg, 
      nang, nfre, rnlcoef, sinth, th, wetail, wp1tail, wp2tail, xkdmin, zpifr, ichnk, 
      nchnk, ij, snonlin_enh, snonlin_xnu, snonlin_sig_th, snonlin_ftemp, snonlin_ad, 
      snonlin_delad, snonlin_delap, snonlin_delam, snonlin_enhfr, peak_ang_mmax, 
      peak_ang_sum0, peak_ang_sum1, peak_ang_sum2, peak_ang_xmax, peak_ang_temp, 
      peak_ang_thmean, peak_ang_sum_s, peak_ang_sum_c);
    
    
    if (lcflx && lwvflx_snl) {
      //     Save source term contributions relevant for the calculation of ocean fluxes
      //!!!!!  SL must only contain contributions contributed to fluxes into the oceans
      //       MODULATE SL BY IMPLICIT FACTOR
      for (m = 1; m <= nfre; m += 1) {
        for (k = 1; k <= nang; k += 1) {
          gtemp1 = fmax((double) (((double) 1.0 - delt5*fld[ij - 1 + kijl*(k - 1 + 
            nang*(m - 1 + nfre*(ichnk - 1)))])), (double) ((double) 1.0));
          ssource[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = 
            sl[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] / gtemp1;
        }
      }
    }
    
    
    
    sdiwbk_c(kijs, kijl, fl1,  (&fld[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), 
       (&sl[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), depth, emaxdpt, 
       (&emean[ + kijl*(ichnk - 1)]),  (&f1mean[ + kijl*(ichnk - 1)]), lbiwbk, nang, 
      nfre, nfre_red, ichnk, nchnk, ij, sdiwbk_sds);
    
    sbottom_c(kijs, kijl, fl1,  (&fld[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), 
       (&sl[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), wavnum, depth, bathymax, gm1, 
      nang, nfre, nfre_red, ichnk, nchnk, ij, sbottom_sbo);
    
    // ----------------------------------------------------------------------
    
    //*    2.4 COMPUTATION OF NEW SPECTRA.
    //         ---------------------------
    
    //     INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE
    //     FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.
    
    for (m = 1; m <= nfre; m += 1) {
      delfl[m - 1] = cofrm4[m - 1]*delt;
    }
    
    usfm[ij - 1 + kijl*(ichnk - 1)] = ufric[ij - 1 + kijl*(ichnk - 1)]*fmax((double) 
      (fmeanws[ij - 1 + kijl*(ichnk - 1)]), (double) (fmean[ij - 1 + kijl*(ichnk - 1)]));
    
    for (m = 1; m <= nfre; m += 1) {
      temp[ij - 1 + kijl*(m - 1 + nfre*(ichnk - 1))] = 
        usfm[ij - 1 + kijl*(ichnk - 1)]*delfl[m - 1];
    }
    
    if (llunstr) {
      for (k = 1; k <= nang; k += 1) {
        for (m = 1; m <= nfre; m += 1) {
          gtemp1 = fmax((double) (((double) 1.0 - delt5*fld[ij - 1 + kijl*(k - 1 + 
            nang*(m - 1 + nfre*(ichnk - 1)))])), (double) ((double) 1.0));
          gtemp2 = 
            delt*sl[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] / gtemp1;
          flhab = fabs((double) (gtemp2));
          flhab = fmin((double) (flhab), (double) (temp[ij - 1 + kijl*(m - 1 + 
            nfre*(ichnk - 1))]));
          fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = fl1[ij - 1 + 
            kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] + iobnd[ij - 1 + kijl*(ichnk 
            - 1)]*copysign((double) (flhab), (double) (gtemp2));
          fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = fmax((double) 
            (iodp[ij - 1 + kijl*(ichnk - 1)]*cireduc[ij - 1 + kijl*(k - 1 + nang*(m - 1 +
             nfre*(ichnk - 1)))]*fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1
            )))]), (double) (flm[ij - 1 + kijl*(k - 1 + nang*(ichnk - 1))]));
          ssource[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = ssource[ij -
             1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] + deltm*fmin((double) 
            (flmax[m - 1] - fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))])
            , (double) ((double) 0.0));
          fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = fmin((double) 
            (fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))]), (double) 
            (flmax[m - 1]));
        }
      }
    } else {
      for (k = 1; k <= nang; k += 1) {
        for (m = 1; m <= nfre; m += 1) {
          gtemp1 = fmax((double) (((double) 1.0 - delt5*fld[ij - 1 + kijl*(k - 1 + 
            nang*(m - 1 + nfre*(ichnk - 1)))])), (double) ((double) 1.0));
          gtemp2 = 
            delt*sl[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] / gtemp1;
          flhab = fabs((double) (gtemp2));
          flhab = fmin((double) (flhab), (double) (temp[ij - 1 + kijl*(m - 1 + 
            nfre*(ichnk - 1))]));
          fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = fl1[ij - 1 + 
            kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] + copysign((double) (flhab), 
            (double) (gtemp2));
          fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = fmax((double) 
            (cireduc[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))]*fl1[ij - 1 
            + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))]), (double) (flm[ij - 1 + 
            kijl*(k - 1 + nang*(ichnk - 1))]));
          ssource[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = ssource[ij -
             1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] + deltm*fmin((double) 
            (flmax[m - 1] - fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))])
            , (double) ((double) 0.0));
          fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))] = fmin((double) 
            (fl1[ij - 1 + kijl*(k - 1 + nang*(m - 1 + nfre*(ichnk - 1)))]), (double) 
            (flmax[m - 1]));
        }
      }
    }
    
    
    if (lcflx) {
      wnfluxes_c(kijs, kijl, mij,  (&rhowgdfth[ + kijl*( + nfre*(ichnk - 1))]), cinv, 
         (&ssource[ + kijl*( + nang*( + nfre*(ichnk - 1)))]), cicover, 
         (&phiwa[ + kijl*(ichnk - 1)]),  (&emean[ + kijl*(ichnk - 1)]), 
         (&f1mean[ + kijl*(ichnk - 1)]), wswave, wdwave, ustra, vstra, ufric, aird, 
        nphieps, ntauoc, nswh, nmwp, nemotaux, nemotauy, nemowswave, nemophif, tauxd, 
        tauyd, tauocxd, tauocyd, tauoc, phiocd, phieps, phiaw, true, afcrv, bfcrv, 
        ciblock, cithrsh, costh, egrcrv, epsu10, epsus, fr, g, licerun, lwamrsetci, 
        lwcouast, lwnemocou, lwnemotauoc, nang, nfre, phiepsmax, phiepsmin, sinth, 
        tauocmax, tauocmin, ichnk, nchnk, ij, wnfluxes_xstress, wnfluxes_ystress, 
        wnfluxes_ustar, wnfluxes_philf, wnfluxes_ooval, wnfluxes_em_oc, wnfluxes_f1_oc, 
        wnfluxes_cmrhowgdfth, wnfluxes_sumt, wnfluxes_sumx, wnfluxes_sumy);
    }
    // ----------------------------------------------------------------------
    
    //*    2.5 REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.
    //         -----------------------------------------------------
    
    fkmean_c(kijs, kijl, fl1, wavnum,  (&emean[ + kijl*(ichnk - 1)]), 
       (&fmean[ + kijl*(ichnk - 1)]),  (&f1mean[ + kijl*(ichnk - 1)]), 
       (&akmean[ + kijl*(ichnk - 1)]),  (&xkmean[ + kijl*(ichnk - 1)]), delth, dfim, 
      dfimfr, dfimofr, epsmin, fr, frtail, g, nang, nfre, wetail, wp1tail, zpi, ichnk, 
      nchnk, ij, fkmean_tempa, fkmean_tempx, fkmean_temp2);
    
    //     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
    femeanws_c(kijs, kijl, fl1, xllws,  (&fmeanws[ + kijl*(ichnk - 1)]), 
       (&emeanws[ + kijl*(ichnk - 1)]), delth, dfim, dfimofr, epsmin, fr, frtail, nang, 
      nfre, wetail, ichnk, nchnk, ij, femeanws_temp2, femeanws_em_loc);
    
    imphftail_c(kijs, kijl, mij,  (&flm[ + kijl*( + nang*(ichnk - 1))]), wavnum, xk2cg, 
      fl1, nang, nfre, ichnk, nchnk, ij, imphftail_temp1, imphftail_temp2);
    
    
    //     UPDATE WINDSEA VARIANCE AND MEAN FREQUENCY IF PASSED TO ATMOSPHERE
    //     ------------------------------------------------------------------
    
    if (lwflux) {
      if (emeanws[ij - 1 + kijl*(ichnk - 1)] < wsemean_min) {
        wsemean[ij - 1 + kijl*(ichnk - 1)] = wsemean_min;
        wsfmean[ij - 1 + kijl*(ichnk - 1)] = (double) 2.*fr[nfre - 1];
      } else {
        wsemean[ij - 1 + kijl*(ichnk - 1)] = emeanws[ij - 1 + kijl*(ichnk - 1)];
        wsfmean[ij - 1 + kijl*(ichnk - 1)] = fmeanws[ij - 1 + kijl*(ichnk - 1)];
      }
    }
    
    
    
    //*    2.6 SET FL1 ON ICE POINTS TO ZERO
    //         -----------------------------
    
    if (licerun && lmaskice) {
      setice_c(kijs, kijl, fl1, cicover,  (&coswdif[ + kijl*( + nang*(ichnk - 1))]), 
        cithrsh, epsmin, flmin, nang, nfre, ichnk, nchnk, ij, setice_cireduc, 
        setice_temp, setice_icefree);
    }
    
    
    //*    2.7 SURFACE STOKES DRIFT AND STRAIN IN SEA ICE
    //         ------------------------------------------
    
    stokestrn_c(kijs, kijl, fl1, wavnum, stokfac, depth, wswave, wdwave, cicover, 
      cithick, ustokes, vstokes, strnms, nemoustokes, nemovstokes, nemostrn, cithrsh, 
      costh, delth, dfim, dfim_sim, flmin, fr, g, licerun, lwamrsetci, lwcou, lwnemocou, 
      lwnemocousend, lwnemocoustk, lwnemocoustrn, nang, nfre, nfre_odd, rowater, sinth, 
      zpi, ichnk, nchnk, ij, cimsstrn_xki, cimsstrn_e, cimsstrn_sume, stokesdrift_stfac);
    
    // ----------------------------------------------------------------------
    
    
  }
}
