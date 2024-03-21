#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "stokestrn_c.h"
#include "setice_c.h"
#include "imphftail_c.h"
#include "femeanws_c.h"
#include "wnfluxes_c.h"
#include "sbottom_c.h"
#include "sdiwbk_c.h"
#include "snonlin_c.h"
#include "sdissip_c.h"
#include "sinflx_c.h"
#include "ciwabr_c.h"
#include "fkmean_c.h"
#include "sdepthlim_c.h"
__global__ void __launch_bounds__(128, 1) implsch_c(int kijs, int kijl, double * fl1, 
  const double * wavnum, const double * cgroup, const double * ciwa, 
  const double * cinv, const double * xk2cg, const double * stokfac, 
  const double * emaxdpt, const int * indep, const double * depth, const int * iobnd, 
  const int * iodp, const double * aird, const double * wdwave, const double * cicover, 
  double * wswave, const double * wstar, double * ufric, double * tauw, 
  double * tauwdir, double * z0m, double * z0b, double * chrnck, const double * cithick, 
  double * nemoustokes, double * nemovstokes, double * nemostrn, double * nphieps, 
  double * ntauoc, double * nswh, double * nmwp, double * nemotaux, double * nemotauy, 
  double * nemowswave, double * nemophif, double * wsemean, double * wsfmean, 
  double * ustokes, double * vstokes, double * strnms, double * tauxd, double * tauyd, 
  double * tauocxd, double * tauocyd, double * tauoc, double * phiocd, double * phieps, 
  double * phiaw, int * mij, double * xllws, double abmax, double abmin, double acd, 
  double acdlin, const double * af11, double afcrv, double alpha, double alphamax, 
  double alphamin, double alphapmax, double ang_gc_a, double ang_gc_b, double ang_gc_c, 
  double bathymax, double bcd, double bcdlin, double betamaxoxkappa2, double bfcrv, 
  double bmaxokap, const double * c2osqrtvg_gc, double cdicwa, double cdis, 
  double cdisvis, double cdmax, double chnkmin_u, double ciblock, double cithrsh, 
  double cithrsh_tail, const double * cm_gc, const double * cofrm4, 
  const double * costh, const double * cumulw, double dal1, double dal2, 
  const double * delkcc_gc_ns, const double * delkcc_omxkm3_gc, double delta_sdis, 
  double delth, const double * dfim, const double * dfimfr, const double * dfimfr2, 
  const double * dfimofr, const double * dfim_sim, double dkmax, double dthrn_a, 
  double dthrn_u, double egrcrv, double eps1, double epsmin, double epsu10, 
  double epsus, const double * fklam, const double * fklam1, const double * fklap, 
  const double * fklap1, const double * flmax, double flmin, double flogsprdm1, 
  const double * fr, const double * fr5, double fratio, double fric, double frtail, 
  double g, double gamnconst, double gm1, int iab, int icode, int icode_cpl, 
  int idamping, int idelt, const int * ikm, const int * ikm1, const int * ikp, 
  const int * ikp1, const int * indicessat, const int * inlcoef, int iphys, int ipsat, 
  int isnonlin, int jtot_tauhf, const int * k11w, const int * k1w, const int * k21w, 
  const int * k2w, int kfrh, int lbiwbk, int lciwabr, int licerun, int llcapchnk, 
  int llgcbz0, int llnormagam, int llunstr, int lmaskice, int lwamrsetci, int lwcou, 
  int lwflux, int lwfluxout, int lwnemocou, int lwnemocousend, int lwnemocoustk, 
  int lwnemocoustrn, int lwnemotauoc, int lwvflx_snl, int mfrstlw, double miche, 
  int mlsthg, int nang, int ndepth, int ndikcumul, int nfre, int nfre_odd, int nfre_red, 
  int nsdsnth, int nwav_gc, const double * om3gmkm_gc, const double * omega_gc, 
  const double * omxkm3_gc, double phiepsmax, double phiepsmin, 
  const double * rhowg_dfim, double rn1_rn, const double * rnlcoef, double rnu, 
  double rnum, double rowater, double rowaterm1, const double * satweights, 
  double sdsbr, const double * sinth, double sqrtgosurft, double ssdsc2, double ssdsc3, 
  double ssdsc4, double ssdsc5, double ssdsc6, double swellf, double swellf2, 
  double swellf3, double swellf4, double swellf5, double swellf6, double swellf7, 
  double swellf7m1, const double * swellft, double tailfactor, double tailfactor_pm, 
  double tauocmax, double tauocmin, double tauwshelter, const double * th, 
  double wetail, double wp1tail, double wp2tail, double wsemean_min, double wspmin, 
  const double * wtauhf, double x0tauhf, double xkappa, double xkdmin, 
  const double * xkmsqrtvgoc2_gc, const double * xkm_gc, const double * xk_gc, 
  double xlogkratiom1_gc, double xnlev, double z0rat, double z0tubmax, double zalp, 
  double zpi, double zpi4gm1, double zpi4gm2, const double * zpifr, int ichnk_start, 
  int ichnk_end, int ichnk_step, int nchnk, int nproma_wam, double * raorw, 
  double * emean, double * fmean, double * halp, double * emeanws, double * fmeanws, 
  double * f1mean, double * akmean, double * xkmean, double * phiwa, double * flm, 
  double * coswdif, double * sinwdif2, double * rhowgdfth, double * fld, double * sl, 
  double * spos, double * cireduc, double * ssource, double * sinflx_rnfac, 
  double * sinflx_tmp_em, double * stresso_xstress, double * stresso_ystress, 
  double * stresso_tauhf, double * stresso_phihf, double * stresso_usdirp, 
  double * stresso_ust, double * snonlin_enh, double * snonlin_xnu, 
  double * snonlin_sig_th);
#include "implsch_c_launch.h"

__global__ void implsch_c(int kijs, int kijl, double * fl1, const double * wavnum, 
  const double * cgroup, const double * ciwa, const double * cinv, const double * xk2cg, 
  const double * stokfac, const double * emaxdpt, const int * indep, 
  const double * depth, const int * iobnd, const int * iodp, const double * aird, 
  const double * wdwave, const double * cicover, double * wswave, const double * wstar, 
  double * ufric, double * tauw, double * tauwdir, double * z0m, double * z0b, 
  double * chrnck, const double * cithick, double * nemoustokes, double * nemovstokes, 
  double * nemostrn, double * nphieps, double * ntauoc, double * nswh, double * nmwp, 
  double * nemotaux, double * nemotauy, double * nemowswave, double * nemophif, 
  double * wsemean, double * wsfmean, double * ustokes, double * vstokes, 
  double * strnms, double * tauxd, double * tauyd, double * tauocxd, double * tauocyd, 
  double * tauoc, double * phiocd, double * phieps, double * phiaw, int * mij, 
  double * xllws, double abmax, double abmin, double acd, double acdlin, 
  const double * af11, double afcrv, double alpha, double alphamax, double alphamin, 
  double alphapmax, double ang_gc_a, double ang_gc_b, double ang_gc_c, double bathymax, 
  double bcd, double bcdlin, double betamaxoxkappa2, double bfcrv, double bmaxokap, 
  const double * c2osqrtvg_gc, double cdicwa, double cdis, double cdisvis, double cdmax, 
  double chnkmin_u, double ciblock, double cithrsh, double cithrsh_tail, 
  const double * cm_gc, const double * cofrm4, const double * costh, 
  const double * cumulw, double dal1, double dal2, const double * delkcc_gc_ns, 
  const double * delkcc_omxkm3_gc, double delta_sdis, double delth, const double * dfim, 
  const double * dfimfr, const double * dfimfr2, const double * dfimofr, 
  const double * dfim_sim, double dkmax, double dthrn_a, double dthrn_u, double egrcrv, 
  double eps1, double epsmin, double epsu10, double epsus, const double * fklam, 
  const double * fklam1, const double * fklap, const double * fklap1, 
  const double * flmax, double flmin, double flogsprdm1, const double * fr, 
  const double * fr5, double fratio, double fric, double frtail, double g, 
  double gamnconst, double gm1, int iab, int icode, int icode_cpl, int idamping, 
  int idelt, const int * ikm, const int * ikm1, const int * ikp, const int * ikp1, 
  const int * indicessat, const int * inlcoef, int iphys, int ipsat, int isnonlin, 
  int jtot_tauhf, const int * k11w, const int * k1w, const int * k21w, const int * k2w, 
  int kfrh, int lbiwbk, int lciwabr, int licerun, int llcapchnk, int llgcbz0, 
  int llnormagam, int llunstr, int lmaskice, int lwamrsetci, int lwcou, int lwflux, 
  int lwfluxout, int lwnemocou, int lwnemocousend, int lwnemocoustk, int lwnemocoustrn, 
  int lwnemotauoc, int lwvflx_snl, int mfrstlw, double miche, int mlsthg, int nang, 
  int ndepth, int ndikcumul, int nfre, int nfre_odd, int nfre_red, int nsdsnth, 
  int nwav_gc, const double * om3gmkm_gc, const double * omega_gc, 
  const double * omxkm3_gc, double phiepsmax, double phiepsmin, 
  const double * rhowg_dfim, double rn1_rn, const double * rnlcoef, double rnu, 
  double rnum, double rowater, double rowaterm1, const double * satweights, 
  double sdsbr, const double * sinth, double sqrtgosurft, double ssdsc2, double ssdsc3, 
  double ssdsc4, double ssdsc5, double ssdsc6, double swellf, double swellf2, 
  double swellf3, double swellf4, double swellf5, double swellf6, double swellf7, 
  double swellf7m1, const double * swellft, double tailfactor, double tailfactor_pm, 
  double tauocmax, double tauocmin, double tauwshelter, const double * th, 
  double wetail, double wp1tail, double wp2tail, double wsemean_min, double wspmin, 
  const double * wtauhf, double x0tauhf, double xkappa, double xkdmin, 
  const double * xkmsqrtvgoc2_gc, const double * xkm_gc, const double * xk_gc, 
  double xlogkratiom1_gc, double xnlev, double z0rat, double z0tubmax, double zalp, 
  double zpi, double zpi4gm1, double zpi4gm2, const double * zpifr, int ichnk_start, 
  int ichnk_end, int ichnk_step, int nchnk, int nproma_wam, double * raorw, 
  double * emean, double * fmean, double * halp, double * emeanws, double * fmeanws, 
  double * f1mean, double * akmean, double * xkmean, double * phiwa, double * flm, 
  double * coswdif, double * sinwdif2, double * rhowgdfth, double * fld, double * sl, 
  double * spos, double * cireduc, double * ssource, double * sinflx_rnfac, 
  double * sinflx_tmp_em, double * stresso_xstress, double * stresso_ystress, 
  double * stresso_tauhf, double * stresso_phihf, double * stresso_usdirp, 
  double * stresso_ust, double * snonlin_enh, double * snonlin_xnu, 
  double * snonlin_sig_th) {
  
  
  













  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  const int nrnl = 25;
  const int ninl = 5;
  
  
  int ij;
  int k;
  int m;
  
  double delt;
  double deltm;
  double ximp;
  double delt5;
  double gtemp1;
  double gtemp2;
  double flhab;
  double usfm;
  
  
  int lcflx;
  int ichnk;
  dim3 blockdim;
  dim3 griddim;
  ij = threadIdx.x;
  ichnk = blockIdx.x;

  if (ichnk <= nchnk && ij <= kijl) {
    // blockdim = dim3(nproma_wam, 1, 1);
    // griddim = dim3(nchnk, 1, 1);

    ichnk++;
    ij++;
    // START of Loki inserted loop ICHNK

    delt = idelt;
    deltm = (double) 1.0 / delt;
    ximp = (double) 1.0;
    delt5 = ximp*delt;
    
    lcflx = lwflux || lwfluxout || lwnemocou;

    raorw[ij - 1 + kijl*(ichnk - 1)] = 
      max((double) (aird[ij - 1 + kijl*(ichnk - 1)]), (double) ((double) 1.0))*rowaterm1;
    
    for (k = 1; k <= nang; k += 1) {
      coswdif[ij - 1 + kijl*(k - 1 + nang_loki_param*(ichnk - 1))] = 
        cos(th[k - 1] - wdwave[ij - 1 + kijl*(ichnk - 1)]);
      sinwdif2[ij - 1 + kijl*(k - 1 + nang_loki_param*(ichnk - 1))] = 
        pow(sin(th[k - 1] - wdwave[ij - 1 + kijl*(ichnk - 1)]), 2);
    }

    if (lbiwbk) {
      sdepthlim_c(kijs, kijl, emaxdpt, fl1, delth, dfim, epsmin, fr, nang, nfre, wetail, 
        ichnk, nchnk, ij);
    }
    fkmean_c(kijs, kijl, fl1, wavnum,  (&emean[ + kijl*(ichnk - 1)]), 
       (&fmean[ + kijl*(ichnk - 1)]),  (&f1mean[ + kijl*(ichnk - 1)]), 
       (&akmean[ + kijl*(ichnk - 1)]),  (&xkmean[ + kijl*(ichnk - 1)]), delth, dfim, 
      dfimfr, dfimofr, epsmin, fr, frtail, g, nang, nfre, wetail, wp1tail, zpi, ichnk, 
      nchnk, ij);
    

    for (k = 1; k <= nang; k += 1) {
      flm[ij - 1 + kijl*(k - 1 + nang_loki_param*(ichnk - 1))] = flmin*(pow(max((double) 
        ((double) 0.0), (double) (coswdif[ij - 1 + kijl*(k - 1 + nang_loki_param*(ichnk -
         1))])), 2));
    }

    if (lciwabr) {
      ciwabr_c(kijs, kijl, cicover, fl1, wavnum, cgroup, 
         (&cireduc[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), 
        cdicwa, dfim, epsmin, idelt, licerun, lmaskice, nang, nfre, ichnk, nchnk, ij);

      for (m = 1; m <= nfre; m += 1) {
        for (k = 1; k <= nang; k += 1) {
          cireduc[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk 
            - 1)))] = ciwa[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))
            ]*cireduc[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
            nfre_loki_param*(ichnk - 1)))];
        }
      }

    } else {

      for (m = 1; m <= nfre; m += 1) {
        for (k = 1; k <= nang; k += 1) {
          cireduc[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk 
            - 1)))] = ciwa[ij - 1 + kijl*(m - 1 + nfre_loki_param*(ichnk - 1))];
        }
      }

    }
    sinflx_c(1, kijs, kijl, true, fl1, wavnum, cinv, xk2cg, wswave, wdwave, aird, 
       (&raorw[ + kijl*(ichnk - 1)]), wstar, cicover, 
       (&coswdif[ + kijl*( + nang_loki_param*(ichnk - 1))]), 
       (&sinwdif2[ + kijl*( + nang_loki_param*(ichnk - 1))]), 
       (&fmean[ + kijl*(ichnk - 1)]),  (&halp[ + kijl*(ichnk - 1)]), 
       (&fmeanws[ + kijl*(ichnk - 1)]), 
       (&flm[ + kijl*( + nang_loki_param*(ichnk - 1))]), ufric, tauw, tauwdir, z0m, z0b, 
      chrnck,  (&phiwa[ + kijl*(ichnk - 1)]), 
       (&fld[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), 
       (&sl[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), 
       (&spos[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), mij, 
       (&rhowgdfth[ + kijl*( + nfre_loki_param*(ichnk - 1))]), xllws, abmax, abmin, acd, 
      acdlin, alpha, alphamax, alphamin, alphapmax, ang_gc_a, ang_gc_b, ang_gc_c, bcd, 
      bcdlin, betamaxoxkappa2, bmaxokap, c2osqrtvg_gc, cdmax, chnkmin_u, cithrsh_tail, 
      cm_gc, costh, delkcc_gc_ns, delkcc_omxkm3_gc, delth, dfim, dfimofr, dthrn_a, 
      dthrn_u, eps1, epsmin, epsus, flogsprdm1, fr, fr5, fric, frtail, g, gamnconst, 
      gm1, iab, icode, icode_cpl, idamping, iphys, jtot_tauhf, llcapchnk, llgcbz0, 
      llnormagam, lwcou, nang, nfre, nwav_gc, om3gmkm_gc, omega_gc, omxkm3_gc, 
      rhowg_dfim, rn1_rn, rnu, rnum, sinth, sqrtgosurft, swellf, swellf2, swellf3, 
      swellf4, swellf5, swellf6, swellf7, swellf7m1, swellft, tailfactor, tailfactor_pm, 
      tauwshelter, th, wetail, wspmin, wtauhf, x0tauhf, xkappa, xkmsqrtvgoc2_gc, xkm_gc, 
      xk_gc, xlogkratiom1_gc, xnlev, z0rat, z0tubmax, zalp, zpi, zpi4gm1, zpi4gm2, 
      zpifr, ichnk, nchnk, ij, sinflx_rnfac, sinflx_tmp_em, stresso_xstress, 
      stresso_ystress, stresso_tauhf, stresso_phihf, stresso_usdirp, stresso_ust);
    sinflx_c(2, kijs, kijl, true, fl1, wavnum, cinv, xk2cg, wswave, wdwave, aird, 
       (&raorw[ + kijl*(ichnk - 1)]), wstar, cicover, 
       (&coswdif[ + kijl*( + nang_loki_param*(ichnk - 1))]), 
       (&sinwdif2[ + kijl*( + nang_loki_param*(ichnk - 1))]), 
       (&fmean[ + kijl*(ichnk - 1)]),  (&halp[ + kijl*(ichnk - 1)]), 
       (&fmeanws[ + kijl*(ichnk - 1)]), 
       (&flm[ + kijl*( + nang_loki_param*(ichnk - 1))]), ufric, tauw, tauwdir, z0m, z0b, 
      chrnck,  (&phiwa[ + kijl*(ichnk - 1)]), 
       (&fld[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), 
       (&sl[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), 
       (&spos[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), mij, 
       (&rhowgdfth[ + kijl*( + nfre_loki_param*(ichnk - 1))]), xllws, abmax, abmin, acd, 
      acdlin, alpha, alphamax, alphamin, alphapmax, ang_gc_a, ang_gc_b, ang_gc_c, bcd, 
      bcdlin, betamaxoxkappa2, bmaxokap, c2osqrtvg_gc, cdmax, chnkmin_u, cithrsh_tail, 
      cm_gc, costh, delkcc_gc_ns, delkcc_omxkm3_gc, delth, dfim, dfimofr, dthrn_a, 
      dthrn_u, eps1, epsmin, epsus, flogsprdm1, fr, fr5, fric, frtail, g, gamnconst, 
      gm1, iab, icode, icode_cpl, idamping, iphys, jtot_tauhf, llcapchnk, llgcbz0, 
      llnormagam, lwcou, nang, nfre, nwav_gc, om3gmkm_gc, omega_gc, omxkm3_gc, 
      rhowg_dfim, rn1_rn, rnu, rnum, sinth, sqrtgosurft, swellf, swellf2, swellf3, 
      swellf4, swellf5, swellf6, swellf7, swellf7m1, swellft, tailfactor, tailfactor_pm, 
      tauwshelter, th, wetail, wspmin, wtauhf, x0tauhf, xkappa, xkmsqrtvgoc2_gc, xkm_gc, 
      xk_gc, xlogkratiom1_gc, xnlev, z0rat, z0tubmax, zalp, zpi, zpi4gm1, zpi4gm2, 
      zpifr, ichnk, nchnk, ij, sinflx_rnfac, sinflx_tmp_em, stresso_xstress, 
      stresso_ystress, stresso_tauhf, stresso_phihf, stresso_usdirp, stresso_ust);
    sdissip_c(kijs, kijl, fl1, 
       (&fld[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), 
       (&sl[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), indep, 
      wavnum, xk2cg,  (&emean[ + kijl*(ichnk - 1)]),  (&f1mean[ + kijl*(ichnk - 1)]), 
       (&xkmean[ + kijl*(ichnk - 1)]), ufric, 
       (&coswdif[ + kijl*( + nang_loki_param*(ichnk - 1))]), 
       (&raorw[ + kijl*(ichnk - 1)]), cdis, cdisvis, cumulw, delta_sdis, g, indicessat, 
      iphys, ipsat, miche, nang, ndepth, ndikcumul, nfre, nsdsnth, rnu, satweights, 
      sdsbr, ssdsc2, ssdsc3, ssdsc4, ssdsc5, ssdsc6, zpi, zpifr, ichnk, nchnk, ij);

    if (lcflx && !lwvflx_snl) {
      for (m = 1; m <= nfre; m += 1) {
        for (k = 1; k <= nang; k += 1) {
          ssource[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk 
            - 1)))] = sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
            nfre_loki_param*(ichnk - 1)))];
        }
      }
    }

    
    snonlin_c(kijs, kijl, fl1, 
       (&fld[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), 
       (&sl[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), wavnum, 
      depth,  (&akmean[ + kijl*(ichnk - 1)]), af11, bathymax, costh, dal1, dal2, delth, 
      dfim, dfimfr, dfimfr2, dkmax, fklam, fklam1, fklap, fklap1, fr, fratio, g, gm1, 
      ikm, ikm1, ikp, ikp1, inlcoef, isnonlin, k11w, k1w, k21w, k2w, kfrh, mfrstlw, 
      mlsthg, nang, nfre, rnlcoef, sinth, th, wetail, wp1tail, wp2tail, xkdmin, zpifr, 
      ichnk, nchnk, ij, snonlin_enh, snonlin_xnu, snonlin_sig_th);
    

    if (lcflx && lwvflx_snl) {
      for (m = 1; m <= nfre; m += 1) {
        for (k = 1; k <= nang; k += 1) {
          gtemp1 = max((double) (((double) 1.0 - delt5*fld[ij - 1 + kijl*(k - 1 + 
            nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))])), (double) ((double)
             1.0));
          ssource[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk 
            - 1)))] = sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
            nfre_loki_param*(ichnk - 1)))] / gtemp1;
        }
      }
    }

    sdiwbk_c(kijs, kijl, fl1, 
       (&fld[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), 
       (&sl[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), depth, 
      emaxdpt,  (&emean[ + kijl*(ichnk - 1)]),  (&f1mean[ + kijl*(ichnk - 1)]), lbiwbk, 
      nang, nfre_red, ichnk, nchnk, ij);
    
    sbottom_c(kijs, kijl, fl1, 
       (&fld[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), 
       (&sl[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), wavnum, 
      depth, bathymax, gm1, nang, nfre_red, ichnk, nchnk, ij);

    usfm = ufric[ij - 1 + kijl*(ichnk - 1)]*max((double) (fmeanws[ij - 1 + kijl*(ichnk - 
      1)]), (double) (fmean[ij - 1 + kijl*(ichnk - 1)]));
    
    if (llunstr) {
      for (k = 1; k <= nang; k += 1) {
        for (m = 1; m <= nfre; m += 1) {
          gtemp1 = max((double) (((double) 1.0 - delt5*fld[ij - 1 + kijl*(k - 1 + 
            nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))])), (double) ((double)
             1.0));
          gtemp2 = delt*sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
            nfre_loki_param*(ichnk - 1)))] / gtemp1;
          flhab = abs((double) (gtemp2));
          flhab = min((double) (flhab), (double) (usfm*cofrm4[m - 1]*delt));
          fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)
            ))] = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
            nfre_loki_param*(ichnk - 1)))] + iobnd[ij - 1 + kijl*(ichnk - 1)
            ]*copysign((double) (flhab), (double) (gtemp2));
          fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)
            ))] = max((double) (iodp[ij - 1 + kijl*(ichnk - 1)]*cireduc[ij - 1 + kijl*(k 
            - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))]*fl1[ij - 1 + 
            kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))]), 
            (double) (flm[ij - 1 + kijl*(k - 1 + nang_loki_param*(ichnk - 1))]));
          ssource[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk 
            - 1)))] = ssource[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
            nfre_loki_param*(ichnk - 1)))] + deltm*min((double) (flmax[m - 1] - fl1[ij - 
            1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))]), 
            (double) ((double) 0.0));
          fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)
            ))] = min((double) (fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
            nfre_loki_param*(ichnk - 1)))]), (double) (flmax[m - 1]));
        }
      }
    } else {
      for (k = 1; k <= nang; k += 1) {
        for (m = 1; m <= nfre; m += 1) {
          gtemp1 = max((double) (((double) 1.0 - delt5*fld[ij - 1 + kijl*(k - 1 + 
            nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))])), (double) ((double)
             1.0));
          gtemp2 = delt*sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
            nfre_loki_param*(ichnk - 1)))] / gtemp1;
          flhab = abs((double) (gtemp2));
          flhab = min((double) (flhab), (double) (usfm*cofrm4[m - 1]*delt));
          fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)
            ))] = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
            nfre_loki_param*(ichnk - 1)))] + copysign((double) (flhab), (double) (gtemp2));
          fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)
            ))] = max((double) (cireduc[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
            nfre_loki_param*(ichnk - 1)))]*fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m 
            - 1 + nfre_loki_param*(ichnk - 1)))]), (double) (flm[ij - 1 + kijl*(k - 1 + 
            nang_loki_param*(ichnk - 1))]));
          ssource[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk 
            - 1)))] = ssource[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
            nfre_loki_param*(ichnk - 1)))] + deltm*min((double) (flmax[m - 1] - fl1[ij - 
            1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)))]), 
            (double) ((double) 0.0));
          fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + nfre_loki_param*(ichnk - 1)
            ))] = min((double) (fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(m - 1 + 
            nfre_loki_param*(ichnk - 1)))]), (double) (flmax[m - 1]));
        }
      }
    }

    
    if (lcflx) {
      wnfluxes_c(kijs, kijl, mij, 
         (&rhowgdfth[ + kijl*( + nfre_loki_param*(ichnk - 1))]), cinv, 
         (&ssource[ + kijl*( + nang_loki_param*( + nfre_loki_param*(ichnk - 1)))]), 
        cicover,  (&phiwa[ + kijl*(ichnk - 1)]),  (&emean[ + kijl*(ichnk - 1)]), 
         (&f1mean[ + kijl*(ichnk - 1)]), wswave, wdwave, ufric, aird, nphieps, ntauoc, 
        nswh, nmwp, nemotaux, nemotauy, nemowswave, nemophif, tauxd, tauyd, tauocxd, 
        tauocyd, tauoc, phiocd, phieps, phiaw, true, afcrv, bfcrv, ciblock, cithrsh, 
        costh, egrcrv, epsu10, epsus, fr, g, licerun, lwamrsetci, lwnemocou, 
        lwnemotauoc, nang, nfre, phiepsmax, phiepsmin, sinth, tauocmax, tauocmin, ichnk, 
        nchnk, ij);
    }
    fkmean_c(kijs, kijl, fl1, wavnum,  (&emean[ + kijl*(ichnk - 1)]), 
       (&fmean[ + kijl*(ichnk - 1)]),  (&f1mean[ + kijl*(ichnk - 1)]), 
       (&akmean[ + kijl*(ichnk - 1)]),  (&xkmean[ + kijl*(ichnk - 1)]), delth, dfim, 
      dfimfr, dfimofr, epsmin, fr, frtail, g, nang, nfre, wetail, wp1tail, zpi, ichnk, 
      nchnk, ij);
    femeanws_c(kijs, kijl, fl1, xllws,  (&fmeanws[ + kijl*(ichnk - 1)]), 
       (&emeanws[ + kijl*(ichnk - 1)]), delth, dfim, dfimofr, epsmin, fr, frtail, nang, 
      nfre, wetail, ichnk, nchnk, ij);
    
    imphftail_c(kijs, kijl, mij,  (&flm[ + kijl*( + nang_loki_param*(ichnk - 1))]), 
      wavnum, xk2cg, fl1, nang, nfre, ichnk, nchnk, ij);

    if (lwflux) {
      if (emeanws[ij - 1 + kijl*(ichnk - 1)] < wsemean_min) {
        wsemean[ij - 1 + kijl*(ichnk - 1)] = wsemean_min;
        wsfmean[ij - 1 + kijl*(ichnk - 1)] = (double) 2.*fr[nfre - 1];
      } else {
        wsemean[ij - 1 + kijl*(ichnk - 1)] = emeanws[ij - 1 + kijl*(ichnk - 1)];
        wsfmean[ij - 1 + kijl*(ichnk - 1)] = fmeanws[ij - 1 + kijl*(ichnk - 1)];
      }
    }

    if (licerun && lmaskice) {
      setice_c(kijs, kijl, fl1, cicover, 
         (&coswdif[ + kijl*( + nang_loki_param*(ichnk - 1))]), cithrsh, epsmin, flmin, 
        nang, nfre, ichnk, nchnk, ij);
    }
    stokestrn_c(kijs, kijl, fl1, wavnum, stokfac, depth, wswave, wdwave, cicover, 
      cithick, ustokes, vstokes, strnms, nemoustokes, nemovstokes, nemostrn, cithrsh, 
      costh, delth, dfim, dfim_sim, flmin, fr, g, licerun, lwamrsetci, lwcou, lwnemocou, 
      lwnemocousend, lwnemocoustk, lwnemocoustrn, nang, nfre, nfre_odd, rowater, sinth, 
      zpi, ichnk, nchnk, ij);
    // END of Loki inserted loop ICHNK
  }
}
