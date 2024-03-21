#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
extern "C" {
void implsch_c_launch(int kijs, int kijl, double * fl1, const double * wavnum, 
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
  // here should be the launcher ....
  griddim = dim3(nchnk, 1, 1);
  blockdim = dim3(nproma_wam, 1, 1);
  implsch_c<<<griddim,blockdim>>>(kijs, kijl, fl1, wavnum, cgroup, ciwa, cinv, xk2cg, 
    stokfac, emaxdpt, indep, depth, iobnd, iodp, aird, wdwave, cicover, wswave, wstar, 
    ufric, tauw, tauwdir, z0m, z0b, chrnck, cithick, nemoustokes, nemovstokes, nemostrn, 
    nphieps, ntauoc, nswh, nmwp, nemotaux, nemotauy, nemowswave, nemophif, wsemean, 
    wsfmean, ustokes, vstokes, strnms, tauxd, tauyd, tauocxd, tauocyd, tauoc, phiocd, 
    phieps, phiaw, mij, xllws, abmax, abmin, acd, acdlin, af11, afcrv, alpha, alphamax, 
    alphamin, alphapmax, ang_gc_a, ang_gc_b, ang_gc_c, bathymax, bcd, bcdlin, 
    betamaxoxkappa2, bfcrv, bmaxokap, c2osqrtvg_gc, cdicwa, cdis, cdisvis, cdmax, 
    chnkmin_u, ciblock, cithrsh, cithrsh_tail, cm_gc, cofrm4, costh, cumulw, dal1, dal2, 
    delkcc_gc_ns, delkcc_omxkm3_gc, delta_sdis, delth, dfim, dfimfr, dfimfr2, dfimofr, 
    dfim_sim, dkmax, dthrn_a, dthrn_u, egrcrv, eps1, epsmin, epsu10, epsus, fklam, 
    fklam1, fklap, fklap1, flmax, flmin, flogsprdm1, fr, fr5, fratio, fric, frtail, g, 
    gamnconst, gm1, iab, icode, icode_cpl, idamping, idelt, ikm, ikm1, ikp, ikp1, 
    indicessat, inlcoef, iphys, ipsat, isnonlin, jtot_tauhf, k11w, k1w, k21w, k2w, kfrh, 
    lbiwbk, lciwabr, licerun, llcapchnk, llgcbz0, llnormagam, llunstr, lmaskice, 
    lwamrsetci, lwcou, lwflux, lwfluxout, lwnemocou, lwnemocousend, lwnemocoustk, 
    lwnemocoustrn, lwnemotauoc, lwvflx_snl, mfrstlw, miche, mlsthg, nang, ndepth, 
    ndikcumul, nfre, nfre_odd, nfre_red, nsdsnth, nwav_gc, om3gmkm_gc, omega_gc, 
    omxkm3_gc, phiepsmax, phiepsmin, rhowg_dfim, rn1_rn, rnlcoef, rnu, rnum, rowater, 
    rowaterm1, satweights, sdsbr, sinth, sqrtgosurft, ssdsc2, ssdsc3, ssdsc4, ssdsc5, 
    ssdsc6, swellf, swellf2, swellf3, swellf4, swellf5, swellf6, swellf7, swellf7m1, 
    swellft, tailfactor, tailfactor_pm, tauocmax, tauocmin, tauwshelter, th, wetail, 
    wp1tail, wp2tail, wsemean_min, wspmin, wtauhf, x0tauhf, xkappa, xkdmin, 
    xkmsqrtvgoc2_gc, xkm_gc, xk_gc, xlogkratiom1_gc, xnlev, z0rat, z0tubmax, zalp, zpi, 
    zpi4gm1, zpi4gm2, zpifr, ichnk_start, ichnk_end, ichnk_step, nchnk, nproma_wam, 
    raorw, emean, fmean, halp, emeanws, fmeanws, f1mean, akmean, xkmean, phiwa, flm, 
    coswdif, sinwdif2, rhowgdfth, fld, sl, spos, cireduc, ssource, sinflx_rnfac, 
    sinflx_tmp_em, stresso_xstress, stresso_ystress, stresso_tauhf, stresso_phihf, 
    stresso_usdirp, stresso_ust, snonlin_enh, snonlin_xnu, snonlin_sig_th);
  cudaDeviceSynchronize();
}
}
