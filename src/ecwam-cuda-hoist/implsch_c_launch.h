

extern "C" {

void implsch_c_launch(int kijs, int kijl, double * __restrict__ fl1, 
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
  // here should be the launcher ....
  griddim = dim3(nchnk, 1, 1);
  blockdim = dim3(nproma_wam, 1, 1);
  implsch_c<<<griddim,blockdim>>>(kijs, kijl, fl1, wavnum, cgroup, ciwa, cinv, xk2cg, 
    stokfac, emaxdpt, depth, iobnd, iodp, aird, wdwave, cicover, wswave, wstar, ustra, 
    vstra, ufric, tauw, tauwdir, z0m, z0b, chrnck, cithick, nemoustokes, nemovstokes, 
    nemostrn, nphieps, ntauoc, nswh, nmwp, nemotaux, nemotauy, nemowswave, nemophif, 
    wsemean, wsfmean, ustokes, vstokes, strnms, tauxd, tauyd, tauocxd, tauocyd, tauoc, 
    phiocd, phieps, phiaw, mij, xllws, abmax, abmin, acd, acdlin, af11, afcrv, alpha, 
    alphamax, alphamin, alphapmax, ang_gc_a, ang_gc_b, ang_gc_c, bathymax, bcd, bcdlin, 
    betamaxoxkappa2, bfcrv, bmaxokap, brkpbcoef, c2osqrtvg_gc, cdicwa, cdis, cdisvis, 
    cdmax, chnkmin_u, ciblock, cithrsh, cithrsh_tail, cm_gc, cofrm4, costh, dal1, dal2, 
    delkcc_gc_ns, delkcc_omxkm3_gc, delta_sdis, delth, dfim, dfimfr, dfimfr2, dfimofr, 
    dfim_sim, dkmax, dthrn_a, dthrn_u, egrcrv, eps1, epsmin, epsu10, epsus, fklam, 
    fklam1, fklap, fklap1, flmax, flmin, flogsprdm1, fr, fr5, fratio, fric, frtail, g, 
    gamnconst, gm1, iab, icode, icode_cpl, idamping, idelt, ikm, ikm1, ikp, ikp1, 
    indicessat, inlcoef, iphys, ipsat, isnonlin, jtot_tauhf, k11w, k1w, k21w, k2w, kfrh, 
    lbiwbk, lciwabr, licerun, llcapchnk, llgcbz0, llnormagam, llunstr, lmaskice, 
    lwamrsetci, lwcou, lwcouast, lwflux, lwfluxout, lwnemocou, lwnemocousend, 
    lwnemocoustk, lwnemocoustrn, lwnemotauoc, lwvflx_snl, mfrstlw, miche, mlsthg, nang, 
    nfre, nfre_odd, nfre_red, nsdsnth, nwav_gc, om3gmkm_gc, omega_gc, omxkm3_gc, 
    phiepsmax, phiepsmin, rhowg_dfim, rn1_rn, rnlcoef, rnu, rnum, rowater, rowaterm1, 
    satweights, sdsbr, sinth, sqrtgosurft, ssdsbrf1, ssdsc2, ssdsc3, ssdsc4, ssdsc5, 
    ssdsc6, swellf, swellf2, swellf3, swellf4, swellf5, swellf6, swellf7, swellf7m1, 
    swellft, tailfactor, tailfactor_pm, tauocmax, tauocmin, tauwshelter, th, wetail, 
    wp1tail, wp2tail, wsemean_min, wspmin, wtauhf, x0tauhf, xkappa, xkdmin, 
    xkmsqrtvgoc2_gc, xkm_gc, xk_gc, xlogkratiom1_gc, xnlev, z0rat, z0tubmax, zalp, zpi, 
    zpi4gm1, zpi4gm2, zpifr, nchnk, nproma_wam, delfl, raorw, emean, fmean, halp, 
    emeanws, fmeanws, usfm, f1mean, akmean, xkmean, phiwa, flm, coswdif, sinwdif2, temp, 
    rhowgdfth, fld, sl, spos, cireduc, ssource, femeanws_temp2, femeanws_em_loc, 
    fkmean_tempa, fkmean_tempx, fkmean_temp2, sbottom_sbo, imphftail_temp1, 
    imphftail_temp2, sdepthlim_em, semean_temp, sdissip_ard_facturb, sdissip_ard_facsat, 
    sdissip_ard_facwtrb, sdissip_ard_temp1, sdissip_ard_bth0, sdissip_ard_c_, 
    sdissip_ard_c_c, sdissip_ard_dsip, sdissip_ard_trpz_dsip, sdissip_ard_bth, 
    sdissip_ard_temp2, sdissip_ard_d, sdissip_ard_scumul, sdissip_ard_renewalfreq, 
    sdissip_ard_wcumul, sdissip_jan_temp1, sdissip_jan_sds, sdissip_jan_x, 
    sdissip_jan_xk2, sdiwbk_sds, setice_cireduc, setice_temp, setice_icefree, 
    sinflx_rnfac, sinflx_tmp_em, taut_z0_alphaog, taut_z0_xmin, taut_z0_w1, 
    taut_z0_tauwact, taut_z0_tauweff, taut_z0_ang_gc, taut_z0_tauunr, taut_z0_llcosdiff, 
    stress_gc_gam_w, z0wave_alphaog, halphap_alphap, halphap_xmss, halphap_em, 
    halphap_fm, halphap_f1d, halphap_wd, halphap_flwd, femean_temp2, meansqs_lf_fd, 
    meansqs_lf_temp1, meansqs_lf_temp2, sinput_ard_constf, sinput_ard_const11, 
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
    sinput_jan_gam0, sinput_jan_lz, stresso_xstress, stresso_ystress, stresso_tauhf, 
    stresso_phihf, stresso_cmrhowgdfth, stresso_us2, stresso_taux, stresso_tauy, 
    stresso_taupx, stresso_taupy, stresso_usdirp, stresso_ust, stresso_sumt, 
    stresso_sumx, stresso_sumy, tau_phi_hf_ns, tau_phi_hf_xks, tau_phi_hf_oms, 
    tau_phi_hf_sqrtz0og, tau_phi_hf_zsup, tau_phi_hf_zinf, tau_phi_hf_delz, 
    tau_phi_hf_taul, tau_phi_hf_xloggz0, tau_phi_hf_sqrtgz0, tau_phi_hf_ustph, 
    tau_phi_hf_const1, tau_phi_hf_const2, tau_phi_hf_consttau, tau_phi_hf_constphi, 
    tau_phi_hf_f1dcos2, tau_phi_hf_f1dcos3, tau_phi_hf_f1d, tau_phi_hf_f1dsin2, 
    snonlin_enh, snonlin_xnu, snonlin_sig_th, snonlin_ftemp, snonlin_ad, snonlin_delad, 
    snonlin_delap, snonlin_delam, snonlin_enhfr, peak_ang_mmax, peak_ang_sum0, 
    peak_ang_sum1, peak_ang_sum2, peak_ang_xmax, peak_ang_temp, peak_ang_thmean, 
    peak_ang_sum_s, peak_ang_sum_c, cimsstrn_xki, cimsstrn_e, cimsstrn_sume, 
    stokesdrift_stfac, wnfluxes_xstress, wnfluxes_ystress, wnfluxes_ustar, 
    wnfluxes_philf, wnfluxes_ooval, wnfluxes_em_oc, wnfluxes_f1_oc, 
    wnfluxes_cmrhowgdfth, wnfluxes_sumt, wnfluxes_sumx, wnfluxes_sumy);
  cudaDeviceSynchronize();
}


} // extern
