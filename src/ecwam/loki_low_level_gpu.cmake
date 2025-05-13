list( APPEND ecwam_srcs_old
  readmdlconf.F90
  mfr.F90
  propconnect.F90
  wsigstar.F90
  omegagc.F90
  sinput_jan.F90
  wgrib_edition.F90
  meansqs_gc.F90
  ctcor.F90
  wamintgr.F90
  abort1.F90
  adjust.F90
  aki.F90
  alphap_tail.F90
  bouinpt.F90
  buildstress.F90
  cal_second_order_spec.F90
  cdustarz0.F90
  check.F90
  checkcfl.F90
  chkoops.F90
  cigetdeac.F90
  cireduce.F90
  ciwaf.F90
  ctuw.F90
  ctuwdrv.F90
  ctuwini.F90
  ctuwupdt.F90
  current2wam.F90
  difdate.F90
  dominant_period.F90
  depthprpt.F90
  expand_string.F90
  femean.F90
  file_transfer.F90
  findb.F90
  fldinter.F90
  fndprt.F90
  get_preset_wgrib_template.F90
  getcurr.F90
  getfrstwnd.F90
  getspec.F90
  getstress.F90
  getwnd.F90
  gradi.F90
  grib2wgrid.F90
  grstname.F90
  gsfile_new.F90
  h_max.F90
  headbc.F90
  incdate.F90
  inisnonlin.F90
  init_fieldg.F90
  init_sdiss_ardh.F90
  init_x0tauhf.F90
  initdpthflds.F90
  initgc.F90
  initialint.F90
  initmdl.F90
  initnemocpl.F90
  iniwcst.F90
  intpol.F90
  intspec.F90
  inwgrib.F90
  iwam_get_unit.F90
  jafu.F90
  jonswap.F90
  kerkei.F90
  kgribsize.F90
  kurtosis.F90
  kzeone.F90
  makegrid.F90
  mblock.F90
  mbounc.F90
  mbounf.F90
  mboxb.F90
  mchunk.F90
  mcout.F90
  means.F90
  meansqs.F90
  meansqs_lf.F90
  mfredir.F90
  mgrid.F90
  micep.F90
  mintf.F90
  mnintw.F90
  mpabort.F90
  mpbcastgrid.F90
  mpbcastintfld.F90
  mpclose_unit.F90
  mpcrtbl.F90
  mpdecomp.F90
  mpdistribfl.F90
  mpdistribscfld.F90
  mpexchng.F90
  mpfldtoifs.F90
  mpgatherbc.F90
  mpgatherfl.F90
  mpgatherscfld.F90
  mpminmaxavg.F90
  mpuserin.F90
  mstart.F90
  mswell.F90
  # mtabs.F90
  mubuf.F90
  mwp1.F90
  mwp2.F90
  newwind.F90
  nlweigt.F90
  notim.F90
  out_onegrdpt.F90
  out_onegrdpt_sp.F90
  outbc.F90
  outbeta.F90
  # outblock.F90
  outbs.F90
  outcom.F90
  outgrid.F90
  outint.F90
  outmdldcp.F90
  outnam.F90
  outpp.F90
  outsetwmask.F90
  outspec.F90
  outstep0.F90
  outwint.F90
  outwnorm.F90
  outwpsp.F90
  outwspec.F90
  packi.F90
  packr.F90
  parkind_wave.F90
  parmean.F90
  peak.F90
  peak_freq.F90
  peakfri.F90
  preset_wgrib_template.F90
  prewind.F90
  proenvhalo.F90
  propag_wam.F90
  propags.F90
  propags1.F90
  propags2.F90
  propdot.F90
  readbou.F90
  readfl.F90
  readpre.F90
  readsta.F90
  readstress.F90
  readwgrib.F90
  readwind.F90
  recvnemofields.F90
  rotspec.F90
  runwam.F90
  savspec.F90
  savstress.F90
  scosfl.F90
  se10mean.F90
  sebtmean.F90
  second_order_lib.F90
  secondhh.F90
  secondhh_gen.F90
  secspom.F90
  semean.F90
  sep3tr.F90
  sepwisw.F90
  set_wflags.F90
  setmarstype.F90
  setwavphys.F90
  skewness.F90
  spectra.F90
  spr.F90
  # stack_mod.F90
  stat_nl.F90
  sthq.F90
  strspec.F90
  tables_2nd.F90
  tabu_swellft.F90
  topoar.F90
  transf_bfi.F90
  transf_r.F90
  uibou.F90
  uiprep.F90
  unsetice.F90
  updnemofields.F90
  updnemostress.F90
  userin.F90
  vmin.F90
  vmin_d.F90
  vplus.F90
  vplus_d.F90
  w_maxh.F90
  w_mode_st.F90
  # wam_init_gpu_mod.F90
  wam_multio_mod.F90
  wam_nproma.F90
  wam_sorti.F90
  wam_sortini.F90
  wam_u2l1cr.F90
  wam_user_clock.F90
  wamadswstar.F90
  wamcur.F90
  wamodel.F90
  wamwnd.F90
  wavemdl.F90
  wdfluxes.F90
  wdirspread.F90
  weflux.F90
  wgrib2fdb.F90
  wgribencode.F90
  wgribencode_model.F90
  wgribenout.F90
  wgribout.F90
  wposnam.F90
  writefl.F90
  writestress.F90
  writsta.F90
  wsmfen.F90
  wstream_strg.F90
  wvalloc.F90
  wvdealloc.F90
  wvfricvelo.F90
  wvwamdecomp.F90
  wvwaminit.F90
  wvwaminit1.F90
  yowabort.F90
  yowassi.F90
  yowcard.F90
  yowcinp.F90
  yowcoer.F90
  yowconst_2nd.F90
  yowcpbo.F90
  yowcurg.F90
  yowcurr.F90
  yowdes.F90
  # yowdrvtype.F90
  yowfpbo.F90
  yowgrib.F90
  yowgrib_handles.F90
  yowgribhd.F90
  yowgrid.F90
  yowgstats.F90
  yowintp.F90
  yowjons.F90
  yowmap.F90
  yowmean.F90
  yowmespas.F90
  yowmpp.F90
  yownemoflds.F90
  yownemoio.F90
  yowprproc.F90
  yowrefd.F90
  yowshal.F90
  yowspec.F90
  yowsphere.F90
  yowtemp.F90
  yowtest.F90
  yowtext.F90
  yowtrains.F90
  yowubuf.F90
  yowunit.F90
  yowwami.F90
  # new
  sdepthlim.F90
  outblock.F90
  setice.F90
  fkmean.F90
  sinflx.F90
  sdissip.F90
  snonlin.F90
  wnfluxes.F90
  femeanws.F90
  stokestrn.F90
  peak_ang.F90
  halphap.F90
  airsea.F90
  sinput.F90
  frcutindex.F90
  stresso.F90
  sdissip_jan.F90
  sdissip_ard.F90
  transf_snl.F90
  transf.F90
  stokesdrift.F90
  cimsstrn.F90
  chnkmin.F90
  ns_gc.F90
  z0wave.F90
  taut_z0.F90
  sinput_ard.F90
  tau_phi_hf.F90
  aki_ice.F90
  stress_gc.F90
)

list(APPEND global_var_mods
     yowaltas.F90
     yowcoup.F90
     yowcout.F90
     yowfred.F90
     yowice.F90
     yowindn.F90
     yowparam.F90
     yowpcons.F90
     yowphys.F90
     yowstat.F90
     yowtabl.F90
     yowwind.F90
     yowwndg.F90
     # yownemoio.F90
)
list(APPEND phys_srcs
     airsea.F90
     aki_ice.F90
     chnkmin.F90
     # cimsstrn.F90
     # ciwabr.F90
     femeanws.F90
     fkmean.F90
     frcutindex.F90
     # halphap.F90
     imphftail.F90
     implsch.F90
     ns_gc.F90
     peak_ang.F90
     sbottom.F90
     sdepthlim.F90
     sdissip.F90
     sdissip_ard.F90
     # sdissip_jan.F90
     sdiwbk.F90
     setice.F90
     sinflx.F90
     sinput.F90
     sinput_ard.F90
     # new
     # sinput_jan.F90
     snonlin.F90
     stokesdrift.F90
     stokestrn.F90
     stress_gc.F90
     stresso.F90
     tau_phi_hf.F90
     taut_z0.F90
     transf.F90
     transf_snl.F90
     wnfluxes.F90
     # z0wave.F90
     # new ...
     omegagc.F90
     # meansqs_lf.F90
     femean.F90
     semean.F90
     wsigstar.F90
     # # new
     # ciwaf.F90
     # outblock.F90
     # # newnew
     # intpol.F90
     # cal_second_order_spec.F90
     # femean.F90
     # dominant_period.F90
     # kurtosis.F90
     # sepwisw.F90
     # sthq.F90
     # outbeta.F90
     # meansqs.F90
     # mwp1.F90
     # mwp2.F90
     # wdirspread.F90
     # se10mean.F90
     # weflux.F90
     # sebtmean.F90
     # w_maxh.F90
     # ctcor.F90
     # outsetwmask.F90
)

list(APPEND inlined_srcs
  # chnkmin.F90
  # ns_gc.F90
  # stress_gc.F90
  # transf_snl.F90
  # transf.F90
  aki_ice.F90
  aki.F90
  peakfri.F90
  femeanws.F90
  frcutindex.F90
  omegagc.F90
  tau_phi_hf.F90
  stresso.F90
  wsigstar.F90
  sinput.F90
  sinput_ard.F90
  taut_z0.F90
  airsea.F90
  femean.F90
  wnfluxes.F90
  sdiwbk.F90
  sbottom.F90
  fkmean.F90
  imphftail.F90
  setice.F90
  stokestrn.F90
  stokesdrift.F90
  semean.F90
  sdepthlim.F90
  sinflx.F90
  sdissip_ard.F90
  sdissip.F90
  peak_ang.F90
)

  # foreach(src ${phys_srcs}) # wamintgr_cuda_mod.F90 ${global_var_mods})
  #    string(REPLACE ".F90" "" fnc ${src})
  #    string(CONCAT fnm "${CMAKE_CURRENT_SOURCE_DIR}/" ${fnc} "_c.c")
  #    list(APPEND wam_scc_cuda_srcs ${fnm})
  # endforeach()

  # list(REMOVE_ITEM ecwam_srcs wamintgr.F90)
  list(REMOVE_ITEM ecwam_srcs_old wamintgr.F90)

  foreach(src ${phys_srcs}) # wamintgr_cuda_mod.F90 ${global_var_mods})
     if (src IN_LIST inlined_srcs)
        message("skipping ${src} since inlined!")
     else()
     string(REPLACE ".F90" "" fnc ${src})
     string(CONCAT fnm "${CMAKE_CURRENT_BINARY_DIR}/loki-hip-hoist/" ${fnc} "_c.c")
     # string(CONCAT fnm "../cuda-ecwam-3/" ${fnc} "_c.c")
     ## string(CONCAT fnm "../cuda-ecwam-1-small-testcase/" ${fnc} "_c.c")
     ## string(CONCAT fnm "../phys-scc-cuda/" ${fnc} "_c.c")
     list(APPEND loki_wam_scc_cuda_srcs ${fnm})
     endif()
  endforeach()

  # foreach(src ${global_var_mods} wamintgr_loki_gpu.F90) # wamintgr_cuda_mod.F90 ${global_var_mods})
  # foreach(src wamintgr_loki_gpu.F90 cireduce_loki_gpu.F90 outbs_loki_gpu.F90 ${global_var_mods}) # wamintgr_cuda_mod.F90 ${global_var_mods})
  foreach(src wamintgr_loki_gpu.F90 ${global_var_mods})
     string(REPLACE ".F90" "" fnc ${src})
     string(CONCAT fnm "${CMAKE_CURRENT_BINARY_DIR}/loki-hip-hoist/" ${fnc} ".hip_hoist.F90")
     # string(CONCAT fnm "../cuda-ecwam-3/" ${fnc} ".cuda_hoist.F90")
     ## string(CONCAT fnm "../cuda-ecwam-1-small-testcase/" ${fnc} ".cuda_hoist.F90")
     ## string(CONCAT fnm "../phys-scc-cuda/" ${fnc} ".c_hoist.F90")
     list(APPEND loki_wam_scc_cuda_srcs_2 ${fnm})
  endforeach()
