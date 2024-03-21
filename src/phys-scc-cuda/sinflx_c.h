#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
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
  double * stresso_usdirp, double * stresso_ust);
