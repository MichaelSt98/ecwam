#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "z0wave_c.h"
#include "taut_z0_c.h"

__device__ void airsea_c(int kijs, int kijl, const double * halp, double * u10, 
  const double * u10dir, const double * tauw, const double * tauwdir, 
  const double * rnfac, double * us, double * z0, double * z0b, double * chrnck, 
  int icode_wnd, int iusfg, double acd, double alpha, double alphamax, double alphamin, 
  double ang_gc_a, double ang_gc_b, double ang_gc_c, double bcd, double betamaxoxkappa2, 
  double bmaxokap, const double * c2osqrtvg_gc, double cdmax, double chnkmin_u, 
  const double * cm_gc, const double * delkcc_gc_ns, const double * delkcc_omxkm3_gc, 
  double eps1, double epsmin, double epsus, double g, double gm1, int llcapchnk, 
  int llgcbz0, int llnormagam, int nwav_gc, const double * om3gmkm_gc, 
  const double * omxkm3_gc, double rn1_rn, double rnu, double rnum, double sqrtgosurft, 
  double wspmin, double xkappa, const double * xkmsqrtvgoc2_gc, const double * xkm_gc, 
  const double * xk_gc, double xlogkratiom1_gc, double xnlev, double zalp, int ichnk, 
  int nchnk, int ij);
