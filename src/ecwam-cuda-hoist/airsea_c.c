#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "airsea_c.h"
#include "z0wave_c.h"
#include "taut_z0_c.h"


__device__ void airsea_c(int kijs, int kijl, const double * __restrict__ halp, 
  double * __restrict__ u10, const double * __restrict__ u10dir, 
  const double * __restrict__ tauw, const double * __restrict__ tauwdir, 
  const double * __restrict__ rnfac, double * __restrict__ us, double * __restrict__ z0, 
  double * __restrict__ z0b, double * __restrict__ chrnck, int icode_wnd, int iusfg, 
  double acd, double alpha, double alphamax, double alphamin, double ang_gc_a, 
  double ang_gc_b, double ang_gc_c, double bcd, double betamaxoxkappa2, double bmaxokap, 
  const double * __restrict__ c2osqrtvg_gc, double cdmax, double chnkmin_u, 
  const double * __restrict__ cm_gc, const double * __restrict__ delkcc_gc_ns, 
  const double * __restrict__ delkcc_omxkm3_gc, double eps1, double epsmin, 
  double epsus, double g, double gm1, int llcapchnk, int llgcbz0, int llnormagam, 
  int nwav_gc, const double * __restrict__ om3gmkm_gc, 
  const double * __restrict__ omxkm3_gc, double rn1_rn, double rnu, double rnum, 
  double sqrtgosurft, double wspmin, double xkappa, 
  const double * __restrict__ xkmsqrtvgoc2_gc, const double * __restrict__ xkm_gc, 
  const double * __restrict__ xk_gc, double xlogkratiom1_gc, double xnlev, double zalp, 
  int ichnk, int nchnk, int ij, double * __restrict__ taut_z0_alphaog, 
  double * __restrict__ taut_z0_xmin, double * __restrict__ taut_z0_w1, 
  double * __restrict__ taut_z0_tauwact, double * __restrict__ taut_z0_tauweff, 
  double * __restrict__ taut_z0_ang_gc, double * __restrict__ taut_z0_tauunr, 
  int * __restrict__ taut_z0_llcosdiff, double * __restrict__ stress_gc_gam_w, 
  double * __restrict__ z0wave_alphaog) {
  
  // needed for Loki
  
  
  // ----------------------------------------------------------------------



  
  int i;
  int j;
  
  double xi;
  double xj;
  double deli1;
  double deli2;
  double delj1;
  double delj2;
  double ust2;
  double arg;
  double sqrtcdm1;
  double xkappad;
  double xloglev;
  double xlev;
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  //*    2. DETERMINE TOTAL STRESS (if needed)
  //        ----------------------------------
  
  if (icode_wnd == 3) {
    
    taut_z0_c(kijs, kijl, iusfg, halp, u10, u10dir, tauw, tauwdir, rnfac, us, z0, z0b, 
      chrnck, acd, alpha, alphamax, alphamin, ang_gc_a, ang_gc_b, ang_gc_c, bcd, 
      betamaxoxkappa2, bmaxokap, c2osqrtvg_gc, cdmax, chnkmin_u, cm_gc, delkcc_gc_ns, 
      delkcc_omxkm3_gc, eps1, epsmin, epsus, g, gm1, llcapchnk, llgcbz0, llnormagam, 
      nwav_gc, om3gmkm_gc, omxkm3_gc, rn1_rn, rnu, rnum, sqrtgosurft, xkappa, 
      xkmsqrtvgoc2_gc, xkm_gc, xk_gc, xlogkratiom1_gc, xnlev, zalp, ichnk, nchnk, ij, 
      taut_z0_alphaog, taut_z0_xmin, taut_z0_w1, taut_z0_tauwact, taut_z0_tauweff, 
      taut_z0_ang_gc, taut_z0_tauunr, taut_z0_llcosdiff, stress_gc_gam_w);
    
  } else if (icode_wnd == 1 || icode_wnd == 2) {
    
    //*    3. DETERMINE ROUGHNESS LENGTH (if needed).
    //        ---------------------------
    
    z0wave_c(kijs, kijl, us, tauw, u10, z0, z0b, chrnck, alpha, alphamin, chnkmin_u, 
      eps1, g, gm1, llcapchnk, ichnk, nchnk, ij, z0wave_alphaog);
    
    //*    3. DETERMINE U10 (if needed).
    //        ---------------------------
    
    xkappad = (double) 1.0 / xkappa;
    xloglev = log(xnlev);
    
    
    u10[ij - 1 + kijl*(ichnk - 1)] = xkappad*us[ij - 1 + kijl*(ichnk - 1)]*(xloglev - 
      log(z0[ij - 1 + kijl*(ichnk - 1)]));
    u10[ij - 1 + kijl*(ichnk - 1)] = 
      fmax((double) (u10[ij - 1 + kijl*(ichnk - 1)]), (double) (wspmin));
    
    
  }
  
  
}
