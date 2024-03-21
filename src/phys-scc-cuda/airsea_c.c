#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "airsea_c.h"
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
  int nchnk, int ij) {
  
  


  
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
  if (icode_wnd == 3) {
    
    taut_z0_c(kijs, kijl, iusfg, halp, u10, u10dir, tauw, tauwdir, rnfac, us, z0, z0b, 
      chrnck, acd, alpha, alphamax, alphamin, ang_gc_a, ang_gc_b, ang_gc_c, bcd, 
      betamaxoxkappa2, bmaxokap, c2osqrtvg_gc, cdmax, chnkmin_u, cm_gc, delkcc_gc_ns, 
      delkcc_omxkm3_gc, eps1, epsmin, epsus, g, gm1, llcapchnk, llgcbz0, llnormagam, 
      nwav_gc, om3gmkm_gc, omxkm3_gc, rn1_rn, rnu, rnum, sqrtgosurft, xkappa, 
      xkmsqrtvgoc2_gc, xkm_gc, xk_gc, xlogkratiom1_gc, xnlev, zalp, ichnk, nchnk, ij);
    
  } else if (icode_wnd == 1 || icode_wnd == 2) {
    z0wave_c(kijs, kijl, us, tauw, u10, z0, z0b, chrnck, alpha, alphamin, chnkmin_u, 
      eps1, g, gm1, llcapchnk, ichnk, nchnk, ij);
    xkappad = (double) 1.0 / xkappa;
    xloglev = log(xnlev);
    

    u10[ij - 1 + kijl*(ichnk - 1)] = xkappad*us[ij - 1 + kijl*(ichnk - 1)]*(xloglev - 
      log(z0[ij - 1 + kijl*(ichnk - 1)]));
    u10[ij - 1 + kijl*(ichnk - 1)] = 
      max((double) (u10[ij - 1 + kijl*(ichnk - 1)]), (double) (wspmin));

    
  }
  
  
}
