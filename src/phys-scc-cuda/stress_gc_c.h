#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "ns_gc_c.h"

__device__ double stress_gc_c(double ang_gc, double ustar, double z0, double z0min, 
  double halp, double rnfac, double betamaxoxkappa2, double bmaxokap, 
  const double * c2osqrtvg_gc, const double * cm_gc, const double * delkcc_gc_ns, 
  const double * delkcc_omxkm3_gc, double epsus, int llnormagam, int nwav_gc, 
  const double * om3gmkm_gc, const double * omxkm3_gc, double rn1_rn, 
  double sqrtgosurft, double xkappa, const double * xkmsqrtvgoc2_gc, 
  const double * xkm_gc, const double * xk_gc, double xlogkratiom1_gc, double zalp);
