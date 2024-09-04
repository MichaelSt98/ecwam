#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "omegagc_c.h"
#include "ns_gc_c.h"


__device__ void omegagc_c(int kijs, int kijl, const double * __restrict__ ust, 
  int * __restrict__ ns, double * __restrict__ xks, double * __restrict__ oms, 
  int nwav_gc, const double * __restrict__ omega_gc, double sqrtgosurft, 
  const double * __restrict__ xkm_gc, const double * __restrict__ xk_gc, 
  double xlogkratiom1_gc, int ichnk, int nchnk, int ij) {
  
  
  
  //----------------------------------------------------------------------
  
  
  
  
  double zhook_handle;

  // ----------------------------------------------------------------------
  
  
  
  ns[ij - 1] = ns_gc_c(ust[ij - 1 + kijl*(ichnk - 1)], nwav_gc, sqrtgosurft, xkm_gc, 
    xlogkratiom1_gc);
  xks[ij - 1] = xk_gc[ns[ij - 1] - 1];
  oms[ij - 1] = omega_gc[ns[ij - 1] - 1];
  
  
  
}
