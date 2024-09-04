#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "ns_gc_c.h"


__device__ int ns_gc_c(double ustar, int nwav_gc, double sqrtgosurft, 
  const double * __restrict__ xkm_gc, double xlogkratiom1_gc) {
  
  
  // ----------------------------------------------------------------------
  
  
  int ns_gc;
  
  double y;
  double xks;
  // <Pragma:: acc routine seq>
  
  // ----------------------------------------------------------------------
  
  
  //!!Y = 1.0_JWRB/(1.48_JWRB+2.05_JWRB*UST)
  //!!Y = (1.0_JWRB + UST**2)/(1.0_JWRB+10.0_JWRB*UST**2)
  
  xks = sqrtgosurft / ((double) 1.48 + (double) 2.05*ustar);
  
  ns_gc = fmin((double) ((int) (log(fmax((double) (xks*xkm_gc[1 - 1]), (double) ((double)
     1.0)))*xlogkratiom1_gc) + 1), (double) (nwav_gc - 1));
  
  
  return ns_gc;
}
