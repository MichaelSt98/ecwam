#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "ns_gc_c.h"

__device__ int ns_gc_c(double ustar, int nwav_gc, double sqrtgosurft, 
  const double * xkm_gc, double xlogkratiom1_gc) {
  
  
  int ns_gc;
  
  double y, xks;
  // <Pragma:: acc routine seq>
  xks = sqrtgosurft / ((double) 1.48 + (double) 2.05*ustar);
  
  ns_gc = min((double) ((int) (log(xks*xkm_gc[1 - 1])*xlogkratiom1_gc) + 1), (double) 
    (nwav_gc - 1));

  return ns_gc;
}
