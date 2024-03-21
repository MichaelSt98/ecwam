#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ int ns_gc_c(double ustar, int nwav_gc, double sqrtgosurft, 
  const double * xkm_gc, double xlogkratiom1_gc);
