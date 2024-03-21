#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "ns_gc_c.h"

__device__ void omegagc_c(double ust, int *ns, double *xks, double *oms, int nwav_gc, 
  const double * omega_gc, double sqrtgosurft, const double * xkm_gc, 
  const double * xk_gc, double xlogkratiom1_gc);
