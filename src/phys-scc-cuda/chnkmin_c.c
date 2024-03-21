#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "chnkmin_c.h"

__device__ double chnkmin_c(double u10, double alpha, double alphamin, double chnkmin_u) {
  
  
  double chnkmin;
  // <Pragma:: acc routine seq>
  
  chnkmin = 
    alphamin + (alpha - alphamin)*(double) 0.5*((double) 1.0 - tanh(u10 - chnkmin_u));
  
  return chnkmin;  
}
