#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "chnkmin_c.h"

__device__ void z0wave_c(int kijs, int kijl, const double * us, const double * tauw, 
  const double * utop, double * z0, double * z0b, double * chrnck, double alpha, 
  double alphamin, double chnkmin_u, double eps1, double g, double gm1, int llcapchnk, 
  int ichnk, int nchnk, int ij);
