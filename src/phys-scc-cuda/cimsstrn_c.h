#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "aki_ice_c.h"

__device__ void cimsstrn_c(int kijs, int kijl, const double * fl1, 
  const double * wavnum, const double * depth, const double * cithick, double * strn, 
  double delth, const double * dfim, double flmin, double g, int nang, int nfre, 
  double rowater, int ichnk, int nchnk, int ij);
