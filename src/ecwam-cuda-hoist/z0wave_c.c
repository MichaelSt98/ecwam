#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "z0wave_c.h"
#include "chnkmin_c.h"


__device__ void z0wave_c(int kijs, int kijl, const double * __restrict__ us, 
  const double * __restrict__ tauw, const double * __restrict__ utop, 
  double * __restrict__ z0, double * __restrict__ z0b, double * __restrict__ chrnck, 
  double alpha, double alphamin, double chnkmin_u, double eps1, double g, double gm1, 
  int llcapchnk, int ichnk, int nchnk, int ij, double * __restrict__ alphaog) {
  
  
  
  // ----------------------------------------------------------------------
  

  
  
  double ust2;
  double ust3;
  double arg;
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  
  
  if (llcapchnk) {
    alphaog[ij - 1 + kijl*(ichnk - 1)] = 
      chnkmin_c(utop[ij - 1 + kijl*(ichnk - 1)], alpha, alphamin, chnkmin_u)*gm1;
  } else {
    alphaog[ij - 1 + kijl*(ichnk - 1)] = alpha*gm1;
  }
  
  ust2 = pow(us[ij - 1 + kijl*(ichnk - 1)], 2);
  ust3 = pow(us[ij - 1 + kijl*(ichnk - 1)], 3);
  arg = fmax((double) (ust2 - tauw[ij - 1 + kijl*(ichnk - 1)]), (double) (eps1));
  z0[ij - 1 + kijl*(ichnk - 1)] = 
    alphaog[ij - 1 + kijl*(ichnk - 1)]*ust3 / sqrt((double) (arg));
  z0b[ij - 1 + kijl*(ichnk - 1)] = alphaog[ij - 1 + kijl*(ichnk - 1)]*ust2;
  chrnck[ij - 1 + kijl*(ichnk - 1)] = g*z0[ij - 1 + kijl*(ichnk - 1)] / ust2;
  
  
  
}
