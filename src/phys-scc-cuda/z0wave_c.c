#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "z0wave_c.h"
#include "chnkmin_c.h"

__device__ void z0wave_c(int kijs, int kijl, const double * us, const double * tauw, 
  const double * utop, double * z0, double * z0b, double * chrnck, double alpha, 
  double alphamin, double chnkmin_u, double eps1, double g, double gm1, int llcapchnk, 
  int ichnk, int nchnk, int ij) {
  
  

  double ust2;
  double ust3;
  double arg;
  double alphaog;
  

  if (llcapchnk) {
    alphaog = chnkmin_c(utop[ij - 1 + kijl*(ichnk - 1)], alpha, alphamin, chnkmin_u)*gm1;
  } else {
    alphaog = alpha*gm1;
  }
  
  ust2 = pow(us[ij - 1 + kijl*(ichnk - 1)], 2);
  ust3 = pow(us[ij - 1 + kijl*(ichnk - 1)], 3);
  arg = max((double) (ust2 - tauw[ij - 1 + kijl*(ichnk - 1)]), (double) (eps1));
  z0[ij - 1 + kijl*(ichnk - 1)] = alphaog*ust3 / sqrt((double) (arg));
  z0b[ij - 1 + kijl*(ichnk - 1)] = alphaog*ust2;
  chrnck[ij - 1 + kijl*(ichnk - 1)] = g*z0[ij - 1 + kijl*(ichnk - 1)] / ust2;

  
  
}
