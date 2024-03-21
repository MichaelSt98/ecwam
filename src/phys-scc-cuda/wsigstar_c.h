#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void wsigstar_c(double wswave, double ufric, double z0m, double wstar, 
  double *sig_n, double acdlin, double alphamax, double alphamin, double bcdlin, 
  double epsus, double g, int llgcbz0, double rnum, double wspmin, double xkappa);
