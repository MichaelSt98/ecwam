#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ double transf_snl_c(double xk0, double d, double xnu, double sig_th, 
  double bathymax, double dkmax, double g, double xkdmin);
