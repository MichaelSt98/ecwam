#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void sdepthlim_c(int kijs, int kijl, const double * emaxdpt, double * fl1, 
  double delth, const double * dfim, double epsmin, const double * fr, int nang, 
  int nfre, double wetail, int ichnk, int nchnk, int ij);
