#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void setice_c(int kijs, int kijl, double * fl1, const double * cicover, 
  const double * coswdif, double cithrsh, double epsmin, double flmin, int nang, 
  int nfre, int ichnk, int nchnk, int ij);
