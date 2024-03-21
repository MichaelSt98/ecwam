#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void imphftail_c(int kijs, int kijl, const int * mij, const double * flm, 
  const double * wavnum, const double * xk2cg, double * fl1, int nang, int nfre, 
  int ichnk, int nchnk, int ij);
