#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void sbottom_c(int kijs, int kijl, const double * fl1, double * fld, 
  double * sl, const double * wavnum, const double * depth, double bathymax, double gm1, 
  int nang, int nfre_red, int ichnk, int nchnk, int ij);
