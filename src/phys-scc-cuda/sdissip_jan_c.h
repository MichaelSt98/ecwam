#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void sdissip_jan_c(int kijs, int kijl, const double * fl1, double * fld, 
  double * sl, const double * wavnum, const double * emean, const double * f1mean, 
  const double * xkmean, double cdis, double cdisvis, double delta_sdis, int nang, 
  int nfre, double rnu, double zpi, int ichnk, int nchnk, int ij);
