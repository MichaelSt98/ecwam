#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void sdiwbk_c(int kijs, int kijl, const double * fl1, double * fld, 
  double * sl, const double * depth, const double * emaxdpt, const double * emean, 
  const double * f1mean, int lbiwbk, int nang, int nfre_red, int ichnk, int nchnk, int ij
  );
