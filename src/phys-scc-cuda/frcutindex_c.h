#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void frcutindex_c(int kijs, int kijl, const double * fm, const double * fmws, 
  const double * ufric, const double * cicover, int * mij, double * rhowgdfth, 
  double cithrsh_tail, double epsmin, double flogsprdm1, const double * fr, double fric, 
  double g, int nfre, const double * rhowg_dfim, double tailfactor, 
  double tailfactor_pm, const double * zpifr, int ichnk, int nchnk, int ij);
