#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void femeanws_c(int kijs, int kijl, const double * fl1, const double * xllws, 
  double * fm, double * em, double delth, const double * dfim, const double * dfimofr, 
  double epsmin, const double * fr, double frtail, int nang, int nfre, double wetail, 
  int ichnk, int nchnk, int ij);
