#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void fkmean_c(int kijs, int kijl, const double * fl1, const double * wavnum, 
  double * em, double * fm1, double * f1, double * ak, double * xk, double delth, 
  const double * dfim, const double * dfimfr, const double * dfimofr, double epsmin, 
  const double * fr, double frtail, double g, int nang, int nfre, double wetail, 
  double wp1tail, double zpi, int ichnk, int nchnk, int ij);
