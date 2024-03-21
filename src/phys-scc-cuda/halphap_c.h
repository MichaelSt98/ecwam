#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void halphap_c(int kijs, int kijl, const double * wavnum, 
  const double * coswdif, const double * fl1, double * halp, double alphapmax, 
  double delth, const double * dfim, const double * dfimofr, double epsmin, 
  const double * fr, const double * fr5, double frtail, int nang, int nfre, 
  double wetail, double zpi4gm2, int ichnk, int nchnk, int ij);
