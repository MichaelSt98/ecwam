#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void peak_ang_c(int kijs, int kijl, const double * fl1, double * xnu, 
  double * sig_th, const double * costh, double delth, const double * dfim, 
  const double * dfimfr, const double * dfimfr2, const double * fr, double fratio, 
  int nang, int nfre, const double * sinth, const double * th, double wetail, 
  double wp1tail, double wp2tail, int ichnk, int nchnk, int ij);
