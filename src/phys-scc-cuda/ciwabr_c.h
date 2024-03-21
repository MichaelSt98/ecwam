#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void ciwabr_c(int kijs, int kijl, const double * cicover, const double * fl1, 
  const double * wavnum, const double * cgroup, double * ciwab, double cdicwa, 
  const double * dfim, double epsmin, int idelt, int licerun, int lmaskice, int nang, 
  int nfre, int ichnk, int nchnk, int ij);
