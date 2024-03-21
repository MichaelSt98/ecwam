#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void stokesdrift_c(int kijs, int kijl, const double * fl1, 
  const double * stokfac, const double * wswave, const double * wdwave, 
  const double * cicover, double * ustokes, double * vstokes, double cithrsh, 
  const double * costh, double delth, const double * dfim_sim, const double * fr, 
  double g, int licerun, int lwamrsetci, int nang, int nfre_odd, const double * sinth, 
  double zpi, int ichnk, int nchnk, int ij);
