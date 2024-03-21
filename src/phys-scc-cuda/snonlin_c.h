#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "transf_snl_c.h"
#include "transf_c.h"
#include "peak_ang_c.h"

__device__ void snonlin_c(int kijs, int kijl, const double * fl1, double * fld, 
  double * sl, const double * wavnum, const double * depth, const double * akmean, 
  const double * af11, double bathymax, const double * costh, double dal1, double dal2, 
  double delth, const double * dfim, const double * dfimfr, const double * dfimfr2, 
  double dkmax, const double * fklam, const double * fklam1, const double * fklap, 
  const double * fklap1, const double * fr, double fratio, double g, double gm1, 
  const int * ikm, const int * ikm1, const int * ikp, const int * ikp1, 
  const int * inlcoef, int isnonlin, const int * k11w, const int * k1w, 
  const int * k21w, const int * k2w, int kfrh, int mfrstlw, int mlsthg, int nang, 
  int nfre, const double * rnlcoef, const double * sinth, const double * th, 
  double wetail, double wp1tail, double wp2tail, double xkdmin, const double * zpifr, 
  int ichnk, int nchnk, int ij, double * enh, double * xnu, double * sig_th);
