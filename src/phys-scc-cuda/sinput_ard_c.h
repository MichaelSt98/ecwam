#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "wsigstar_c.h"

__device__ void sinput_ard_c(int ngst, int llsneg, int kijs, int kijl, 
  const double * fl1, const double * wavnum, const double * cinv, const double * xk2cg, 
  const double * wdwave, const double * wswave, const double * ufric, 
  const double * z0m, const double * coswdif, const double * sinwdif2, 
  const double * raorw, const double * wstar, const double * rnfac, double * fld, 
  double * sl, double * spos, double * xllws, double abmax, double abmin, double acdlin, 
  double alphamax, double alphamin, double bcdlin, double betamaxoxkappa2, 
  const double * costh, double delth, const double * dfim, double epsmin, double epsus, 
  double g, int iab, int llgcbz0, int llnormagam, int nang, int nfre, double rnu, 
  double rnum, const double * sinth, double swellf, double swellf2, double swellf3, 
  double swellf4, double swellf5, double swellf6, double swellf7, double swellf7m1, 
  const double * swellft, double tauwshelter, const double * th, double wspmin, 
  double xkappa, double z0rat, double z0tubmax, double zalp, double zpi, 
  const double * zpifr, int ichnk, int nchnk, int ij);
