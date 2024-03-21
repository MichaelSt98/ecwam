#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "wsigstar_c.h"

__device__ void sinput_jan_c(int ngst, int llsneg, int kijs, int kijl, 
  const double * fl1, const double * wavnum, const double * cinv, const double * xk2cg, 
  const double * wswave, const double * ufric, const double * z0m, 
  const double * coswdif, const double * sinwdif2, const double * raorw, 
  const double * wstar, const double * rnfac, double * fld, double * sl, double * spos, 
  double * xllws, double acdlin, double alphamax, double alphamin, double bcdlin, 
  double betamaxoxkappa2, double delth, double epsus, double g, int idamping, 
  int llgcbz0, int llnormagam, int nang, int nfre, double rnum, double wspmin, 
  double xkappa, double zalp, double zpi, const double * zpifr, int ichnk, int nchnk, 
  int ij);
