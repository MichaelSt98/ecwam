#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sinput_c.h"
#include "sinput_ard_c.h"
#include "sinput_jan_c.c"

__device__ void sinput_c(int ngst, int llsneg, int kijs, int kijl, const double * fl1, 
  const double * wavnum, const double * cinv, const double * xk2cg, 
  const double * wdwave, const double * wswave, const double * ufric, 
  const double * z0m, const double * coswdif, const double * sinwdif2, 
  const double * raorw, const double * wstar, const double * rnfac, double * fld, 
  double * sl, double * spos, double * xllws, double abmax, double abmin, double acdlin, 
  double alphamax, double alphamin, double bcdlin, double betamaxoxkappa2, 
  const double * costh, double delth, const double * dfim, double epsmin, double epsus, 
  double g, int iab, int idamping, int iphys, int llgcbz0, int llnormagam, int nang, 
  int nfre, double rnu, double rnum, const double * sinth, double swellf, 
  double swellf2, double swellf3, double swellf4, double swellf5, double swellf6, 
  double swellf7, double swellf7m1, const double * swellft, double tauwshelter, 
  const double * th, double wspmin, double xkappa, double z0rat, double z0tubmax, 
  double zalp, double zpi, const double * zpifr, int ichnk, int nchnk, int ij) {
  
  
  
  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  
  
  
  
  switch (iphys) {
  case 0:
    sinput_jan_c(ngst, llsneg, kijs, kijl, fl1, wavnum, cinv, xk2cg, wswave, ufric, z0m, 
      coswdif, sinwdif2, raorw, wstar, rnfac, fld, sl, spos, xllws, acdlin, alphamax, 
      alphamin, bcdlin, betamaxoxkappa2, delth, epsus, g, idamping, llgcbz0, llnormagam, 
      nang, nfre, rnum, wspmin, xkappa, zalp, zpi, zpifr, ichnk, nchnk, ij);
  break;
  case 1:
    sinput_ard_c(ngst, llsneg, kijs, kijl, fl1, wavnum, cinv, xk2cg, wdwave, wswave, 
      ufric, z0m, coswdif, sinwdif2, raorw, wstar, rnfac, fld, sl, spos, xllws, abmax, 
      abmin, acdlin, alphamax, alphamin, bcdlin, betamaxoxkappa2, costh, delth, dfim, 
      epsmin, epsus, g, iab, llgcbz0, llnormagam, nang, nfre, rnu, rnum, sinth, swellf, 
      swellf2, swellf3, swellf4, swellf5, swellf6, swellf7, swellf7m1, swellft, 
      tauwshelter, th, wspmin, xkappa, z0rat, z0tubmax, zalp, zpi, zpifr, ichnk, nchnk, 
      ij);
  break;
  }
  
  
}
