#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sdissip_ard_c.h"
#include "sdissip_jan_c.h"

__device__ void sdissip_c(int kijs, int kijl, const double * fl1, double * fld, 
  double * sl, const int * indep, const double * wavnum, const double * xk2cg, 
  const double * emean, const double * f1mean, const double * xkmean, 
  const double * ufric, const double * coswdif, const double * raorw, double cdis, 
  double cdisvis, const double * cumulw, double delta_sdis, double g, 
  const int * indicessat, int iphys, int ipsat, double miche, int nang, int ndepth, 
  int ndikcumul, int nfre, int nsdsnth, double rnu, const double * satweights, 
  double sdsbr, double ssdsc2, double ssdsc3, double ssdsc4, double ssdsc5, 
  double ssdsc6, double zpi, const double * zpifr, int ichnk, int nchnk, int ij);
