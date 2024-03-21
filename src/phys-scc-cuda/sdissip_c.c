#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sdissip_c.h"
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
  double ssdsc6, double zpi, const double * zpifr, int ichnk, int nchnk, int ij) {
  
  


  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  
  
  switch (iphys) {
  case 0:
    sdissip_jan_c(kijs, kijl, fl1, fld, sl, wavnum, emean, f1mean, xkmean, cdis, 
      cdisvis, delta_sdis, nang, nfre, rnu, zpi, ichnk, nchnk, ij);
    
  break;
  case 1:
    sdissip_ard_c(kijs, kijl, fl1, fld, sl, indep, wavnum, xk2cg, ufric, coswdif, raorw, 
      cumulw, g, indicessat, ipsat, miche, nang, ndepth, ndikcumul, nfre, nsdsnth, 
      satweights, sdsbr, ssdsc2, ssdsc3, ssdsc4, ssdsc5, ssdsc6, zpi, zpifr, ichnk, 
      nchnk, ij);
  break;
  }
  
  
}
