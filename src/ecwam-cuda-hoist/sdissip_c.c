#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sdissip_c.h"
#include "sdissip_jan_c.h"
#include "sdissip_ard_c.h"


__device__ void sdissip_c(int kijs, int kijl, const double * __restrict__ fl1, 
  double * __restrict__ fld, double * __restrict__ sl, 
  const double * __restrict__ wavnum, const double * __restrict__ cgroup, 
  const double * __restrict__ xk2cg, const double * __restrict__ emean, 
  const double * __restrict__ f1mean, const double * __restrict__ xkmean, 
  const double * __restrict__ ufric, const double * __restrict__ coswdif, 
  const double * __restrict__ raorw, double brkpbcoef, double cdis, double cdisvis, 
  double delta_sdis, double delth, double fratio, double g, 
  const int * __restrict__ indicessat, int iphys, int ipsat, double miche, int nang, 
  int nfre, int nsdsnth, double rnu, const double * __restrict__ satweights, 
  double sdsbr, double ssdsbrf1, double ssdsc2, double ssdsc3, double ssdsc4, 
  double ssdsc5, double ssdsc6, double zpi, const double * __restrict__ zpifr, 
  int ichnk, int nchnk, int ij, double * __restrict__ sdissip_ard_facturb, 
  double * __restrict__ sdissip_ard_facsat, double * __restrict__ sdissip_ard_facwtrb, 
  double * __restrict__ sdissip_ard_temp1, double * __restrict__ sdissip_ard_bth0, 
  double * __restrict__ sdissip_ard_c_, double * __restrict__ sdissip_ard_c_c, 
  double * __restrict__ sdissip_ard_dsip, double * __restrict__ sdissip_ard_trpz_dsip, 
  double * __restrict__ sdissip_ard_bth, double * __restrict__ sdissip_ard_temp2, 
  double * __restrict__ sdissip_ard_d, double * __restrict__ sdissip_ard_scumul, 
  double * __restrict__ sdissip_ard_renewalfreq, 
  double * __restrict__ sdissip_ard_wcumul, double * __restrict__ sdissip_jan_temp1, 
  double * __restrict__ sdissip_jan_sds, double * __restrict__ sdissip_jan_x, 
  double * __restrict__ sdissip_jan_xk2) {
  
  
  
  // ----------------------------------------------------------------------
  


  
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  
  switch (iphys) {
  case 0:
    sdissip_jan_c(kijs, kijl, fl1, fld, sl, wavnum, emean, f1mean, xkmean, cdis, 
      cdisvis, delta_sdis, nang, nfre, rnu, zpi, ichnk, nchnk, ij, sdissip_jan_temp1, 
      sdissip_jan_sds, sdissip_jan_x, sdissip_jan_xk2);
    
  break;
  case 1:
    sdissip_ard_c(kijs, kijl, fl1, fld, sl, wavnum, cgroup, xk2cg, ufric, coswdif, 
      raorw, brkpbcoef, delth, fratio, g, indicessat, ipsat, miche, nang, nfre, nsdsnth, 
      satweights, sdsbr, ssdsbrf1, ssdsc2, ssdsc3, ssdsc4, ssdsc5, ssdsc6, zpi, zpifr, 
      ichnk, nchnk, ij, sdissip_ard_facturb, sdissip_ard_facsat, sdissip_ard_facwtrb, 
      sdissip_ard_temp1, sdissip_ard_bth0, sdissip_ard_c_, sdissip_ard_c_c, 
      sdissip_ard_dsip, sdissip_ard_trpz_dsip, sdissip_ard_bth, sdissip_ard_temp2, 
      sdissip_ard_d, sdissip_ard_scumul, sdissip_ard_renewalfreq, sdissip_ard_wcumul);
  break;
  }
  
  
}
