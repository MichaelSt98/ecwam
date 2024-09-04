#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "wsigstar_c.h"


__device__ void wsigstar_c(int kijs, int kijl, const double * __restrict__ wswave, 
  const double * __restrict__ ufric, const double * __restrict__ z0m, 
  const double * __restrict__ wstar, double * __restrict__ sig_n, double acdlin, 
  double alphamax, double alphamin, double bcdlin, double epsus, double g, int llgcbz0, 
  double rnum, double wspmin, double xkappa, int ichnk, int nchnk, int ij) {
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  
  double bg_gust = (double) 0.0;  // NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
  double onethird = (double) 1.0 / (double) 3.0;
  double sig_nmax = (double) 0.9;  // MAX OF RELATIVE STANDARD DEVIATION OF USTAR
  
  double log10 = log((double) 10.0);
  double c1 = (double) 1.03E-3;
  double c2 = (double) 0.04E-3;
  double p1 = (double) 1.48;
  double p2 = -(double) 0.21;
  
  double zchar;
  double c_d;
  double dc_ddu;
  double sig_conv;
  double xkappad;
  double u10;
  double c2u10p1;
  double u10p2;
  double bcd;
  double u10m1;
  double zn;
  double z0vis;
  double zhook_handle;
  
  
  // ----------------------------------------------------------------------
  
  
  
  if (llgcbz0) {
    zn = rnum;
    
    u10m1 = (double) 1.0 / fmax((double) (wswave[ij - 1 + kijl*(ichnk - 1)]), (double) 
      (wspmin));
    // CHARNOCK:
    z0vis = zn / fmax((double) (ufric[ij - 1 + kijl*(ichnk - 1)]), (double) (epsus));
    zchar = g*(z0m[ij - 1 + kijl*(ichnk - 1)] - z0vis) / fmax((double) (pow(ufric[ij - 1 
      + kijl*(ichnk - 1)], 2)), (double) (epsus));
    zchar = 
      fmax((double) (fmin((double) (zchar), (double) (alphamax))), (double) (alphamin));
    
    bcd = bcdlin*sqrt((double) (zchar));
    c_d = acdlin + bcd*wswave[ij - 1 + kijl*(ichnk - 1)];
    dc_ddu = bcd;
    sig_conv = (double) 1.0 + (double) 0.5*wswave[ij - 1 + kijl*(ichnk - 1)] / c_d*dc_ddu
      ;
    sig_n[ij - 1] = fmin((double) (sig_nmax), (double) 
      (sig_conv*u10m1*(pow((bg_gust*(pow(ufric[ij - 1 + kijl*(ichnk - 1)], 3)) + (double)
       0.5*xkappa*(pow(wstar[ij - 1 + kijl*(ichnk - 1)], 3))), onethird))));
    
  } else {
    zn = (double) 0.0;
    
    
    //!! for consistency I have kept the old method, even though the new method above could be used,
    //!! but until LLGCBZ0 is the default, keep the old scheme whe it is not...
    //
    //       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
    //       BASED ON U*
    //
    xkappad = (double) 1.0 / xkappa;
    u10 = ufric[ij - 1 + kijl*(ichnk - 1)]*xkappad*(log10 - log(z0m[ij - 1 + kijl*(ichnk 
      - 1)]));
    u10 = fmax((double) (u10), (double) (wspmin));
    u10m1 = (double) 1.0 / u10;
    c2u10p1 = c2*(pow(u10, p1));
    u10p2 = pow(u10, p2);
    c_d = (c1 + c2u10p1)*u10p2;
    dc_ddu = (p2*c1 + (p1 + p2)*c2u10p1)*u10p2*u10m1;
    sig_conv = (double) 1.0 + (double) 0.5*u10 / c_d*dc_ddu;
    sig_n[ij - 1] = fmin((double) (sig_nmax), (double) 
      (sig_conv*u10m1*(pow((bg_gust*(pow(ufric[ij - 1 + kijl*(ichnk - 1)], 3)) + (double)
       0.5*xkappa*(pow(wstar[ij - 1 + kijl*(ichnk - 1)], 3))), onethird))));
    
  }
  
  
  
}
