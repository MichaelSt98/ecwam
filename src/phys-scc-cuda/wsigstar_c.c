#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "wsigstar_c.h"

__device__ void wsigstar_c(double wswave, double ufric, double z0m, double wstar, 
  double *sig_n, double acdlin, double alphamax, double alphamin, double bcdlin, 
  double epsus, double g, int llgcbz0, double rnum, double wspmin, double xkappa) {
  
  
  
  
  double bg_gust = (double) 0.0;  // NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
  double onethird = (double) 1.0 / (double) 3.0;
  double sig_nmax = (double) 0.9;  // MAX OF RELATIVE STANDARD DEVIATION OF USTAR
  
  double log10 = log((double) 10.0);
  double c1 = (double) 1.03E-3;
  double c2 = (double) 0.04E-3;
  double p1 = (double) 1.48;
  double p2 = -(double) 0.21;
  double zchar, c_d, dc_ddu, sig_conv;
  double xkappad, u10, c2u10p1, u10p2;
  double bcd, u10m1, zn, z0vis;
  // <Pragma:: acc routine seq>
  
  
  if (llgcbz0) {
    zn = rnum;
    
    u10m1 = (double) 1.0 / max((double) (wswave), (double) (wspmin));
    // CHARNOCK:
    z0vis = zn / max((double) (ufric), (double) (epsus));
    zchar = g*(z0m - z0vis) / max((double) (pow(ufric, 2)), (double) (epsus));
    zchar = 
      max((double) (min((double) (zchar), (double) (alphamax))), (double) (alphamin));
    
    bcd = bcdlin*sqrt((double) (zchar));
    c_d = acdlin + bcd*wswave;
    dc_ddu = bcd;
    sig_conv = (double) 1.0 + (double) 0.5*wswave / c_d*dc_ddu;
     (*sig_n) = min((double) (sig_nmax), (double) 
      (sig_conv*u10m1*(pow((bg_gust*(pow(ufric, 3)) + (double) 0.5*xkappa*(pow(wstar, 3))
      ), onethird))));
  } else {
    zn = (double) 0.0;
    xkappad = (double) 1.0 / xkappa;
    u10 = ufric*xkappad*(log10 - log(z0m));
    u10 = max((double) (u10), (double) (wspmin));
    u10m1 = (double) 1.0 / u10;
    c2u10p1 = c2*(pow(u10, p1));
    u10p2 = pow(u10, p2);
    c_d = (c1 + c2u10p1)*u10p2;
    dc_ddu = (p2*c1 + (p1 + p2)*c2u10p1)*u10p2*u10m1;
    sig_conv = (double) 1.0 + (double) 0.5*u10 / c_d*dc_ddu;
     (*sig_n) = min((double) (sig_nmax), (double) 
      (sig_conv*u10m1*(pow((bg_gust*(pow(ufric, 3)) + (double) 0.5*xkappa*(pow(wstar, 3))
      ), onethird))));
  }
  
  
}
