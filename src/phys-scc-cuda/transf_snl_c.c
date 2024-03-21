#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "transf_snl_c.h"

__device__ double transf_snl_c(double xk0, double d, double xnu, double sig_th, 
  double bathymax, double dkmax, double g, double xkdmin) {
  
  
  
  double transf_snl;
  
  double eps = (double) 0.0001;
  double transf_snl_min = (double) 0.1;
  double transf_snl_max = (double) 10.;
  
  double x, xk, t_0, t_0_sq, om, c_0, v_g, v_g_sq, dv_g;
  double xnl_1, xnl_2, xnl_3, xnl_4, xnl;
  double c_s_sq, alp, zfac;
  // <Pragma:: acc routine seq>
  if (d < bathymax && d > (double) 0.) {
    x = xk0*d;
    if (x > dkmax) {
      transf_snl = (double) 1.;
    } else {
      xk = max((double) (xk0), (double) (xkdmin / d));
      x = xk*d;
      t_0 = tanh(x);
      t_0_sq = pow(t_0, 2);
      om = sqrt((double) (g*xk*t_0));
      c_0 = om / xk;
      c_s_sq = g*d;
      if (x < eps) {
        v_g = c_0;
      } else {
        v_g = (double) 0.5*c_0*((double) 1. + (double) 2.*x / sinh((double) 2.*x));
      }
      v_g_sq = pow(v_g, 2);
      dv_g = (pow((t_0 - x*(1. - t_0_sq)), 2)) + (double) 4.*(pow(x, 2))*t_0_sq*((double)
         1. - t_0_sq);
      
      xnl_1 = ((double) 9.*(pow(t_0_sq, 2)) - (double) 10.*t_0_sq + (double) 9.) / 
        ((double) 8.*t_0_sq*t_0);
      xnl_2 = ((pow(((double) 2.*v_g - (double) 0.5*c_0), 2)) / (g*d - v_g_sq) + (double)
         1.) / x;
      xnl_4 = 1. / ((double) 4.*t_0)*(pow(((double) 2.*c_0 + v_g*((double) 1. - t_0_sq)),
         2)) / (c_s_sq - v_g_sq);
      alp = (1. - v_g_sq / c_s_sq)*(pow(c_0, 2)) / v_g_sq;
      zfac = (pow(sig_th, 2)) / ((pow(sig_th, 2)) + alp*(pow(xnu, 2)));
      xnl_3 = zfac*xnl_4;
      
      xnl = xnl_1 - xnl_2 + xnl_3;
      transf_snl = (pow(xnl, 2)) / (dv_g*(pow(t_0_sq, 4)));
      transf_snl = max((double) (min((double) (transf_snl_max), (double) (transf_snl))), 
        (double) (transf_snl_min));
    }
  } else {
    transf_snl = (double) 1.;
  }
  return transf_snl; 
}
