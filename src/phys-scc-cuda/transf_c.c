#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "transf_c.h"

__device__ double transf_c(double xk, double d, double dkmax, double g) {
  
  
  double transf;
  double eps, x, t_0, om, c_0, v_g, dv_g, xnl_1, xnl_2, xnl;
  // <Pragma:: acc routine seq>
  
  eps = (double) 0.0001;
  if (d < (double) 999.0 && d > (double) 0.0) {
    x = xk*d;
    if (x > dkmax) {
      transf = (double) 1.0;
    } else {
      t_0 = tanh(x);
      om = sqrt((double) (g*xk*t_0));
      c_0 = om / xk;
      if (x < eps) {
        v_g = (double) 0.5*c_0;
        v_g = c_0;
      } else {
        v_g = (double) 0.5*c_0*((double) 1.0 + (double) 2.0*x / sinh((double) 2.0*x));
      }
      dv_g = (pow((t_0 - x*((double) 1.0 - (pow(t_0, 2)))), 2)) + (double) 4.0*(pow(x, 2)
        )*(pow(t_0, 2))*((double) 1.0 - (pow(t_0, 2)));
      
      xnl_1 = ((double) 9.0*(pow(t_0, 4)) - (double) 10.0*(pow(t_0, 2)) + (double) 9.0) /
         ((double) 8.0*(pow(t_0, 3)));
      xnl_2 = ((pow(((double) 2.0*v_g - (double) 0.5*c_0), 2)) / (g*d - (pow(v_g, 2))) + 
        (double) 1.0) / x;
      
      xnl = xnl_1 - xnl_2;
      transf = (pow(xnl, 2)) / (dv_g*(pow(t_0, 8)));
    }
  } else {
    transf = (double) 1.0;
  }
  return transf;
  //
}
