#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "aki_ice_c.h"

__device__ double aki_ice_c(double g, double xk, double depth, double rhow, double cith) {
  
  
  double aki_ice;
  double ymice = (double) 5.5E+9;  // typical value of Young modulus of sea ice
  double rmuice = (double) 0.3;  // Poisson's ratio of sea ice
  double rhoi = (double) 922.5;  // typical value of the sea ice density
  double ebs = (double) 0.000001;
  //     MAXIMUM WAVE NUMBER
  double aki_max = (double) 20.0;
  
  double ficstf, rdh;
  double om2, aki, akiold, f, fprime, akid;
  // <Pragma:: acc routine seq>
  if (cith <= (double) 0.0) {
    aki = xk;
  } else {
    //       BENDING STIFFNESS / WATER DENSITY
    ficstf = (ymice*(pow(cith, 3)) / (12*(1 - (pow(rmuice, 2))))) / rhow;
    rdh = (rhoi / rhow)*cith;
    om2 = g*xk*tanh(xk*depth);
    akiold = (double) 0.0;
    aki = min((double) (xk), (double) (pow((om2 / max((double) (ficstf), (double) 
      ((double) 1.0))), (double) 0.2)));
    
    while (abs((double) (aki - akiold)) > ebs*akiold && aki < aki_max) {
      akiold = aki;
      akid = min((double) (depth*aki), (double) ((double) 50.0));
      f = ficstf*(pow(aki, 5)) + g*aki - om2*(rdh*aki + 1. / tanh(akid));
      fprime = 
        (double) 5.*ficstf*(pow(aki, 4)) + g - om2*(rdh - depth / (pow(sinh(akid), 2)));
      aki = aki - f / fprime;
      //         in case of overshoot because it is trying to find a very large wave number
      if (aki <= (double) 0.0) {
        aki = aki_max;
      }
    }
    
  }
  
  aki_ice = aki;
  return aki_ice;  
}
