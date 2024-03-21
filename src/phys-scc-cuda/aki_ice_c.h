#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ double aki_ice_c(double g, double xk, double depth, double rhow, double cith);
