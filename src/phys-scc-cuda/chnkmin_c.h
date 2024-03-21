#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ double chnkmin_c(double u10, double alpha, double alphamin, double chnkmin_u);
