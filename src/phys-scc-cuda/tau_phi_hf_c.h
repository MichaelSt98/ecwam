#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "omegagc_c.h"

__device__ void tau_phi_hf_c(int kijs, int kijl, const int * mij, int ltauwshelter, 
  const double * ufric, const double * z0m, const double * fl1, const double * aird, 
  const double * rnfac, const double * coswdif, const double * sinwdif2, double * ust, 
  double * tauhf, double * phihf, int llphihf, double delth, const double * fr5, 
  double g, double gamnconst, double gm1, int jtot_tauhf, int llgcbz0, int llnormagam, 
  int nang, int nwav_gc, const double * omega_gc, double sqrtgosurft, 
  double tauwshelter, const double * wtauhf, double x0tauhf, double xkappa, 
  const double * xkm_gc, const double * xk_gc, double xlogkratiom1_gc, double zalp, 
  double zpi4gm1, double zpi4gm2, const double * zpifr, int ichnk, int nchnk, int ij);
