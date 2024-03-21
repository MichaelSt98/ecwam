#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "tau_phi_hf_c.h"

__device__ void stresso_c(int kijs, int kijl, const int * mij, const double * rhowgdfth, 
  const double * fl1, const double * sl, const double * spos, const double * cinv, 
  const double * wdwave, const double * ufric, const double * z0m, const double * aird, 
  const double * rnfac, const double * coswdif, const double * sinwdif2, double * tauw, 
  double * tauwdir, double * phiwa, int llphiwa, const double * costh, double delth, 
  double eps1, const double * fr5, double g, double gamnconst, double gm1, int iphys, 
  int jtot_tauhf, int llgcbz0, int llnormagam, int nang, int nfre, int nwav_gc, 
  const double * omega_gc, const double * rhowg_dfim, const double * sinth, 
  double sqrtgosurft, double tauwshelter, const double * wtauhf, double x0tauhf, 
  double xkappa, const double * xkm_gc, const double * xk_gc, double xlogkratiom1_gc, 
  double zalp, double zpi4gm1, double zpi4gm2, const double * zpifr, int ichnk, 
  int nchnk, int ij, double * xstress, double * ystress, double * tauhf, double * phihf, 
  double * usdirp, double * ust);
