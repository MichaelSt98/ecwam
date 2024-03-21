#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


__device__ void wnfluxes_c(int kijs, int kijl, const int * mij, 
  const double * rhowgdfth, const double * cinv, const double * ssurf, 
  const double * cicover, const double * phiwa, const double * em, const double * f1, 
  const double * wswave, const double * wdwave, const double * ufric, 
  const double * aird, double * nphieps, double * ntauoc, double * nswh, double * nmwp, 
  double * nemotaux, double * nemotauy, double * nemowswave, double * nemophif, 
  double * tauxd, double * tauyd, double * tauocxd, double * tauocyd, double * tauoc, 
  double * phiocd, double * phieps, double * phiaw, int lnupd, double afcrv, 
  double bfcrv, double ciblock, double cithrsh, const double * costh, double egrcrv, 
  double epsu10, double epsus, const double * fr, double g, int licerun, int lwamrsetci, 
  int lwnemocou, int lwnemotauoc, int nang, int nfre, double phiepsmax, 
  double phiepsmin, const double * sinth, double tauocmax, double tauocmin, int ichnk, 
  int nchnk, int ij);
