#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "cimsstrn_c.h"
#include "stokesdrift_c.h"

__device__ void stokestrn_c(int kijs, int kijl, const double * fl1, 
  const double * wavnum, const double * stokfac, const double * depth, 
  const double * wswave, const double * wdwave, const double * cicover, 
  const double * cithick, double * ustokes, double * vstokes, double * strnms, 
  double * nemoustokes, double * nemovstokes, double * nemostrn, double cithrsh, 
  const double * costh, double delth, const double * dfim, const double * dfim_sim, 
  double flmin, const double * fr, double g, int licerun, int lwamrsetci, int lwcou, 
  int lwnemocou, int lwnemocousend, int lwnemocoustk, int lwnemocoustrn, int nang, 
  int nfre, int nfre_odd, double rowater, const double * sinth, double zpi, int ichnk, 
  int nchnk, int ij);
