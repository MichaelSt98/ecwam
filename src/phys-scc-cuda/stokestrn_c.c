#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "stokestrn_c.h"
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
  int nchnk, int ij) {
  
  


  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  
  
  stokesdrift_c(kijs, kijl, fl1, stokfac, wswave, wdwave, cicover, ustokes, vstokes, 
    cithrsh, costh, delth, dfim_sim, fr, g, licerun, lwamrsetci, nang, nfre_odd, sinth, 
    zpi, ichnk, nchnk, ij);
  
  if (lwnemocoustrn) {
    cimsstrn_c(kijs, kijl, fl1, wavnum, depth, cithick, strnms, delth, dfim, flmin, g, 
      nang, nfre, rowater, ichnk, nchnk, ij);
  }

  if (lwnemocou && (lwnemocousend && lwcou || !lwcou)) {
    if (lwnemocoustk) {
      nemoustokes[ij - 1 + kijl*(ichnk - 1)] = ustokes[ij - 1 + kijl*(ichnk - 1)];
      nemovstokes[ij - 1 + kijl*(ichnk - 1)] = vstokes[ij - 1 + kijl*(ichnk - 1)];
    } else {
      nemoustokes[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
      nemovstokes[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
    }
    
    if (lwnemocoustrn) {
      nemostrn[ij - 1 + kijl*(ichnk - 1)] = strnms[ij - 1 + kijl*(ichnk - 1)];
    }
  }

  
}
