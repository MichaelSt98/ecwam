#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "stokestrn_c.h"
#include "stokesdrift_c.h"
#include "cimsstrn_c.h"


__device__ void stokestrn_c(int kijs, int kijl, const double * __restrict__ fl1, 
  const double * __restrict__ wavnum, const double * __restrict__ stokfac, 
  const double * __restrict__ depth, const double * __restrict__ wswave, 
  const double * __restrict__ wdwave, const double * __restrict__ cicover, 
  const double * __restrict__ cithick, double * __restrict__ ustokes, 
  double * __restrict__ vstokes, double * __restrict__ strnms, 
  double * __restrict__ nemoustokes, double * __restrict__ nemovstokes, 
  double * __restrict__ nemostrn, double cithrsh, const double * __restrict__ costh, 
  double delth, const double * __restrict__ dfim, const double * __restrict__ dfim_sim, 
  double flmin, const double * __restrict__ fr, double g, int licerun, int lwamrsetci, 
  int lwcou, int lwnemocou, int lwnemocousend, int lwnemocoustk, int lwnemocoustrn, 
  int nang, int nfre, int nfre_odd, double rowater, const double * __restrict__ sinth, 
  double zpi, int ichnk, int nchnk, int ij, double * __restrict__ cimsstrn_xki, 
  double * __restrict__ cimsstrn_e, double * __restrict__ cimsstrn_sume, 
  double * __restrict__ stokesdrift_stfac) {
  
  
  
  // ----------------------------------------------------------------------
  


  
  
  
  double zhook_handle;
  
  // ----------------------------------------------------------------------
  
  
  stokesdrift_c(kijs, kijl, fl1, stokfac, wswave, wdwave, cicover, ustokes, vstokes, 
    cithrsh, costh, delth, dfim_sim, fr, g, licerun, lwamrsetci, nang, nfre, nfre_odd, 
    sinth, zpi, ichnk, nchnk, ij, stokesdrift_stfac);
  
  if (lwnemocoustrn) {
    cimsstrn_c(kijs, kijl, fl1, wavnum, depth, cithick, strnms, delth, dfim, flmin, g, 
      nang, nfre, rowater, ichnk, nchnk, ij, cimsstrn_xki, cimsstrn_e, cimsstrn_sume);
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
  
  
  
  // ----------------------------------------------------------------------
  
}
