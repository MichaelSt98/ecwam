

__device__ void cimsstrn_c(int kijs, int kijl, const double * __restrict__ fl1, 
  const double * __restrict__ wavnum, const double * __restrict__ depth, 
  const double * __restrict__ cithick, double * __restrict__ strn, double delth, 
  const double * __restrict__ dfim, double flmin, double g, int nang, int nfre, 
  double rowater, int ichnk, int nchnk, int ij, double * __restrict__ xki, 
  double * __restrict__ e, double * __restrict__ sume);
