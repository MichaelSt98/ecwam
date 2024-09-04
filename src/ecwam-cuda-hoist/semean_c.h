

__device__ void semean_c(const double * __restrict__ fl1, int kijs, int kijl, 
  double * __restrict__ em, int llepsmin, double delth, 
  const double * __restrict__ dfim, double epsmin, const double * __restrict__ fr, 
  int nang, int nfre, double wetail, int ichnk, int nchnk, int ij, 
  double * __restrict__ temp);
