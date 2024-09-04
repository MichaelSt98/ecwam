

__device__ void sdepthlim_c(int kijs, int kijl, const double * __restrict__ emaxdpt, 
  double * __restrict__ fl1, double delth, const double * __restrict__ dfim, 
  double epsmin, const double * __restrict__ fr, int nang, int nfre, double wetail, 
  int ichnk, int nchnk, int ij, double * __restrict__ em, 
  double * __restrict__ semean_temp);
