

__device__ void setice_c(int kijs, int kijl, double * __restrict__ fl1, 
  const double * __restrict__ cicover, const double * __restrict__ coswdif, 
  double cithrsh, double epsmin, double flmin, int nang, int nfre, int ichnk, int nchnk, 
  int ij, double * __restrict__ cireduc, double * __restrict__ temp, 
  double * __restrict__ icefree);
