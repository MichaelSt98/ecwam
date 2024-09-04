

__device__ void wsigstar_c(int kijs, int kijl, const double * __restrict__ wswave, 
  const double * __restrict__ ufric, const double * __restrict__ z0m, 
  const double * __restrict__ wstar, double * __restrict__ sig_n, double acdlin, 
  double alphamax, double alphamin, double bcdlin, double epsus, double g, int llgcbz0, 
  double rnum, double wspmin, double xkappa, int ichnk, int nchnk, int ij);
