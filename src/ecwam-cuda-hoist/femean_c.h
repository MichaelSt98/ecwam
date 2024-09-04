

__device__ void femean_c(int kijs, int kijl, const double * __restrict__ f, 
  double * __restrict__ em, double * __restrict__ fm, double delth, 
  const double * __restrict__ dfim, const double * __restrict__ dfimofr, double epsmin, 
  const double * __restrict__ fr, double frtail, int nang, int nfre, double wetail, 
  int ichnk, int nchnk, int ij, double * __restrict__ temp2);
