

__device__ void femeanws_c(int kijs, int kijl, const double * __restrict__ fl1, 
  const double * __restrict__ xllws, double * __restrict__ fm, double * __restrict__ em, 
  double delth, const double * __restrict__ dfim, const double * __restrict__ dfimofr, 
  double epsmin, const double * __restrict__ fr, double frtail, int nang, int nfre, 
  double wetail, int ichnk, int nchnk, int ij, double * __restrict__ temp2, 
  double * __restrict__ em_loc);
