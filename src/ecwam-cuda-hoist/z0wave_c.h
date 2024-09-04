

__device__ void z0wave_c(int kijs, int kijl, const double * __restrict__ us, 
  const double * __restrict__ tauw, const double * __restrict__ utop, 
  double * __restrict__ z0, double * __restrict__ z0b, double * __restrict__ chrnck, 
  double alpha, double alphamin, double chnkmin_u, double eps1, double g, double gm1, 
  int llcapchnk, int ichnk, int nchnk, int ij, double * __restrict__ alphaog);
