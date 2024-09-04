

__device__ void omegagc_c(int kijs, int kijl, const double * __restrict__ ust, 
  int * __restrict__ ns, double * __restrict__ xks, double * __restrict__ oms, 
  int nwav_gc, const double * __restrict__ omega_gc, double sqrtgosurft, 
  const double * __restrict__ xkm_gc, const double * __restrict__ xk_gc, 
  double xlogkratiom1_gc, int ichnk, int nchnk, int ij);
