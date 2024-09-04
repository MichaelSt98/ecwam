

__device__ void fkmean_c(int kijs, int kijl, const double * __restrict__ fl1, 
  const double * __restrict__ wavnum, double * __restrict__ em, 
  double * __restrict__ fm1, double * __restrict__ f1, double * __restrict__ ak, 
  double * __restrict__ xk, double delth, const double * __restrict__ dfim, 
  const double * __restrict__ dfimfr, const double * __restrict__ dfimofr, 
  double epsmin, const double * __restrict__ fr, double frtail, double g, int nang, 
  int nfre, double wetail, double wp1tail, double zpi, int ichnk, int nchnk, int ij, 
  double * __restrict__ tempa, double * __restrict__ tempx, double * __restrict__ temp2);
