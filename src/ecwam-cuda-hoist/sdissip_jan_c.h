

__device__ void sdissip_jan_c(int kijs, int kijl, const double * __restrict__ fl1, 
  double * __restrict__ fld, double * __restrict__ sl, 
  const double * __restrict__ wavnum, const double * __restrict__ emean, 
  const double * __restrict__ f1mean, const double * __restrict__ xkmean, double cdis, 
  double cdisvis, double delta_sdis, int nang, int nfre, double rnu, double zpi, 
  int ichnk, int nchnk, int ij, double * __restrict__ temp1, double * __restrict__ sds, 
  double * __restrict__ x, double * __restrict__ xk2);
