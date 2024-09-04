

__device__ void sbottom_c(int kijs, int kijl, const double * __restrict__ fl1, 
  double * __restrict__ fld, double * __restrict__ sl, 
  const double * __restrict__ wavnum, const double * __restrict__ depth, 
  double bathymax, double gm1, int nang, int nfre, int nfre_red, int ichnk, int nchnk, 
  int ij, double * __restrict__ sbo);
