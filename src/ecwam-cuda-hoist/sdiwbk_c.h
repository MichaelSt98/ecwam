

__device__ void sdiwbk_c(int kijs, int kijl, const double * __restrict__ fl1, 
  double * __restrict__ fld, double * __restrict__ sl, 
  const double * __restrict__ depth, const double * __restrict__ emaxdpt, 
  const double * __restrict__ emean, const double * __restrict__ f1mean, int lbiwbk, 
  int nang, int nfre, int nfre_red, int ichnk, int nchnk, int ij, 
  double * __restrict__ sds);
