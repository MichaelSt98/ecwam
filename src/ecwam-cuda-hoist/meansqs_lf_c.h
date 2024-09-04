

__device__ void meansqs_lf_c(int nfre_eff, int kijs, int kijl, 
  const double * __restrict__ f, const double * __restrict__ wavnum, 
  double * __restrict__ xmss, const double * __restrict__ dfim, int nang, int nfre, 
  int ichnk, int nchnk, int ij, double * __restrict__ fd, double * __restrict__ temp1, 
  double * __restrict__ temp2);
