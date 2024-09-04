

__device__ void peak_ang_c(int kijs, int kijl, const double * __restrict__ fl1, 
  double * __restrict__ xnu, double * __restrict__ sig_th, 
  const double * __restrict__ costh, double delth, const double * __restrict__ dfim, 
  const double * __restrict__ dfimfr, const double * __restrict__ dfimfr2, 
  const double * __restrict__ fr, double fratio, int nang, int nfre, 
  const double * __restrict__ sinth, const double * __restrict__ th, double wetail, 
  double wp1tail, double wp2tail, int ichnk, int nchnk, int ij, int * __restrict__ mmax, 
  double * __restrict__ sum0, double * __restrict__ sum1, double * __restrict__ sum2, 
  double * __restrict__ xmax, double * __restrict__ temp, double * __restrict__ thmean, 
  double * __restrict__ sum_s, double * __restrict__ sum_c);
