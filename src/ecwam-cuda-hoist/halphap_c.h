

__device__ void halphap_c(int kijs, int kijl, const double * __restrict__ wavnum, 
  const double * __restrict__ coswdif, const double * __restrict__ fl1, 
  double * __restrict__ halp, double alphapmax, double delth, 
  const double * __restrict__ dfim, const double * __restrict__ dfimofr, double epsmin, 
  const double * __restrict__ fr, const double * __restrict__ fr5, double frtail, 
  int nang, int nfre, double wetail, double zpi4gm2, int ichnk, int nchnk, int ij, 
  double * __restrict__ alphap, double * __restrict__ xmss, double * __restrict__ em, 
  double * __restrict__ fm, double * __restrict__ f1d, double * __restrict__ wd, 
  double * __restrict__ flwd, double * __restrict__ femean_temp2, 
  double * __restrict__ meansqs_lf_fd, double * __restrict__ meansqs_lf_temp1, 
  double * __restrict__ meansqs_lf_temp2);
