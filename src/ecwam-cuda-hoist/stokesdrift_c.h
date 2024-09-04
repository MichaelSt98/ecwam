

__device__ void stokesdrift_c(int kijs, int kijl, const double * __restrict__ fl1, 
  const double * __restrict__ stokfac, const double * __restrict__ wswave, 
  const double * __restrict__ wdwave, const double * __restrict__ cicover, 
  double * __restrict__ ustokes, double * __restrict__ vstokes, double cithrsh, 
  const double * __restrict__ costh, double delth, const double * __restrict__ dfim_sim, 
  const double * __restrict__ fr, double g, int licerun, int lwamrsetci, int nang, 
  int nfre, int nfre_odd, const double * __restrict__ sinth, double zpi, int ichnk, 
  int nchnk, int ij, double * __restrict__ stfac);
