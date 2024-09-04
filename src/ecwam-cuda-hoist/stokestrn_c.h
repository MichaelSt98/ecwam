

__device__ void stokestrn_c(int kijs, int kijl, const double * __restrict__ fl1, 
  const double * __restrict__ wavnum, const double * __restrict__ stokfac, 
  const double * __restrict__ depth, const double * __restrict__ wswave, 
  const double * __restrict__ wdwave, const double * __restrict__ cicover, 
  const double * __restrict__ cithick, double * __restrict__ ustokes, 
  double * __restrict__ vstokes, double * __restrict__ strnms, 
  double * __restrict__ nemoustokes, double * __restrict__ nemovstokes, 
  double * __restrict__ nemostrn, double cithrsh, const double * __restrict__ costh, 
  double delth, const double * __restrict__ dfim, const double * __restrict__ dfim_sim, 
  double flmin, const double * __restrict__ fr, double g, int licerun, int lwamrsetci, 
  int lwcou, int lwnemocou, int lwnemocousend, int lwnemocoustk, int lwnemocoustrn, 
  int nang, int nfre, int nfre_odd, double rowater, const double * __restrict__ sinth, 
  double zpi, int ichnk, int nchnk, int ij, double * __restrict__ cimsstrn_xki, 
  double * __restrict__ cimsstrn_e, double * __restrict__ cimsstrn_sume, 
  double * __restrict__ stokesdrift_stfac);
