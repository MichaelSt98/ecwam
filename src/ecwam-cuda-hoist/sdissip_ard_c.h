

__device__ void sdissip_ard_c(int kijs, int kijl, const double * __restrict__ fl1, 
  double * __restrict__ fld, double * __restrict__ sl, 
  const double * __restrict__ wavnum, const double * __restrict__ cgroup, 
  const double * __restrict__ xk2cg, const double * __restrict__ ufric, 
  const double * __restrict__ coswdif, const double * __restrict__ raorw, 
  double brkpbcoef, double delth, double fratio, double g, 
  const int * __restrict__ indicessat, int ipsat, double miche, int nang, int nfre, 
  int nsdsnth, const double * __restrict__ satweights, double sdsbr, double ssdsbrf1, 
  double ssdsc2, double ssdsc3, double ssdsc4, double ssdsc5, double ssdsc6, double zpi, 
  const double * __restrict__ zpifr, int ichnk, int nchnk, int ij, 
  double * __restrict__ facturb, double * __restrict__ facsat, 
  double * __restrict__ facwtrb, double * __restrict__ temp1, 
  double * __restrict__ bth0, double * __restrict__ c_, double * __restrict__ c_c, 
  double * __restrict__ dsip, double * __restrict__ trpz_dsip, 
  double * __restrict__ bth, double * __restrict__ temp2, double * __restrict__ d, 
  double * __restrict__ scumul, double * __restrict__ renewalfreq, 
  double * __restrict__ wcumul);
