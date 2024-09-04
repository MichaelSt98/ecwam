

__device__ void frcutindex_c(int kijs, int kijl, const double * __restrict__ fm, 
  const double * __restrict__ fmws, const double * __restrict__ ufric, 
  const double * __restrict__ cicover, int * __restrict__ mij, 
  double * __restrict__ rhowgdfth, double cithrsh_tail, double epsmin, 
  double flogsprdm1, const double * __restrict__ fr, double fric, double g, int nfre, 
  const double * __restrict__ rhowg_dfim, double tailfactor, double tailfactor_pm, 
  const double * __restrict__ zpifr, int ichnk, int nchnk, int ij);
