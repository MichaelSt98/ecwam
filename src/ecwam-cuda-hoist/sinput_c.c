#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "sinput_c.h"
#include "sinput_jan_c.h"
#include "sinput_ard_c.h"


__device__ void sinput_c(int ngst, int llsneg, int kijs, int kijl, 
  const double * __restrict__ fl1, const double * __restrict__ wavnum, 
  const double * __restrict__ cinv, const double * __restrict__ xk2cg, 
  const double * __restrict__ wdwave, const double * __restrict__ wswave, 
  const double * __restrict__ ufric, const double * __restrict__ z0m, 
  const double * __restrict__ coswdif, const double * __restrict__ sinwdif2, 
  const double * __restrict__ raorw, const double * __restrict__ wstar, 
  const double * __restrict__ rnfac, double * __restrict__ fld, 
  double * __restrict__ sl, double * __restrict__ spos, double * __restrict__ xllws, 
  double abmax, double abmin, double acdlin, double alphamax, double alphamin, 
  double bcdlin, double betamaxoxkappa2, const double * __restrict__ costh, 
  double delth, const double * __restrict__ dfim, double epsmin, double epsus, double g, 
  int iab, int idamping, int iphys, int llgcbz0, int llnormagam, int nang, int nfre, 
  double rnu, double rnum, const double * __restrict__ sinth, double swellf, 
  double swellf2, double swellf3, double swellf4, double swellf5, double swellf6, 
  double swellf7, double swellf7m1, const double * __restrict__ swellft, 
  double tauwshelter, const double * __restrict__ th, double wspmin, double xkappa, 
  double z0rat, double z0tubmax, double zalp, double zpi, 
  const double * __restrict__ zpifr, int ichnk, int nchnk, int ij, 
  double * __restrict__ sinput_ard_constf, double * __restrict__ sinput_ard_const11, 
  double * __restrict__ sinput_ard_const22, double * __restrict__ sinput_ard_z0vis, 
  double * __restrict__ sinput_ard_z0noz, double * __restrict__ sinput_ard_fww, 
  double * __restrict__ sinput_ard_pvisc, double * __restrict__ sinput_ard_pturb, 
  double * __restrict__ sinput_ard_zcn, double * __restrict__ sinput_ard_sig_n, 
  double * __restrict__ sinput_ard_uorbt, double * __restrict__ sinput_ard_aorb, 
  double * __restrict__ sinput_ard_temp, double * __restrict__ sinput_ard_re, 
  double * __restrict__ sinput_ard_re_c, double * __restrict__ sinput_ard_zorb, 
  double * __restrict__ sinput_ard_cnsn, double * __restrict__ sinput_ard_sumf, 
  double * __restrict__ sinput_ard_sumfsin2, double * __restrict__ sinput_ard_cstrnfac, 
  double * __restrict__ sinput_ard_flp_avg, double * __restrict__ sinput_ard_slp_avg, 
  double * __restrict__ sinput_ard_rogoroair, 
  double * __restrict__ sinput_ard_aird_pvisc, double * __restrict__ sinput_ard_xstress, 
  double * __restrict__ sinput_ard_ystress, double * __restrict__ sinput_ard_flp, 
  double * __restrict__ sinput_ard_slp, double * __restrict__ sinput_ard_usg2, 
  double * __restrict__ sinput_ard_taux, double * __restrict__ sinput_ard_tauy, 
  double * __restrict__ sinput_ard_ustp, double * __restrict__ sinput_ard_ustpm1, 
  double * __restrict__ sinput_ard_usdirp, double * __restrict__ sinput_ard_ucn, 
  double * __restrict__ sinput_ard_ucnzalpd, 
  double * __restrict__ sinput_ard_xngamconst, 
  double * __restrict__ sinput_ard_gamnorma, double * __restrict__ sinput_ard_dstab1, 
  double * __restrict__ sinput_ard_temp1, double * __restrict__ sinput_ard_temp2, 
  double * __restrict__ sinput_ard_gam0, double * __restrict__ sinput_ard_dstab, 
  double * __restrict__ sinput_ard_coslp, double * __restrict__ sinput_jan_ztanhkd, 
  double * __restrict__ sinput_jan_sig_n, double * __restrict__ sinput_jan_cnsn, 
  double * __restrict__ sinput_jan_sumf, double * __restrict__ sinput_jan_sumfsin2, 
  double * __restrict__ sinput_jan_cstrnfac, 
  double * __restrict__ sinput_jan_xngamconst, 
  double * __restrict__ sinput_jan_gamnorma, double * __restrict__ sinput_jan_sigdev, 
  double * __restrict__ sinput_jan_us, double * __restrict__ sinput_jan_z0, 
  double * __restrict__ sinput_jan_ucn, double * __restrict__ sinput_jan_zcn, 
  double * __restrict__ sinput_jan_ustpm1, double * __restrict__ sinput_jan_xvd, 
  double * __restrict__ sinput_jan_ucnd, double * __restrict__ sinput_jan_const3_ucn2, 
  double * __restrict__ sinput_jan_ufac1, double * __restrict__ sinput_jan_ufac2, 
  double * __restrict__ sinput_jan_tempd, double * __restrict__ sinput_jan_gam0, 
  int * __restrict__ sinput_jan_lz) {
  
  
  
  // ----------------------------------------------------------------------
  


  
  
  
  
  
  double zhook_handle;
  // ----------------------------------------------------------------------
  
  
  switch (iphys) {
  case 0:
    sinput_jan_c(ngst, llsneg, kijs, kijl, fl1, wavnum, cinv, xk2cg, wswave, ufric, z0m, 
      coswdif, sinwdif2, raorw, wstar, rnfac, fld, sl, spos, xllws, acdlin, alphamax, 
      alphamin, bcdlin, betamaxoxkappa2, delth, epsus, g, idamping, llgcbz0, llnormagam, 
      nang, nfre, rnum, wspmin, xkappa, zalp, zpi, zpifr, ichnk, nchnk, ij, 
      sinput_jan_ztanhkd, sinput_jan_sig_n, sinput_jan_cnsn, sinput_jan_sumf, 
      sinput_jan_sumfsin2, sinput_jan_cstrnfac, sinput_jan_xngamconst, 
      sinput_jan_gamnorma, sinput_jan_sigdev, sinput_jan_us, sinput_jan_z0, 
      sinput_jan_ucn, sinput_jan_zcn, sinput_jan_ustpm1, sinput_jan_xvd, 
      sinput_jan_ucnd, sinput_jan_const3_ucn2, sinput_jan_ufac1, sinput_jan_ufac2, 
      sinput_jan_tempd, sinput_jan_gam0, sinput_jan_lz);
  break;
  case 1:
    sinput_ard_c(ngst, llsneg, kijs, kijl, fl1, wavnum, cinv, xk2cg, wdwave, wswave, 
      ufric, z0m, coswdif, sinwdif2, raorw, wstar, rnfac, fld, sl, spos, xllws, abmax, 
      abmin, acdlin, alphamax, alphamin, bcdlin, betamaxoxkappa2, costh, delth, dfim, 
      epsmin, epsus, g, iab, llgcbz0, llnormagam, nang, nfre, rnu, rnum, sinth, swellf, 
      swellf2, swellf3, swellf4, swellf5, swellf6, swellf7, swellf7m1, swellft, 
      tauwshelter, th, wspmin, xkappa, z0rat, z0tubmax, zalp, zpi, zpifr, ichnk, nchnk, 
      ij, sinput_ard_constf, sinput_ard_const11, sinput_ard_const22, sinput_ard_z0vis, 
      sinput_ard_z0noz, sinput_ard_fww, sinput_ard_pvisc, sinput_ard_pturb, 
      sinput_ard_zcn, sinput_ard_sig_n, sinput_ard_uorbt, sinput_ard_aorb, 
      sinput_ard_temp, sinput_ard_re, sinput_ard_re_c, sinput_ard_zorb, sinput_ard_cnsn, 
      sinput_ard_sumf, sinput_ard_sumfsin2, sinput_ard_cstrnfac, sinput_ard_flp_avg, 
      sinput_ard_slp_avg, sinput_ard_rogoroair, sinput_ard_aird_pvisc, 
      sinput_ard_xstress, sinput_ard_ystress, sinput_ard_flp, sinput_ard_slp, 
      sinput_ard_usg2, sinput_ard_taux, sinput_ard_tauy, sinput_ard_ustp, 
      sinput_ard_ustpm1, sinput_ard_usdirp, sinput_ard_ucn, sinput_ard_ucnzalpd, 
      sinput_ard_xngamconst, sinput_ard_gamnorma, sinput_ard_dstab1, sinput_ard_temp1, 
      sinput_ard_temp2, sinput_ard_gam0, sinput_ard_dstab, sinput_ard_coslp);
  break;
  }
  
  
}
