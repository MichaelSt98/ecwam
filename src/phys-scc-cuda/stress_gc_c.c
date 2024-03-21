#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "stress_gc_c.h"
#include "ns_gc_c.h"

__device__ double stress_gc_c(double ang_gc, double ustar, double z0, double z0min, 
  double halp, double rnfac, double betamaxoxkappa2, double bmaxokap, 
  const double * c2osqrtvg_gc, const double * cm_gc, const double * delkcc_gc_ns, 
  const double * delkcc_omxkm3_gc, double epsus, int llnormagam, int nwav_gc, 
  const double * om3gmkm_gc, const double * omxkm3_gc, double rn1_rn, 
  double sqrtgosurft, double xkappa, const double * xkmsqrtvgoc2_gc, 
  const double * xkm_gc, const double * xk_gc, double xlogkratiom1_gc, double zalp) {
  
  

  
  double stress_gc;
  int ns;
  int i;
  
  double xlambda;  // Correction factor in the wave growth for gravity-capillary waves
  // XLAMBDA = 1.0_JWRB + XLAMA * TANH(XLAMB * USTAR**NLAM)
  double xlama = (double) 0.25;
  double xlamb = (double) 4.0;
  const int nlam = 4;
  
  double tauwcg_min;
  double tauwcg;
  double zabhrc;
  double x, xlog, zlog, zlog2x;
  double const_var, zn;
  double gamnorma;  // RENORMALISATION FACTOR OF THE GROWTH RATE
  double gam_w;
  // <Pragma:: acc routine seq>
  ns = ns_gc_c(ustar, nwav_gc, sqrtgosurft,  xkm_gc, xlogkratiom1_gc);
  
  tauwcg_min = pow((ustar*(z0min / z0)), 2);
  
  xlambda = (double) 1.0 + xlama*tanh(xlamb*(pow(ustar, nlam)));
  
  zabhrc = ang_gc*betamaxoxkappa2*halp*c2osqrtvg_gc[ns - 1];
  if (llnormagam) {
    const_var = 
      rnfac*bmaxokap*halp*c2osqrtvg_gc[ns - 1] / max((double) (ustar), (double) (epsus));
  } else {
    const_var = (double) 0.0;
  }
  
  for (i = ns; i <= nwav_gc; i += 1) {
    x = ustar*cm_gc[i - 1];
    xlog = log(xk_gc[i - 1]*z0) + xkappa / (x + zalp);
    zlog = xlog - log(xlambda);
    zlog = min((double) (zlog), (double) ((double) 0.0));
    zlog2x = zlog*zlog*x;
  }
  
  gam_w = zlog2x*zlog2x*exp((double) (xlog))*om3gmkm_gc[ns - 1];
  zn = const_var*xkmsqrtvgoc2_gc[ns - 1]*gam_w;
  gamnorma = ((double) 1.0 + rn1_rn*zn) / ((double) 1.0 + zn);
  tauwcg = gam_w*delkcc_gc_ns[ns - 1]*omxkm3_gc[ns - 1]*gamnorma;
  for (i = ns + 1; i <= nwav_gc; i += 1) {
    gam_w = zlog2x*zlog2x*exp((double) (xlog))*om3gmkm_gc[i - 1];
    zn = const_var*xkmsqrtvgoc2_gc[i - 1]*gam_w;
    gamnorma = ((double) 1.0 + rn1_rn*zn) / ((double) 1.0 + zn);
    tauwcg = tauwcg + gam_w*delkcc_omxkm3_gc[i - 1]*gamnorma;
  }
  stress_gc = max((double) (zabhrc*tauwcg), (double) (tauwcg_min));
  return stress_gc;  
  
}
