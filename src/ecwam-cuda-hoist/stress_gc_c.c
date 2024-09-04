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
  const double * __restrict__ c2osqrtvg_gc, const double * __restrict__ cm_gc, 
  const double * __restrict__ delkcc_gc_ns, 
  const double * __restrict__ delkcc_omxkm3_gc, double epsus, int llnormagam, 
  int nwav_gc, const double * __restrict__ om3gmkm_gc, 
  const double * __restrict__ omxkm3_gc, double rn1_rn, double sqrtgosurft, 
  double xkappa, const double * __restrict__ xkmsqrtvgoc2_gc, 
  const double * __restrict__ xkm_gc, const double * __restrict__ xk_gc, 
  double xlogkratiom1_gc, double zalp, double * __restrict__ gam_w) {
  
  
  
  //----------------------------------------------------------------------
  
  
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
  double x;
  double xlog;
  double zlog;
  double zlog2x;
  double const_var;
  double zn;
  double gamnorma;  // RENORMALISATION FACTOR OF THE GROWTH RATE
  double zhook_handle;
  
  //     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS

  // <Pragma:: acc routine seq>
  // ----------------------------------------------------------------------
  
  
  //*    1.0  DETERMINE GRAV_CAP SPECTRUM, TAUWHF.
  //          ------------------------------------
  
  //     FIND NS:
  ns = ns_gc_c(ustar, nwav_gc, sqrtgosurft, xkm_gc, xlogkratiom1_gc);
  
  tauwcg_min = pow((ustar*(z0min / z0)), 2);
  
  xlambda = (double) 1.0 + xlama*tanh(xlamb*(pow(ustar, nlam)));
  
  zabhrc = ang_gc*betamaxoxkappa2*halp*c2osqrtvg_gc[ns - 1];
  if (llnormagam) {
    const_var = 
      rnfac*bmaxokap*halp*c2osqrtvg_gc[ns - 1] / fmax((double) (ustar), (double) (epsus))
      ;
  } else {
    const_var = (double) 0.0;
  }
  
  for (i = ns; i <= nwav_gc; i += 1) {
    //       GROWTHRATE BY WIND WITHOUT the multiplicative factor representing the ratio of air density to water density (eps)
    //       and BETAMAXOXKAPPA2
    x = ustar*cm_gc[i - 1];
    xlog = log(xk_gc[i - 1]*z0) + xkappa / (x + zalp);
    zlog = xlog - log(xlambda);
    zlog = fmin((double) (zlog), (double) ((double) 0.0));
    zlog2x = zlog*zlog*x;
    gam_w[i - 1] = zlog2x*zlog2x*exp((double) (xlog))*om3gmkm_gc[i - 1];
  }
  
  zn = const_var*xkmsqrtvgoc2_gc[ns - 1]*gam_w[ns - 1];
  gamnorma = ((double) 1.0 + rn1_rn*zn) / ((double) 1.0 + zn);
  tauwcg = gam_w[ns - 1]*delkcc_gc_ns[ns - 1]*omxkm3_gc[ns - 1]*gamnorma;
  for (i = ns + 1; i <= nwav_gc; i += 1) {
    //       ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
    //       BB = HALP * C2OSQRTVG_GC(NS)*SQRT(VG_GC(I))/C_GC(I)**2
    //       Tauwcg : (rhow * g /rhoa) * integral of (1/c) * gammma * F(k)  k dk
    //       with omega=g*k and omega=k*c,  then
    //       Tauwcg : (rhow /rhoa) * integral of omega * gammma * F(k)  k dk
    //       but gamma is computed wihtout the rhoa/rhow factor so
    //       Tauwcg : integral of omega * gammma_wam * F(k)  k dk
    //       It should be done in vector form with actual directional spreading information
    //       It simplified here by using the ANG_GC factor.
    zn = const_var*xkmsqrtvgoc2_gc[i - 1]*gam_w[i - 1];
    gamnorma = ((double) 1.0 + rn1_rn*zn) / ((double) 1.0 + zn);
    tauwcg = tauwcg + gam_w[i - 1]*delkcc_omxkm3_gc[i - 1]*gamnorma;
  }
  stress_gc = fmax((double) (zabhrc*tauwcg), (double) (tauwcg_min));
  
  
  return stress_gc;
}
