

__device__ double stress_gc_c(double ang_gc, double ustar, double z0, double z0min, 
  double halp, double rnfac, double betamaxoxkappa2, double bmaxokap, 
  const double * __restrict__ c2osqrtvg_gc, const double * __restrict__ cm_gc, 
  const double * __restrict__ delkcc_gc_ns, 
  const double * __restrict__ delkcc_omxkm3_gc, double epsus, int llnormagam, 
  int nwav_gc, const double * __restrict__ om3gmkm_gc, 
  const double * __restrict__ omxkm3_gc, double rn1_rn, double sqrtgosurft, 
  double xkappa, const double * __restrict__ xkmsqrtvgoc2_gc, 
  const double * __restrict__ xkm_gc, const double * __restrict__ xk_gc, 
  double xlogkratiom1_gc, double zalp, double * __restrict__ gam_w);
