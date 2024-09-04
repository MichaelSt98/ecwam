

__device__ void tau_phi_hf_c(int kijs, int kijl, const int * __restrict__ mij, 
  int ltauwshelter, const double * __restrict__ ufric, const double * __restrict__ z0m, 
  const double * __restrict__ fl1, const double * __restrict__ aird, 
  const double * __restrict__ rnfac, const double * __restrict__ coswdif, 
  const double * __restrict__ sinwdif2, double * __restrict__ ust, 
  double * __restrict__ tauhf, double * __restrict__ phihf, int llphihf, double delth, 
  const double * __restrict__ fr5, double g, double gamnconst, double gm1, 
  int jtot_tauhf, int llgcbz0, int llnormagam, int nang, int nfre, int nwav_gc, 
  const double * __restrict__ omega_gc, double sqrtgosurft, double tauwshelter, 
  const double * __restrict__ wtauhf, double x0tauhf, double xkappa, 
  const double * __restrict__ xkm_gc, const double * __restrict__ xk_gc, 
  double xlogkratiom1_gc, double zalp, double zpi4gm1, double zpi4gm2, 
  const double * __restrict__ zpifr, int ichnk, int nchnk, int ij, 
  int * __restrict__ ns, double * __restrict__ xks, double * __restrict__ oms, 
  double * __restrict__ sqrtz0og, double * __restrict__ zsup, 
  double * __restrict__ zinf, double * __restrict__ delz, double * __restrict__ taul, 
  double * __restrict__ xloggz0, double * __restrict__ sqrtgz0, 
  double * __restrict__ ustph, double * __restrict__ const1, 
  double * __restrict__ const2, double * __restrict__ consttau, 
  double * __restrict__ constphi, double * __restrict__ f1dcos2, 
  double * __restrict__ f1dcos3, double * __restrict__ f1d, double * __restrict__ f1dsin2
  );
