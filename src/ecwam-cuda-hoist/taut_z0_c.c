#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "taut_z0_c.h"
#include "stress_gc_c.h"
#include "chnkmin_c.h"


__device__ void taut_z0_c(int kijs, int kijl, int iusfg, 
  const double * __restrict__ halp, const double * __restrict__ utop, 
  const double * __restrict__ udir, const double * __restrict__ tauw, 
  const double * __restrict__ tauwdir, const double * __restrict__ rnfac, 
  double * __restrict__ ustar, double * __restrict__ z0, double * __restrict__ z0b, 
  double * __restrict__ chrnck, double acd, double alpha, double alphamax, 
  double alphamin, double ang_gc_a, double ang_gc_b, double ang_gc_c, double bcd, 
  double betamaxoxkappa2, double bmaxokap, const double * __restrict__ c2osqrtvg_gc, 
  double cdmax, double chnkmin_u, const double * __restrict__ cm_gc, 
  const double * __restrict__ delkcc_gc_ns, 
  const double * __restrict__ delkcc_omxkm3_gc, double eps1, double epsmin, 
  double epsus, double g, double gm1, int llcapchnk, int llgcbz0, int llnormagam, 
  int nwav_gc, const double * __restrict__ om3gmkm_gc, 
  const double * __restrict__ omxkm3_gc, double rn1_rn, double rnu, double rnum, 
  double sqrtgosurft, double xkappa, const double * __restrict__ xkmsqrtvgoc2_gc, 
  const double * __restrict__ xkm_gc, const double * __restrict__ xk_gc, 
  double xlogkratiom1_gc, double xnlev, double zalp, int ichnk, int nchnk, int ij, 
  double * __restrict__ alphaog, double * __restrict__ xmin, double * __restrict__ w1, 
  double * __restrict__ tauwact, double * __restrict__ tauweff, 
  double * __restrict__ ang_gc, double * __restrict__ tauunr, 
  int * __restrict__ llcosdiff, double * __restrict__ stress_gc_gam_w) {
  
  // needed for Loki
  
  
  // ----------------------------------------------------------------------
  


  
  
  const int niter = 17;
  
  double twoxmp1 = (double) 3.0;
  
  int iter;
  int ifrph;
  
  // Cd and Z0 from Hersbach 2010, ECMWF Tech Memo (without the viscous part)
  //     CD = ACDLIN + BCDLIN*SQRT(PCHAR) * U10
  double acdlin = (double) 0.0008;
  double bcdlin = (double) 0.00047;
  double alphagm1;
  
  double z0min = (double) 0.000001;
  double pce_gc;
  double z0minrst;
  double charnock_min;
  double cosdiff;
  double zchar;
  double us2totauw;
  double usmax;
  double xlogxl;
  double xkutop;
  double xologz0;
  double ustold;
  double ustnew;
  double tauold;
  double taunew;
  double x;
  double f;
  double delf;
  double cdfg;
  double ustm1;
  double z0tot;
  double z0ch;
  double z0vis;
  double hz0viso1mx;
  double zz;
  double const_var;
  double tauv;
  double del;
  double rnueff;
  double rnukappam1;
  double zhook_handle;
  
  
  // ----------------------------------------------------------------------
  
  //     INLINE FUNCTION.
  //     ----------------
  
  //     Simple empirical fit to model drag coefficient
  double cdm;
  double u10;
  
  // CDM(U10) = MAX(MIN(0.0006_JWRB+0.00008_JWRB*U10, 0.001_JWRB+0.0018_JWRB*EXP(-0.05_JWRB*(U10-33._JWRB))),0.001_JWRB)
  // ----------------------------------------------------------------------
  
  
  xlogxl = log(xnlev);
  us2totauw = (double) 1.0 + eps1;
  
  //     ONLY take the contribution of TAUW that is in the wind direction
  
  cosdiff = cos(udir[ij - 1 + kijl*(ichnk - 1)] - tauwdir[ij - 1 + kijl*(ichnk - 1)]);
  tauwact[ij - 1 + kijl*(ichnk - 1)] = 
    fmax((double) (tauw[ij - 1 + kijl*(ichnk - 1)]*cosdiff), (double) (epsmin));
  llcosdiff[ij - 1 + kijl*(ichnk - 1)] = cosdiff > (double) 0.9;
  
  //  USING THE CG MODEL:
  if (llgcbz0) {
    
    if (llcapchnk) {
      charnock_min = 
        chnkmin_c(utop[ij - 1 + kijl*(ichnk - 1)], alpha, alphamin, chnkmin_u);
      alphaog[ij - 1 + kijl*(ichnk - 1)] = charnock_min*gm1;
    } else {
      alphaog[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
    }
    
    usmax = fmax((double) (-(double) 0.21339 + (double) 0.093698*utop[ij - 1 + 
      kijl*(ichnk - 1)] - (double) 0.0020944*(pow(utop[ij - 1 + kijl*(ichnk - 1)], 2)) + 
      (double) 5.5091E-5*(pow(utop[ij - 1 + kijl*(ichnk - 1)], 3))), (double) ((double) 
      0.03));
    tauweff[ij - 1 + kijl*(ichnk - 1)] = fmin((double) (tauwact[ij - 1 + kijl*(ichnk - 1)
      ]*us2totauw), (double) (pow(usmax, 2)));
    
    rnueff = (double) 0.04*rnu;
    
    rnukappam1 = rnueff / xkappa;
    
    pce_gc = (double) 0.001*iusfg + (1 - iusfg)*(double) 0.005;
    
    if (iusfg == 0) {
      alphagm1 = alpha*gm1;
      if (utop[ij - 1 + kijl*(ichnk - 1)] < (double) 1.0) {
        cdfg = (double) 0.002;
      } else if (llcosdiff[ij - 1 + kijl*(ichnk - 1)]) {
        x = fmin((double) (tauwact[ij - 1 + kijl*(ichnk - 1)] / (pow(fmax((double) 
          (ustar[ij - 1 + kijl*(ichnk - 1)]), (double) (epsus)), 2))), (double) ((double)
           0.99));
        zchar = fmin((double) (alphagm1*(pow(ustar[ij - 1 + kijl*(ichnk - 1)], 2)) / 
          sqrt((double) ((double) 1.0 - x))), (double) ((double) 0.05*exp((double) 
          (-(double) 0.05*(utop[ij - 1 + kijl*(ichnk - 1)] - (double) 35.)))));
        zchar = fmin((double) (zchar), (double) (alphamax));
        cdfg = acdlin + bcdlin*sqrt((double) (zchar))*utop[ij - 1 + kijl*(ichnk - 1)];
      } else {
        // CDFG = CDM(UTOP(IJ))
        cdfg = fmax((double) (fmin((double) ((double) 0.0006 + (double) 0.00008*utop[ij -
           1 + kijl*(ichnk - 1)]), (double) ((double) 0.001 + (double) 0.0018*exp((double
          ) (-(double) 0.05*(utop[ij - 1 + kijl*(ichnk - 1)] - (double) 33.)))))), 
          (double) ((double) 0.001));
      }
      ustar[ij - 1 + kijl*(ichnk - 1)] = 
        utop[ij - 1 + kijl*(ichnk - 1)]*sqrt((double) (cdfg));
    }
    
    w1[ij - 1 + kijl*(ichnk - 1)] = (double) 0.85 - (double) 0.05*(tanh((double) 
      10.0*(utop[ij - 1 + kijl*(ichnk - 1)] - (double) 5.0)) + (double) 1.0);
    
    xkutop = xkappa*utop[ij - 1 + kijl*(ichnk - 1)];
    
    ustold = ustar[ij - 1 + kijl*(ichnk - 1)];
    tauold = pow(ustold, 2);
    
    for (iter = 1; iter <= niter; iter += 1) {
      //         Z0 IS DERIVED FROM THE NEUTRAL LOG PROFILE: UTOP = (USTAR/XKAPPA)*LOG((XNLEV+Z0)/Z0)
      z0[ij - 1 + kijl*(ichnk - 1)] = fmax((double) (xnlev / (exp((double) (fmin((double)
         (xkutop / ustold), (double) ((double) 50.0)))) - (double) 1.0)), (double) (z0min
        ));
      // Viscous kinematic stress nu_air * dU/dz at z=0 of the neutral log profile reduced by factor 25 (0.04)
      tauv = rnukappam1*ustold / z0[ij - 1 + kijl*(ichnk - 1)];
      
      ang_gc[ij - 1 + kijl*(ichnk - 1)] = ang_gc_a + ang_gc_b*tanh(ang_gc_c*tauold);
      
      tauunr[ij - 1 + kijl*(ichnk - 1)] = stress_gc_c(ang_gc[ij - 1 + kijl*(ichnk - 1)], 
        ustar[ij - 1 + kijl*(ichnk - 1)], z0[ij - 1 + kijl*(ichnk - 1)], z0min, halp[ij -
         1], rnfac[ij - 1], betamaxoxkappa2, bmaxokap, c2osqrtvg_gc, cm_gc, delkcc_gc_ns,
         delkcc_omxkm3_gc, epsus, llnormagam, nwav_gc, om3gmkm_gc, omxkm3_gc, rn1_rn, 
        sqrtgosurft, xkappa, xkmsqrtvgoc2_gc, xkm_gc, xk_gc, xlogkratiom1_gc, zalp, 
        stress_gc_gam_w);
      
      //         TOTAL kinematic STRESS:
      taunew = 
        tauweff[ij - 1 + kijl*(ichnk - 1)] + tauv + tauunr[ij - 1 + kijl*(ichnk - 1)];
      ustnew = sqrt((double) (taunew));
      ustar[ij - 1 + kijl*(ichnk - 1)] = w1[ij - 1 + kijl*(ichnk - 1)]*ustold + ((double)
         1.0 - w1[ij - 1 + kijl*(ichnk - 1)])*ustnew;
      
      //         CONVERGENCE ?
      del = ustar[ij - 1 + kijl*(ichnk - 1)] - ustold;
      // IF (ABS(DEL) < PCE_GC*USTAR(IJ)) EXIT
      tauold = pow(ustar[ij - 1 + kijl*(ichnk - 1)], 2);
      ustold = ustar[ij - 1 + kijl*(ichnk - 1)];
    }
    // protection just in case there is no convergence
    if (iter > niter) {
      // CDFG = CDM(UTOP(IJ))
      cdfg = fmax((double) (fmin((double) ((double) 0.0006 + (double) 0.00008*utop[ij - 1
         + kijl*(ichnk - 1)]), (double) ((double) 0.001 + (double) 0.0018*exp((double) 
        (-(double) 0.05*(utop[ij - 1 + kijl*(ichnk - 1)] - (double) 33.)))))), (double) 
        ((double) 0.001));
      ustar[ij - 1 + kijl*(ichnk - 1)] = 
        utop[ij - 1 + kijl*(ichnk - 1)]*sqrt((double) (cdfg));
      z0minrst = (pow(ustar[ij - 1 + kijl*(ichnk - 1)], 2))*alpha*gm1;
      z0[ij - 1 + kijl*(ichnk - 1)] = fmax((double) (xnlev / (exp((double) (xkutop / 
        ustar[ij - 1 + kijl*(ichnk - 1)])) - (double) 1.0)), (double) (z0minrst));
      z0b[ij - 1 + kijl*(ichnk - 1)] = z0minrst;
    } else {
      z0[ij - 1 + kijl*(ichnk - 1)] = fmax((double) (xnlev / (exp((double) (xkutop / 
        ustar[ij - 1 + kijl*(ichnk - 1)])) - (double) 1.0)), (double) (z0min));
      z0b[ij - 1 + kijl*(ichnk - 1)] = z0[ij - 1 + kijl*(ichnk - 1)]*sqrt((double) 
        (tauunr[ij - 1 + kijl*(ichnk - 1)] / tauold));
    }
    
    //       Refine solution
    x = tauweff[ij - 1 + kijl*(ichnk - 1)] / tauold;
    
    if (x < (double) 0.99) {
      ustold = ustar[ij - 1 + kijl*(ichnk - 1)];
      tauold = 
        fmax((double) (pow(ustold, 2)), (double) (tauweff[ij - 1 + kijl*(ichnk - 1)]));
      
      for (iter = 1; iter <= niter; iter += 1) {
        x = fmin((double) (tauweff[ij - 1 + kijl*(ichnk - 1)] / tauold), (double) 
          ((double) 0.99));
        ustm1 = (double) 1.0 / fmax((double) (ustold), (double) (epsus));
        //!!! Limit how small z0 could become
        //!!! This is a bit of a compromise to limit very low Charnock for intermediate high winds (15 -25 m/s)
        //!!! It is not ideal !!!
        z0[ij - 1 + kijl*(ichnk - 1)] = fmax((double) (xnlev / (exp((double) 
          (fmin((double) (xkutop / ustold), (double) ((double) 50.0)))) - (double) 1.0)),
           (double) (z0min));
        
        tauunr[ij - 1 + kijl*(ichnk - 1)] = stress_gc_c(ang_gc[ij - 1 + kijl*(ichnk - 1)
          ], ustold, z0[ij - 1 + kijl*(ichnk - 1)], z0min, halp[ij - 1], rnfac[ij - 1], 
          betamaxoxkappa2, bmaxokap, c2osqrtvg_gc, cm_gc, delkcc_gc_ns, delkcc_omxkm3_gc,
           epsus, llnormagam, nwav_gc, om3gmkm_gc, omxkm3_gc, rn1_rn, sqrtgosurft, 
          xkappa, xkmsqrtvgoc2_gc, xkm_gc, xk_gc, xlogkratiom1_gc, zalp, stress_gc_gam_w)
          ;
        
        z0b[ij - 1 + kijl*(ichnk - 1)] = fmax((double) (z0[ij - 1 + kijl*(ichnk - 1)
          ]*sqrt((double) (tauunr[ij - 1 + kijl*(ichnk - 1)] / tauold))), (double) 
          (alphaog[ij - 1 + kijl*(ichnk - 1)]*tauold));
        z0vis = rnum*ustm1;
        hz0viso1mx = (double) 0.5*z0vis / ((double) 1.0 - x);
        z0[ij - 1 + kijl*(ichnk - 1)] = hz0viso1mx + sqrt((double) ((pow(hz0viso1mx, 2)) 
          + (pow(z0b[ij - 1 + kijl*(ichnk - 1)], 2)) / ((double) 1.0 - x)));
        
        xologz0 = (double) 1.0 / (xlogxl - log(z0[ij - 1 + kijl*(ichnk - 1)]));
        f = ustold - xkutop*xologz0;
        zz = (double) 2.0*ustm1*((double) 3.0*(pow(z0b[ij - 1 + kijl*(ichnk - 1)], 2)) + 
          (double) 0.5*z0vis*z0[ij - 1 + kijl*(ichnk - 1)] - (pow(z0[ij - 1 + kijl*(ichnk
           - 1)], 2))) / ((double) 2.0*(pow(z0[ij - 1 + kijl*(ichnk - 1)], 2))*((double) 
          1.0 - x) - z0vis*z0[ij - 1 + kijl*(ichnk - 1)]);
        
        delf = (double) 1.0 - xkutop*(pow(xologz0, 2))*zz;
        if (delf != (double) 0.0) {
          ustar[ij - 1 + kijl*(ichnk - 1)] = ustold - f / delf;
        }
        
        //           CONVERGENCE ?
        del = ustar[ij - 1 + kijl*(ichnk - 1)] - ustold;
        
        // IF (ABS(DEL) < PCE_GC*USTAR(IJ)) EXIT
        ustold = ustar[ij - 1 + kijl*(ichnk - 1)];
        tauold = 
          fmax((double) (pow(ustold, 2)), (double) (tauweff[ij - 1 + kijl*(ichnk - 1)]));
      }
      // protection just in case there is no convergence
      if (iter > niter) {
        // CDFG = CDM(UTOP(IJ))
        cdfg = fmax((double) (fmin((double) ((double) 0.0006 + (double) 0.00008*utop[ij -
           1 + kijl*(ichnk - 1)]), (double) ((double) 0.001 + (double) 0.0018*exp((double
          ) (-(double) 0.05*(utop[ij - 1 + kijl*(ichnk - 1)] - (double) 33.)))))), 
          (double) ((double) 0.001));
        ustar[ij - 1 + kijl*(ichnk - 1)] = 
          utop[ij - 1 + kijl*(ichnk - 1)]*sqrt((double) (cdfg));
        z0minrst = (pow(ustar[ij - 1 + kijl*(ichnk - 1)], 2))*alpha*gm1;
        z0[ij - 1 + kijl*(ichnk - 1)] = fmax((double) (xnlev / (exp((double) (xkutop / 
          ustar[ij - 1 + kijl*(ichnk - 1)])) - (double) 1.0)), (double) (z0minrst));
        z0b[ij - 1 + kijl*(ichnk - 1)] = z0minrst;
        chrnck[ij - 1 + kijl*(ichnk - 1)] = fmax((double) (g*z0[ij - 1 + kijl*(ichnk - 1)
          ] / (pow(ustar[ij - 1 + kijl*(ichnk - 1)], 2))), (double) (alphamin));
      } else {
        chrnck[ij - 1 + kijl*(ichnk - 1)] = fmax((double) (g*(z0b[ij - 1 + kijl*(ichnk - 
          1)] / sqrt((double) ((double) 1.0 - x))) / (pow(fmax((double) (ustar[ij - 1 + 
          kijl*(ichnk - 1)]), (double) (epsus)), 2))), (double) (alphamin));
      }
      
    } else {
      ustm1 = (double) 1.0 / fmax((double) (ustar[ij - 1 + kijl*(ichnk - 1)]), (double) 
        (epsus));
      z0vis = rnum*ustm1;
      chrnck[ij - 1 + kijl*(ichnk - 1)] = fmax((double) (g*(z0[ij - 1 + kijl*(ichnk - 1)]
         - z0vis)*(pow(ustm1, 2))), (double) (alphamin));
    }
    
    
    
  } else {
    
    tauweff[ij - 1 + kijl*(ichnk - 1)] = tauwact[ij - 1 + kijl*(ichnk - 1)]*us2totauw;
    
    if (llcapchnk) {
      charnock_min = 
        chnkmin_c(utop[ij - 1 + kijl*(ichnk - 1)], alpha, alphamin, chnkmin_u);
      xmin[ij - 1 + kijl*(ichnk - 1)] = (double) 0.15*(alpha - charnock_min);
      alphaog[ij - 1 + kijl*(ichnk - 1)] = charnock_min*gm1;
    } else {
      xmin[ij - 1 + kijl*(ichnk - 1)] = (double) 0.0;
      alphaog[ij - 1 + kijl*(ichnk - 1)] = alpha*gm1;
    }
    
    xkutop = xkappa*utop[ij - 1 + kijl*(ichnk - 1)];
    
    ustold = (1 - iusfg)*utop[ij - 1 + kijl*(ichnk - 1)]*sqrt((double) (fmin((double) 
      (acd + bcd*utop[ij - 1 + kijl*(ichnk - 1)]), (double) (cdmax)))) + iusfg*ustar[ij -
       1 + kijl*(ichnk - 1)];
    tauold = 
      fmax((double) (pow(ustold, 2)), (double) (tauweff[ij - 1 + kijl*(ichnk - 1)]));
    ustar[ij - 1 + kijl*(ichnk - 1)] = sqrt((double) (tauold));
    ustm1 = 
      (double) 1.0 / fmax((double) (ustar[ij - 1 + kijl*(ichnk - 1)]), (double) (epsus));
    
    for (iter = 1; iter <= niter; iter += 1) {
      x = fmax((double) (tauwact[ij - 1 + kijl*(ichnk - 1)] / tauold), (double) (xmin[ij 
        - 1 + kijl*(ichnk - 1)]));
      z0ch = 
        alphaog[ij - 1 + kijl*(ichnk - 1)]*tauold / sqrt((double) ((double) 1.0 - x));
      z0vis = rnum*ustm1;
      z0tot = z0ch + z0vis;
      
      xologz0 = (double) 1.0 / (xlogxl - log(z0tot));
      f = ustar[ij - 1 + kijl*(ichnk - 1)] - xkutop*xologz0;
      zz = ustm1*(z0ch*((double) 2.0 - twoxmp1*x) / ((double) 1.0 - x) - z0vis) / z0tot;
      delf = (double) 1.0 - xkutop*(pow(xologz0, 2))*zz;
      
      if (delf != (double) 0.0) {
        ustar[ij - 1 + kijl*(ichnk - 1)] = ustar[ij - 1 + kijl*(ichnk - 1)] - f / delf;
      }
      taunew = fmax((double) (pow(ustar[ij - 1 + kijl*(ichnk - 1)], 2)), (double) 
        (tauweff[ij - 1 + kijl*(ichnk - 1)]));
      ustar[ij - 1 + kijl*(ichnk - 1)] = sqrt((double) (taunew));
      // IF (TAUNEW == TAUOLD) EXIT
      ustm1 = (double) 1.0 / fmax((double) (ustar[ij - 1 + kijl*(ichnk - 1)]), (double) 
        (epsus));
      tauold = taunew;
    }
    
    z0[ij - 1 + kijl*(ichnk - 1)] = z0ch;
    z0b[ij - 1 + kijl*(ichnk - 1)] = alphaog[ij - 1 + kijl*(ichnk - 1)]*tauold;
    chrnck[ij - 1 + kijl*(ichnk - 1)] = fmax((double) (g*z0[ij - 1 + kijl*(ichnk - 1)
      ]*(pow(ustm1, 2))), (double) (alphamin));
    
    
  }
  
  
  
}
