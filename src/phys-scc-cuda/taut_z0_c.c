#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "taut_z0_c.h"
#include "stress_gc_c.h"
#include "chnkmin_c.h"

__device__ void taut_z0_c(int kijs, int kijl, int iusfg, const double * halp, 
  const double * utop, const double * udir, const double * tauw, const double * tauwdir, 
  const double * rnfac, double * ustar, double * z0, double * z0b, double * chrnck, 
  double acd, double alpha, double alphamax, double alphamin, double ang_gc_a, 
  double ang_gc_b, double ang_gc_c, double bcd, double betamaxoxkappa2, double bmaxokap, 
  const double * c2osqrtvg_gc, double cdmax, double chnkmin_u, const double * cm_gc, 
  const double * delkcc_gc_ns, const double * delkcc_omxkm3_gc, double eps1, 
  double epsmin, double epsus, double g, double gm1, int llcapchnk, int llgcbz0, 
  int llnormagam, int nwav_gc, const double * om3gmkm_gc, const double * omxkm3_gc, 
  double rn1_rn, double rnu, double rnum, double sqrtgosurft, double xkappa, 
  const double * xkmsqrtvgoc2_gc, const double * xkm_gc, const double * xk_gc, 
  double xlogkratiom1_gc, double xnlev, double zalp, int ichnk, int nchnk, int ij) {
  
  


  const int niter = 17;
  
  double twoxmp1 = (double) 3.0;
  
  int iter;
  int ifrph;
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
  double alphaog;
  double xmin;
  double w1;
  double tauwact;
  double tauweff;
  double ang_gc;
  double tauunr;
  
  int llcosdiff;
  
  xlogxl = log(xnlev);
  us2totauw = (double) 1.0 + eps1;

  cosdiff = cos(udir[ij - 1 + kijl*(ichnk - 1)] - tauwdir[ij - 1 + kijl*(ichnk - 1)]);
  tauwact = max((double) (tauw[ij - 1 + kijl*(ichnk - 1)]*cosdiff), (double) (epsmin));
  llcosdiff = cosdiff > (double) 0.9;
  if (llgcbz0) {
    
    if (llcapchnk) {
      charnock_min = 
        chnkmin_c(utop[ij - 1 + kijl*(ichnk - 1)], alpha, alphamin, chnkmin_u);
      alphaog = charnock_min*gm1;
    } else {
      alphaog = (double) 0.0;
    }
    
    usmax = max((double) (-(double) 0.21339 + (double) 0.093698*utop[ij - 1 + kijl*(ichnk
       - 1)] - (double) 0.0020944*(pow(utop[ij - 1 + kijl*(ichnk - 1)], 2)) + (double) 
      5.5091E-5*(pow(utop[ij - 1 + kijl*(ichnk - 1)], 3))), (double) ((double) 0.03));
    tauweff = min((double) (tauwact*us2totauw), (double) (pow(usmax, 2)));
    
    rnueff = (double) 0.04*rnu;
    
    rnukappam1 = rnueff / xkappa;
    
    pce_gc = (double) 0.001*iusfg + (1 - iusfg)*(double) 0.005;
    
    if (iusfg == 0) {
      alphagm1 = alpha*gm1;
      if (utop[ij - 1 + kijl*(ichnk - 1)] < (double) 1.0) {
        cdfg = (double) 0.002;
      } else if (llcosdiff) {
        x = min((double) (tauwact / (pow(max((double) (ustar[ij - 1 + kijl*(ichnk - 1)]),
           (double) (epsus)), 2))), (double) ((double) 0.99));
        zchar = min((double) (alphagm1*(pow(ustar[ij - 1 + kijl*(ichnk - 1)], 2)) / 
          sqrt((double) ((double) 1.0 - x))), (double) ((double) 0.05*exp((double) 
          (-(double) 0.05*(utop[ij - 1 + kijl*(ichnk - 1)] - (double) 35.)))));
        zchar = min((double) (zchar), (double) (alphamax));
        cdfg = acdlin + bcdlin*sqrt((double) (zchar))*utop[ij - 1 + kijl*(ichnk - 1)];
      } else {
        // CDFG = CDM(UTOP(IJ)) ! TODO: revert and automate
        cdfg = max((double) (min((double) ((double) 0.0006 + (double) 0.00008*utop[ij - 1
           + kijl*(ichnk - 1)]), (double) ((double) 0.001 + (double) 0.0018*exp((double) 
          (-(double) 0.05*(utop[ij - 1 + kijl*(ichnk - 1)] - (double) 33.)))))), (double)
           ((double) 0.001));
      }
      ustar[ij - 1 + kijl*(ichnk - 1)] = 
        utop[ij - 1 + kijl*(ichnk - 1)]*sqrt((double) (cdfg));
    }
    
    w1 = (double) 0.85 - (double) 0.05*(tanh((double) 10.0*(utop[ij - 1 + kijl*(ichnk - 1
      )] - (double) 5.0)) + (double) 1.0);
    
    xkutop = xkappa*utop[ij - 1 + kijl*(ichnk - 1)];
    
    ustold = ustar[ij - 1 + kijl*(ichnk - 1)];
    tauold = pow(ustold, 2);
    
    for (iter = 1; iter <= niter; iter += 1) {
      //         Z0 IS DERIVED FROM THE NEUTRAL LOG PROFILE: UTOP = (USTAR/XKAPPA)*LOG((XNLEV+Z0)/Z0)
      z0[ij - 1 + kijl*(ichnk - 1)] = max((double) (xnlev / (exp((double) (min((double) 
        (xkutop / ustold), (double) ((double) 50.0)))) - (double) 1.0)), (double) (z0min)
        );
      // Viscous kinematic stress nu_air * dU/dz at z=0 of the neutral log profile reduced by factor 25 (0.04)
      tauv = rnukappam1*ustold / z0[ij - 1 + kijl*(ichnk - 1)];
      
      ang_gc = ang_gc_a + ang_gc_b*tanh(ang_gc_c*tauold);
      
      tauunr = stress_gc_c(ang_gc, ustar[ij - 1 + kijl*(ichnk - 1)], z0[ij - 1 + 
        kijl*(ichnk - 1)], z0min, halp[ij - 1], rnfac[ij - 1], betamaxoxkappa2, bmaxokap,
          c2osqrtvg_gc,  cm_gc,  delkcc_gc_ns,  delkcc_omxkm3_gc, epsus, 
        llnormagam, nwav_gc,  om3gmkm_gc,  omxkm3_gc, rn1_rn, sqrtgosurft, xkappa, 
         xkmsqrtvgoc2_gc,  xkm_gc,  xk_gc, xlogkratiom1_gc, zalp);
      taunew = tauweff + tauv + tauunr;
      ustnew = sqrt((double) (taunew));
      ustar[ij - 1 + kijl*(ichnk - 1)] = w1*ustold + ((double) 1.0 - w1)*ustnew;
      del = ustar[ij - 1 + kijl*(ichnk - 1)] - ustold;
      if (abs((double) (del)) < pce_gc*ustar[ij - 1 + kijl*(ichnk - 1)]) {
        // EXIT
      }
      tauold = pow(ustar[ij - 1 + kijl*(ichnk - 1)], 2);
      ustold = ustar[ij - 1 + kijl*(ichnk - 1)];
    }
    // protection just in case there is no convergence
    if (iter > niter) {
      // CDFG = CDM(UTOP(IJ))
      cdfg = max((double) (min((double) ((double) 0.0006 + (double) 0.00008*utop[ij - 1 +
         kijl*(ichnk - 1)]), (double) ((double) 0.001 + (double) 0.0018*exp((double) 
        (-(double) 0.05*(utop[ij - 1 + kijl*(ichnk - 1)] - (double) 33.)))))), (double) 
        ((double) 0.001));
      ustar[ij - 1 + kijl*(ichnk - 1)] = 
        utop[ij - 1 + kijl*(ichnk - 1)]*sqrt((double) (cdfg));
      z0minrst = (pow(ustar[ij - 1 + kijl*(ichnk - 1)], 2))*alpha*gm1;
      z0[ij - 1 + kijl*(ichnk - 1)] = max((double) (xnlev / (exp((double) (xkutop / 
        ustar[ij - 1 + kijl*(ichnk - 1)])) - (double) 1.0)), (double) (z0minrst));
      z0b[ij - 1 + kijl*(ichnk - 1)] = z0minrst;
    } else {
      z0[ij - 1 + kijl*(ichnk - 1)] = max((double) (xnlev / (exp((double) (xkutop / 
        ustar[ij - 1 + kijl*(ichnk - 1)])) - (double) 1.0)), (double) (z0min));
      z0b[ij - 1 + kijl*(ichnk - 1)] = 
        z0[ij - 1 + kijl*(ichnk - 1)]*sqrt((double) (tauunr / tauold));
    }
    x = tauweff / tauold;
    
    if (x < (double) 0.99) {
      ustold = ustar[ij - 1 + kijl*(ichnk - 1)];
      tauold = max((double) (pow(ustold, 2)), (double) (tauweff));
      
      for (iter = 1; iter <= niter; iter += 1) {
        x = min((double) (tauweff / tauold), (double) ((double) 0.99));
        ustm1 = (double) 1.0 / max((double) (ustold), (double) (epsus));
        z0[ij - 1 + kijl*(ichnk - 1)] = max((double) (xnlev / (exp((double) (min((double)
           (xkutop / ustold), (double) ((double) 50.0)))) - (double) 1.0)), (double) 
          (z0min));
        
        tauunr = stress_gc_c(ang_gc, ustold, z0[ij - 1 + kijl*(ichnk - 1)], z0min, 
          halp[ij - 1], rnfac[ij - 1], betamaxoxkappa2, bmaxokap,  c2osqrtvg_gc,  
          cm_gc,  delkcc_gc_ns,  delkcc_omxkm3_gc, epsus, llnormagam, nwav_gc,  
          om3gmkm_gc,  omxkm3_gc, rn1_rn, sqrtgosurft, xkappa,  xkmsqrtvgoc2_gc,
            xkm_gc,  xk_gc, xlogkratiom1_gc, zalp);
        
        z0b[ij - 1 + kijl*(ichnk - 1)] = max((double) (z0[ij - 1 + kijl*(ichnk - 1)
          ]*sqrt((double) (tauunr / tauold))), (double) (alphaog*tauold));
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
        del = ustar[ij - 1 + kijl*(ichnk - 1)] - ustold;
        
        if (abs((double) (del)) < pce_gc*ustar[ij - 1 + kijl*(ichnk - 1)]) {
          // EXIT
        }
        ustold = ustar[ij - 1 + kijl*(ichnk - 1)];
        tauold = max((double) (pow(ustold, 2)), (double) (tauweff));
      }
      // protection just in case there is no convergence
      if (iter > niter) {
        // CDFG = CDM(UTOP(IJ))
        cdfg = max((double) (min((double) ((double) 0.0006 + (double) 0.00008*utop[ij - 1
           + kijl*(ichnk - 1)]), (double) ((double) 0.001 + (double) 0.0018*exp((double) 
          (-(double) 0.05*(utop[ij - 1 + kijl*(ichnk - 1)] - (double) 33.)))))), (double)
           ((double) 0.001));
        ustar[ij - 1 + kijl*(ichnk - 1)] = 
          utop[ij - 1 + kijl*(ichnk - 1)]*sqrt((double) (cdfg));
        z0minrst = (pow(ustar[ij - 1 + kijl*(ichnk - 1)], 2))*alpha*gm1;
        z0[ij - 1 + kijl*(ichnk - 1)] = max((double) (xnlev / (exp((double) (xkutop / 
          ustar[ij - 1 + kijl*(ichnk - 1)])) - (double) 1.0)), (double) (z0minrst));
        z0b[ij - 1 + kijl*(ichnk - 1)] = z0minrst;
        chrnck[ij - 1 + kijl*(ichnk - 1)] = max((double) (g*z0[ij - 1 + kijl*(ichnk - 1)]
           / (pow(ustar[ij - 1 + kijl*(ichnk - 1)], 2))), (double) (alphamin));
      } else {
        chrnck[ij - 1 + kijl*(ichnk - 1)] = max((double) (g*(z0b[ij - 1 + kijl*(ichnk - 1
          )] / sqrt((double) ((double) 1.0 - x))) / (pow(max((double) (ustar[ij - 1 + 
          kijl*(ichnk - 1)]), (double) (epsus)), 2))), (double) (alphamin));
      }
      
    } else {
      ustm1 = 
        (double) 1.0 / max((double) (ustar[ij - 1 + kijl*(ichnk - 1)]), (double) (epsus))
        ;
      z0vis = rnum*ustm1;
      chrnck[ij - 1 + kijl*(ichnk - 1)] = max((double) (g*(z0[ij - 1 + kijl*(ichnk - 1)] 
        - z0vis)*(pow(ustm1, 2))), (double) (alphamin));
    }
    
  } else {
    
    tauweff = tauwact*us2totauw;
    
    if (llcapchnk) {
      charnock_min = 
        chnkmin_c(utop[ij - 1 + kijl*(ichnk - 1)], alpha, alphamin, chnkmin_u);
      xmin = (double) 0.15*(alpha - charnock_min);
      alphaog = charnock_min*gm1;
    } else {
      xmin = (double) 0.0;
      alphaog = alpha*gm1;
    }
    
    xkutop = xkappa*utop[ij - 1 + kijl*(ichnk - 1)];
    
    ustold = (1 - iusfg)*utop[ij - 1 + kijl*(ichnk - 1)]*sqrt((double) (min((double) (acd
       + bcd*utop[ij - 1 + kijl*(ichnk - 1)]), (double) (cdmax)))) + iusfg*ustar[ij - 1 +
       kijl*(ichnk - 1)];
    tauold = max((double) (pow(ustold, 2)), (double) (tauweff));
    ustar[ij - 1 + kijl*(ichnk - 1)] = sqrt((double) (tauold));
    ustm1 = 
      (double) 1.0 / max((double) (ustar[ij - 1 + kijl*(ichnk - 1)]), (double) (epsus));
    
    for (iter = 1; iter <= niter; iter += 1) {
      x = max((double) (tauwact / tauold), (double) (xmin));
      z0ch = alphaog*tauold / sqrt((double) ((double) 1.0 - x));
      z0vis = rnum*ustm1;
      z0tot = z0ch + z0vis;
      
      xologz0 = (double) 1.0 / (xlogxl - log(z0tot));
      f = ustar[ij - 1 + kijl*(ichnk - 1)] - xkutop*xologz0;
      zz = ustm1*(z0ch*((double) 2.0 - twoxmp1*x) / ((double) 1.0 - x) - z0vis) / z0tot;
      delf = (double) 1.0 - xkutop*(pow(xologz0, 2))*zz;
      
      if (delf != (double) 0.0) {
        ustar[ij - 1 + kijl*(ichnk - 1)] = ustar[ij - 1 + kijl*(ichnk - 1)] - f / delf;
      }
      taunew = 
        max((double) (pow(ustar[ij - 1 + kijl*(ichnk - 1)], 2)), (double) (tauweff));
      ustar[ij - 1 + kijl*(ichnk - 1)] = sqrt((double) (taunew));
      if (taunew == tauold) {
        // EXIT
      }
      ustm1 = 
        (double) 1.0 / max((double) (ustar[ij - 1 + kijl*(ichnk - 1)]), (double) (epsus))
        ;
      tauold = taunew;
    }
    
    z0[ij - 1 + kijl*(ichnk - 1)] = z0ch;
    z0b[ij - 1 + kijl*(ichnk - 1)] = alphaog*tauold;
    chrnck[ij - 1 + kijl*(ichnk - 1)] = max((double) (g*z0[ij - 1 + kijl*(ichnk - 1)
      ]*(pow(ustm1, 2))), (double) (alphamin));
    
    
  }

  
}
