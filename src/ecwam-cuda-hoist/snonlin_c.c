#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "snonlin_c.h"
#include "transf_snl_c.h"
#include "transf_c.h"
#include "peak_ang_c.h"


__device__ void snonlin_c(int kijs, int kijl, const double * __restrict__ fl1, 
  double * __restrict__ fld, double * __restrict__ sl, 
  const double * __restrict__ wavnum, const double * __restrict__ depth, 
  const double * __restrict__ akmean, const double * __restrict__ af11, double bathymax, 
  const double * __restrict__ costh, double dal1, double dal2, double delth, 
  const double * __restrict__ dfim, const double * __restrict__ dfimfr, 
  const double * __restrict__ dfimfr2, double dkmax, const double * __restrict__ fklam, 
  const double * __restrict__ fklam1, const double * __restrict__ fklap, 
  const double * __restrict__ fklap1, const double * __restrict__ fr, double fratio, 
  double g, double gm1, const int * __restrict__ ikm, const int * __restrict__ ikm1, 
  const int * __restrict__ ikp, const int * __restrict__ ikp1, 
  const int * __restrict__ inlcoef, int isnonlin, const int * __restrict__ k11w, 
  const int * __restrict__ k1w, const int * __restrict__ k21w, 
  const int * __restrict__ k2w, int kfrh, int mfrstlw, int mlsthg, int nang, int nfre, 
  const double * __restrict__ rnlcoef, const double * __restrict__ sinth, 
  const double * __restrict__ th, double wetail, double wp1tail, double wp2tail, 
  double xkdmin, const double * __restrict__ zpifr, int ichnk, int nchnk, int ij, 
  double * __restrict__ enh, double * __restrict__ xnu, double * __restrict__ sig_th, 
  double * __restrict__ ftemp, double * __restrict__ ad, double * __restrict__ delad, 
  double * __restrict__ delap, double * __restrict__ delam, double * __restrict__ enhfr, 
  int * __restrict__ peak_ang_mmax, double * __restrict__ peak_ang_sum0, 
  double * __restrict__ peak_ang_sum1, double * __restrict__ peak_ang_sum2, 
  double * __restrict__ peak_ang_xmax, double * __restrict__ peak_ang_temp, 
  double * __restrict__ peak_ang_thmean, double * __restrict__ peak_ang_sum_s, 
  double * __restrict__ peak_ang_sum_c) {
  
  
  // ----------------------------------------------------------------------
  



  
  
  int k;
  int m;
  int mc;
  int kh;
  int k1;
  int k2;
  int k11;
  int k21;
  int mp;
  int mp1;
  int mm;
  int mm1;
  int ic;
  int ip;
  int ip1;
  int im;
  int im1;
  int mfr1stfr;
  int mfrlstfr;
  
  double enh_max = (double) 10.0;
  double enh_min = (double) 0.1;  // to prevent ENH to become too small
  double xk;
  
  double zhook_handle;
  double ftail;
  double fklamp;
  double gw1;
  double gw2;
  double gw3;
  double gw4;
  double fklampa;
  double fklampb;
  double fklamp2;
  double fklamp1;
  double fklapa2;
  double fklapb2;
  double fklap12;
  double fklap22;
  double fklamm;
  double fklamm1;
  double gw5;
  double gw6;
  double gw7;
  double gw8;
  double fklamma;
  double fklammb;
  double fklamm2;
  double fklama2;
  double fklamb2;
  double fklam12;
  double fklam22;
  double sap;
  double sam;
  double fij;
  double fad1;
  double fad2;
  double fcen;
  
  
  
  // ----------------------------------------------------------------------
  
  
  
  //*    1. SHALLOW WATER SCALING
  //        ---------------------
  
  switch (isnonlin) {
  case 0:
    
    enhfr[ij - 1 + kijl*(ichnk - 1)] = fmax((double) ((double) 0.75*depth[ij - 1 + 
      kijl*(ichnk - 1)]*akmean[ij - 1]), (double) ((double) 0.5));
    enhfr[ij - 1 + kijl*(ichnk - 1)] = (double) 1.0 + ((double) 5.5 / enhfr[ij - 1 + 
      kijl*(ichnk - 1)])*((double) 1.0 - (double) .833*enhfr[ij - 1 + kijl*(ichnk - 1)])
      *exp((double) (-(double) 1.25*enhfr[ij - 1 + kijl*(ichnk - 1)]));
    for (mc = 1; mc <= mlsthg; mc += 1) {
      enh[ij - 1 + kijl*(mc - 1 + mlsthg*(ichnk - 1))] = enhfr[ij - 1 + kijl*(ichnk - 1)]
        ;
    }
    
    
  break;
  case 1:
    
    for (mc = 1; mc <= nfre; mc += 1) {
      enh[ij - 1 + kijl*(mc - 1 + mlsthg*(ichnk - 1))] = fmax((double) (fmin((double) 
        (enh_max), (double) (transf_c(wavnum[ij - 1 + kijl*(mc - 1 + nfre*(ichnk - 1))], 
        depth[ij - 1 + kijl*(ichnk - 1)], dkmax, g)))), (double) (enh_min));
    }
    for (mc = nfre + 1; mc <= mlsthg; mc += 1) {
      xk = gm1*(pow((zpifr[nfre - 1]*(pow(fratio, (mc - nfre)))), 2));
      enh[ij - 1 + kijl*(mc - 1 + mlsthg*(ichnk - 1))] = fmax((double) (fmin((double) 
        (enh_max), (double) (transf_c(xk, depth[ij - 1 + kijl*(ichnk - 1)], dkmax, g)))),
         (double) (enh_min));
    }
    
    
  break;
  case 2:
    peak_ang_c(kijs, kijl, fl1,  (&xnu[ + kijl*(ichnk - 1)]), 
       (&sig_th[ + kijl*(ichnk - 1)]), costh, delth, dfim, dfimfr, dfimfr2, fr, fratio, 
      nang, nfre, sinth, th, wetail, wp1tail, wp2tail, ichnk, nchnk, ij, peak_ang_mmax, 
      peak_ang_sum0, peak_ang_sum1, peak_ang_sum2, peak_ang_xmax, peak_ang_temp, 
      peak_ang_thmean, peak_ang_sum_s, peak_ang_sum_c);
    
    for (mc = 1; mc <= nfre; mc += 1) {
      enh[ij - 1 + kijl*(mc - 1 + mlsthg*(ichnk - 1))] = transf_snl_c(wavnum[ij - 1 + 
        kijl*(mc - 1 + nfre*(ichnk - 1))], depth[ij - 1 + kijl*(ichnk - 1)], xnu[ij - 1 +
         kijl*(ichnk - 1)], sig_th[ij - 1 + kijl*(ichnk - 1)], bathymax, dkmax, g, xkdmin
        );
    }
    for (mc = nfre + 1; mc <= mlsthg; mc += 1) {
      xk = gm1*(pow((zpifr[nfre - 1]*(pow(fratio, (mc - nfre)))), 2));
      enh[ij - 1 + kijl*(mc - 1 + mlsthg*(ichnk - 1))] = transf_snl_c(xk, depth[ij - 1 + 
        kijl*(ichnk - 1)], xnu[ij - 1 + kijl*(ichnk - 1)], sig_th[ij - 1 + kijl*(ichnk - 
        1)], bathymax, dkmax, g, xkdmin);
    }
    
  break;
  }
  
  
  //*    2. FREQUENCY LOOP.
  //        ---------------
  
  mfr1stfr = -mfrstlw + 1;
  mfrlstfr = nfre - kfrh + mfr1stfr;
  
  
  for (mc = 1; mc <= mlsthg; mc += 1) {
    mp = ikp[1 + mc - mfrstlw - 1];
    mp1 = ikp1[1 + mc - mfrstlw - 1];
    mm = ikm[1 + mc - mfrstlw - 1];
    mm1 = ikm1[1 + mc - mfrstlw - 1];
    ic = inlcoef[1 - 1 + 5*(mc - 1)];
    ip = inlcoef[2 - 1 + 5*(mc - 1)];
    ip1 = inlcoef[3 - 1 + 5*(mc - 1)];
    im = inlcoef[4 - 1 + 5*(mc - 1)];
    im1 = inlcoef[5 - 1 + 5*(mc - 1)];
    
    ftail = rnlcoef[1 - 1 + 25*(mc - 1)];
    
    fklamp = fklap[1 + mc - mfrstlw - 1];
    fklamp1 = fklap1[1 + mc - mfrstlw - 1];
    gw1 = rnlcoef[2 - 1 + 25*(mc - 1)];
    gw2 = rnlcoef[3 - 1 + 25*(mc - 1)];
    gw3 = rnlcoef[4 - 1 + 25*(mc - 1)];
    gw4 = rnlcoef[5 - 1 + 25*(mc - 1)];
    fklampa = rnlcoef[6 - 1 + 25*(mc - 1)];
    fklampb = rnlcoef[7 - 1 + 25*(mc - 1)];
    fklamp2 = rnlcoef[8 - 1 + 25*(mc - 1)];
    fklamp1 = rnlcoef[9 - 1 + 25*(mc - 1)];
    fklapa2 = rnlcoef[10 - 1 + 25*(mc - 1)];
    fklapb2 = rnlcoef[11 - 1 + 25*(mc - 1)];
    fklap12 = rnlcoef[12 - 1 + 25*(mc - 1)];
    fklap22 = rnlcoef[13 - 1 + 25*(mc - 1)];
    
    fklamm = fklam[1 + mc - mfrstlw - 1];
    fklamm1 = fklam1[1 + mc - mfrstlw - 1];
    gw5 = rnlcoef[14 - 1 + 25*(mc - 1)];
    gw6 = rnlcoef[15 - 1 + 25*(mc - 1)];
    gw7 = rnlcoef[16 - 1 + 25*(mc - 1)];
    gw8 = rnlcoef[17 - 1 + 25*(mc - 1)];
    fklamma = rnlcoef[18 - 1 + 25*(mc - 1)];
    fklammb = rnlcoef[19 - 1 + 25*(mc - 1)];
    fklamm2 = rnlcoef[20 - 1 + 25*(mc - 1)];
    fklamm1 = rnlcoef[21 - 1 + 25*(mc - 1)];
    fklama2 = rnlcoef[22 - 1 + 25*(mc - 1)];
    fklamb2 = rnlcoef[23 - 1 + 25*(mc - 1)];
    fklam12 = rnlcoef[24 - 1 + 25*(mc - 1)];
    fklam22 = rnlcoef[25 - 1 + 25*(mc - 1)];
    
    ftemp[ij - 1 + kijl*(ichnk - 1)] = 
      af11[1 + mc - mfrstlw - 1]*enh[ij - 1 + kijl*(mc - 1 + mlsthg*(ichnk - 1))];
    
    
    if (mc > mfr1stfr && mc < mfrlstfr) {
      //       the interactions for MC are all within the fully resolved spectral domain
      
      for (kh = 1; kh <= 2; kh += 1) {
        for (k = 1; k <= nang; k += 1) {
          k1 = k1w[k - 1 + nang*(kh - 1)];
          k2 = k2w[k - 1 + nang*(kh - 1)];
          k11 = k11w[k - 1 + nang*(kh - 1)];
          k21 = k21w[k - 1 + nang*(kh - 1)];
          
          //*    2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND
          //*            DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
          //             ----------------------------------------------
          sap = gw1*fl1[ij - 1 + kijl*(k1 - 1 + nang*(ip - 1 + nfre*(ichnk - 1)))] + 
            gw2*fl1[ij - 1 + kijl*(k11 - 1 + nang*(ip - 1 + nfre*(ichnk - 1)))] + 
            gw3*fl1[ij - 1 + kijl*(k1 - 1 + nang*(ip1 - 1 + nfre*(ichnk - 1)))] + 
            gw4*fl1[ij - 1 + kijl*(k11 - 1 + nang*(ip1 - 1 + nfre*(ichnk - 1)))];
          sam = gw5*fl1[ij - 1 + kijl*(k2 - 1 + nang*(im - 1 + nfre*(ichnk - 1)))] + 
            gw6*fl1[ij - 1 + kijl*(k21 - 1 + nang*(im - 1 + nfre*(ichnk - 1)))] + 
            gw7*fl1[ij - 1 + kijl*(k2 - 1 + nang*(im1 - 1 + nfre*(ichnk - 1)))] + 
            gw8*fl1[ij - 1 + kijl*(k21 - 1 + nang*(im1 - 1 + nfre*(ichnk - 1)))];
          //!!! not needed ftail always=1.                FIJ = FL1(IJ,K  ,IC )*FTAIL
          fij = fl1[ij - 1 + kijl*(k - 1 + nang*(ic - 1 + nfre*(ichnk - 1)))];
          fad1 = fij*(sap + sam);
          fad2 = fad1 - (double) 2.0*sap*sam;
          fad1 = fad1 + fad2;
          fcen = ftemp[ij - 1 + kijl*(ichnk - 1)]*fij;
          ad[ij - 1 + kijl*(ichnk - 1)] = fad2*fcen;
          delad[ij - 1 + kijl*(ichnk - 1)] = fad1*ftemp[ij - 1 + kijl*(ichnk - 1)];
          delap[ij - 1 + kijl*(ichnk - 1)] = (fij - (double) 2.0*sam)*dal1*fcen;
          delam[ij - 1 + kijl*(ichnk - 1)] = (fij - (double) 2.0*sap)*dal2*fcen;
          
          sl[ij - 1 + kijl*(k - 1 + nang*(mc - 1))] = sl[ij - 1 + kijl*(k - 1 + nang*(mc 
            - 1))] - (double) 2.0*ad[ij - 1 + kijl*(ichnk - 1)];
          fld[ij - 1 + kijl*(k - 1 + nang*(mc - 1))] = fld[ij - 1 + kijl*(k - 1 + 
            nang*(mc - 1))] - (double) 2.0*delad[ij - 1 + kijl*(ichnk - 1)];
          sl[ij - 1 + kijl*(k2 - 1 + nang*(mm - 1))] = sl[ij - 1 + kijl*(k2 - 1 + 
            nang*(mm - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklamm1;
          fld[ij - 1 + kijl*(k2 - 1 + nang*(mm - 1))] = fld[ij - 1 + kijl*(k2 - 1 + 
            nang*(mm - 1))] + delam[ij - 1 + kijl*(ichnk - 1)]*fklam12;
          sl[ij - 1 + kijl*(k21 - 1 + nang*(mm - 1))] = sl[ij - 1 + kijl*(k21 - 1 + 
            nang*(mm - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklamm2;
          fld[ij - 1 + kijl*(k21 - 1 + nang*(mm - 1))] = fld[ij - 1 + kijl*(k21 - 1 + 
            nang*(mm - 1))] + delam[ij - 1 + kijl*(ichnk - 1)]*fklam22;
          sl[ij - 1 + kijl*(k2 - 1 + nang*(mm1 - 1))] = sl[ij - 1 + kijl*(k2 - 1 + 
            nang*(mm1 - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklamma;
          fld[ij - 1 + kijl*(k2 - 1 + nang*(mm1 - 1))] = fld[ij - 1 + kijl*(k2 - 1 + 
            nang*(mm1 - 1))] + delam[ij - 1 + kijl*(ichnk - 1)]*fklama2;
          sl[ij - 1 + kijl*(k21 - 1 + nang*(mm1 - 1))] = sl[ij - 1 + kijl*(k21 - 1 + 
            nang*(mm1 - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklammb;
          fld[ij - 1 + kijl*(k21 - 1 + nang*(mm1 - 1))] = fld[ij - 1 + kijl*(k21 - 1 + 
            nang*(mm1 - 1))] + delam[ij - 1 + kijl*(ichnk - 1)]*fklamb2;
          sl[ij - 1 + kijl*(k1 - 1 + nang*(mp - 1))] = sl[ij - 1 + kijl*(k1 - 1 + 
            nang*(mp - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklamp1;
          fld[ij - 1 + kijl*(k1 - 1 + nang*(mp - 1))] = fld[ij - 1 + kijl*(k1 - 1 + 
            nang*(mp - 1))] + delap[ij - 1 + kijl*(ichnk - 1)]*fklap12;
          sl[ij - 1 + kijl*(k11 - 1 + nang*(mp - 1))] = sl[ij - 1 + kijl*(k11 - 1 + 
            nang*(mp - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklamp2;
          fld[ij - 1 + kijl*(k11 - 1 + nang*(mp - 1))] = fld[ij - 1 + kijl*(k11 - 1 + 
            nang*(mp - 1))] + delap[ij - 1 + kijl*(ichnk - 1)]*fklap22;
          sl[ij - 1 + kijl*(k1 - 1 + nang*(mp1 - 1))] = sl[ij - 1 + kijl*(k1 - 1 + 
            nang*(mp1 - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklampa;
          fld[ij - 1 + kijl*(k1 - 1 + nang*(mp1 - 1))] = fld[ij - 1 + kijl*(k1 - 1 + 
            nang*(mp1 - 1))] + delap[ij - 1 + kijl*(ichnk - 1)]*fklapa2;
          sl[ij - 1 + kijl*(k11 - 1 + nang*(mp1 - 1))] = sl[ij - 1 + kijl*(k11 - 1 + 
            nang*(mp1 - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklampb;
          fld[ij - 1 + kijl*(k11 - 1 + nang*(mp1 - 1))] = fld[ij - 1 + kijl*(k11 - 1 + 
            nang*(mp1 - 1))] + delap[ij - 1 + kijl*(ichnk - 1)]*fklapb2;
        }
      }
      
    } else if (mc >= mfrlstfr) {
      for (kh = 1; kh <= 2; kh += 1) {
        for (k = 1; k <= nang; k += 1) {
          k1 = k1w[k - 1 + nang*(kh - 1)];
          k2 = k2w[k - 1 + nang*(kh - 1)];
          k11 = k11w[k - 1 + nang*(kh - 1)];
          k21 = k21w[k - 1 + nang*(kh - 1)];
          
          sap = gw1*fl1[ij - 1 + kijl*(k1 - 1 + nang*(ip - 1 + nfre*(ichnk - 1)))] + 
            gw2*fl1[ij - 1 + kijl*(k11 - 1 + nang*(ip - 1 + nfre*(ichnk - 1)))] + 
            gw3*fl1[ij - 1 + kijl*(k1 - 1 + nang*(ip1 - 1 + nfre*(ichnk - 1)))] + 
            gw4*fl1[ij - 1 + kijl*(k11 - 1 + nang*(ip1 - 1 + nfre*(ichnk - 1)))];
          sam = gw5*fl1[ij - 1 + kijl*(k2 - 1 + nang*(im - 1 + nfre*(ichnk - 1)))] + 
            gw6*fl1[ij - 1 + kijl*(k21 - 1 + nang*(im - 1 + nfre*(ichnk - 1)))] + 
            gw7*fl1[ij - 1 + kijl*(k2 - 1 + nang*(im1 - 1 + nfre*(ichnk - 1)))] + 
            gw8*fl1[ij - 1 + kijl*(k21 - 1 + nang*(im1 - 1 + nfre*(ichnk - 1)))];
          fij = fl1[ij - 1 + kijl*(k - 1 + nang*(ic - 1 + nfre*(ichnk - 1)))]*ftail;
          fad1 = fij*(sap + sam);
          fad2 = fad1 - (double) 2.0*sap*sam;
          fad1 = fad1 + fad2;
          fcen = ftemp[ij - 1 + kijl*(ichnk - 1)]*fij;
          ad[ij - 1 + kijl*(ichnk - 1)] = fad2*fcen;
          delad[ij - 1 + kijl*(ichnk - 1)] = fad1*ftemp[ij - 1 + kijl*(ichnk - 1)];
          delap[ij - 1 + kijl*(ichnk - 1)] = (fij - (double) 2.0*sam)*dal1*fcen;
          delam[ij - 1 + kijl*(ichnk - 1)] = (fij - (double) 2.0*sap)*dal2*fcen;
          
          sl[ij - 1 + kijl*(k2 - 1 + nang*(mm - 1))] = sl[ij - 1 + kijl*(k2 - 1 + 
            nang*(mm - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklamm1;
          fld[ij - 1 + kijl*(k2 - 1 + nang*(mm - 1))] = fld[ij - 1 + kijl*(k2 - 1 + 
            nang*(mm - 1))] + delam[ij - 1 + kijl*(ichnk - 1)]*fklam12;
          sl[ij - 1 + kijl*(k21 - 1 + nang*(mm - 1))] = sl[ij - 1 + kijl*(k21 - 1 + 
            nang*(mm - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklamm2;
          fld[ij - 1 + kijl*(k21 - 1 + nang*(mm - 1))] = fld[ij - 1 + kijl*(k21 - 1 + 
            nang*(mm - 1))] + delam[ij - 1 + kijl*(ichnk - 1)]*fklam22;
          
          if (mm1 <= nfre) {
            sl[ij - 1 + kijl*(k2 - 1 + nang*(mm1 - 1))] = sl[ij - 1 + kijl*(k2 - 1 + 
              nang*(mm1 - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklamma;
            fld[ij - 1 + kijl*(k2 - 1 + nang*(mm1 - 1))] = fld[ij - 1 + kijl*(k2 - 1 + 
              nang*(mm1 - 1))] + delam[ij - 1 + kijl*(ichnk - 1)]*fklama2;
            sl[ij - 1 + kijl*(k21 - 1 + nang*(mm1 - 1))] = sl[ij - 1 + kijl*(k21 - 1 + 
              nang*(mm1 - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklammb;
            fld[ij - 1 + kijl*(k21 - 1 + nang*(mm1 - 1))] = fld[ij - 1 + kijl*(k21 - 1 + 
              nang*(mm1 - 1))] + delam[ij - 1 + kijl*(ichnk - 1)]*fklamb2;
            
            if (mc <= nfre) {
              sl[ij - 1 + kijl*(k - 1 + nang*(mc - 1))] = sl[ij - 1 + kijl*(k - 1 + 
                nang*(mc - 1))] - (double) 2.0*ad[ij - 1 + kijl*(ichnk - 1)];
              fld[ij - 1 + kijl*(k - 1 + nang*(mc - 1))] = fld[ij - 1 + kijl*(k - 1 + 
                nang*(mc - 1))] - (double) 2.0*delad[ij - 1 + kijl*(ichnk - 1)];
              
              if (mp <= nfre) {
                sl[ij - 1 + kijl*(k1 - 1 + nang*(mp - 1))] = sl[ij - 1 + kijl*(k1 - 1 + 
                  nang*(mp - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklamp1;
                fld[ij - 1 + kijl*(k1 - 1 + nang*(mp - 1))] = fld[ij - 1 + kijl*(k1 - 1 +
                   nang*(mp - 1))] + delap[ij - 1 + kijl*(ichnk - 1)]*fklap12;
                sl[ij - 1 + kijl*(k11 - 1 + nang*(mp - 1))] = sl[ij - 1 + kijl*(k11 - 1 +
                   nang*(mp - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklamp2;
                fld[ij - 1 + kijl*(k11 - 1 + nang*(mp - 1))] = fld[ij - 1 + kijl*(k11 - 1
                   + nang*(mp - 1))] + delap[ij - 1 + kijl*(ichnk - 1)]*fklap22;
                
                if (mp1 <= nfre) {
                  sl[ij - 1 + kijl*(k1 - 1 + nang*(mp1 - 1))] = sl[ij - 1 + kijl*(k1 - 1 
                    + nang*(mp1 - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklampa;
                  fld[ij - 1 + kijl*(k1 - 1 + nang*(mp1 - 1))] = fld[ij - 1 + kijl*(k1 - 
                    1 + nang*(mp1 - 1))] + delap[ij - 1 + kijl*(ichnk - 1)]*fklapa2;
                  sl[ij - 1 + kijl*(k11 - 1 + nang*(mp1 - 1))] = sl[ij - 1 + kijl*(k11 - 
                    1 + nang*(mp1 - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklampb;
                  fld[ij - 1 + kijl*(k11 - 1 + nang*(mp1 - 1))] = fld[ij - 1 + kijl*(k11 
                    - 1 + nang*(mp1 - 1))] + delap[ij - 1 + kijl*(ichnk - 1)]*fklapb2;
                }
              }
            }
          }
        }
      }
      
    } else {
      
      for (kh = 1; kh <= 2; kh += 1) {
        for (k = 1; k <= nang; k += 1) {
          k1 = k1w[k - 1 + nang*(kh - 1)];
          k2 = k2w[k - 1 + nang*(kh - 1)];
          k11 = k11w[k - 1 + nang*(kh - 1)];
          k21 = k21w[k - 1 + nang*(kh - 1)];
          
          sap = gw1*fl1[ij - 1 + kijl*(k1 - 1 + nang*(ip - 1 + nfre*(ichnk - 1)))] + 
            gw2*fl1[ij - 1 + kijl*(k11 - 1 + nang*(ip - 1 + nfre*(ichnk - 1)))] + 
            gw3*fl1[ij - 1 + kijl*(k1 - 1 + nang*(ip1 - 1 + nfre*(ichnk - 1)))] + 
            gw4*fl1[ij - 1 + kijl*(k11 - 1 + nang*(ip1 - 1 + nfre*(ichnk - 1)))];
          sam = gw5*fl1[ij - 1 + kijl*(k2 - 1 + nang*(im - 1 + nfre*(ichnk - 1)))] + 
            gw6*fl1[ij - 1 + kijl*(k21 - 1 + nang*(im - 1 + nfre*(ichnk - 1)))] + 
            gw7*fl1[ij - 1 + kijl*(k2 - 1 + nang*(im1 - 1 + nfre*(ichnk - 1)))] + 
            gw8*fl1[ij - 1 + kijl*(k21 - 1 + nang*(im1 - 1 + nfre*(ichnk - 1)))];
          fij = fl1[ij - 1 + kijl*(k - 1 + nang*(ic - 1 + nfre*(ichnk - 1)))]*ftail;
          fad1 = fij*(sap + sam);
          fad2 = fad1 - (double) 2.0*sap*sam;
          fad1 = fad1 + fad2;
          fcen = ftemp[ij - 1 + kijl*(ichnk - 1)]*fij;
          ad[ij - 1 + kijl*(ichnk - 1)] = fad2*fcen;
          delad[ij - 1 + kijl*(ichnk - 1)] = fad1*ftemp[ij - 1 + kijl*(ichnk - 1)];
          delap[ij - 1 + kijl*(ichnk - 1)] = (fij - (double) 2.0*sam)*dal1*fcen;
          delam[ij - 1 + kijl*(ichnk - 1)] = (fij - (double) 2.0*sap)*dal2*fcen;
          
          if (mm1 >= 1) {
            sl[ij - 1 + kijl*(k2 - 1 + nang*(mm1 - 1))] = sl[ij - 1 + kijl*(k2 - 1 + 
              nang*(mm1 - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklamma;
            fld[ij - 1 + kijl*(k2 - 1 + nang*(mm1 - 1))] = fld[ij - 1 + kijl*(k2 - 1 + 
              nang*(mm1 - 1))] + delam[ij - 1 + kijl*(ichnk - 1)]*fklama2;
            sl[ij - 1 + kijl*(k21 - 1 + nang*(mm1 - 1))] = sl[ij - 1 + kijl*(k21 - 1 + 
              nang*(mm1 - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklammb;
            fld[ij - 1 + kijl*(k21 - 1 + nang*(mm1 - 1))] = fld[ij - 1 + kijl*(k21 - 1 + 
              nang*(mm1 - 1))] + delam[ij - 1 + kijl*(ichnk - 1)]*fklamb2;
          }
          
          sl[ij - 1 + kijl*(k - 1 + nang*(mc - 1))] = sl[ij - 1 + kijl*(k - 1 + nang*(mc 
            - 1))] - (double) 2.0*ad[ij - 1 + kijl*(ichnk - 1)];
          fld[ij - 1 + kijl*(k - 1 + nang*(mc - 1))] = fld[ij - 1 + kijl*(k - 1 + 
            nang*(mc - 1))] - (double) 2.0*delad[ij - 1 + kijl*(ichnk - 1)];
          sl[ij - 1 + kijl*(k1 - 1 + nang*(mp - 1))] = sl[ij - 1 + kijl*(k1 - 1 + 
            nang*(mp - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklamp1;
          fld[ij - 1 + kijl*(k1 - 1 + nang*(mp - 1))] = fld[ij - 1 + kijl*(k1 - 1 + 
            nang*(mp - 1))] + delap[ij - 1 + kijl*(ichnk - 1)]*fklap12;
          sl[ij - 1 + kijl*(k11 - 1 + nang*(mp - 1))] = sl[ij - 1 + kijl*(k11 - 1 + 
            nang*(mp - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklamp2;
          fld[ij - 1 + kijl*(k11 - 1 + nang*(mp - 1))] = fld[ij - 1 + kijl*(k11 - 1 + 
            nang*(mp - 1))] + delap[ij - 1 + kijl*(ichnk - 1)]*fklap22;
          sl[ij - 1 + kijl*(k1 - 1 + nang*(mp1 - 1))] = sl[ij - 1 + kijl*(k1 - 1 + 
            nang*(mp1 - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklampa;
          fld[ij - 1 + kijl*(k1 - 1 + nang*(mp1 - 1))] = fld[ij - 1 + kijl*(k1 - 1 + 
            nang*(mp1 - 1))] + delap[ij - 1 + kijl*(ichnk - 1)]*fklapa2;
          sl[ij - 1 + kijl*(k11 - 1 + nang*(mp1 - 1))] = sl[ij - 1 + kijl*(k11 - 1 + 
            nang*(mp1 - 1))] + ad[ij - 1 + kijl*(ichnk - 1)]*fklampb;
          fld[ij - 1 + kijl*(k11 - 1 + nang*(mp1 - 1))] = fld[ij - 1 + kijl*(k11 - 1 + 
            nang*(mp1 - 1))] + delap[ij - 1 + kijl*(ichnk - 1)]*fklapb2;
        }
      }
      
    }
    
    //*    BRANCH BACK TO 2. FOR NEXT FREQUENCY.
    
  }
  
  
  
}
