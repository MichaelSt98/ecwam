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

__device__ void snonlin_c(int kijs, int kijl, const double * fl1, double * fld, 
  double * sl, const double * wavnum, const double * depth, const double * akmean, 
  const double * af11, double bathymax, const double * costh, double dal1, double dal2, 
  double delth, const double * dfim, const double * dfimfr, const double * dfimfr2, 
  double dkmax, const double * fklam, const double * fklam1, const double * fklap, 
  const double * fklap1, const double * fr, double fratio, double g, double gm1, 
  const int * ikm, const int * ikm1, const int * ikp, const int * ikp1, 
  const int * inlcoef, int isnonlin, const int * k11w, const int * k1w, 
  const int * k21w, const int * k2w, int kfrh, int mfrstlw, int mlsthg, int nang, 
  int nfre, const double * rnlcoef, const double * sinth, const double * th, 
  double wetail, double wp1tail, double wp2tail, double xkdmin, const double * zpifr, 
  int ichnk, int nchnk, int ij, double * enh, double * xnu, double * sig_th) {
  
  



  const int nang_loki_param = 24;
  const int nfre_loki_param = 36;
  const int nrnl = 25;
  const int ninl = 5;
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
  
  double ftemp;
  double ad;
  double delad;
  double delap;
  double delam;
  double enhfr;
  switch (isnonlin) {
  case 0:

    enhfr = max((double) ((double) 0.75*depth[ij - 1 + kijl*(ichnk - 1)]*akmean[ij - 1]),
       (double) ((double) 0.5));
    enhfr = (double) 1.0 + ((double) 5.5 / enhfr)*((double) 1.0 - (double) .833*enhfr)
      *exp((double) (-(double) 1.25*enhfr));
    for (mc = 1; mc <= mlsthg; mc += 1) {
      enh[ij - 1 + kijl*(mc - 1 + mlsthg*(ichnk - 1))] = enhfr;
    }

    
  break;
  case 1:

    for (mc = 1; mc <= nfre; mc += 1) {
      enh[ij - 1 + kijl*(mc - 1 + mlsthg*(ichnk - 1))] = max((double) (min((double) 
        (enh_max), (double) (transf_c(wavnum[ij - 1 + kijl*(mc - 1 + 
        nfre_loki_param*(ichnk - 1))], depth[ij - 1 + kijl*(ichnk - 1)], dkmax, g)))), 
        (double) (enh_min));
    }
    for (mc = nfre + 1; mc <= mlsthg; mc += 1) {
      xk = gm1*(pow((zpifr[nfre - 1]*(pow(fratio, (mc - nfre)))), 2));
      enh[ij - 1 + kijl*(mc - 1 + mlsthg*(ichnk - 1))] = max((double) (min((double) 
        (enh_max), (double) (transf_c(xk, depth[ij - 1 + kijl*(ichnk - 1)], dkmax, g)))),
         (double) (enh_min));
    }

    
  break;
  case 2:
    peak_ang_c(kijs, kijl, fl1,  (&xnu[ + kijl*(ichnk - 1)]), 
       (&sig_th[ + kijl*(ichnk - 1)]), costh, delth, dfim, dfimfr, dfimfr2, fr, fratio, 
      nang, nfre, sinth, th, wetail, wp1tail, wp2tail, ichnk, nchnk, ij);

    for (mc = 1; mc <= nfre; mc += 1) {
      enh[ij - 1 + kijl*(mc - 1 + mlsthg*(ichnk - 1))] = transf_snl_c(wavnum[ij - 1 + 
        kijl*(mc - 1 + nfre_loki_param*(ichnk - 1))], depth[ij - 1 + kijl*(ichnk - 1)], 
        xnu[ij - 1 + kijl*(ichnk - 1)], sig_th[ij - 1 + kijl*(ichnk - 1)], bathymax, 
        dkmax, g, xkdmin);
    }
    for (mc = nfre + 1; mc <= mlsthg; mc += 1) {
      xk = gm1*(pow((zpifr[nfre - 1]*(pow(fratio, (mc - nfre)))), 2));
      enh[ij - 1 + kijl*(mc - 1 + mlsthg*(ichnk - 1))] = transf_snl_c(xk, depth[ij - 1 + 
        kijl*(ichnk - 1)], xnu[ij - 1 + kijl*(ichnk - 1)], sig_th[ij - 1 + kijl*(ichnk - 
        1)], bathymax, dkmax, g, xkdmin);
    }

  break;
  }
  mfr1stfr = -mfrstlw + 1;
  mfrlstfr = nfre - kfrh + mfr1stfr;
  

  for (mc = 1; mc <= mlsthg; mc += 1) {
    mp = ikp[1 + mc - mfrstlw - 1];
    mp1 = ikp1[1 + mc - mfrstlw - 1];
    mm = ikm[1 + mc - mfrstlw - 1];
    mm1 = ikm1[1 + mc - mfrstlw - 1];
    ic = inlcoef[1 - 1 + ninl*(mc - 1)];
    ip = inlcoef[2 - 1 + ninl*(mc - 1)];
    ip1 = inlcoef[3 - 1 + ninl*(mc - 1)];
    im = inlcoef[4 - 1 + ninl*(mc - 1)];
    im1 = inlcoef[5 - 1 + ninl*(mc - 1)];
    
    ftail = rnlcoef[1 - 1 + nrnl*(mc - 1)];
    
    fklamp = fklap[1 + mc - mfrstlw - 1];
    fklamp1 = fklap1[1 + mc - mfrstlw - 1];
    gw1 = rnlcoef[2 - 1 + nrnl*(mc - 1)];
    gw2 = rnlcoef[3 - 1 + nrnl*(mc - 1)];
    gw3 = rnlcoef[4 - 1 + nrnl*(mc - 1)];
    gw4 = rnlcoef[5 - 1 + nrnl*(mc - 1)];
    fklampa = rnlcoef[6 - 1 + nrnl*(mc - 1)];
    fklampb = rnlcoef[7 - 1 + nrnl*(mc - 1)];
    fklamp2 = rnlcoef[8 - 1 + nrnl*(mc - 1)];
    fklamp1 = rnlcoef[9 - 1 + nrnl*(mc - 1)];
    fklapa2 = rnlcoef[10 - 1 + nrnl*(mc - 1)];
    fklapb2 = rnlcoef[11 - 1 + nrnl*(mc - 1)];
    fklap12 = rnlcoef[12 - 1 + nrnl*(mc - 1)];
    fklap22 = rnlcoef[13 - 1 + nrnl*(mc - 1)];
    
    fklamm = fklam[1 + mc - mfrstlw - 1];
    fklamm1 = fklam1[1 + mc - mfrstlw - 1];
    gw5 = rnlcoef[14 - 1 + nrnl*(mc - 1)];
    gw6 = rnlcoef[15 - 1 + nrnl*(mc - 1)];
    gw7 = rnlcoef[16 - 1 + nrnl*(mc - 1)];
    gw8 = rnlcoef[17 - 1 + nrnl*(mc - 1)];
    fklamma = rnlcoef[18 - 1 + nrnl*(mc - 1)];
    fklammb = rnlcoef[19 - 1 + nrnl*(mc - 1)];
    fklamm2 = rnlcoef[20 - 1 + nrnl*(mc - 1)];
    fklamm1 = rnlcoef[21 - 1 + nrnl*(mc - 1)];
    fklama2 = rnlcoef[22 - 1 + nrnl*(mc - 1)];
    fklamb2 = rnlcoef[23 - 1 + nrnl*(mc - 1)];
    fklam12 = rnlcoef[24 - 1 + nrnl*(mc - 1)];
    fklam22 = rnlcoef[25 - 1 + nrnl*(mc - 1)];
    
    ftemp = af11[1 + mc - mfrstlw - 1]*enh[ij - 1 + kijl*(mc - 1 + mlsthg*(ichnk - 1))];
    if (mc > mfr1stfr && mc < mfrlstfr) {
      for (kh = 1; kh <= 2; kh += 1) {
        for (k = 1; k <= nang; k += 1) {
          k1 = k1w[k - 1 + nang_loki_param*(kh - 1)];
          k2 = k2w[k - 1 + nang_loki_param*(kh - 1)];
          k11 = k11w[k - 1 + nang_loki_param*(kh - 1)];
          k21 = k21w[k - 1 + nang_loki_param*(kh - 1)];
          sap = gw1*fl1[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(ip - 1 + 
            nfre_loki_param*(ichnk - 1)))] + gw2*fl1[ij - 1 + kijl*(k11 - 1 + 
            nang_loki_param*(ip - 1 + nfre_loki_param*(ichnk - 1)))] + gw3*fl1[ij - 1 + 
            kijl*(k1 - 1 + nang_loki_param*(ip1 - 1 + nfre_loki_param*(ichnk - 1)))] + 
            gw4*fl1[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(ip1 - 1 + 
            nfre_loki_param*(ichnk - 1)))];
          sam = gw5*fl1[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(im - 1 + 
            nfre_loki_param*(ichnk - 1)))] + gw6*fl1[ij - 1 + kijl*(k21 - 1 + 
            nang_loki_param*(im - 1 + nfre_loki_param*(ichnk - 1)))] + gw7*fl1[ij - 1 + 
            kijl*(k2 - 1 + nang_loki_param*(im1 - 1 + nfre_loki_param*(ichnk - 1)))] + 
            gw8*fl1[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(im1 - 1 + 
            nfre_loki_param*(ichnk - 1)))];
          //!!! not needed ftail always=1.                FIJ = FL1(IJ,K  ,IC )*FTAIL
          fij = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(ic - 1 + 
            nfre_loki_param*(ichnk - 1)))];
          fad1 = fij*(sap + sam);
          fad2 = fad1 - (double) 2.0*sap*sam;
          fad1 = fad1 + fad2;
          fcen = ftemp*fij;
          ad = fad2*fcen;
          delad = fad1*ftemp;
          delap = (fij - (double) 2.0*sam)*dal1*fcen;
          delam = (fij - (double) 2.0*sap)*dal2*fcen;
          
          sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(mc - 1))] = 
            sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(mc - 1))] - (double) 2.0*ad;
          fld[ij - 1 + kijl*(k - 1 + nang_loki_param*(mc - 1))] = 
            fld[ij - 1 + kijl*(k - 1 + nang_loki_param*(mc - 1))] - (double) 2.0*delad;
          sl[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm - 1))] = 
            sl[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm - 1))] + ad*fklamm1;
          fld[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm - 1))] = 
            fld[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm - 1))] + delam*fklam12;
          sl[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm - 1))] = 
            sl[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm - 1))] + ad*fklamm2;
          fld[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm - 1))] = 
            fld[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm - 1))] + delam*fklam22;
          sl[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm1 - 1))] = 
            sl[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm1 - 1))] + ad*fklamma;
          fld[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm1 - 1))] = 
            fld[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm1 - 1))] + delam*fklama2;
          sl[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm1 - 1))] = 
            sl[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm1 - 1))] + ad*fklammb;
          fld[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm1 - 1))] = 
            fld[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm1 - 1))] + delam*fklamb2;
          sl[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp - 1))] = 
            sl[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp - 1))] + ad*fklamp1;
          fld[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp - 1))] = 
            fld[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp - 1))] + delap*fklap12;
          sl[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp - 1))] = 
            sl[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp - 1))] + ad*fklamp2;
          fld[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp - 1))] = 
            fld[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp - 1))] + delap*fklap22;
          sl[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp1 - 1))] = 
            sl[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp1 - 1))] + ad*fklampa;
          fld[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp1 - 1))] = 
            fld[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp1 - 1))] + delap*fklapa2;
          sl[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp1 - 1))] = 
            sl[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp1 - 1))] + ad*fklampb;
          fld[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp1 - 1))] = 
            fld[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp1 - 1))] + delap*fklapb2;
        }
      }
      
    } else if (mc >= mfrlstfr) {
      for (kh = 1; kh <= 2; kh += 1) {
        for (k = 1; k <= nang; k += 1) {
          k1 = k1w[k - 1 + nang_loki_param*(kh - 1)];
          k2 = k2w[k - 1 + nang_loki_param*(kh - 1)];
          k11 = k11w[k - 1 + nang_loki_param*(kh - 1)];
          k21 = k21w[k - 1 + nang_loki_param*(kh - 1)];
          
          sap = gw1*fl1[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(ip - 1 + 
            nfre_loki_param*(ichnk - 1)))] + gw2*fl1[ij - 1 + kijl*(k11 - 1 + 
            nang_loki_param*(ip - 1 + nfre_loki_param*(ichnk - 1)))] + gw3*fl1[ij - 1 + 
            kijl*(k1 - 1 + nang_loki_param*(ip1 - 1 + nfre_loki_param*(ichnk - 1)))] + 
            gw4*fl1[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(ip1 - 1 + 
            nfre_loki_param*(ichnk - 1)))];
          sam = gw5*fl1[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(im - 1 + 
            nfre_loki_param*(ichnk - 1)))] + gw6*fl1[ij - 1 + kijl*(k21 - 1 + 
            nang_loki_param*(im - 1 + nfre_loki_param*(ichnk - 1)))] + gw7*fl1[ij - 1 + 
            kijl*(k2 - 1 + nang_loki_param*(im1 - 1 + nfre_loki_param*(ichnk - 1)))] + 
            gw8*fl1[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(im1 - 1 + 
            nfre_loki_param*(ichnk - 1)))];
          fij = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(ic - 1 + 
            nfre_loki_param*(ichnk - 1)))]*ftail;
          fad1 = fij*(sap + sam);
          fad2 = fad1 - (double) 2.0*sap*sam;
          fad1 = fad1 + fad2;
          fcen = ftemp*fij;
          ad = fad2*fcen;
          delad = fad1*ftemp;
          delap = (fij - (double) 2.0*sam)*dal1*fcen;
          delam = (fij - (double) 2.0*sap)*dal2*fcen;
          
          sl[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm - 1))] = 
            sl[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm - 1))] + ad*fklamm1;
          fld[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm - 1))] = 
            fld[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm - 1))] + delam*fklam12;
          sl[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm - 1))] = 
            sl[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm - 1))] + ad*fklamm2;
          fld[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm - 1))] = 
            fld[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm - 1))] + delam*fklam22;
          
          if (mm1 <= nfre) {
            sl[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm1 - 1))] = 
              sl[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm1 - 1))] + ad*fklamma;
            fld[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm1 - 1))] = 
              fld[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm1 - 1))] + delam*fklama2;
            sl[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm1 - 1))] = 
              sl[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm1 - 1))] + ad*fklammb;
            fld[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm1 - 1))] = 
              fld[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm1 - 1))] + delam*fklamb2;
            
            if (mc <= nfre) {
              sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(mc - 1))] = 
                sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(mc - 1))] - (double) 2.0*ad;
              fld[ij - 1 + kijl*(k - 1 + nang_loki_param*(mc - 1))] = fld[ij - 1 + 
                kijl*(k - 1 + nang_loki_param*(mc - 1))] - (double) 2.0*delad;
              
              if (mp <= nfre) {
                sl[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp - 1))] = 
                  sl[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp - 1))] + ad*fklamp1;
                fld[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp - 1))] = 
                  fld[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp - 1))] + delap*fklap12;
                sl[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp - 1))] = 
                  sl[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp - 1))] + ad*fklamp2;
                fld[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp - 1))] = 
                  fld[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp - 1))] + delap*fklap22
                  ;
                
                if (mp1 <= nfre) {
                  sl[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp1 - 1))] = 
                    sl[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp1 - 1))] + ad*fklampa;
                  fld[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp1 - 1))] = fld[ij - 1 + 
                    kijl*(k1 - 1 + nang_loki_param*(mp1 - 1))] + delap*fklapa2;
                  sl[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp1 - 1))] = 
                    sl[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp1 - 1))] + ad*fklampb;
                  fld[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp1 - 1))] = fld[ij - 1 +
                     kijl*(k11 - 1 + nang_loki_param*(mp1 - 1))] + delap*fklapb2;
                }
              }
            }
          }
        }
      }
      
    } else {
      
      for (kh = 1; kh <= 2; kh += 1) {
        for (k = 1; k <= nang; k += 1) {
          k1 = k1w[k - 1 + nang_loki_param*(kh - 1)];
          k2 = k2w[k - 1 + nang_loki_param*(kh - 1)];
          k11 = k11w[k - 1 + nang_loki_param*(kh - 1)];
          k21 = k21w[k - 1 + nang_loki_param*(kh - 1)];
          
          sap = gw1*fl1[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(ip - 1 + 
            nfre_loki_param*(ichnk - 1)))] + gw2*fl1[ij - 1 + kijl*(k11 - 1 + 
            nang_loki_param*(ip - 1 + nfre_loki_param*(ichnk - 1)))] + gw3*fl1[ij - 1 + 
            kijl*(k1 - 1 + nang_loki_param*(ip1 - 1 + nfre_loki_param*(ichnk - 1)))] + 
            gw4*fl1[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(ip1 - 1 + 
            nfre_loki_param*(ichnk - 1)))];
          sam = gw5*fl1[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(im - 1 + 
            nfre_loki_param*(ichnk - 1)))] + gw6*fl1[ij - 1 + kijl*(k21 - 1 + 
            nang_loki_param*(im - 1 + nfre_loki_param*(ichnk - 1)))] + gw7*fl1[ij - 1 + 
            kijl*(k2 - 1 + nang_loki_param*(im1 - 1 + nfre_loki_param*(ichnk - 1)))] + 
            gw8*fl1[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(im1 - 1 + 
            nfre_loki_param*(ichnk - 1)))];
          fij = fl1[ij - 1 + kijl*(k - 1 + nang_loki_param*(ic - 1 + 
            nfre_loki_param*(ichnk - 1)))]*ftail;
          fad1 = fij*(sap + sam);
          fad2 = fad1 - (double) 2.0*sap*sam;
          fad1 = fad1 + fad2;
          fcen = ftemp*fij;
          ad = fad2*fcen;
          delad = fad1*ftemp;
          delap = (fij - (double) 2.0*sam)*dal1*fcen;
          delam = (fij - (double) 2.0*sap)*dal2*fcen;
          
          if (mm1 >= 1) {
            sl[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm1 - 1))] = 
              sl[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm1 - 1))] + ad*fklamma;
            fld[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm1 - 1))] = 
              fld[ij - 1 + kijl*(k2 - 1 + nang_loki_param*(mm1 - 1))] + delam*fklama2;
            sl[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm1 - 1))] = 
              sl[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm1 - 1))] + ad*fklammb;
            fld[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm1 - 1))] = 
              fld[ij - 1 + kijl*(k21 - 1 + nang_loki_param*(mm1 - 1))] + delam*fklamb2;
          }
          
          sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(mc - 1))] = 
            sl[ij - 1 + kijl*(k - 1 + nang_loki_param*(mc - 1))] - (double) 2.0*ad;
          fld[ij - 1 + kijl*(k - 1 + nang_loki_param*(mc - 1))] = 
            fld[ij - 1 + kijl*(k - 1 + nang_loki_param*(mc - 1))] - (double) 2.0*delad;
          sl[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp - 1))] = 
            sl[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp - 1))] + ad*fklamp1;
          fld[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp - 1))] = 
            fld[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp - 1))] + delap*fklap12;
          sl[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp - 1))] = 
            sl[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp - 1))] + ad*fklamp2;
          fld[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp - 1))] = 
            fld[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp - 1))] + delap*fklap22;
          sl[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp1 - 1))] = 
            sl[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp1 - 1))] + ad*fklampa;
          fld[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp1 - 1))] = 
            fld[ij - 1 + kijl*(k1 - 1 + nang_loki_param*(mp1 - 1))] + delap*fklapa2;
          sl[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp1 - 1))] = 
            sl[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp1 - 1))] + ad*fklampb;
          fld[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp1 - 1))] = 
            fld[ij - 1 + kijl*(k11 - 1 + nang_loki_param*(mp1 - 1))] + delap*fklapb2;
        }
      }
      
    }
  }

  
  
}
