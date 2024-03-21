! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SNONLIN_CUF_HOIST_NEW_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE SNONLIN_CUF_HOIST_NEW (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN, AF11, BATHYMAX, COSTH,  &
  & DAL1, DAL2, DELTH, DFIM, DFIMFR, DFIMFR2, DKMAX, FKLAM, FKLAM1, FKLAP, FKLAP1, FR, FRATIO, G, GM1, IKM, IKM1, IKP, IKP1,  &
  & INLCOEF, ISNONLIN, K11W, K1W, K21W, K2W, KFRH, MFRSTLW, MLSTHG, NANG, NFRE, RNLCOEF, SINTH, TH, WETAIL, WP1TAIL, WP2TAIL,  &
  & XKDMIN, ZPIFR, ICHNK, NCHNK, IJ, XNU, SIG_TH, ENH)
    
    ! ----------------------------------------------------------------------
    
    !**** *SNONLIN* - COMPUTATION OF NONLINEAR TRANSFER RATE AND ITS
    !****             FUNCTIONAL DERIVATIVE (DIAGONAL TERMS ONLY) AND
    !****             ADDITION TO CORRESPONDING NET EXPRESSIONS.
    
    !     S.D. HASSELMANN.  MPI
    
    !     G. KOMEN, P. JANSSEN   KNMI             MODIFIED TO SHALLOW WATER
    !     H. GUENTHER, L. ZAMBRESKY               OPTIMIZED
    !     H. GUENTHER       GKSS/ECMWF  JUNE 1991 INTERACTIONS BETWEEN DIAG-
    !                                             AND PROGNOSTIC PART.
    !     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
    !     P. JANSSEN  ECMWF  JUNE 2005       IMPROVED SCALING IN SHALLOW
    !                                        WATER
    !     J. BIDLOT   ECMWF  AUGUST 2006     KEEP THE OLD FORMULATION
    !                                        UNDER A SWITCH (ISNONLIN = 0 for OLD
    !                                                                 = 1 for NEW
    !                                        BE AWARE THAT THE OLD FORMULATION
    !                                        REQUIRES THE MEAN WAVE NUMBER AKMEAN.
    !     J. BIDLOT   ECMWF  JANUARY 2012    ADD EXTENSION TO LOW FREQUENCIES
    !                                        OPTIMISATION FOR IBM.
    
    !*    PURPOSE.
    !     --------
    
    !       SEE ABOVE.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *SNONLIN (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)*
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *FL1*    - SPECTRUM.
    !          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*     - TOTAL SOURCE FUNCTION ARRAY.
    !          *WAVNUM* - WAVE NUMBER.
    !          *DEPTH*  - WATER DEPTH.
    !          *AKMEAN* - MEAN WAVE NUMBER  BASED ON sqrt(1/k)*F INTGRATION
    
    !     METHOD.
    !     -------
    
    !       NONE.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       NONE.
    
    ! ----------------------------------------------------------------------
    USE TRANSF_SNL_CUF_HOIST_NEW_MOD, ONLY: TRANSF_SNL_CUF_HOIST_NEW
    USE TRANSF_CUF_HOIST_NEW_MOD, ONLY: TRANSF_CUF_HOIST_NEW
    USE PEAK_ANG_CUF_HOIST_NEW_MOD, ONLY: PEAK_ANG_CUF_HOIST_NEW
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
    INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
    INTEGER(KIND=JWIM), PARAMETER :: NRNL = 25
    INTEGER(KIND=JWIM), PARAMETER :: NINL = 5
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    REAL(KIND=JWRB), INTENT(IN) :: FL1(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(INOUT) :: FLD(KIJL, NANG_loki_param, NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(INOUT) :: SL(KIJL, NANG_loki_param, NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(IN) :: WAVNUM(KIJL, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: DEPTH(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: AKMEAN(KIJL)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: K
    INTEGER(KIND=JWIM) :: M
    INTEGER(KIND=JWIM) :: MC
    INTEGER(KIND=JWIM) :: KH
    INTEGER(KIND=JWIM) :: K1
    INTEGER(KIND=JWIM) :: K2
    INTEGER(KIND=JWIM) :: K11
    INTEGER(KIND=JWIM) :: K21
    INTEGER(KIND=JWIM) :: MP
    INTEGER(KIND=JWIM) :: MP1
    INTEGER(KIND=JWIM) :: MM
    INTEGER(KIND=JWIM) :: MM1
    INTEGER(KIND=JWIM) :: IC
    INTEGER(KIND=JWIM) :: IP
    INTEGER(KIND=JWIM) :: IP1
    INTEGER(KIND=JWIM) :: IM
    INTEGER(KIND=JWIM) :: IM1
    INTEGER(KIND=JWIM) :: MFR1STFR
    INTEGER(KIND=JWIM) :: MFRLSTFR
    
    REAL(KIND=JWRB), PARAMETER :: ENH_MAX = 10.0_JWRB
    REAL(KIND=JWRB), PARAMETER :: ENH_MIN = 0.1_JWRB    ! to prevent ENH to become too small
    REAL(KIND=JWRB) :: XK
    ! REAL(KIND=JWRB) :: ENH(MLSTHG)
    REAL(KIND=JWRB), INTENT(INOUT), DEVICE, DIMENSION(KIJL, MLSTHG) :: ENH
    
    REAL(KIND=JWRB) :: FTAIL
    REAL(KIND=JWRB) :: FKLAMP
    REAL(KIND=JWRB) :: GW1
    REAL(KIND=JWRB) :: GW2
    REAL(KIND=JWRB) :: GW3
    REAL(KIND=JWRB) :: GW4
    REAL(KIND=JWRB) :: FKLAMPA
    REAL(KIND=JWRB) :: FKLAMPB
    REAL(KIND=JWRB) :: FKLAMP2
    REAL(KIND=JWRB) :: FKLAMP1
    REAL(KIND=JWRB) :: FKLAPA2
    REAL(KIND=JWRB) :: FKLAPB2
    REAL(KIND=JWRB) :: FKLAP12
    REAL(KIND=JWRB) :: FKLAP22
    REAL(KIND=JWRB) :: FKLAMM
    REAL(KIND=JWRB) :: FKLAMM1
    REAL(KIND=JWRB) :: GW5
    REAL(KIND=JWRB) :: GW6
    REAL(KIND=JWRB) :: GW7
    REAL(KIND=JWRB) :: GW8
    REAL(KIND=JWRB) :: FKLAMMA
    REAL(KIND=JWRB) :: FKLAMMB
    REAL(KIND=JWRB) :: FKLAMM2
    REAL(KIND=JWRB) :: FKLAMA2
    REAL(KIND=JWRB) :: FKLAMB2
    REAL(KIND=JWRB) :: FKLAM12
    REAL(KIND=JWRB) :: FKLAM22
    REAL(KIND=JWRB) :: SAP
    REAL(KIND=JWRB) :: SAM
    REAL(KIND=JWRB) :: FIJ
    REAL(KIND=JWRB) :: FAD1
    REAL(KIND=JWRB) :: FAD2
    REAL(KIND=JWRB) :: FCEN
    
    REAL(KIND=JWRB), INTENT(INOUT) :: XNU(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(INOUT) :: SIG_TH(KIJL, NCHNK)
    REAL(KIND=JWRB) :: FTEMP
    REAL(KIND=JWRB) :: AD
    REAL(KIND=JWRB) :: DELAD
    REAL(KIND=JWRB) :: DELAP
    REAL(KIND=JWRB) :: DELAM
    REAL(KIND=JWRB) :: ENHFR
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: AF11(MFRSTLW:MLSTHG)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: BATHYMAX
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: COSTH(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DAL1
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DAL2
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: DFIM(NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: DFIMFR(NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: DFIMFR2(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DKMAX
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: FKLAM(MFRSTLW:MLSTHG)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: FKLAM1(MFRSTLW:MLSTHG)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: FKLAP(MFRSTLW:MLSTHG)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: FKLAP1(MFRSTLW:MLSTHG)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: FR(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: FRATIO
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: GM1
    INTEGER(KIND=JWIM), INTENT(IN), DEVICE :: IKM(MFRSTLW:MLSTHG)
    INTEGER(KIND=JWIM), INTENT(IN), DEVICE :: IKM1(MFRSTLW:MLSTHG)
    INTEGER(KIND=JWIM), INTENT(IN), DEVICE :: IKP(MFRSTLW:MLSTHG)
    INTEGER(KIND=JWIM), INTENT(IN), DEVICE :: IKP1(MFRSTLW:MLSTHG)
    INTEGER(KIND=JWIM), INTENT(IN), DEVICE :: INLCOEF(NINL, 1:MLSTHG)
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ISNONLIN
    INTEGER(KIND=JWIM), INTENT(IN), DEVICE :: K11W(NANG_loki_param, 2)
    INTEGER(KIND=JWIM), INTENT(IN), DEVICE :: K1W(NANG_loki_param, 2)
    INTEGER(KIND=JWIM), INTENT(IN), DEVICE :: K21W(NANG_loki_param, 2)
    INTEGER(KIND=JWIM), INTENT(IN), DEVICE :: K2W(NANG_loki_param, 2)
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KFRH
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: MFRSTLW
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: MLSTHG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: RNLCOEF(NRNL, 1:MLSTHG)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: SINTH(NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: TH(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WETAIL
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WP1TAIL
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WP2TAIL
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: XKDMIN
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: ZPIFR(NFRE_loki_param)
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    
    
    ! ----------------------------------------------------------------------
    
    
    
    !*    1. SHALLOW WATER SCALING
    !        ---------------------
    
    SELECT CASE (ISNONLIN)
    CASE (0)
      
      ENHFR = MAX(0.75_JWRB*DEPTH(IJ, ICHNK)*AKMEAN(IJ), 0.5_JWRB)
      ENHFR = 1.0_JWRB + (5.5_JWRB / ENHFR)*(1.0_JWRB - .833_JWRB*ENHFR)*EXP(-1.25_JWRB*ENHFR)
      DO MC=1,MLSTHG
        ENH(IJ, MC) = ENHFR
      END DO
      
      
    CASE (1)
      
      DO MC=1,NFRE
        ENH(IJ, MC) = MAX(MIN(ENH_MAX, TRANSF_CUF_HOIST_NEW(WAVNUM(IJ, MC, ICHNK), DEPTH(IJ, ICHNK), DKMAX, G)), ENH_MIN)
      END DO
      DO MC=NFRE + 1,MLSTHG
        XK = GM1*(ZPIFR(NFRE)*FRATIO**(MC - NFRE))**2
        ENH(IJ, MC) = MAX(MIN(ENH_MAX, TRANSF_CUF_HOIST_NEW(XK, DEPTH(IJ, ICHNK), DKMAX, G)), ENH_MIN)
      END DO
      
      
    CASE (2)
      CALL PEAK_ANG_CUF_HOIST_NEW(KIJS, KIJL, FL1(:, :, :, :), XNU(:, ICHNK), SIG_TH(:, ICHNK), COSTH(:), DELTH, DFIM(:),  &
      & DFIMFR(:), DFIMFR2(:), FR(:), FRATIO, NANG, NFRE, SINTH(:), TH(:), WETAIL, WP1TAIL, WP2TAIL, ICHNK, NCHNK, IJ)
      
      DO MC=1,NFRE
        ENH(IJ, MC) = TRANSF_SNL_CUF_HOIST_NEW(WAVNUM(IJ, MC, ICHNK), DEPTH(IJ, ICHNK), XNU(IJ, ICHNK), SIG_TH(IJ, ICHNK), BATHYMAX,  &
        & DKMAX, G, XKDMIN)
      END DO
      DO MC=NFRE + 1,MLSTHG
        XK = GM1*(ZPIFR(NFRE)*FRATIO**(MC - NFRE))**2
        ENH(IJ, MC) = TRANSF_SNL_CUF_HOIST_NEW(XK, DEPTH(IJ, ICHNK), XNU(IJ, ICHNK), SIG_TH(IJ, ICHNK), BATHYMAX, DKMAX, G, XKDMIN)
      END DO
      
    END SELECT
    
    
    !*    2. FREQUENCY LOOP.
    !        ---------------
    
    MFR1STFR = -MFRSTLW + 1
    MFRLSTFR = NFRE - KFRH + MFR1STFR
    
    
    DO MC=1,MLSTHG
      MP = IKP(MC)
      MP1 = IKP1(MC)
      MM = IKM(MC)
      MM1 = IKM1(MC)
      IC = INLCOEF(1, MC)
      IP = INLCOEF(2, MC)
      IP1 = INLCOEF(3, MC)
      IM = INLCOEF(4, MC)
      IM1 = INLCOEF(5, MC)
      
      FTAIL = RNLCOEF(1, MC)
      
      FKLAMP = FKLAP(MC)
      FKLAMP1 = FKLAP1(MC)
      GW1 = RNLCOEF(2, MC)
      GW2 = RNLCOEF(3, MC)
      GW3 = RNLCOEF(4, MC)
      GW4 = RNLCOEF(5, MC)
      FKLAMPA = RNLCOEF(6, MC)
      FKLAMPB = RNLCOEF(7, MC)
      FKLAMP2 = RNLCOEF(8, MC)
      FKLAMP1 = RNLCOEF(9, MC)
      FKLAPA2 = RNLCOEF(10, MC)
      FKLAPB2 = RNLCOEF(11, MC)
      FKLAP12 = RNLCOEF(12, MC)
      FKLAP22 = RNLCOEF(13, MC)
      
      FKLAMM = FKLAM(MC)
      FKLAMM1 = FKLAM1(MC)
      GW5 = RNLCOEF(14, MC)
      GW6 = RNLCOEF(15, MC)
      GW7 = RNLCOEF(16, MC)
      GW8 = RNLCOEF(17, MC)
      FKLAMMA = RNLCOEF(18, MC)
      FKLAMMB = RNLCOEF(19, MC)
      FKLAMM2 = RNLCOEF(20, MC)
      FKLAMM1 = RNLCOEF(21, MC)
      FKLAMA2 = RNLCOEF(22, MC)
      FKLAMB2 = RNLCOEF(23, MC)
      FKLAM12 = RNLCOEF(24, MC)
      FKLAM22 = RNLCOEF(25, MC)
      
      FTEMP = AF11(MC)*ENH(IJ, MC)
      
      
      IF (MC > MFR1STFR .and. MC < MFRLSTFR) THEN
        !       the interactions for MC are all within the fully resolved spectral domain
        
        DO KH=1,2
          DO K=1,NANG
            K1 = K1W(K, KH)
            K2 = K2W(K, KH)
            K11 = K11W(K, KH)
            K21 = K21W(K, KH)
            
            !*    2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND
            !*            DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
            !             ----------------------------------------------
            SAP = GW1*FL1(IJ, K1, IP, ICHNK) + GW2*FL1(IJ, K11, IP, ICHNK) + GW3*FL1(IJ, K1, IP1, ICHNK) + GW4*FL1(IJ, K11, IP1,  &
            & ICHNK)
            SAM = GW5*FL1(IJ, K2, IM, ICHNK) + GW6*FL1(IJ, K21, IM, ICHNK) + GW7*FL1(IJ, K2, IM1, ICHNK) + GW8*FL1(IJ, K21, IM1,  &
            & ICHNK)
            !!!! not needed ftail always=1.                FIJ = FL1(IJ,K  ,IC )*FTAIL
            FIJ = FL1(IJ, K, IC, ICHNK)
            FAD1 = FIJ*(SAP + SAM)
            FAD2 = FAD1 - 2.0_JWRB*SAP*SAM
            FAD1 = FAD1 + FAD2
            FCEN = FTEMP*FIJ
            AD = FAD2*FCEN
            DELAD = FAD1*FTEMP
            DELAP = (FIJ - 2.0_JWRB*SAM)*DAL1*FCEN
            DELAM = (FIJ - 2.0_JWRB*SAP)*DAL2*FCEN
            
            SL(IJ, K, MC) = SL(IJ, K, MC) - 2.0_JWRB*AD
            FLD(IJ, K, MC) = FLD(IJ, K, MC) - 2.0_JWRB*DELAD
            SL(IJ, K2, MM) = SL(IJ, K2, MM) + AD*FKLAMM1
            FLD(IJ, K2, MM) = FLD(IJ, K2, MM) + DELAM*FKLAM12
            SL(IJ, K21, MM) = SL(IJ, K21, MM) + AD*FKLAMM2
            FLD(IJ, K21, MM) = FLD(IJ, K21, MM) + DELAM*FKLAM22
            SL(IJ, K2, MM1) = SL(IJ, K2, MM1) + AD*FKLAMMA
            FLD(IJ, K2, MM1) = FLD(IJ, K2, MM1) + DELAM*FKLAMA2
            SL(IJ, K21, MM1) = SL(IJ, K21, MM1) + AD*FKLAMMB
            FLD(IJ, K21, MM1) = FLD(IJ, K21, MM1) + DELAM*FKLAMB2
            SL(IJ, K1, MP) = SL(IJ, K1, MP) + AD*FKLAMP1
            FLD(IJ, K1, MP) = FLD(IJ, K1, MP) + DELAP*FKLAP12
            SL(IJ, K11, MP) = SL(IJ, K11, MP) + AD*FKLAMP2
            FLD(IJ, K11, MP) = FLD(IJ, K11, MP) + DELAP*FKLAP22
            SL(IJ, K1, MP1) = SL(IJ, K1, MP1) + AD*FKLAMPA
            FLD(IJ, K1, MP1) = FLD(IJ, K1, MP1) + DELAP*FKLAPA2
            SL(IJ, K11, MP1) = SL(IJ, K11, MP1) + AD*FKLAMPB
            FLD(IJ, K11, MP1) = FLD(IJ, K11, MP1) + DELAP*FKLAPB2
          END DO
        END DO
        
      ELSE IF (MC >= MFRLSTFR) THEN
        DO KH=1,2
          DO K=1,NANG
            K1 = K1W(K, KH)
            K2 = K2W(K, KH)
            K11 = K11W(K, KH)
            K21 = K21W(K, KH)
            
            SAP = GW1*FL1(IJ, K1, IP, ICHNK) + GW2*FL1(IJ, K11, IP, ICHNK) + GW3*FL1(IJ, K1, IP1, ICHNK) + GW4*FL1(IJ, K11, IP1,  &
            & ICHNK)
            SAM = GW5*FL1(IJ, K2, IM, ICHNK) + GW6*FL1(IJ, K21, IM, ICHNK) + GW7*FL1(IJ, K2, IM1, ICHNK) + GW8*FL1(IJ, K21, IM1,  &
            & ICHNK)
            FIJ = FL1(IJ, K, IC, ICHNK)*FTAIL
            FAD1 = FIJ*(SAP + SAM)
            FAD2 = FAD1 - 2.0_JWRB*SAP*SAM
            FAD1 = FAD1 + FAD2
            FCEN = FTEMP*FIJ
            AD = FAD2*FCEN
            DELAD = FAD1*FTEMP
            DELAP = (FIJ - 2.0_JWRB*SAM)*DAL1*FCEN
            DELAM = (FIJ - 2.0_JWRB*SAP)*DAL2*FCEN
            
            SL(IJ, K2, MM) = SL(IJ, K2, MM) + AD*FKLAMM1
            FLD(IJ, K2, MM) = FLD(IJ, K2, MM) + DELAM*FKLAM12
            SL(IJ, K21, MM) = SL(IJ, K21, MM) + AD*FKLAMM2
            FLD(IJ, K21, MM) = FLD(IJ, K21, MM) + DELAM*FKLAM22
            
            IF (MM1 <= NFRE) THEN
              SL(IJ, K2, MM1) = SL(IJ, K2, MM1) + AD*FKLAMMA
              FLD(IJ, K2, MM1) = FLD(IJ, K2, MM1) + DELAM*FKLAMA2
              SL(IJ, K21, MM1) = SL(IJ, K21, MM1) + AD*FKLAMMB
              FLD(IJ, K21, MM1) = FLD(IJ, K21, MM1) + DELAM*FKLAMB2
              
              IF (MC <= NFRE) THEN
                SL(IJ, K, MC) = SL(IJ, K, MC) - 2.0_JWRB*AD
                FLD(IJ, K, MC) = FLD(IJ, K, MC) - 2.0_JWRB*DELAD
                
                IF (MP <= NFRE) THEN
                  SL(IJ, K1, MP) = SL(IJ, K1, MP) + AD*FKLAMP1
                  FLD(IJ, K1, MP) = FLD(IJ, K1, MP) + DELAP*FKLAP12
                  SL(IJ, K11, MP) = SL(IJ, K11, MP) + AD*FKLAMP2
                  FLD(IJ, K11, MP) = FLD(IJ, K11, MP) + DELAP*FKLAP22
                  
                  IF (MP1 <= NFRE) THEN
                    SL(IJ, K1, MP1) = SL(IJ, K1, MP1) + AD*FKLAMPA
                    FLD(IJ, K1, MP1) = FLD(IJ, K1, MP1) + DELAP*FKLAPA2
                    SL(IJ, K11, MP1) = SL(IJ, K11, MP1) + AD*FKLAMPB
                    FLD(IJ, K11, MP1) = FLD(IJ, K11, MP1) + DELAP*FKLAPB2
                  END IF
                END IF
              END IF
            END IF
          END DO
        END DO
        
      ELSE
        
        DO KH=1,2
          DO K=1,NANG
            K1 = K1W(K, KH)
            K2 = K2W(K, KH)
            K11 = K11W(K, KH)
            K21 = K21W(K, KH)
            
            SAP = GW1*FL1(IJ, K1, IP, ICHNK) + GW2*FL1(IJ, K11, IP, ICHNK) + GW3*FL1(IJ, K1, IP1, ICHNK) + GW4*FL1(IJ, K11, IP1,  &
            & ICHNK)
            SAM = GW5*FL1(IJ, K2, IM, ICHNK) + GW6*FL1(IJ, K21, IM, ICHNK) + GW7*FL1(IJ, K2, IM1, ICHNK) + GW8*FL1(IJ, K21, IM1,  &
            & ICHNK)
            FIJ = FL1(IJ, K, IC, ICHNK)*FTAIL
            FAD1 = FIJ*(SAP + SAM)
            FAD2 = FAD1 - 2.0_JWRB*SAP*SAM
            FAD1 = FAD1 + FAD2
            FCEN = FTEMP*FIJ
            AD = FAD2*FCEN
            DELAD = FAD1*FTEMP
            DELAP = (FIJ - 2.0_JWRB*SAM)*DAL1*FCEN
            DELAM = (FIJ - 2.0_JWRB*SAP)*DAL2*FCEN
            
            IF (MM1 >= 1) THEN
              SL(IJ, K2, MM1) = SL(IJ, K2, MM1) + AD*FKLAMMA
              FLD(IJ, K2, MM1) = FLD(IJ, K2, MM1) + DELAM*FKLAMA2
              SL(IJ, K21, MM1) = SL(IJ, K21, MM1) + AD*FKLAMMB
              FLD(IJ, K21, MM1) = FLD(IJ, K21, MM1) + DELAM*FKLAMB2
            END IF
            
            SL(IJ, K, MC) = SL(IJ, K, MC) - 2.0_JWRB*AD
            FLD(IJ, K, MC) = FLD(IJ, K, MC) - 2.0_JWRB*DELAD
            SL(IJ, K1, MP) = SL(IJ, K1, MP) + AD*FKLAMP1
            FLD(IJ, K1, MP) = FLD(IJ, K1, MP) + DELAP*FKLAP12
            SL(IJ, K11, MP) = SL(IJ, K11, MP) + AD*FKLAMP2
            FLD(IJ, K11, MP) = FLD(IJ, K11, MP) + DELAP*FKLAP22
            SL(IJ, K1, MP1) = SL(IJ, K1, MP1) + AD*FKLAMPA
            FLD(IJ, K1, MP1) = FLD(IJ, K1, MP1) + DELAP*FKLAPA2
            SL(IJ, K11, MP1) = SL(IJ, K11, MP1) + AD*FKLAMPB
            FLD(IJ, K11, MP1) = FLD(IJ, K11, MP1) + DELAP*FKLAPB2
          END DO
        END DO
        
      END IF
      
      !*    BRANCH BACK TO 2. FOR NEXT FREQUENCY.
      
    END DO
    
    
    
  END SUBROUTINE SNONLIN_CUF_HOIST_NEW
END MODULE SNONLIN_CUF_HOIST_NEW_MOD
