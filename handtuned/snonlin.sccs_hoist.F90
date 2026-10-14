! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SNONLIN_LOKI_MOD
  CONTAINS
  SUBROUTINE SNONLIN_LOKI (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN, IJ, ENH)

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
    USE TRANSF_SNL_LOKI_MOD, ONLY: TRANSF_SNL_LOKI
    USE TRANSF_LOKI_MOD, ONLY: TRANSF_LOKI
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWFRED, ONLY: FR, FRATIO, ZPIFR, DFIM, DFIMFR, DFIMOFR, DFIMFR2, DELTH, TH, SINTH, COSTH, WETAIL, WP1TAIL, WP2TAIL
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWINDN, ONLY: IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W, AF11, FKLAP, FKLAP1, FKLAM, FKLAM1, DAL1, DAL2, MFRSTLW,  &
    & MLSTHG, KFRH, INLCOEF, RNLCOEF
    USE YOWSTAT, ONLY: ISNONLIN
    USE YOWPCONS, ONLY: GM1

    USE YOMHOOK, ONLY: LHOOK, JPHOOK
    ! ----------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: DEPTH, AKMEAN


    INTEGER(KIND=JWIM) :: K, M, MC, KH, K1, K2, K11, K21
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: MP, MP1, MM, MM1, IC, IP, IP1, IM, IM1
    INTEGER(KIND=JWIM) :: MFR1STFR, MFRLSTFR

    REAL(KIND=JWRB), PARAMETER :: ENH_MAX = 10.0_JWRB
    REAL(KIND=JWRB), PARAMETER :: ENH_MIN = 0.1_JWRB    ! to prevent ENH to become too small
    REAL(KIND=JWRB) :: XK
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, MLSTHG) :: ENH

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    REAL(KIND=JWRB) :: FTAIL, FKLAMP, GW1, GW2, GW3, GW4
    REAL(KIND=JWRB) :: FKLAMPA, FKLAMPB, FKLAMP2, FKLAMP1, FKLAPA2
    REAL(KIND=JWRB) :: FKLAPB2, FKLAP12, FKLAP22, FKLAMM, FKLAMM1
    REAL(KIND=JWRB) :: GW5, GW6, GW7, GW8, FKLAMMA, FKLAMMB, FKLAMM2
    REAL(KIND=JWRB) :: FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
    REAL(KIND=JWRB) :: SAP, SAM, FIJ, FAD1, FAD2, FCEN

    REAL(KIND=JWRB) :: XNU, SIG_TH
    REAL(KIND=JWRB) :: FTEMP, AD, DELAD, DELAP, DELAM, ENHFR
    REAL(KIND=JWRB) :: loki_temp_30
    REAL(KIND=JWRB) :: loki_temp_14
    REAL(KIND=JWRB) :: loki_temp_16
    REAL(KIND=JWRB) :: loki_temp_28
    REAL(KIND=JWRB) :: loki_temp_59
    REAL(KIND=JWRB) :: loki_temp_56
    REAL(KIND=JWRB) :: loki_temp_10
    REAL(KIND=JWRB) :: loki_temp_17
    REAL(KIND=JWRB) :: loki_temp_15
    REAL(KIND=JWRB) :: loki_temp_6
    REAL(KIND=JWRB) :: loki_temp_13
    REAL(KIND=JWRB) :: loki_temp_57
    REAL(KIND=JWRB) :: loki_temp_21
    REAL(KIND=JWRB) :: loki_temp_3
    REAL(KIND=JWRB) :: loki_temp_32
    REAL(KIND=JWRB) :: loki_temp_58
    REAL(KIND=JWRB) :: loki_temp_55
    REAL(KIND=JWRB) :: loki_temp_62
    REAL(KIND=JWRB) :: loki_temp_19
    REAL(KIND=JWRB) :: loki_temp_7
    REAL(KIND=JWRB) :: loki_temp_66
    REAL(KIND=JWRB) :: loki_temp_23
    REAL(KIND=JWRB) :: loki_temp_31
    REAL(KIND=JWRB) :: loki_temp_35
    REAL(KIND=JWRB) :: loki_temp_34
    REAL(KIND=JWRB) :: loki_temp_64
    REAL(KIND=JWRB) :: loki_temp_5
    REAL(KIND=JWRB) :: loki_temp_9
    REAL(KIND=JWRB) :: loki_temp_29
    REAL(KIND=JWRB) :: loki_temp_63
    REAL(KIND=JWRB) :: loki_temp_26
    REAL(KIND=JWRB) :: loki_temp_61
    REAL(KIND=JWRB) :: loki_temp_20
    REAL(KIND=JWRB) :: loki_temp_27
    REAL(KIND=JWRB) :: loki_temp_12
    REAL(KIND=JWRB) :: loki_temp_4
    REAL(KIND=JWRB) :: loki_temp_0
    REAL(KIND=JWRB) :: loki_temp_8
    REAL(KIND=JWRB) :: loki_temp_33
    REAL(KIND=JWRB) :: loki_temp_11
    REAL(KIND=JWRB) :: loki_temp_65
    REAL(KIND=JWRB) :: loki_temp_22
    REAL(KIND=JWRB) :: loki_temp_60
    REAL(KIND=JWRB) :: loki_temp_24
    REAL(KIND=JWRB) :: loki_temp_54
    REAL(KIND=JWRB) :: loki_temp_1
    REAL(KIND=JWRB) :: loki_temp_25
    REAL(KIND=JWRB) :: loki_temp_18
    REAL(KIND=JWRB) :: loki_temp_67
    REAL(KIND=JWRB) :: loki_temp_2
    INTEGER(KIND=JWIM) :: NSH
    INTEGER(KIND=JWIM) :: PEAK_ANG_M, PEAK_ANG_K
    INTEGER(KIND=JWIM) :: MMAX
    INTEGER(KIND=JWIM) :: MMSTART, MMSTOP
    REAL(KIND=JWRB), PARAMETER :: CONST_SIG = 1.0_JWRB
    REAL(KIND=JWRB) :: R1
    REAL(KIND=JWRB) :: DELT25, COEF_FR, COEF_FR2
    REAL(KIND=JPHOOK) :: PEAK_ANG_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: ZEPSILON
    REAL(KIND=JWRB) :: SUM0, SUM1, SUM2, XMAX, TEMP
    REAL(KIND=JWRB) :: THMEAN, SUM_S, SUM_C
!$acc routine seq
!$acc data present( FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN, ENH )


    ! ----------------------------------------------------------------------



    !*    1. SHALLOW WATER SCALING
    !        ---------------------

    SELECT CASE (ISNONLIN)
    CASE (0)
      ENHFR = MAX(0.75_JWRB*DEPTH(IJ)*AKMEAN(IJ), 0.5_JWRB)
      ENHFR = 1.0_JWRB + (5.5_JWRB / ENHFR)*(1.0_JWRB - .833_JWRB*ENHFR)*EXP(-1.25_JWRB*ENHFR)
!$acc loop seq
      DO MC=1,MLSTHG
        ENH(IJ, MC) = ENHFR
      END DO

    CASE (1)
!$acc loop seq
      DO MC=1,NFRE
        ENH(IJ, MC) = MAX(MIN(ENH_MAX, TRANSF_LOKI(WAVNUM(IJ, MC), DEPTH(IJ))), ENH_MIN)
      END DO
!$acc loop seq
      DO MC=NFRE + 1,MLSTHG
        XK = GM1*(ZPIFR(NFRE)*FRATIO**(MC - NFRE))**2
        ENH(IJ, MC) = MAX(MIN(ENH_MAX, TRANSF_LOKI(XK, DEPTH(IJ))), ENH_MIN)
      END DO

    CASE (2)
      ! [Loki] inlined child subroutine: PEAK_ANG
      ! =========================================

      ! ----------------------------------------------------------------------

      !***  1. DETERMINE L-H SPECTRAL WIDTH OF THE 2-D SPECTRUM.
      !     ---------------------------------------------------

      ZEPSILON = 10._JWRB*EPSILON(ZEPSILON)

      NSH = 1 + INT(LOG(1.5_JWRB) / LOG(FRATIO))

      SUM0 = ZEPSILON
      SUM1 = 0._JWRB
      SUM2 = 0._JWRB

!$acc loop seq
      DO PEAK_ANG_M=1,NFRE
        PEAK_ANG_K = 1
        TEMP = FL1(IJ, PEAK_ANG_K, PEAK_ANG_M)
!$acc loop seq
        DO PEAK_ANG_K=2,NANG
          TEMP = TEMP + FL1(IJ, PEAK_ANG_K, PEAK_ANG_M)
        END DO
        SUM0 = SUM0 + TEMP*DFIM(PEAK_ANG_M)
        SUM1 = SUM1 + TEMP*DFIMFR(PEAK_ANG_M)
        SUM2 = SUM2 + TEMP*DFIMFR2(PEAK_ANG_M)
      END DO

      !     ADD TAIL CORRECTIONS
      DELT25 = WETAIL*FR(NFRE)*DELTH
      COEF_FR = WP1TAIL*DELTH*FR(NFRE)**2
      COEF_FR2 = WP2TAIL*DELTH*FR(NFRE)**3
      SUM0 = SUM0 + DELT25*TEMP
      SUM1 = SUM1 + COEF_FR*TEMP
      SUM2 = SUM2 + COEF_FR2*TEMP

      IF (SUM0 > ZEPSILON) THEN
        XNU = SQRT(MAX(ZEPSILON, SUM2*SUM0 / SUM1**2 - 1._JWRB))
      ELSE
        XNU = ZEPSILON
      END IF

      !***  2. DETERMINE ANGULAR WIDTH OF THE 2-D SPECTRUM.
      !     ----------------------------------------------

      !     MAX OF 2-D SPECTRUM
      XMAX = 0._JWRB
      MMAX = 2

!$acc loop seq
      DO PEAK_ANG_M=2,NFRE - 1
!$acc loop seq
        DO PEAK_ANG_K=1,NANG
          IF (FL1(IJ, PEAK_ANG_K, PEAK_ANG_M) > XMAX) THEN
            MMAX = PEAK_ANG_M
            XMAX = FL1(IJ, PEAK_ANG_K, PEAK_ANG_M)
          END IF
        END DO
      END DO

      SUM1 = ZEPSILON
      SUM2 = 0._JWRB

      MMSTART = MAX(1, MMAX - NSH)
      MMSTOP = MIN(NFRE, MMAX + NSH)

      SUM_S = 0._JWRB
      SUM_C = ZEPSILON
!$acc loop seq
      DO PEAK_ANG_M=MMSTART,MMSTOP
!$acc loop seq
        DO PEAK_ANG_K=1,NANG
          SUM_S = SUM_S + SINTH(PEAK_ANG_K)*FL1(IJ, PEAK_ANG_K, PEAK_ANG_M)
          SUM_C = SUM_C + COSTH(PEAK_ANG_K)*FL1(IJ, PEAK_ANG_K, PEAK_ANG_M)
        END DO
        THMEAN = ATAN2(SUM_S, SUM_C)
!$acc loop seq
        DO PEAK_ANG_K=1,NANG
          SUM1 = SUM1 + FL1(IJ, PEAK_ANG_K, PEAK_ANG_M)*DFIM(PEAK_ANG_M)
          SUM2 = SUM2 + COS(TH(PEAK_ANG_K) - THMEAN)*FL1(IJ, PEAK_ANG_K, PEAK_ANG_M)*DFIM(PEAK_ANG_M)
        END DO
      END DO

      IF (SUM1 > ZEPSILON) THEN
        R1 = SUM2 / SUM1
        SIG_TH = CONST_SIG*SQRT(2._JWRB*(1._JWRB - R1))
      ELSE
        SIG_TH = 0._JWRB
      END IF


      ! =========================================
!$acc loop seq
      DO MC=1,NFRE
        ENH(IJ, MC) = TRANSF_SNL_LOKI(WAVNUM(IJ, MC), DEPTH(IJ), XNU, SIG_TH)
      END DO
!$acc loop seq
      DO MC=NFRE + 1,MLSTHG
        XK = GM1*(ZPIFR(NFRE)*FRATIO**(MC - NFRE))**2
        ENH(IJ, MC) = TRANSF_SNL_LOKI(XK, DEPTH(IJ), XNU, SIG_TH)
      END DO
    END SELECT


    !*    2. FREQUENCY LOOP.
    !        ---------------

    MFR1STFR = -MFRSTLW + 1
    MFRLSTFR = NFRE - KFRH + MFR1STFR

!$acc loop seq
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

        ! Loki loop-interchange (kh, k <--> k, kh)
!$acc loop seq
        DO K=1,NANG
!$acc loop seq
          DO KH=1,2
            K1 = K1W(K, KH)
            K2 = K2W(K, KH)
            K11 = K11W(K, KH)
            K21 = K21W(K, KH)

            !*    2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND
            !*            DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
            !             ----------------------------------------------
            SAP = GW1*FL1(IJ, K1, IP) + GW2*FL1(IJ, K11, IP) + GW3*FL1(IJ, K1, IP1) + GW4*FL1(IJ, K11, IP1)
            SAM = GW5*FL1(IJ, K2, IM) + GW6*FL1(IJ, K21, IM) + GW7*FL1(IJ, K2, IM1) + GW8*FL1(IJ, K21, IM1)
            !!!! not needed ftail always=1.                FIJ = FL1(IJ,K  ,IC )*FTAIL
            FIJ = FL1(IJ, K, IC)
            FAD1 = FIJ*(SAP + SAM)
            FAD2 = FAD1 - 2.0_JWRB*SAP*SAM
            FAD1 = FAD1 + FAD2
            FCEN = FTEMP*FIJ
            AD = FAD2*FCEN
            DELAD = FAD1*FTEMP
            DELAP = (FIJ - 2.0_JWRB*SAM)*DAL1*FCEN
            DELAM = (FIJ - 2.0_JWRB*SAP)*DAL2*FCEN

!$loki split-read-write
            loki_temp_0 = SL(IJ, K, MC) - 2.0_JWRB*AD
            loki_temp_1 = FLD(IJ, K, MC) - 2.0_JWRB*DELAD
            loki_temp_2 = SL(IJ, K2, MM) + AD*FKLAMM1
            loki_temp_3 = FLD(IJ, K2, MM) + DELAM*FKLAM12
            loki_temp_4 = SL(IJ, K21, MM) + AD*FKLAMM2
            loki_temp_5 = FLD(IJ, K21, MM) + DELAM*FKLAM22
            loki_temp_6 = SL(IJ, K2, MM1) + AD*FKLAMMA
            loki_temp_7 = FLD(IJ, K2, MM1) + DELAM*FKLAMA2
            loki_temp_8 = SL(IJ, K21, MM1) + AD*FKLAMMB
            loki_temp_9 = FLD(IJ, K21, MM1) + DELAM*FKLAMB2
            loki_temp_10 = SL(IJ, K1, MP) + AD*FKLAMP1
            loki_temp_11 = FLD(IJ, K1, MP) + DELAP*FKLAP12
            loki_temp_12 = SL(IJ, K11, MP) + AD*FKLAMP2
            loki_temp_13 = FLD(IJ, K11, MP) + DELAP*FKLAP22
            loki_temp_14 = SL(IJ, K1, MP1) + AD*FKLAMPA
            loki_temp_15 = FLD(IJ, K1, MP1) + DELAP*FKLAPA2
            loki_temp_16 = SL(IJ, K11, MP1) + AD*FKLAMPB
            loki_temp_17 = FLD(IJ, K11, MP1) + DELAP*FKLAPB2
            SL(IJ, K, MC) = loki_temp_0
            FLD(IJ, K, MC) = loki_temp_1
            SL(IJ, K2, MM) = loki_temp_2
            FLD(IJ, K2, MM) = loki_temp_3
            SL(IJ, K21, MM) = loki_temp_4
            FLD(IJ, K21, MM) = loki_temp_5
            SL(IJ, K2, MM1) = loki_temp_6
            FLD(IJ, K2, MM1) = loki_temp_7
            SL(IJ, K21, MM1) = loki_temp_8
            FLD(IJ, K21, MM1) = loki_temp_9
            SL(IJ, K1, MP) = loki_temp_10
            FLD(IJ, K1, MP) = loki_temp_11
            SL(IJ, K11, MP) = loki_temp_12
            FLD(IJ, K11, MP) = loki_temp_13
            SL(IJ, K1, MP1) = loki_temp_14
            FLD(IJ, K1, MP1) = loki_temp_15
            SL(IJ, K11, MP1) = loki_temp_16
            FLD(IJ, K11, MP1) = loki_temp_17
!$loki end split-read-write
          END DO
        END DO

      ELSE IF (MC >= MFRLSTFR) THEN
        ! Loki loop-interchange (kh, k <--> k, kh)
!$acc loop seq
        DO K=1,NANG
!$acc loop seq
          DO KH=1,2
            K1 = K1W(K, KH)
            K2 = K2W(K, KH)
            K11 = K11W(K, KH)
            K21 = K21W(K, KH)

            SAP = GW1*FL1(IJ, K1, IP) + GW2*FL1(IJ, K11, IP) + GW3*FL1(IJ, K1, IP1) + GW4*FL1(IJ, K11, IP1)
            SAM = GW5*FL1(IJ, K2, IM) + GW6*FL1(IJ, K21, IM) + GW7*FL1(IJ, K2, IM1) + GW8*FL1(IJ, K21, IM1)
            FIJ = FL1(IJ, K, IC)*FTAIL
            FAD1 = FIJ*(SAP + SAM)
            FAD2 = FAD1 - 2.0_JWRB*SAP*SAM
            FAD1 = FAD1 + FAD2
            FCEN = FTEMP*FIJ
            AD = FAD2*FCEN
            DELAD = FAD1*FTEMP
            DELAP = (FIJ - 2.0_JWRB*SAM)*DAL1*FCEN
            DELAM = (FIJ - 2.0_JWRB*SAP)*DAL2*FCEN

!$loki split-read-write
            loki_temp_18 = SL(IJ, K2, MM) + AD*FKLAMM1
            loki_temp_19 = FLD(IJ, K2, MM) + DELAM*FKLAM12
            loki_temp_20 = SL(IJ, K21, MM) + AD*FKLAMM2
            loki_temp_21 = FLD(IJ, K21, MM) + DELAM*FKLAM22

            IF (MM1 <= NFRE) THEN
              loki_temp_22 = SL(IJ, K2, MM1) + AD*FKLAMMA
              loki_temp_23 = FLD(IJ, K2, MM1) + DELAM*FKLAMA2
              loki_temp_24 = SL(IJ, K21, MM1) + AD*FKLAMMB
              loki_temp_25 = FLD(IJ, K21, MM1) + DELAM*FKLAMB2

              IF (MC <= NFRE) THEN
                loki_temp_26 = SL(IJ, K, MC) - 2.0_JWRB*AD
                loki_temp_27 = FLD(IJ, K, MC) - 2.0_JWRB*DELAD

                IF (MP <= NFRE) THEN
                  loki_temp_28 = SL(IJ, K1, MP) + AD*FKLAMP1
                  loki_temp_29 = FLD(IJ, K1, MP) + DELAP*FKLAP12
                  loki_temp_30 = SL(IJ, K11, MP) + AD*FKLAMP2
                  loki_temp_31 = FLD(IJ, K11, MP) + DELAP*FKLAP22

                  IF (MP1 <= NFRE) THEN
                    loki_temp_32 = SL(IJ, K1, MP1) + AD*FKLAMPA
                    loki_temp_33 = FLD(IJ, K1, MP1) + DELAP*FKLAPA2
                    loki_temp_34 = SL(IJ, K11, MP1) + AD*FKLAMPB
                    loki_temp_35 = FLD(IJ, K11, MP1) + DELAP*FKLAPB2
                  END IF
                END IF
              END IF
            END IF
            SL(IJ, K2, MM) = loki_temp_18
            FLD(IJ, K2, MM) = loki_temp_19
            SL(IJ, K21, MM) = loki_temp_20
            FLD(IJ, K21, MM) = loki_temp_21
            IF (MM1 <= NFRE) THEN
              SL(IJ, K2, MM1) = loki_temp_22
              FLD(IJ, K2, MM1) = loki_temp_23
              SL(IJ, K21, MM1) = loki_temp_24
              FLD(IJ, K21, MM1) = loki_temp_25
              IF (MC <= NFRE) THEN
                SL(IJ, K, MC) = loki_temp_26
                FLD(IJ, K, MC) = loki_temp_27
                IF (MP <= NFRE) THEN
                  SL(IJ, K1, MP) = loki_temp_28
                  FLD(IJ, K1, MP) = loki_temp_29
                  SL(IJ, K11, MP) = loki_temp_30
                  FLD(IJ, K11, MP) = loki_temp_31
                  IF (MP1 <= NFRE) THEN
                    SL(IJ, K1, MP1) = loki_temp_32
                    FLD(IJ, K1, MP1) = loki_temp_33
                    SL(IJ, K11, MP1) = loki_temp_34
                    FLD(IJ, K11, MP1) = loki_temp_35
                  END IF
                END IF
              END IF
            END IF
!$loki end split-read-write
          END DO
        END DO

      ELSE

        ! Loki loop-interchange (kh, k <--> k, kh)
!$acc loop seq
        DO K=1,NANG
!$acc loop seq
          DO KH=1,2
            K1 = K1W(K, KH)
            K2 = K2W(K, KH)
            K11 = K11W(K, KH)
            K21 = K21W(K, KH)

            SAP = GW1*FL1(IJ, K1, IP) + GW2*FL1(IJ, K11, IP) + GW3*FL1(IJ, K1, IP1) + GW4*FL1(IJ, K11, IP1)
            SAM = GW5*FL1(IJ, K2, IM) + GW6*FL1(IJ, K21, IM) + GW7*FL1(IJ, K2, IM1) + GW8*FL1(IJ, K21, IM1)
            FIJ = FL1(IJ, K, IC)*FTAIL
            FAD1 = FIJ*(SAP + SAM)
            FAD2 = FAD1 - 2.0_JWRB*SAP*SAM
            FAD1 = FAD1 + FAD2
            FCEN = FTEMP*FIJ
            AD = FAD2*FCEN
            DELAD = FAD1*FTEMP
            DELAP = (FIJ - 2.0_JWRB*SAM)*DAL1*FCEN
            DELAM = (FIJ - 2.0_JWRB*SAP)*DAL2*FCEN

!$loki split-read-write
            IF (MM1 >= 1) THEN
              loki_temp_54 = SL(IJ, K2, MM1) + AD*FKLAMMA
              loki_temp_55 = FLD(IJ, K2, MM1) + DELAM*FKLAMA2
              loki_temp_56 = SL(IJ, K21, MM1) + AD*FKLAMMB
              loki_temp_57 = FLD(IJ, K21, MM1) + DELAM*FKLAMB2
            END IF

            loki_temp_58 = SL(IJ, K, MC) - 2.0_JWRB*AD
            loki_temp_59 = FLD(IJ, K, MC) - 2.0_JWRB*DELAD
            loki_temp_60 = SL(IJ, K1, MP) + AD*FKLAMP1
            loki_temp_61 = FLD(IJ, K1, MP) + DELAP*FKLAP12
            loki_temp_62 = SL(IJ, K11, MP) + AD*FKLAMP2
            loki_temp_63 = FLD(IJ, K11, MP) + DELAP*FKLAP22
            loki_temp_64 = SL(IJ, K1, MP1) + AD*FKLAMPA
            loki_temp_65 = FLD(IJ, K1, MP1) + DELAP*FKLAPA2
            loki_temp_66 = SL(IJ, K11, MP1) + AD*FKLAMPB
            loki_temp_67 = FLD(IJ, K11, MP1) + DELAP*FKLAPB2
            IF (MM1 >= 1) THEN
              SL(IJ, K2, MM1) = loki_temp_54
              FLD(IJ, K2, MM1) = loki_temp_55
              SL(IJ, K21, MM1) = loki_temp_56
              FLD(IJ, K21, MM1) = loki_temp_57
            END IF
            SL(IJ, K, MC) = loki_temp_58
            FLD(IJ, K, MC) = loki_temp_59
            SL(IJ, K1, MP) = loki_temp_60
            FLD(IJ, K1, MP) = loki_temp_61
            SL(IJ, K11, MP) = loki_temp_62
            FLD(IJ, K11, MP) = loki_temp_63
            SL(IJ, K1, MP1) = loki_temp_64
            FLD(IJ, K1, MP1) = loki_temp_65
            SL(IJ, K11, MP1) = loki_temp_66
            FLD(IJ, K11, MP1) = loki_temp_67
!$loki end split-read-write
          END DO
        END DO

      END IF

      !*    BRANCH BACK TO 2. FOR NEXT FREQUENCY.

    END DO



!$acc end data
  END SUBROUTINE SNONLIN_LOKI
END MODULE SNONLIN_LOKI_MOD
