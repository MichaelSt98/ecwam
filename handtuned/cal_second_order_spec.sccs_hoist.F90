! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE CAL_SECOND_ORDER_SPEC_LOKI_MOD
  CONTAINS
  SUBROUTINE CAL_SECOND_ORDER_SPEC_LOKI (KIJS, KIJL, F1, WAVNUM, DEPTH, SIG, IJ, EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN, F3, PF1,  &
  & PF3, SECSPOM_LL, SECSPOM_OMEGA_EXT, SECSPOM_F2)

    !***  *CAL_SEC_ORDER_SPEC*   DETERMINES SECOND_ORDER SPECTRUM

    !     PETER JANSSEN JULY 2010

    !     PURPOSE.
    !     --------

    !              DETERMINATION OF SECOND-ORDER SPECTRUM

    !     INTERFACE.
    !     ----------
    !              *CALL*  *CAL_SEC_ORDER_SPEC(KIJS,KIJL,F1,WAVNUM,DEPTH,SIG)*

    !                       INPUT:
    !                            *KIJS*   - FIRST GRIDPOINT
    !                            *KIJL*   - LAST GRIDPOINT
    !                            *F1*     - 2-D FREE WAVE SPECTRUM (at input)
    !                            *WAVNUM* - WAVE NUMBER
    !                            *DEPTH*  - DEPTH ARRAY
    !                            *SIG*    - DIRECTION OF MAPPING:
    !                                       FORWARD: SIG = +1
    !                                       INVERSE: SIG = -1

    !                       OUTPUT:
    !                            *F1*    - 2-D SPECTRUM INCLUDING SECOND-ORDER
    !                                      CORRECTION (at output)

    !     METHOD.
    !     -------
    !              IS DESCRIBED IN JANSSEN (2009), JFM, 637, 1-44.

    !     EXTERNALS.
    !     ----------
    !              SECSPOM

    !-----------------------------------------------------------------------

    USE SECSPOM_LOKI_MOD, ONLY: SECSPOM_LOKI
    USE FKMEAN_LOKI_MOD, ONLY: FKMEAN_LOKI
    ![Loki::HoistVariablesTransformation] -------- Added hoisted temporary size imports -------------------------------
    USE YOWTABL, ONLY: NMAX
    ![Loki::HoistVariablesTransformation] ---------------------------------------
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWFRED, ONLY: FR, DFIM, DELTH, TH, FRATIO
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: G, PI, ZPI
    USE YOWSHAL, ONLY: NDEPTH, DEPTHA, DEPTHD
    USE YOWTABL, ONLY: MR, XMR, MA, XMA, NFREH, NANGH, OMEGA, DFDTH, THH, DELTHH, IM_P, IM_M, TA, TB, TC_QL, TT_4M, TT_4P
    USE YOWTEST, ONLY: IU06

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    !-----------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: F1
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: DEPTH
    REAL(KIND=JWRB), INTENT(IN) :: SIG

    INTEGER(KIND=JWIM) :: M, K, K0, M0, MP, KP, MM, KM, KL, KLL, ML
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ

    REAL(KIND=JWRB) :: FRAC, CO1, DEL, DELF, D1, D2, D3, D4, C1
    REAL(KIND=JWRB) :: C2, XM, XK, OMSTART, AREA, SUM, SUM1, SUM3, GAM_B_J, ZFAC
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: EMEAN
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: FMEAN
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: F1MEAN
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: AKMEAN
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: XKMEAN
    REAL(KIND=JWRB) :: EMAXL
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: F3
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANGH, NFREH) :: PF1
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANGH, NFREH) :: PF3
!$acc routine seq
    INTEGER(KIND=JWIM), INTENT(INOUT) :: SECSPOM_LL(NANGH, NANGH)
    REAL(KIND=JWRB), INTENT(INOUT) :: SECSPOM_OMEGA_EXT(NMAX)
    REAL(KIND=JWRB), INTENT(INOUT) :: SECSPOM_F2(KIJL, NANGH, NMAX)
!$acc data present( F1, WAVNUM, DEPTH, EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN, F3, PF1, PF3, SECSPOM_LL, SECSPOM_OMEGA_EXT,  &
!$acc & SECSPOM_F2 )

    !-----------------------------------------------------------------------


    !***  1. DETERMINE SECOND ORDER CORRECTION TO THE SPECTRUM
    !     ----------------------------------------------------

    GAM_B_J = 0.6_JWRB
    ZFAC = GAM_B_J**2 / 16.0_JWRB
    FRAC = FRATIO - 1.0_JWRB
    OMSTART = ZPI*FR(1)

    CALL FKMEAN_LOKI(KIJS, KIJL, F1, WAVNUM, EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN, IJ=IJ)

    !***  1.1 INTERPOLATION OR NOT.
    !     ------------------------

    IF (MR == 1 .and. MA == 1) THEN

      !***     1.11 NO INTERPOLATION.
      !        ----------------------

      CALL SECSPOM_LOKI(F1, F3, KIJS, KIJL, NFRE, NANG, NDEPTH, DEPTHA, DEPTHD, OMSTART, FRAC, MR, OMEGA, DEPTH, AKMEAN, TA, TB,  &
      & TC_QL, TT_4M, TT_4P, IM_P, IM_M, IJ=IJ, LL=SECSPOM_LL, OMEGA_EXT=SECSPOM_OMEGA_EXT, F2=SECSPOM_F2)
!$acc loop seq
      DO M=1,NFRE
!$acc loop seq
        DO K=1,NANG
          DELF = F3(IJ, K, M)
          F1(IJ, K, M) = MAX(MIN(0.000001_JWRB, F1(IJ, K, M)), F1(IJ, K, M) + SIG*DELF)
        END DO
      END DO

    ELSE


      !***     1.12 ENERGY CONSERVING INTERPOLATION SCHEME
      !        -------------------------------------------

!$acc loop seq
      DO M=1,NFREH
!$acc loop seq
        DO K=1,NANGH
          PF1(IJ, K, M) = 0._JWRB
        END DO
      END DO

!$acc loop seq
      DO M=1,NFREH
        M0 = MR*M
!$acc loop seq
        DO K=1,NANGH
          K0 = MA*K + 1
          IF (K0 > NANG) K0 = K0 - NANG
          IF (K0 < 1) K0 = K0 + NANG
          PF1(IJ, K, M) = PF1(IJ, K, M) + F1(IJ, K0, M0)
        END DO
      END DO

      !***     1.13 DETERMINE SECOND-ORDER SPEC
      !        --------------------------------

      CALL SECSPOM_LOKI(PF1, PF3, KIJS, KIJL, NFREH, NANGH, NDEPTH, DEPTHA, DEPTHD, OMSTART, FRAC, MR, OMEGA, DEPTH, AKMEAN, TA,  &
      & TB, TC_QL, TT_4M, TT_4P, IM_P, IM_M, IJ=IJ, LL=SECSPOM_LL, OMEGA_EXT=SECSPOM_OMEGA_EXT, F2=SECSPOM_F2)

      !***     2.24 INTERPOLATE TOWARDS HIGH-RES GRID
      !        --------------------------------------

      IF (EMEAN(IJ) <= ZFAC*DEPTH(IJ)**2) THEN
        EMAXL = 1._JWRB
      ELSE
        EMAXL = 0._JWRB
      END IF

!$acc loop seq
      DO M=1,NFRE
        XM = REAL(M / MR)
        !!!            D1 = REAL(M)/REAL(MR)-XM
        M0 = INT(XM)
        IF (M0 < 1) THEN
          M0 = 1
          MP = M0 + 1
          D1 = 1._JWRB
        ELSE IF (M0 < NFREH) THEN
          MP = M0 + 1
          D1 = (FR(M) - FR(MR*M0)) / (FR(MR*MP) - FR(MR*M0))
        ELSE
          M0 = NFREH
          MP = M0
          D1 = 0._JWRB
        END IF
        D2 = 1.0_JWRB - D1

!$acc loop seq
        DO K=1,NANG
          XK = REAL((K - 1) / MA)
          K0 = INT(XK)
          D3 = REAL(K - 1) / REAL(MA) - XK
          D4 = 1. - D3
          IF (K0 < 1) K0 = K0 + NANGH
          KP = K0 + 1
          IF (KP > NANGH) KP = KP - NANGH

          C1 = PF3(IJ, K0, M0)*D4 + PF3(IJ, KP, M0)*D3
          C2 = PF3(IJ, KP, MP)*D3 + PF3(IJ, K0, MP)*D4
          DELF = C1*D2 + C2*D1

          F1(IJ, K, M) = MAX(MIN(0.000001_JWRB, F1(IJ, K, M)), F1(IJ, K, M) + EMAXL*SIG*DELF)

        END DO
      END DO
    END IF



!$acc end data
  END SUBROUTINE CAL_SECOND_ORDER_SPEC_LOKI
END MODULE CAL_SECOND_ORDER_SPEC_LOKI_MOD
