! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE HALPHAP_LOKI_MOD
  CONTAINS
  SUBROUTINE HALPHAP_LOKI (KIJS, KIJL, WAVNUM, COSWDIF, FL1, HALP, IJ, WD, FLWD, FD)

    ! ----------------------------------------------------------------------

    !**** *HALPHAP* - COMPUTATION OF 1/2 PHILLIPS PARAMETER


    !**   INTERFACE.
    !     ----------

    !       *CALL* *HALPHAP(KIJS, KIJL, WAVNUM, UDIR, FL1, HALP)
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *WAVNUM* - WAVE NUMBER
    !          *COSWDIF*- COSINE ( WIND SPEED DIRECTION - WAVE DIRECTIONS)
    !          *FL1*    - SPECTRA
    !          *HALP*   - 1/2 PHILLIPS PARAMETER

    !     METHOD.
    !     -------

    ! ----------------------------------------------------------------------

    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWFRED, ONLY: FR, TH, FR5, DELTH, ZPIFR, DFIM, DFIMOFR, WETAIL, FRTAIL, COSTH, SINTH, FRATIO, FLOGSPRDM1, RHOWG_DFIM,  &
    & FRIC, OMEGA_GC, XK_GC
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: G, ZPI, ZPI4GM2, GM1, EPSMIN, EPSUS, ACD, BCD, ACDLIN, BCDLIN, CDMAX, ZPI4GM1
    USE YOWPHYS, ONLY: ALPHAPMAX

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: HALP

    INTEGER(KIND=JWIM) :: HALPHAP_K, M
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ

    REAL(KIND=JWRB) :: ZLNFRNFRE
    REAL(KIND=JPHOOK) :: HALPHAP_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: ALPHAP
    REAL(KIND=JWRB) :: XMSS, EM, FM, F1D
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG) :: WD
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FLWD
    INTEGER(KIND=JWIM) :: MEANSQS_LF_M, MEANSQS_LF_K, KFRE
    REAL(KIND=JPHOOK) :: MEANSQS_LF_ZHOOK_HANDLE
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NFRE) :: FD
    REAL(KIND=JWRB) :: TEMP1, TEMP2
    INTEGER(KIND=JWIM) :: FEMEAN_M, FEMEAN_K
    REAL(KIND=JWRB) :: DELT25, DELT2, DEL2
    REAL(KIND=JPHOOK) :: FEMEAN_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: FEMEAN_TEMP2
!$acc routine seq
!$acc data present( WAVNUM, COSWDIF, FL1, HALP, WD, FLWD, FD )

    ! ----------------------------------------------------------------------


    ZLNFRNFRE = LOG(FR(NFRE))

    ! Find spectrum in wind direction
!$acc loop seq
    DO HALPHAP_K=1,NANG
      WD(IJ, HALPHAP_K) = 0.5_JWRB + 0.5_JWRB*SIGN(1.0_JWRB, COSWDIF(IJ, HALPHAP_K))
    END DO

!$acc loop seq
    DO M=1,NFRE
!$acc loop seq
      DO HALPHAP_K=1,NANG
        FLWD(IJ, HALPHAP_K, M) = FL1(IJ, HALPHAP_K, M)*WD(IJ, HALPHAP_K)
      END DO
    END DO

    ! [Loki] inlined child subroutine: MEANSQS_LF
    ! =========================================

    ! ----------------------------------------------------------------------

    !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
    !        ------------------------------------------

    KFRE = MIN(NFRE, NFRE)

    XMSS = 0.0_JWRB

    !*    2.2 SHALLOW WATER INTEGRATION.
    !         --------------------------

!$acc loop seq
    DO MEANSQS_LF_M=1,KFRE
      TEMP1 = DFIM(MEANSQS_LF_M)*WAVNUM(IJ, MEANSQS_LF_M)**2
      TEMP2 = 0.0_JWRB
!$acc loop seq
      DO MEANSQS_LF_K=1,NANG
        TEMP2 = TEMP2 + FLWD(IJ, MEANSQS_LF_K, MEANSQS_LF_M)
      END DO
      XMSS = XMSS + TEMP1*TEMP2
    END DO


    ! =========================================

    ! [Loki] inlined child subroutine: FEMEAN
    ! =========================================

    ! ----------------------------------------------------------------------


    !*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
    !        ------------------------------------------------

    EM = 0.0_JWRB
    FM = 0.0_JWRB

    DELT25 = WETAIL*FR(NFRE)*DELTH
    DELT2 = FRTAIL*DELTH

    !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
    !        ------------------------------------------

!$acc loop seq
    DO FEMEAN_M=1,NFRE
      FEMEAN_K = 1
      FEMEAN_TEMP2 = MAX(FLWD(IJ, FEMEAN_K, FEMEAN_M), EPSMIN)
!$acc loop seq
      DO FEMEAN_K=2,NANG
        FEMEAN_TEMP2 = FEMEAN_TEMP2 + MAX(FLWD(IJ, FEMEAN_K, FEMEAN_M), EPSMIN)
      END DO
      EM = EM + FEMEAN_TEMP2*DFIM(FEMEAN_M)
      FM = FM + DFIMOFR(FEMEAN_M)*FEMEAN_TEMP2
    END DO


    !*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
    !*       NORMALIZE WITH TOTAL ENERGY.
    !        ------------------------------------------

    EM = EM + DELT25*FEMEAN_TEMP2
    FM = FM + DELT2*FEMEAN_TEMP2
    FM = EM / FM
    FM = MAX(FM, FR(1))


    ! =========================================

    IF (EM > 0.0_JWRB .and. FM < FR(-2 + NFRE)) THEN
      ALPHAP = XMSS / (ZLNFRNFRE - LOG(FM))
      IF (ALPHAP > ALPHAPMAX) THEN
        ! some odd cases, revert to tail value
        F1D = 0.0_JWRB
!$acc loop seq
        DO HALPHAP_K=1,NANG
          F1D = F1D + FLWD(IJ, HALPHAP_K, NFRE)*DELTH
        END DO
        ALPHAP = ZPI4GM2*FR5(NFRE)*F1D
      END IF
    ELSE
      F1D = 0.0_JWRB
!$acc loop seq
      DO HALPHAP_K=1,NANG
        F1D = F1D + FLWD(IJ, HALPHAP_K, NFRE)*DELTH
      END DO
      ALPHAP = ZPI4GM2*FR5(NFRE)*F1D
    END IF

    !     1/2 ALPHAP:
    HALP(IJ) = 0.5_JWRB*MIN(ALPHAP, ALPHAPMAX)



!$acc end data
  END SUBROUTINE HALPHAP_LOKI
END MODULE HALPHAP_LOKI_MOD
