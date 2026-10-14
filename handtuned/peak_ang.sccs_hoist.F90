! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE PEAK_ANG_LOKI_MOD
  CONTAINS
  SUBROUTINE PEAK_ANG_LOKI (KIJS, KIJL, FL1, XNU, SIG_TH, IJ)

    !***  *PEAK_ANG*   DETERMINES ANGULAR WIDTH NEAR PEAK OF SPECTRUM

    !     PETER JANSSEN

    !     PURPOSE.
    !     --------

    !              DETERMINATION OF PEAK PARAMETERS

    !     INTERFACE.
    !     ----------
    !              *CALL*  *PEAK_ANG(KIJS,KIJL,FL1,XNU,SIG_TH)*

    !               INPUT:
    !                  *KIJS*   - FIRST GRIDPOINT
    !                  *KIJL*   - LAST GRIDPOINT
    !                  *FL1*    - SPECTRUM
    !               OUTPUT:
    !                  *XNU*    - RELATIVE SPECTRAL WIDTH
    !                  *SIG_TH* - RELATIVE WIDTH IN DIRECTION

    !     METHOD.
    !     -------
    !              NONE

    !     EXTERNALS.
    !     ----------
    !              NONE

    !-----------------------------------------------------------------------

    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWFRED, ONLY: FR, DFIM, DFIMFR, DFIMOFR, DFIMFR2, DELTH, TH, SINTH, COSTH, WETAIL, WP1TAIL, WP2TAIL, FRATIO
    USE YOWPARAM, ONLY: NANG, NFRE

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: XNU, SIG_TH


    INTEGER(KIND=JWIM) :: NSH
    INTEGER(KIND=JWIM) :: PEAK_ANG_M, PEAK_ANG_K
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
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
!$acc data present( FL1, XNU, SIG_TH )

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
      XNU(IJ) = SQRT(MAX(ZEPSILON, SUM2*SUM0 / SUM1**2 - 1._JWRB))
    ELSE
      XNU(IJ) = ZEPSILON
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
      SIG_TH(IJ) = CONST_SIG*SQRT(2._JWRB*(1._JWRB - R1))
    ELSE
      SIG_TH(IJ) = 0._JWRB
    END IF



!$acc end data
  END SUBROUTINE PEAK_ANG_LOKI
END MODULE PEAK_ANG_LOKI_MOD
