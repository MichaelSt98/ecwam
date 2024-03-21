! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE PEAK_ANG_CUF_HOIST_NEW_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE PEAK_ANG_CUF_HOIST_NEW (KIJS, KIJL, FL1, XNU, SIG_TH, COSTH, DELTH, DFIM, DFIMFR, DFIMFR2, FR,  &
  & FRATIO, NANG, NFRE, SINTH, TH, WETAIL, WP1TAIL, WP2TAIL, ICHNK, NCHNK, IJ)
    
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
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    USE YOWFRED, ONLY: DFIMOFR
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
    INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    REAL(KIND=JWRB), INTENT(IN) :: FL1(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(OUT) :: XNU(KIJL)
    REAL(KIND=JWRB), INTENT(OUT) :: SIG_TH(KIJL)
    
    
    INTEGER(KIND=JWIM) :: NSH
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: M
    INTEGER(KIND=JWIM) :: K
    INTEGER(KIND=JWIM) :: MMAX
    INTEGER(KIND=JWIM) :: MMSTART
    INTEGER(KIND=JWIM) :: MMSTOP
    REAL(KIND=JWRB), PARAMETER :: CONST_SIG = 1.0_JWRB
    REAL(KIND=JWRB) :: R1
    REAL(KIND=JWRB) :: DELT25
    REAL(KIND=JWRB) :: COEF_FR
    REAL(KIND=JWRB) :: COEF_FR2
    REAL(KIND=JWRB) :: ZEPSILON
    REAL(KIND=JWRB) :: SUM0
    REAL(KIND=JWRB) :: SUM1
    REAL(KIND=JWRB) :: SUM2
    REAL(KIND=JWRB) :: XMAX
    REAL(KIND=JWRB) :: TEMP
    REAL(KIND=JWRB) :: THMEAN
    REAL(KIND=JWRB) :: SUM_S
    REAL(KIND=JWRB) :: SUM_C
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: COSTH(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: DFIM(NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: DFIMFR(NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: DFIMFR2(NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: FR(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: FRATIO
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: SINTH(NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: TH(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WETAIL
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WP1TAIL
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WP2TAIL
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    
    ! ----------------------------------------------------------------------
    
    !***  1. DETERMINE L-H SPECTRAL WIDTH OF THE 2-D SPECTRUM.
    !     ---------------------------------------------------
    
    ZEPSILON = 10._JWRB*EPSILON(ZEPSILON)
    NSH = 1 + INT(LOG(1.5_JWRB) / LOG(FRATIO))
    
    
    SUM0 = ZEPSILON
    SUM1 = 0._JWRB
    SUM2 = 0._JWRB
    
    DO M=1,NFRE
      K = 1
      TEMP = FL1(IJ, K, M, ICHNK)
      DO K=2,NANG
        TEMP = TEMP + FL1(IJ, K, M, ICHNK)
      END DO
      SUM0 = SUM0 + TEMP*DFIM(M)
      SUM1 = SUM1 + TEMP*DFIMFR(M)
      SUM2 = SUM2 + TEMP*DFIMFR2(M)
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
    
    DO M=2,NFRE - 1
      DO K=1,NANG
        IF (FL1(IJ, K, M, ICHNK) > XMAX) THEN
          MMAX = M
          XMAX = FL1(IJ, K, M, ICHNK)
        END IF
      END DO
    END DO
    
    SUM1 = ZEPSILON
    SUM2 = 0._JWRB
    
    MMSTART = MAX(1, MMAX - NSH)
    MMSTOP = MIN(NFRE, MMAX + NSH)
    DO M=MMSTART,MMSTOP
      SUM_S = 0._JWRB
      SUM_C = ZEPSILON
      DO K=1,NANG
        SUM_S = SUM_S + SINTH(K)*FL1(IJ, K, M, ICHNK)
        SUM_C = SUM_C + COSTH(K)*FL1(IJ, K, M, ICHNK)
      END DO
      THMEAN = ATAN2(SUM_S, SUM_C)
      DO K=1,NANG
        SUM1 = SUM1 + FL1(IJ, K, M, ICHNK)*DFIM(M)
        SUM2 = SUM2 + COS(TH(K) - THMEAN)*FL1(IJ, K, M, ICHNK)*DFIM(M)
      END DO
    END DO
    
    IF (SUM1 > ZEPSILON) THEN
      R1 = SUM2 / SUM1
      SIG_TH(IJ) = CONST_SIG*SQRT(2._JWRB*(1._JWRB - R1))
    ELSE
      SIG_TH(IJ) = 0._JWRB
    END IF
    
    
    
  END SUBROUTINE PEAK_ANG_CUF_HOIST_NEW
END MODULE PEAK_ANG_CUF_HOIST_NEW_MOD
