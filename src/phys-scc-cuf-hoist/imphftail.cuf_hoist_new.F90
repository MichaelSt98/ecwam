! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE IMPHFTAIL_CUF_HOIST_NEW_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE IMPHFTAIL_CUF_HOIST_NEW (KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1, NANG, NFRE, ICHNK, NCHNK, IJ)
    ! ----------------------------------------------------------------------
    
    !**** *IMPHFTAIL* - IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM
    
    
    !*    PURPOSE.
    !     --------
    
    !     IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM ABOVE FREQUENCY INDEX MIJ
    
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *IMPHFTAIL (KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1)
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT
    !          *KIJL*    - INDEX OF LAST GRIDPOINT
    !          *MIJ*     - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
    !          *FLM*     - SPECTAL DENSITY MINIMUM VALUE
    !          *WAVNUM*  - WAVENUMBER
    !          *XK2CG*   - (WAVNUM)**2 * GROUP SPEED
    !          *FL1*     - SPECTRUM (INPUT AND OUTPUT).
    
    !     METHOD.
    !     -------
    
    !     EXTERNALS.
    !     ---------
    
    !     REFERENCE.
    !     ----------
    
    ! ----------------------------------------------------------------------
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
    INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    INTEGER(KIND=JWIM), INTENT(IN) :: MIJ(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: FLM(KIJL, NANG_loki_param)
    REAL(KIND=JWRB), INTENT(IN) :: WAVNUM(KIJL, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: XK2CG(KIJL, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(INOUT) :: FL1(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: K
    INTEGER(KIND=JWIM) :: M
    
    REAL(KIND=JWRB) :: AKM1
    REAL(KIND=JWRB) :: TFAC
    REAL(KIND=JWRB) :: TEMP1
    REAL(KIND=JWRB) :: TEMP2
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    
    ! ----------------------------------------------------------------------
    
    
    !*    DIAGNOSTIC TAIL.
    !     ----------------
    
    
    TEMP1 = 1.0_JWRB / XK2CG(IJ, MIJ(IJ, ICHNK), ICHNK) / WAVNUM(IJ, MIJ(IJ, ICHNK), ICHNK)
    
    DO M=MIJ(IJ, ICHNK) + 1,NFRE
      TEMP2 = 1.0_JWRB / XK2CG(IJ, M, ICHNK) / WAVNUM(IJ, M, ICHNK)
      TEMP2 = TEMP2 / TEMP1
      
      !*    MERGE TAIL INTO SPECTRA.
      !     ------------------------
      DO K=1,NANG
        TFAC = FL1(IJ, K, MIJ(IJ, ICHNK), ICHNK)
        FL1(IJ, K, M, ICHNK) = MAX(TEMP2*TFAC, FLM(IJ, K))
      END DO
    END DO
    
    
    ! ----------------------------------------------------------------------
    
    
  END SUBROUTINE IMPHFTAIL_CUF_HOIST_NEW
END MODULE IMPHFTAIL_CUF_HOIST_NEW_MOD
