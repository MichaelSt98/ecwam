! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
ATTRIBUTES(DEVICE) SUBROUTINE IMPHFTAIL_FC (KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1, NANG, NFRE, ICHNK, NCHNK, IJ)
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
  USE PARKIND_WAVE, ONLY: JWRU, JWIM, JWRB
  
  
  ! ----------------------------------------------------------------------
  
  IMPLICIT NONE
  
  INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
  INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
  INTEGER(KIND=JWIM), TARGET, INTENT(IN) :: MIJ(:, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: FLM(:, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: WAVNUM(:, :, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: XK2CG(:, :, :)
  REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: FL1(:, :, :, :)
  
  
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
  
  
END SUBROUTINE IMPHFTAIL_FC
