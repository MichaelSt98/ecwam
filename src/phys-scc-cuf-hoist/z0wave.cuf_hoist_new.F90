! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE Z0WAVE_CUF_HOIST_NEW_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE Z0WAVE_CUF_HOIST_NEW (KIJS, KIJL, US, TAUW, UTOP, Z0, Z0B, CHRNCK, ALPHA, ALPHAMIN, CHNKMIN_U,  &
  & EPS1, G, GM1, LLCAPCHNK, ICHNK, NCHNK, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *Z0WAVE* - DETERMINE THE SEA STATE DEPENDENT ROUGHNESS LENGTH.
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTE ROUGHNESS LENGTH.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *Z0WAVE (KIJS, KIJL, US, TAUW, UTOP, Z0, Z0B, CHRNCK)
    !          *KIJS* - INDEX OF FIRST GRIDPOINT.
    !          *KIJL* - INDEX OF LAST GRIDPOINT.
    !          *US*   - OUTPUT BLOCK OF SURFACE STRESSES.
    !          *TAUW* - INPUT BLOCK OF WAVE STRESSES.
    !          *UTOP* - WIND SPEED.
    !          *Z0*   - OUTPUT BLOCK OF ROUGHNESS LENGTH.
    !          *Z0B*  - BACKGROUND ROUGHNESS LENGTH.
    !          *CHRNCK- CHARNOCK COEFFICIENT
    
    !     METHOD.
    !     -------
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ---------
    
    !       NONE.
    
    ! ----------------------------------------------------------------------
    
    USE CHNKMIN_CUF_HOIST_NEW_MOD, ONLY: CHNKMIN_CUF_HOIST_NEW
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    REAL(KIND=JWRB), INTENT(IN) :: US(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: TAUW(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: UTOP(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(OUT) :: Z0(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(OUT) :: Z0B(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(OUT) :: CHRNCK(KIJL, NCHNK)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    REAL(KIND=JWRB) :: UST2
    REAL(KIND=JWRB) :: UST3
    REAL(KIND=JWRB) :: ARG
    REAL(KIND=JWRB) :: ALPHAOG
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ALPHA
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ALPHAMIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: CHNKMIN_U
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPS1
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: GM1
    LOGICAL, VALUE, INTENT(IN) :: LLCAPCHNK
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    
    ! ----------------------------------------------------------------------
    
    
    
    IF (LLCAPCHNK) THEN
      ALPHAOG = CHNKMIN_CUF_HOIST_NEW(UTOP(IJ, ICHNK), ALPHA, ALPHAMIN, CHNKMIN_U)*GM1
    ELSE
      ALPHAOG = ALPHA*GM1
    END IF
    
    UST2 = US(IJ, ICHNK)**2
    UST3 = US(IJ, ICHNK)**3
    ARG = MAX(UST2 - TAUW(IJ, ICHNK), EPS1)
    Z0(IJ, ICHNK) = ALPHAOG*UST3 / SQRT(ARG)
    Z0B(IJ, ICHNK) = ALPHAOG*UST2
    CHRNCK(IJ, ICHNK) = G*Z0(IJ, ICHNK) / UST2
    
    
    
  END SUBROUTINE Z0WAVE_CUF_HOIST_NEW
END MODULE Z0WAVE_CUF_HOIST_NEW_MOD
