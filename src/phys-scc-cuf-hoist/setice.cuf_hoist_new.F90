! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!-----------------------------------------------------------------------
MODULE SETICE_CUF_HOIST_NEW_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE SETICE_CUF_HOIST_NEW (KIJS, KIJL, FL1, CICOVER, COSWDIF, CITHRSH, EPSMIN, FLMIN, NANG, NFRE,  &
  & ICHNK, NCHNK, IJ)
    
    !-----------------------------------------------------------------------
    
    !**** *SETICE* ROUTINE TO SET SPECTRA ON ICE TO NOISE LEVEL.
    
    !     R.PORTZ      MPI         OKT.1992
    !     J. BIDLOT    ECMWF       JUNE 1996  MESSAGE PASSING
    
    !     PURPOSE.
    !     -------
    
    !          *SETICE* SET ICE SPECTRA (FL1) TO NOISE LEVEL
    
    !**   INTERFACE.
    !     ----------
    
    !         *CALL* *SETICE(KIJS, KIJL, FL1, CICOVER, WSWAVE, COSWDIF)*
    !          *KIJS*    - LOCAL INDEX OF FIRST GRIDPOINT
    !          *KIJL*    - LOCAL  INDEX OF LAST GRIDPOINT
    !          *FL1*     - SPECTRA
    !          *CICOVER* - SEA ICE COVER
    !          *WSWAVE*  - WIND SPEED.
    !          *COSWDIF* - COS(TH(K)-WDWAVE(IJ))
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
    INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    REAL(KIND=JWRB), INTENT(INOUT) :: FL1(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: CICOVER(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: COSWDIF(KIJL, NANG_loki_param)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: M
    INTEGER(KIND=JWIM) :: K
    
    REAL(KIND=JWRB) :: CIREDUC
    REAL(KIND=JWRB) :: TEMP
    REAL(KIND=JWRB) :: ICEFREE
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: CITHRSH
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPSMIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: FLMIN
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    ! ----------------------------------------------------------------------
    
    
    !*    1. SET SPECTRA TO NOISE LEVEL OVER ICE POINTS.
    !     ----------------------------------------------
    
    
    IF (CICOVER(IJ, ICHNK) > CITHRSH) THEN
      CIREDUC = MAX(EPSMIN, (1.0_JWRB - CICOVER(IJ, ICHNK)))
      ICEFREE = 0.0_JWRB
    ELSE
      CIREDUC = 0.0_JWRB
      ICEFREE = 1.0_JWRB
    END IF
    
    TEMP = CIREDUC*FLMIN
    DO M=1,NFRE
      DO K=1,NANG
        FL1(IJ, K, M, ICHNK) = FL1(IJ, K, M, ICHNK)*ICEFREE + TEMP*MAX(0.0_JWRB, COSWDIF(IJ, K))**2
      END DO
    END DO
    
    
    
  END SUBROUTINE SETICE_CUF_HOIST_NEW
END MODULE SETICE_CUF_HOIST_NEW_MOD
