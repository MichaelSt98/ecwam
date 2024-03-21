! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE CIMSSTRN_CUF_HOIST_NEW_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE CIMSSTRN_CUF_HOIST_NEW (KIJS, KIJL, FL1, WAVNUM, DEPTH, CITHICK, STRN, DELTH, DFIM, FLMIN, G,  &
  & NANG, NFRE, ROWATER, ICHNK, NCHNK, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *CIMSSTRN* - COMPUTATION OF THE MEAN SQUARE WAVE STRAIN IN SEA ICE.
    
    !     J. BIDLOT  ECMWF  JANUARY 2013.
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTES MEAN SQUARE WAVE STRAIN AT EACH GRID POINT.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *CIMSSTRN (KIJS, KIJL, FL1, WAVNUM, DEPTH, CITHICK, STRN)*
    !              *KIJS*    - INDEX OF FIRST GRIDPOINT
    !              *KIJL*    - INDEX OF LAST GRIDPOINT
    !              *FL1*     - SPECTRUM.
    !              *WAVNUM*  - OPEN WATER WAVE NUMBER
    !              *DEPTH*   - WATER DEPTH
    !              *CITHICK* - SEA ICE THICKNESS
    !              *STRN*    - MEAN SQUARE WAVE STRAIN IN ICE (OUTPUT).
    
    !     METHOD.
    !     -------
    
    !      !!! IT ASSUMES SO DEFAULT SETTING FOR THE MECHANICAL PROPERTIES OF
    !          THE SEA ICE (SEE AKI_ICE) !!!!!!!
    
    !       NONE.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       NONE.
    
    ! ----------------------------------------------------------------------
    
    USE AKI_ICE_CUF_HOIST_NEW_MOD, ONLY: AKI_ICE_CUF_HOIST_NEW
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    USE YOWFRED, ONLY: FR
    USE YOWPCONS, ONLY: ZPI
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
    INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    REAL(KIND=JWRB), INTENT(IN) :: FL1(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: WAVNUM(KIJL, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: DEPTH(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: CITHICK(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(OUT) :: STRN(KIJL, NCHNK)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: M
    INTEGER(KIND=JWIM) :: K
    REAL(KIND=JWRB) :: F1LIM
    REAL(KIND=JWRB) :: XKI
    REAL(KIND=JWRB) :: E
    REAL(KIND=JWRB) :: SUME
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: DFIM(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: FLMIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ROWATER
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    
    ! ----------------------------------------------------------------------
    
    
    !*    1. INITIALISE
    !        ----------
    
    F1LIM = FLMIN / DELTH
    
    
    STRN(IJ, ICHNK) = 0.0_JWRB
    
    ! ----------------------------------------------------------------------
    
    !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
    !        ------------------------------------------
    
    DO M=1,NFRE
      XKI = AKI_ICE_CUF_HOIST_NEW(G, WAVNUM(IJ, M, ICHNK), DEPTH(IJ, ICHNK), ROWATER, CITHICK(IJ, ICHNK))
      E = 0.5_JWRB*CITHICK(IJ, ICHNK)*XKI**3 / WAVNUM(IJ, M, ICHNK)
      
      SUME = 0.0_JWRB
      DO K=1,NANG
        SUME = SUME + FL1(IJ, K, M, ICHNK)
      END DO
      
      IF (SUME > F1LIM) THEN
        STRN(IJ, ICHNK) = STRN(IJ, ICHNK) + E**2*SUME*DFIM(M)
      END IF
      
    END DO
    
    
    
  END SUBROUTINE CIMSSTRN_CUF_HOIST_NEW
END MODULE CIMSSTRN_CUF_HOIST_NEW_MOD
