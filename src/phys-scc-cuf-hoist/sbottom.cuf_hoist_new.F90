! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SBOTTOM_CUF_HOIST_NEW_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE SBOTTOM_CUF_HOIST_NEW (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, BATHYMAX, GM1, NANG, NFRE_RED,  &
  & ICHNK, NCHNK, IJ)
    
    !SHALLOW
    ! ----------------------------------------------------------------------
    
    !**** *SBOTTOM* - COMPUTATION OF BOTTOM FRICTION.
    
    !     G.J.KOMEN AND Q.D.GAO
    !     OPTIMIZED BY L.F. ZAMBRESKY
    !     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTATION OF BOTTOM FRICTION DISSIPATION
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *SBOTTOM (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH)
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT
    !          *KIJL*    - INDEX OF LAST GRIDPOINT
    !          *FL1*     - SPECTRUM.
    !          *FLD*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*      - TOTAL SOURCE FUNCTION ARRAY
    !          *WAVNUM*  - WAVE NUMBER
    !          *DEPTH*   - WATER DEPTH
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCES.
    
    !     REFERENCES.
    !     -----------
    
    !       HASSELMANN ET AL, D. HYDR. Z SUPPL A12(1973) (JONSWAP)
    !       BOUWS AND KOMEN, JPO 13(1983)1653-1658
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    USE YOWPARAM, ONLY: NFRE
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
    INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    REAL(KIND=JWRB), INTENT(IN) :: FL1(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(INOUT) :: FLD(KIJL, NANG_loki_param, NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(INOUT) :: SL(KIJL, NANG_loki_param, NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(IN) :: WAVNUM(KIJL, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: DEPTH(KIJL, NCHNK)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: K
    INTEGER(KIND=JWIM) :: M
    REAL(KIND=JWRB) :: CONST
    REAL(KIND=JWRB) :: ARG
    REAL(KIND=JWRB) :: SBO
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: BATHYMAX
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: GM1
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE_RED
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    
    ! ----------------------------------------------------------------------
    
    
    CONST = -2.0_JWRB*0.038_JWRB*GM1
    
    DO M=1,NFRE_RED
      IF (DEPTH(IJ, ICHNK) < BATHYMAX) THEN
        ARG = 2.0_JWRB*DEPTH(IJ, ICHNK)*WAVNUM(IJ, M, ICHNK)
        ARG = MIN(ARG, 50.0_JWRB)
        SBO = CONST*WAVNUM(IJ, M, ICHNK) / SINH(ARG)
      ELSE
        SBO = 0.0_JWRB
      END IF
      
      DO K=1,NANG
        SL(IJ, K, M) = SL(IJ, K, M) + SBO*FL1(IJ, K, M, ICHNK)
        FLD(IJ, K, M) = FLD(IJ, K, M) + SBO
      END DO
    END DO
    
    
    
  END SUBROUTINE SBOTTOM_CUF_HOIST_NEW
END MODULE SBOTTOM_CUF_HOIST_NEW_MOD
