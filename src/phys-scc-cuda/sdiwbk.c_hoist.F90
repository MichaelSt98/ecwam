! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
ATTRIBUTES(DEVICE) SUBROUTINE SDIWBK_FC (KIJS, KIJL, FL1, FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN, LBIWBK, NANG, NFRE_RED, ICHNK,  &
& NCHNK, IJ)
  
  ! ----------------------------------------------------------------------
  
  !**** *SDIWBK* - COMPUTATION OF BOTTOM-INDUCED WAVE BREAKING DISSIPATION
  
  
  !*    PURPOSE.
  !     --------
  !       COMPUTE BOTTOM-INDUCED DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
  !       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
  !       OF DISSIPATION SOURCE FUNCTION.
  
  !**   INTERFACE.
  !     ----------
  
  !       *CALL* *SDIWBK (KIJS, KIJL, FL1, FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN)*
  !          *KIJS*    - INDEX OF FIRST GRIDPOINT
  !          *KIJL*    - INDEX OF LAST GRIDPOINT
  !          *FL1*     - SPECTRUM.
  !          *FLD*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
  !          *SL*      - TOTAL SOURCE FUNCTION ARRAY
  !          *DEPTH*   - WATER DEPTH
  !          *EMAXDPT* - MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
  !          *EMEAN*   - MEAN ENERGY DENSITY
  !          *F1MEAN*  - MEAN FREQUENCY BASED ON 1st MOMENT.
  
  !     METHOD.
  !     -------
  
  !       SEE REFERENCES.
  
  !     EXTERNALS.
  !     ----------
  
  !       NONE.
  
  !     REFERENCE.
  !     ----------
  
  ! ----------------------------------------------------------------------
  
  USE PARKIND_WAVE, ONLY: JWRU, JWIM, JWRB
  
  USE YOWPARAM, ONLY: NFRE
  
  
  ! ----------------------------------------------------------------------
  
  IMPLICIT NONE
  
  INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
  INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: FL1(:, :, :, :)
  REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: FLD(:, :, :)
  REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: SL(:, :, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: DEPTH(:, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: EMAXDPT(:, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: EMEAN(:)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: F1MEAN(:)
  
  
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
  INTEGER(KIND=JWIM) :: K
  INTEGER(KIND=JWIM) :: M
  INTEGER(KIND=JWIM) :: IC
  REAL(KIND=JWRB) :: ALPH
  REAL(KIND=JWRB) :: ARG
  REAL(KIND=JWRB) :: Q
  REAL(KIND=JWRB) :: Q_OLD
  REAL(KIND=JWRB) :: REL_ERR
  REAL(KIND=JWRB) :: EXPQ
  REAL(KIND=JWRB) :: SDS
  
  REAL, PARAMETER :: ALPH_B_J = 1.0_JWRB
  REAL, PARAMETER :: COEF_B_J = 2*ALPH_B_J
  REAL, PARAMETER :: DEPTHTRS = 50.0_JWRB
  LOGICAL, VALUE, INTENT(IN) :: LBIWBK
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE_RED
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
  INTEGER, VALUE, INTENT(IN) :: NCHNK
  
  ! ----------------------------------------------------------------------
  
  
  !*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
  !*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
  !        --------------------------------------------------------------
  
  
  IF (LBIWBK) THEN
    !       (FOLLOWING BATTJES-JANSSEN AND BEJI)
    IF (DEPTH(IJ, ICHNK) < DEPTHTRS) THEN
      ALPH = 2.0_JWRB*EMAXDPT(IJ, ICHNK) / EMEAN(IJ)
      ARG = MIN(ALPH, 50.0_JWRB)
      Q_OLD = EXP(-ARG)
      !            USE NEWTON-RAPHSON METHOD
      DO IC=1,15
        EXPQ = EXP(-ARG*(1.0_JWRB - Q_OLD))
        Q = Q_OLD - (EXPQ - Q_OLD) / (ARG*EXPQ - 1.0_JWRB)
        REL_ERR = ABS(Q - Q_OLD) / Q_OLD
        IF (REL_ERR < 0.00001_JWRB) EXIT
        Q_OLD = Q
      END DO
      Q = MIN(Q, 1.0_JWRB)
      SDS = COEF_B_J*ALPH*Q*F1MEAN(IJ)
    END IF
    
    DO M=1,NFRE_RED
      DO K=1,NANG
        IF (DEPTH(IJ, ICHNK) < DEPTHTRS) THEN
          SL(IJ, K, M) = SL(IJ, K, M) - SDS*FL1(IJ, K, M, ICHNK)
          FLD(IJ, K, M) = FLD(IJ, K, M) - SDS
        END IF
      END DO
    END DO
    
  END IF
  
  
  
END SUBROUTINE SDIWBK_FC
