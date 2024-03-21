! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDISSIP_JAN_CUF_HOIST_NEW_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE SDISSIP_JAN_CUF_HOIST_NEW (KIJS, KIJL, FL1, FLD, SL, WAVNUM, EMEAN, F1MEAN, XKMEAN, CDIS,  &
  & CDISVIS, DELTA_SDIS, NANG, NFRE, RNU, ZPI, ICHNK, NCHNK, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *SDISSIP_JAN* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.
    
    !     S.D.HASSELMANN.
    !     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
    !     OPTIMIZATION : L. ZAMBRESKY
    !     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
    !     J. BIDLOT   ECMWF  NOVEMBER 2004  REFORMULATION BASED ON XKMEAN
    !                                       AND F1MEAN.
    !                        AUGUST 2020 Added small viscous dissipation term
    
    !*    PURPOSE.
    !     --------
    !       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
    !       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
    !       OF DISSIPATION SOURCE FUNCTION.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *SDISSIP_JAN (KIJS, KIJ, FL1, FLD, SL,
    !                            WAVNUM,
    !                            EMEAN,F1MEAN, XKMEAN,)*
    !          *FL1*    - SPECTRUM.
    !          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*     - TOTAL SOURCE FUNCTION ARRAY
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *WAVNUM* - WAVE NUMBER
    !          *EMEAN*  - MEAN ENERGY DENSITY
    !          *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
    !          *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.
    
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCES.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       G.KOMEN, S. HASSELMANN AND K. HASSELMANN, ON THE EXISTENCE
    !          OF A FULLY DEVELOPED WINDSEA SPECTRUM, JGR, 1984.
    
    ! ---------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    USE YOWFRED, ONLY: FR, DFIM, DELTH, FRATIO
    USE YOWPCONS, ONLY: ZPI4GM2, G
    
    
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
    REAL(KIND=JWRB), INTENT(IN) :: EMEAN(KIJL)
    REAL(KIND=JWRB), INTENT(IN) :: F1MEAN(KIJL)
    REAL(KIND=JWRB), INTENT(IN) :: XKMEAN(KIJL)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: K
    INTEGER(KIND=JWIM) :: M
    
    REAL(KIND=JWRB) :: SCDFM
    REAL(KIND=JWRB) :: CONSD
    REAL(KIND=JWRB) :: CONSS
    REAL(KIND=JWRB) :: DELTA_SDISM1
    REAL(KIND=JWRB) :: CVIS
    REAL(KIND=JWRB) :: TEMP1
    REAL(KIND=JWRB) :: SDS
    REAL(KIND=JWRB) :: X
    REAL(KIND=JWRB) :: XK2
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: CDIS
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: CDISVIS
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTA_SDIS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: RNU
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZPI
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    
    ! ----------------------------------------------------------------------
    
    
    !*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
    !*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
    !        --------------------------------------------------------------
    
    DELTA_SDISM1 = 1.0_JWRB - DELTA_SDIS
    
    CONSS = CDIS*ZPI
    
    SDS = CONSS*F1MEAN(IJ)*EMEAN(IJ)**2*XKMEAN(IJ)**4
    
    DO M=1,NFRE
      X = WAVNUM(IJ, M, ICHNK) / XKMEAN(IJ)
      XK2 = WAVNUM(IJ, M, ICHNK)**2
      
      CVIS = RNU*CDISVIS
      TEMP1 = SDS*X*(DELTA_SDISM1 + DELTA_SDIS*X) + CVIS*XK2
      
      DO K=1,NANG
        FLD(IJ, K, M) = FLD(IJ, K, M) + TEMP1
        SL(IJ, K, M) = SL(IJ, K, M) + TEMP1*FL1(IJ, K, M, ICHNK)
      END DO
      
    END DO
    
    
    
  END SUBROUTINE SDISSIP_JAN_CUF_HOIST_NEW
END MODULE SDISSIP_JAN_CUF_HOIST_NEW_MOD
