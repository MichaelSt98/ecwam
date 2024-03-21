! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE FEMEANWS_CUF_HOIST_NEW_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE FEMEANWS_CUF_HOIST_NEW (KIJS, KIJL, FL1, XLLWS, FM, EM, DELTH, DFIM, DFIMOFR, EPSMIN, FR,  &
  & FRTAIL, NANG, NFRE, WETAIL, ICHNK, NCHNK, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *FEMEANWS* - COMPUTATION OF MEAN ENERGY, MEAN FREQUENCY
    !                  FOR WINDSEA PART OF THE SPECTRUM AS DETERMINED
    !                  BY XLLWS
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTE MEAN FREQUENCY AT EACH GRID POINT FOR PART OF THE
    !       SPECTRUM WHERE XLLWS IS NON ZERO.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *FEMEANWS (KIJS, KIJL, FL1, XLLWS, EM, FM)*
    !              *KIJS*   - INDEX OF FIRST GRIDPOINT
    !              *KIJL*   - INDEX OF LAST GRIDPOINT
    !              *FL1*    - SPECTRUM.
    !              *XLLWS* - TOTAL WINDSEA MASK FROM INPUT SOURCE TERM
    !              *EM*     - MEAN WAVE ENERGY (OUTPUT)
    !              *FM*     - MEAN WAVE FREQUENCY (OUTPUT)
    
    !     METHOD.
    !     -------
    
    !       NONE.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       NONE.
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
    INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    REAL(KIND=JWRB), INTENT(IN) :: FL1(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: XLLWS(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(OUT) :: FM(KIJL)
    REAL(KIND=JWRB), OPTIONAL, INTENT(OUT) :: EM(KIJL)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: M
    INTEGER(KIND=JWIM) :: K
    
    REAL(KIND=JWRB) :: DELT25
    REAL(KIND=JWRB) :: DELT2
    REAL(KIND=JWRB) :: TEMP2
    REAL(KIND=JWRB) :: EM_LOC
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: DFIM(NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: DFIMOFR(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPSMIN
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: FR(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: FRTAIL
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WETAIL
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    
    ! ----------------------------------------------------------------------
    
    
    !*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
    !        ------------------------------------------------
    
    
    EM_LOC = EPSMIN
    FM(IJ) = EPSMIN
    
    DELT25 = WETAIL*FR(NFRE)*DELTH
    DELT2 = FRTAIL*DELTH
    
    
    !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
    !        ------------------------------------------
    
    DO M=1,NFRE
      TEMP2 = 0.0_JWRB
      DO K=1,NANG
        TEMP2 = TEMP2 + XLLWS(IJ, K, M, ICHNK)*FL1(IJ, K, M, ICHNK)
      END DO
      EM_LOC = EM_LOC + DFIM(M)*TEMP2
      FM(IJ) = FM(IJ) + DFIMOFR(M)*TEMP2
    END DO
    
    !*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
    !*       NORMALIZE WITH TOTAL ENERGY.
    !        ------------------------------------------
    
    EM_LOC = EM_LOC + DELT25*TEMP2
    FM(IJ) = FM(IJ) + DELT2*TEMP2
    FM(IJ) = EM_LOC / FM(IJ)
    
    IF (PRESENT(EM)) THEN
      EM(IJ) = EM_LOC
    END IF
    
    
    
  END SUBROUTINE FEMEANWS_CUF_HOIST_NEW
END MODULE FEMEANWS_CUF_HOIST_NEW_MOD
