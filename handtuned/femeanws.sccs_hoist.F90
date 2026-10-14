! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE FEMEANWS_LOKI_MOD
  CONTAINS
  SUBROUTINE FEMEANWS_LOKI (KIJS, KIJL, FL1, XLLWS, FM, EM, IJ)

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

    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWFRED, ONLY: FR, DFIM, DFIMOFR, DELTH, WETAIL, FRTAIL
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: EPSMIN

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1, XLLWS
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: FM
    REAL(KIND=JWRB), OPTIONAL, INTENT(OUT), DIMENSION(KIJL) :: EM


    INTEGER(KIND=JWIM) :: FEMEANWS_FEMEANWS_M, FEMEANWS_FEMEANWS_K
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ

    REAL(KIND=JWRB) :: FEMEANWS_FEMEANWS_DELT25, FEMEANWS_FEMEANWS_DELT2
    REAL(KIND=JPHOOK) :: FEMEANWS_FEMEANWS_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: FEMEANWS_FEMEANWS_TEMP2, FEMEANWS_EM_LOC
!$acc routine seq
!$acc data present( FL1, XLLWS, FM, EM )

    ! ----------------------------------------------------------------------


    !*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
    !        ------------------------------------------------

    FEMEANWS_EM_LOC = EPSMIN
    FM(IJ) = EPSMIN

    FEMEANWS_FEMEANWS_DELT25 = WETAIL*FR(NFRE)*DELTH
    FEMEANWS_FEMEANWS_DELT2 = FRTAIL*DELTH


    !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
    !        ------------------------------------------

!$acc loop seq
    DO FEMEANWS_FEMEANWS_M=1,NFRE
      FEMEANWS_FEMEANWS_TEMP2 = 0.0_JWRB
!$acc loop seq
      DO FEMEANWS_FEMEANWS_K=1,NANG
        FEMEANWS_FEMEANWS_TEMP2 = FEMEANWS_FEMEANWS_TEMP2 + XLLWS(IJ, FEMEANWS_FEMEANWS_K, FEMEANWS_FEMEANWS_M)*FL1(IJ,  &
        & FEMEANWS_FEMEANWS_K, FEMEANWS_FEMEANWS_M)
      END DO
      FEMEANWS_EM_LOC = FEMEANWS_EM_LOC + DFIM(FEMEANWS_FEMEANWS_M)*FEMEANWS_FEMEANWS_TEMP2
      FM(IJ) = FM(IJ) + DFIMOFR(FEMEANWS_FEMEANWS_M)*FEMEANWS_FEMEANWS_TEMP2
    END DO

    !*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
    !*       NORMALIZE WITH TOTAL ENERGY.
    !        ------------------------------------------

    FEMEANWS_EM_LOC = FEMEANWS_EM_LOC + FEMEANWS_FEMEANWS_DELT25*FEMEANWS_FEMEANWS_TEMP2
    FM(IJ) = FM(IJ) + FEMEANWS_FEMEANWS_DELT2*FEMEANWS_FEMEANWS_TEMP2
    FM(IJ) = FEMEANWS_EM_LOC / FM(IJ)

    IF (PRESENT(EM)) THEN
      EM(IJ) = FEMEANWS_EM_LOC
    END IF



!$acc end data
  END SUBROUTINE FEMEANWS_LOKI
END MODULE FEMEANWS_LOKI_MOD
