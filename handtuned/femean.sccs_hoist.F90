! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE FEMEAN_LOKI_MOD
  CONTAINS
  SUBROUTINE FEMEAN_LOKI (KIJS, KIJL, F, EM, FM, IJ)

    ! ----------------------------------------------------------------------

    !**** *FEMEAN* - COMPUTATION OF VARIANCE and MEAN FREQUENCY AT EACH GRID POINT

    !     S.D. HASSELMANN
    !     MODIFIED : P.JANSSEN (INTEGRATION OF F**-4 TAIL)
    !     OPTIMIZED BY : L. ZAMBRESKY AND H. GUENTHER


    !*    PURPOSE.
    !     --------

    !       COMPUTE MEAN FREQUENCY AT EACH GRID POINT.

    !**   INTERFACE.
    !     ----------

    !       *CALL* *FEMEAN (KIJS, KIJL, F, EM, FM)*
    !              *KIJS*    - INDEX OF FIRST GRIDPOINT
    !              *KIJL*    - INDEX OF LAST GRIDPOINT
    !              *F*       - SPECTRUM.
    !              *EM*      - MEAN WAVE ENERGY (OUTPUT)
    !              *FM*      - MEAN WAVE FREQUENCY (OUTPUT)

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
    USE YOWPCONS, ONLY: G, ZPI, EPSMIN

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE

    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: F

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: EM, FM


    INTEGER(KIND=JWIM) :: FEMEAN_M, FEMEAN_K
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    REAL(KIND=JWRB) :: DELT25, DELT2, DEL2
    REAL(KIND=JPHOOK) :: FEMEAN_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: FEMEAN_TEMP2
!$acc routine seq
!$acc data present( F, EM, FM )

    ! ----------------------------------------------------------------------


    !*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
    !        ------------------------------------------------

    EM(IJ) = 0.0_JWRB
    FM(IJ) = 0.0_JWRB

    DELT25 = WETAIL*FR(NFRE)*DELTH
    DELT2 = FRTAIL*DELTH

    !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
    !        ------------------------------------------

!$acc loop seq
    DO FEMEAN_M=1,NFRE
      FEMEAN_K = 1
      FEMEAN_TEMP2 = MAX(F(IJ, FEMEAN_K, FEMEAN_M), EPSMIN)
!$acc loop seq
      DO FEMEAN_K=2,NANG
        FEMEAN_TEMP2 = FEMEAN_TEMP2 + MAX(F(IJ, FEMEAN_K, FEMEAN_M), EPSMIN)
      END DO
      EM(IJ) = EM(IJ) + FEMEAN_TEMP2*DFIM(FEMEAN_M)
      FM(IJ) = FM(IJ) + DFIMOFR(FEMEAN_M)*FEMEAN_TEMP2
    END DO


    !*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
    !*       NORMALIZE WITH TOTAL ENERGY.
    !        ------------------------------------------

    EM(IJ) = EM(IJ) + DELT25*FEMEAN_TEMP2
    FM(IJ) = FM(IJ) + DELT2*FEMEAN_TEMP2
    FM(IJ) = EM(IJ) / FM(IJ)
    FM(IJ) = MAX(FM(IJ), FR(1))



!$acc end data
  END SUBROUTINE FEMEAN_LOKI
END MODULE FEMEAN_LOKI_MOD
