! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SEMEAN_LOKI_MOD
  CONTAINS
  SUBROUTINE SEMEAN_LOKI (FL1, KIJS, KIJL, EM, LLEPSMIN, IJ)

    ! ----------------------------------------------------------------------

    !**** *SEMEAN* - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.

    !     S.D. HASSELMANN.
    !     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER

    !*    PURPOSE.
    !     --------

    !       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.

    !**   INTERFACE.
    !     ----------

    !       *CALL* *SEMEAN(FL1, KIJS, KIJL, EM, LLEPSMIN)*
    !          *FL1*  - SPECTRUM.
    !          *KIJS* - LOCAL INDEX OF FIRST GRIDPOINT
    !          *KIJL* - LOCAL  INDEX OF LAST GRIDPOIN
    !          *EM*   - MEAN VARIANCE
    !          *LLEPSMIN* - TRUE IF THE WAVE ENERGY IS AT LEAST = EPSMIN

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

    USE YOWFRED, ONLY: FR, DFIM, DELTH, WETAIL
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: EPSMIN

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE

    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: EM
    LOGICAL, INTENT(IN) :: LLEPSMIN


    INTEGER(KIND=JWIM) :: SEMEAN_M, SEMEAN_K
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ

    REAL(KIND=JWRB) :: DELT25
    REAL(KIND=JPHOOK) :: SEMEAN_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: TEMP
!$acc routine seq
!$acc data present( FL1, EM )

    ! ----------------------------------------------------------------------


    !*    1. INITIALISE ENERGY ARRAY.
    !        ------------------------

    IF (LLEPSMIN) THEN
      EM(IJ) = EPSMIN
    ELSE
      EM(IJ) = 0.0_JWRB
    END IF

    ! ----------------------------------------------------------------------

    !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
    !        -----------------------------------------

!$acc loop seq
    DO SEMEAN_M=1,NFRE
      SEMEAN_K = 1
      TEMP = FL1(IJ, SEMEAN_K, SEMEAN_M)
!$acc loop seq
      DO SEMEAN_K=2,NANG
        TEMP = TEMP + FL1(IJ, SEMEAN_K, SEMEAN_M)
      END DO
      EM(IJ) = EM(IJ) + DFIM(SEMEAN_M)*TEMP
    END DO

    ! ----------------------------------------------------------------------

    !*    3. ADD TAIL ENERGY.
    !        ----------------

    DELT25 = WETAIL*FR(NFRE)*DELTH
    EM(IJ) = EM(IJ) + DELT25*TEMP



!$acc end data
  END SUBROUTINE SEMEAN_LOKI
END MODULE SEMEAN_LOKI_MOD
