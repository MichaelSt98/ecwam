! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE STHQ_LOKI_MOD
  CONTAINS
  SUBROUTINE STHQ_LOKI (KIJS, KIJL, FL1, THQ, IJ)

    ! ----------------------------------------------------------------------

    !**** *STHQ* - COMPUTATION OF MEAN WAVE DIRECTION AT EACH GRID POINT.

    !     S.D. HASSELMANN
    !     OPTIMIZED BY L. ZAMBRESKY

    !*    PURPOSE.
    !     --------

    !       TO COMPUTE MEAN WAVE DIRECTION AT EACH GRID POINT.

    !**   INTERFACE.
    !     ----------

    !       *CALL* *STHQ (KIJS, KIJL, FL1, THQ)*
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *FL1*    - SPECTRUM.
    !          *THQ*    - MEAN WAVE DIRECTION

    !     METHOD.
    !     -------

    !       INTEGRATION OF SPECTRUM TIMES SIN AND COS OVER DIRECTION.

    !     EXTERNALS.
    !     ----------

    !       NONE.

    !     REFERENCE.
    !     ----------

    !       NONE.

    ! ----------------------------------------------------------------------

    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: ZPI, EPSMIN
    USE YOWFRED, ONLY: DFIM, COSTH, SINTH

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: THQ

    INTEGER(KIND=JWIM) :: M, K
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    REAL(KIND=JWRB) :: TEMP, SI, CI
!$acc routine seq
!$acc data present( FL1, THQ )

    ! ----------------------------------------------------------------------


    !*    1. INITIALISE SIN AND COS ARRAYS.
    !        ------------------------------

    SI = 0.0_JWRB
    CI = 0.0_JWRB

    ! ----------------------------------------------------------------------

    !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
    !        ------------------------------------------

!$acc loop seq
    DO K=1,NANG
      TEMP = 0.0_JWRB
!$acc loop seq
      DO M=1,NFRE
        TEMP = TEMP + FL1(IJ, K, M)*DFIM(M)
      END DO
      SI = SI + SINTH(K)*TEMP
      CI = CI + COSTH(K)*TEMP
    END DO

    ! ----------------------------------------------------------------------

    !*    3. COMPUTE MEAN DIRECTION.
    !        -----------------------

    IF (CI == 0.0_JWRB) CI = EPSMIN
    THQ(IJ) = ATAN2(SI, CI)
    IF (THQ(IJ) < 0.0_JWRB) THQ(IJ) = THQ(IJ) + ZPI



!$acc end data
  END SUBROUTINE STHQ_LOKI
END MODULE STHQ_LOKI_MOD
