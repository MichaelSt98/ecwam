! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SCOSFL_LOKI_MOD
  CONTAINS
  SUBROUTINE SCOSFL_LOKI (KIJS, KIJL, F, MM, MEANCOSFL, IJ)

    ! ----------------------------------------------------------------------

    !**** *SCOSFL* - COMPUTATION OF THE CENTERED FOURIER COEFFICIENT FOR
    !                A GIVEN FREQUENCY AT EACH GRID POINT.

    !     J.-R. BIDLOT        ECMWF       MARCH 2000

    !*    PURPOSE.
    !     --------

    !     TO COMPUTE THE CENTERED FOURIER COEFFICIENT FOR A GIVEN FREQUENCY

    !**   INTERFACE.
    !     ----------

    !       *CALL* *SCOSFL (KIJS, KIJL, F, MM, MEANCOSFL)*
    !          *KIJS* - INDEX OF FIRST GRIDPOINT
    !          *KIJL* - INDEX OF LAST GRIDPOINT
    !          *F*    - SPECTRUM.
    !          *MM*   - FREQUENCY INDEX (array)
    !          *MEANCOSFL* - CENTERED FOURIER COEFFICIENT

    !     METHOD.
    !     -------

    !       INTEGRATION OF SPECTRUM TIMES COS(THETA-MEAN THETA)
    !       OVER DIRECTION.

    !     EXTERNALS.
    !     ----------

    !       NONE.

    !     REFERENCE.
    !     ----------

    !       NONE.

    ! ----------------------------------------------------------------------

    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWFRED, ONLY: DELTH, TH, COSTH, SINTH

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: F
    INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(KIJL) :: MM
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: MEANCOSFL


    INTEGER(KIND=JWIM) :: SCOSFL_K
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    REAL(KIND=JPHOOK) :: SCOSFL_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: MEANDIR, SI, CI
!$acc routine seq
!$acc data present( F, MM, MEANCOSFL )

    ! ----------------------------------------------------------------------


    !*    1. INITIALISE ARRAYS
    !        -----------------

    SI = 0.0_JWRB
    CI = 0.0_JWRB
    MEANCOSFL(IJ) = 0.0_JWRB

    ! ----------------------------------------------------------------------

    !*    2. FIND MEAN DIRECTION
    !        -------------------

!$acc loop seq
    DO SCOSFL_K=1,NANG
      SI = SI + SINTH(SCOSFL_K)*F(IJ, SCOSFL_K, MM(IJ))
      CI = CI + COSTH(SCOSFL_K)*F(IJ, SCOSFL_K, MM(IJ))
    END DO

    IF (CI == 0.0_JWRB .and. SI == 0.0_JWRB) THEN
      MEANDIR = 0.0_JWRB
    ELSE
      MEANDIR = ATAN2(SI, CI)
    END IF

    !*    3. COMPUTE FOURRIER COEFFICIENT
    !        -----------------------------

!$acc loop seq
    DO SCOSFL_K=1,NANG
      MEANCOSFL(IJ) = MEANCOSFL(IJ) + COS(TH(SCOSFL_K) - MEANDIR)*F(IJ, SCOSFL_K, MM(IJ))
    END DO
    MEANCOSFL(IJ) = DELTH*MEANCOSFL(IJ)



!$acc end data
  END SUBROUTINE SCOSFL_LOKI
END MODULE SCOSFL_LOKI_MOD
