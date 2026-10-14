! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE DOMINANT_PERIOD_LOKI_MOD
  CONTAINS
  SUBROUTINE DOMINANT_PERIOD_LOKI (KIJS, KIJL, FL1, DP, IJ, F1D4)

    ! ----------------------------------------------------------------------

    !**** *DOMINANT_PERIOD* - COMPUTATION OF THE DOMINANT PERIOD

    !*    PURPOSE.
    !     --------

    !       COMPUTE DOMINANT PERIOD AT EACH GRID POINT.

    !**   INTERFACE.
    !     ----------

    !       *CALL* *DOMINANT_PERIOD (KIJS, KIJL, FL1, DP)*
    !              *KIJS*    - INDEX OF FIRST GRIDPOINT
    !              *KIJL*    - INDEX OF LAST GRIDPOINT
    !              *FL1*     - 2D-SPECTRUM.
    !              *DP*      - DOMINANT PERIOD

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

    USE YOWFRED, ONLY: FR, DFIM, DFIMFR, DELTH
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: EPSMIN

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------
    IMPLICIT NONE


    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: DP

    REAL(KIND=JWRB), PARAMETER :: FLTHRS = 0.1_JWRB

    INTEGER(KIND=JWIM) :: K, M
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    REAL(KIND=JWRB) :: TEMP, EM, FCROP
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: F1D4

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!$acc routine seq
!$acc data present( FL1, DP, F1D4 )

    ! ----------------------------------------------------------------------


!$acc loop seq
    DO M=1,NFRE
      F1D4(IJ, M) = 0.0_JWRB
    END DO

    EM = 0.0_JWRB
    DP(IJ) = 0.0_JWRB
    FCROP = 0.0_JWRB

!$acc loop seq
    DO M=1,NFRE
!$acc loop seq
      DO K=1,NANG
        IF (FL1(IJ, K, M) > FCROP) THEN
          FCROP = FL1(IJ, K, M)
        END IF
      END DO
    END DO
    FCROP = FLTHRS*FCROP

!$acc loop seq
    DO M=1,NFRE
!$acc loop seq
      DO K=1,NANG
        IF (FL1(IJ, K, M) > FCROP) F1D4(IJ, M) = F1D4(IJ, M) + FL1(IJ, K, M)*DELTH
      END DO
    END DO

!$acc loop seq
    DO M=1,NFRE
      F1D4(IJ, M) = F1D4(IJ, M)**4
      EM = EM + DFIM(M)*F1D4(IJ, M)
      DP(IJ) = DP(IJ) + DFIMFR(M)*F1D4(IJ, M)
    END DO

    IF (EM > 0.0_JWRB .and. DP(IJ) > EPSMIN) THEN
      DP(IJ) = EM / DP(IJ)
    ELSE
      DP(IJ) = 0.0_JWRB
    END IF



!$acc end data
  END SUBROUTINE DOMINANT_PERIOD_LOKI
END MODULE DOMINANT_PERIOD_LOKI_MOD
