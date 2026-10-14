! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE OUTSETWMASK_LOKI_MOD
  CONTAINS
  SUBROUTINE OUTSETWMASK_LOKI (KIJS, KIJL, IODP, CICVR, BOUT, IJ)
    ! ----------------------------------------------------------------------

    !**** *OUTSETWMASK* -


    !*    PURPOSE.
    !     --------

    !    IMPOSE THE ICE MASK AND LAND MASK ON OUTPUT PARAMETERS


    !**   INTERFACE.
    !     ----------

    !        *CALL* *OUTSETWMASK (KIJS, KIJL, IODP, CICVR, BOUT)
    !         *KIJS*   - INDEX OF FIRST LOCAL GRIDPOINT.
    !         *KIJL*   - INDEX OF LAST LOCAL GRIDPOINT.
    !         *CICVR   - SEA ICE COVER FIELD.
    !         *IODP*   - LAND MASK IF IODP(IJ)=0
    !         *BOUT*   - OUTPUT PARAMETERS

    ! ----------------------------------------------------------------------
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWCOUT, ONLY: JPPFLAG, IPRMINFO, NIPRMOUT, ITOBOUT
    USE YOWICE, ONLY: LICERUN, CITHRSH
    USE YOWPCONS, ONLY: ZMISS
    USE YOWSTAT, ONLY: LLSOURCE

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(KIJL) :: IODP
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CICVR
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NIPRMOUT) :: BOUT

    INTEGER(KIND=JWIM) :: ITG, IR
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!$acc routine seq
!$acc data present( IODP, CICVR, BOUT )

    !----------------------------------------------------------------------

!$acc loop seq
    DO IR=1,JPPFLAG
      ITG = ITOBOUT(IR)
      IF (ITG > 0) THEN
        IF (LICERUN .and. LLSOURCE .and. IPRMINFO(IR, 6) == 1) THEN
          !         SEA ICE MASK IS APPLIED
          IF (CICVR(IJ) > CITHRSH) BOUT(IJ, ITG) = ZMISS
        END IF

        IF (IPRMINFO(IR, 7) == 1) THEN
          !           SEA MASK IS APPLIED
          BOUT(IJ, ITG) = BOUT(IJ, ITG)*IODP(IJ) + (1 - IODP(IJ))*ZMISS
        END IF

      END IF
    END DO



!$acc end data
  END SUBROUTINE OUTSETWMASK_LOKI
END MODULE OUTSETWMASK_LOKI_MOD
