! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE PARMEAN_LOKI_MOD
  CONTAINS
  SUBROUTINE PARMEAN_LOKI (KIJS, KIJL, NPMAX, NPEAK, SPEC, ENE, DIR, PER, IJ, F1D, TEMP, SI, CI, EM, FM, THQ)

    ! ----------------------------------------------------------------------

    !**** *PARMEAN* - COMPUTATION OF TOTAL ENERGY, MEAN DIRECTION,
    !                  AND MEAN PERIOD BASED ON THE FIRST MOMENT

    !      J-R BIDLOT    ECMWF     MARCH 2000
    !      D PETTENUZZO  MAY 2012 MERGED 3 ROUTINES FOR PAR COMPUTATION

    !*    PURPOSE.
    !     --------

    !       COMPUTE TOTAL ENERGY, MEAN DIRECTION AND MEAN PERIOD.

    !**   INTERFACE.
    !     ----------

    !       *CALL* *PARMEAN (KIJS, KIJL, NPMAX, NPEAK,
    !                        SPEC,
    !                        ENE, DIR, PER)
    !              *KIJS*    - INDEX OF FIRST GRIDPOINT
    !              *KIJL*    - INDEX OF LAST GRIDPOINT
    !              *NPMAX*   - MAXIMUM NUMBER OF PARTITIONS
    !              *NPEAK*   - NUMBER OF PEAKS
    !              *SPEC*    - SPECTRUM.
    !              *ENE*     - MEAN WAVE ENERGY FOR THE WAVE SYSTEM
    !              *DIR*     - MEAN WAVE DIRECTION
    !              *PER*  - MEAN PERIOD BASED ON -1 MOMENT.

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

    USE YOWFRED, ONLY: FR, DFIM, DFIMOFR, COSTH, SINTH
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: EPSMIN, ZPI

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL, NPMAX
    INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(KIJL) :: NPEAK
    REAL(KIND=JWRB), INTENT(IN) :: SPEC(NANG, NFRE, NPMAX, KIJL)
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, 0:NPMAX) :: ENE, DIR, PER


    INTEGER(KIND=JWIM) :: IPK, K, M, IP
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    REAL(KIND=JWRB), INTENT(INOUT) :: F1D(NFRE)
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NPMAX) :: TEMP
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NPMAX) :: SI
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NPMAX) :: CI
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NPMAX) :: EM
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NPMAX) :: FM
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NPMAX) :: THQ
!$acc routine seq
!$acc data present( NPEAK, SPEC, ENE, DIR, PER, F1D, TEMP, SI, CI, EM, FM, THQ )

    ! ----------------------------------------------------------------------



    !     ENE must be initialised and index 0 must stay equal to 0
!$acc loop seq
    DO IP=0,NPMAX
      ENE(IJ, IP) = 0.0_JWRB
      DIR(IJ, IP) = 0.0_JWRB
      PER(IJ, IP) = 0.0_JWRB
    END DO

    !     LOOP OVER GRID POINT
!$acc loop seq
    DO IPK=1,NPEAK(IJ)
      EM(IPK) = EPSMIN
      FM(IPK) = EPSMIN
!$acc loop seq
      DO M=1,NFRE
        F1D(M) = SPEC(1, M, IPK, IJ)
!$acc loop seq
        DO K=2,NANG
          F1D(M) = F1D(M) + SPEC(K, M, IPK, IJ)
        END DO
        EM(IPK) = EM(IPK) + F1D(M)*DFIM(M)
        FM(IPK) = FM(IPK) + F1D(M)*DFIMOFR(M)
      END DO

      SI(IPK) = 0.0_JWRB
      CI(IPK) = 0.0_JWRB
!$acc loop seq
      DO K=1,NANG
        TEMP(IPK) = SPEC(K, 1, IPK, IJ)*DFIM(1)
!$acc loop seq
        DO M=2,NFRE
          TEMP(IPK) = TEMP(IPK) + SPEC(K, M, IPK, IJ)*DFIM(M)
        END DO
        SI(IPK) = SI(IPK) + SINTH(K)*TEMP(IPK)
        CI(IPK) = CI(IPK) + COSTH(K)*TEMP(IPK)
        IF (CI(IPK) == 0.0_JWRB) CI(IPK) = EPSMIN
        THQ(IPK) = ATAN2(SI(IPK), CI(IPK))
        IF (THQ(IPK) < 0.0_JWRB) THQ(IPK) = THQ(IPK) + ZPI
      END DO
    END DO

!$acc loop seq
    DO IPK=1,NPEAK(IJ)
      IF (EM(IPK) > EPSMIN) THEN
        ENE(IJ, IPK) = EM(IPK)
        PER(IJ, IPK) = FM(IPK) / EM(IPK)
        DIR(IJ, IPK) = THQ(IPK)
      END IF
    END DO




!$acc end data
  END SUBROUTINE PARMEAN_LOKI
END MODULE PARMEAN_LOKI_MOD
