! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE OUTBETA_LOKI_MOD
  CONTAINS
  SUBROUTINE OUTBETA_LOKI (KIJS, KIJL, U10, USTAR, Z0B, CHRNCK, BETAM, BETAHQ, CD, IJ)

    ! ----------------------------------------------------------------------

    !**** *OUTBETA* - DETERMINES THE CHARNOCK PARAMETER

    !     P.JANSSEN      KNMI/ECMWF  JANUARY 1992
    !     J.BIDLOT       ECMWF       FEBRUARY 1996  MESSAGE PASSING
    !     J.BIDLOT       ECMWF       AUGUST 2008  REMOVE MESSAGE PASSING
    !     J.BIDLOT       ECMWF       MARCH 2014  USE MODEL CHARNOCK FOR ALL
    !                                SEA POINTS

    !*    PURPOSE.
    !     --------

    !       COMPUTES THE CHARNOCK PARAMETER.

    !**   INTERFACE.
    !     ----------

    !       *CALL* *OUTBETA (KIJS, KIJL,
    !                        U10, USTAR, Z0M, Z0B,
    !                        CHRNCK , BETAHQ, CD)
    !         *KIJS*    - INDEX OF FIRST GRIDPOINT.
    !         *KIJL*    - INDEX OF LAST GRIDPOINT.
    !         *U10*     - WIND SPEED.
    !         *USTAR*   - FRICTION VELOCITY.
    !         *Z0M*     - ROUGHNESS LENGTH.
    !         *Z0B*     - BACKGROUND ROUGHNESS LENGTH.
    !         *CHRNCK*  - MODEL CHARNOCK
    !         *BETA*    - OUTOUT CHARNOCK FOR ATMOSPHERE MODEL
    !         *BETAHQ*  - EQUIVALENT CHARNOCK FIELD FOR HEAT AND MOISTURE
    !                    (i.e. Charnock with the background roughness removed)
    !         *CD*      - OPTIONAL : DRAG COEFFICIENT THAT IS CONSISTENT WITH THE DERIVED CHARNOCK
    !                     AND WHETHER OR NOT VISCOSITY EFFECT WAS ACCOUNTED FOR IN Z0M
    !                     ASSUMING A LOGARTHMIC PROFILE FOR THE WIND.
    !

    !     EXTERNALS.
    !     ----------

    !       NONE.

    !     METHOD.
    !     -------

    !       NONE.

    !     REFERENCE.
    !     ----------

    !       NONE.

    ! ----------------------------------------------------------------------

    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWCOUP, ONLY: LLGCBZ0
    USE YOWPCONS, ONLY: G, GM1, EPSUS
    USE YOWPHYS, ONLY: XKAPPA, XNLEV, RNUM, PRCHAR, ALPHAMIN, ALPHAMAX, ALPHA

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: U10
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: USTAR
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: Z0B
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CHRNCK
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: BETAM
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: BETAHQ
    REAL(KIND=JWRB), OPTIONAL, INTENT(OUT), DIMENSION(KIJL) :: CD


    REAL(KIND=JWRB), PARAMETER :: AMAX = 0.02_JWRB
    REAL(KIND=JWRB), PARAMETER :: BMAX = 0.01_JWRB

    !     !!!!! BETAHQ_REDUCE was introduced to reduce the impact of the sea state dependent
    !     !!!!! heat and moisture surface flux. It is just a dirty tuning knob !!

    REAL(KIND=JWRB), PARAMETER :: BETAHQ_REDUCE = 0.25_JWRB

    INTEGER(KIND=JWIM), INTENT(IN) :: IJ

    REAL(KIND=JWRB) :: GUSM2, Z0ATM
    REAL(KIND=JWRB) :: USM
    REAL(KIND=JWRB) :: ALPHAMAXU10
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!$acc routine seq
!$acc data present( U10, USTAR, Z0B, CHRNCK, BETAM, BETAHQ, CD )

    ! ----------------------------------------------------------------------



    !*    COMPUTE CHARNOCK BETAM
    !     ----------------------

    IF (LLGCBZ0) THEN
      ALPHAMAXU10 = ALPHAMAX
    ELSE
      ALPHAMAXU10 = MIN(ALPHAMAX, AMAX + BMAX*U10(IJ))
    END IF


    USM = 1.0_JWRB / MAX(USTAR(IJ), EPSUS)
    GUSM2 = G*USM**2
    BETAM(IJ) = MAX(MIN(CHRNCK(IJ), ALPHAMAXU10), ALPHAMIN)
    BETAHQ(IJ) = MAX(BETAHQ_REDUCE*(BETAM(IJ) - MAX(Z0B(IJ)*GUSM2, ALPHA)), ALPHAMIN)

    IF (PRESENT(CD)) THEN
      !!!     we are assuming here that z0 = RNUM/USTAR + Charnock USTAR**2/g
      !!!     in order to fit with what is used in the IFS.
      Z0ATM = RNUM*USM + GM1*BETAM(IJ)*USTAR(IJ)**2
      CD(IJ) = (XKAPPA / LOG(1.0_JWRB + XNLEV / Z0ATM))**2
    END IF

    ! ----------------------------------------------------------------------

!$acc end data
  END SUBROUTINE OUTBETA_LOKI
END MODULE OUTBETA_LOKI_MOD
