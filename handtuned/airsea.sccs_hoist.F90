! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE AIRSEA_LOKI_MOD
  CONTAINS
  SUBROUTINE AIRSEA_LOKI (KIJS, KIJL, HALP, U10, U10DIR, TAUW, TAUWDIR, RNFAC, US, Z0, Z0B, CHRNCK, ICODE_WND, IUSFG, IJ)

    ! ----------------------------------------------------------------------

    !**** *AIRSEA* - DETERMINE TOTAL STRESS IN SURFACE LAYER.

    !     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
    !     JEAN BIDLOT         ECMWF     FEBRUARY 1999 : TAUT is already
    !                                                   SQRT(TAUT)
    !     JEAN BIDLOT         ECMWF     OCTOBER 2004: QUADRATIC STEP FOR
    !                                                 TAUW

    !*    PURPOSE.
    !     --------

    !       COMPUTE TOTAL STRESS.

    !**   INTERFACE.
    !     ----------

    !       *CALL* *AIRSEA (KIJS, KIJL, FL1, WAVNUM,
    !                       HALP, U10, U10DIR, TAUW, TAUWDIR, RNFAC,
    !                       US, Z0, Z0B, CHRNCK, ICODE_WND, IUSFG)*

    !          *KIJS*    - INDEX OF FIRST GRIDPOINT.
    !          *KIJL*    - INDEX OF LAST GRIDPOINT.
    !          *FL1*     - SPECTRA
    !          *WAVNUM*  - WAVE NUMBER
    !          *HALP*    - 1/2 PHILLIPS PARAMETER
    !          *U10*     - WINDSPEED U10.
    !          *U10DIR*  - WINDSPEED DIRECTION.
    !          *TAUW*    - WAVE STRESS.
    !          *TAUWDIR* - WAVE STRESS DIRECTION.
    !          *RNFAC*   - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !          *US*      - OUTPUT OR OUTPUT BLOCK OF FRICTION VELOCITY.
    !          *Z0*      - OUTPUT BLOCK OF ROUGHNESS LENGTH.
    !          *Z0B*     - BACKGROUND ROUGHNESS LENGTH.
    !          *CHRNCK*  - CHARNOCK COEFFICIENT
    !          *ICODE_WND* SPECIFIES WHICH OF U10 OR US HAS BEEN FILED UPDATED:
    !                     U10: ICODE_WND=3 --> US will be updated
    !                     US:  ICODE_WND=1 OR 2 --> U10 will be updated
    !          *IUSFG*   - IF = 1 THEN USE THE FRICTION VELOCITY (US) AS FIRST GUESS in TAUT_Z0
    !                           0 DO NOT USE THE FIELD US


    ! ----------------------------------------------------------------------

    USE STRESS_GC_LOKI_MOD, ONLY: STRESS_GC_LOKI
    USE CHNKMIN_LOKI_MOD, ONLY: CHNKMIN_LOKI
    USE YOWCOUP, ONLY: LLCAPCHNK, LLGCBZ0
    USE YOWPCONS, ONLY: G, GM1, EPSUS, EPSMIN, ACD, BCD, ACDLIN, BCDLIN, CDMAX
    USE YOWTABL, ONLY: EPS1, IAB, SWELLFT
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPHYS, ONLY: XKAPPA, XNLEV, RNU, RNUM, ALPHA, ALPHAMIN, ALPHAMAX, ANG_GC_A, ANG_GC_B, ANG_GC_C
    USE YOWTEST, ONLY: IU06
    USE YOWWIND, ONLY: WSPMIN

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    INTERFACE
      SUBROUTINE ABORT1 ()
      END SUBROUTINE ABORT1
    END INTERFACE
    INTERFACE
      SUBROUTINE Z0WAVE (KIJS, KIJL, US, TAUW, UTOP, Z0, Z0B, CHRNCK)
        USE parkind_wave, ONLY: jwim, jwrb
        INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
        REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: US, TAUW, UTOP
        REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: Z0, Z0B, CHRNCK
      END SUBROUTINE Z0WAVE
    END INTERFACE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL, ICODE_WND, IUSFG
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: HALP, U10DIR, TAUW, TAUWDIR, RNFAC
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: U10, US
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: Z0, Z0B, CHRNCK

    INTEGER(KIND=JWIM) :: I, J
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ

    REAL(KIND=JWRB) :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, UST2, ARG, SQRTCDM1
    REAL(KIND=JWRB) :: XKAPPAD, XLOGLEV
    REAL(KIND=JWRB) :: XLEV
    REAL(KIND=JPHOOK) :: AIRSEA_ZHOOK_HANDLE
    INTEGER(KIND=JWIM), PARAMETER :: NITER = 18
    REAL(KIND=JWRB), PARAMETER :: TWOXMP1 = 3.0_JWRB
    REAL(KIND=JWRB), PARAMETER :: PMAX = 0.99_JWRB
    INTEGER(KIND=JWIM) :: ITER
    INTEGER(KIND=JWIM) :: IFRPH
    REAL(KIND=JWRB) :: ALPHAGM1
    REAL(KIND=JWRB), PARAMETER :: Z0MIN = 0.000001_JWRB
    REAL(KIND=JWRB) :: PCE_GC
    REAL(KIND=JWRB) :: Z0MINRST
    REAL(KIND=JWRB) :: CHARNOCK_MIN
    REAL(KIND=JWRB) :: COSDIFF
    REAL(KIND=JWRB) :: ZCHAR
    REAL(KIND=JWRB) :: US2TOTAUW, USMAX
    REAL(KIND=JWRB) :: XLOGXL, XKUTOP, XOLOGZ0
    REAL(KIND=JWRB) :: USTOLD, USTNEW, TAUOLD, TAUNEW, X, F, DELF, CDFG
    REAL(KIND=JWRB) :: USNRF, Z0NRF, Z0BNRF, ALPOG
    REAL(KIND=JWRB) :: USTM1, Z0TOT, Z0CH, Z0VIS, HZ0VISO1MX, ZZ
    REAL(KIND=JWRB) :: CONST, TAUV, DEL
    REAL(KIND=JWRB) :: RNUEFF, RNUKAPPAM1
    REAL(KIND=JPHOOK) :: TAUT_Z0_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: ALPHAOG, XMIN
    REAL(KIND=JWRB) :: W1
    REAL(KIND=JWRB) :: TAUWACT, TAUWEFF
    REAL(KIND=JWRB) :: ANG_GC, TAUUNR
    LOGICAL :: LLCOSDIFF
#include "cdm.func.h"
!$acc routine seq
!$acc data present( HALP, U10, U10DIR, TAUW, TAUWDIR, RNFAC, US, Z0, Z0B, CHRNCK )

    ! ----------------------------------------------------------------------

    !*    2. DETERMINE TOTAL STRESS (if needed)
    !        ----------------------------------

    IF (ICODE_WND == 3) THEN

      ! [Loki] inlined child subroutine: TAUT_Z0
      ! =========================================

      ! ----------------------------------------------------------------------


      XLOGXL = LOG(XNLEV)
      US2TOTAUW = 1.0_JWRB + EPS1

      RNUEFF = 0.04_JWRB*RNU
      RNUKAPPAM1 = RNUEFF / XKAPPA

      PCE_GC = 0.001_JWRB*IUSFG + (1 - IUSFG)*0.005_JWRB

      !     ONLY take the contribution of TAUW that is in the wind direction
      COSDIFF = COS(U10DIR(IJ) - TAUWDIR(IJ))
      TAUWACT = MAX(TAUW(IJ)*COSDIFF, EPSMIN)
      LLCOSDIFF = COSDIFF > 0.9_JWRB

      !  USING THE CG MODEL:
      IF (LLGCBZ0) THEN

        IF (LLCAPCHNK) THEN
          CHARNOCK_MIN = CHNKMIN_LOKI(U10(IJ))
          ALPHAOG = CHARNOCK_MIN*GM1
        ELSE
          ALPHAOG = 0.0_JWRB
        END IF

        USMAX = MAX(-0.21339_JWRB + 0.093698_JWRB*U10(IJ) - 0.0020944_JWRB*U10(IJ)**2 + 5.5091E-5_JWRB*U10(IJ)**3, 0.03_JWRB)
        TAUWEFF = MIN(TAUWACT*US2TOTAUW, USMAX**2)


        IF (IUSFG == 0) THEN
          ALPHAGM1 = ALPHA*GM1
          IF (U10(IJ) < 1.0_JWRB) THEN
            CDFG = 0.002_JWRB
          ELSE IF (LLCOSDIFF) THEN
            X = MIN(TAUWACT / MAX(US(IJ), EPSUS)**2, PMAX)
            ZCHAR = MIN(ALPHAGM1*US(IJ)**2 / SQRT(1.0_JWRB - X), 0.05_JWRB*EXP(-0.05_JWRB*(U10(IJ) - 35._JWRB)))
            ZCHAR = MIN(ZCHAR, ALPHAMAX)
            CDFG = ACDLIN + BCDLIN*SQRT(ZCHAR)*U10(IJ)
          ELSE
            CDFG = CDM(U10(IJ))
          END IF
          US(IJ) = U10(IJ)*SQRT(CDFG)
        END IF

        W1 = 0.85_JWRB - 0.05_JWRB*(TANH(10.0_JWRB*(U10(IJ) - 5.0_JWRB)) + 1.0_JWRB)

        XKUTOP = XKAPPA*U10(IJ)

        USTOLD = US(IJ)
        TAUOLD = USTOLD**2

!$acc loop seq
        DO ITER=1,NITER
          !         Z0 IS DERIVED FROM THE NEUTRAL LOG PROFILE: UTOP = (USTAR/XKAPPA)*LOG((XNLEV+Z0)/Z0)
          Z0(IJ) = MAX(XNLEV / (EXP(MIN(XKUTOP / USTOLD, 50.0_JWRB)) - 1.0_JWRB), Z0MIN)
          ! Viscous kinematic stress nu_air * dU/dz at z=0 of the neutral log profile reduced by factor 25 (0.04)
          TAUV = RNUKAPPAM1*USTOLD / Z0(IJ)

          ANG_GC = ANG_GC_A + ANG_GC_B*TANH(ANG_GC_C*TAUOLD)

          TAUUNR = STRESS_GC_LOKI(ANG_GC, US(IJ), Z0(IJ), Z0MIN, HALP(IJ), RNFAC(IJ))

          !         TOTAL kinematic STRESS:
          TAUNEW = TAUWEFF + TAUV + TAUUNR
          USTNEW = SQRT(TAUNEW)
          US(IJ) = W1*USTOLD + (1.0_JWRB - W1)*USTNEW

          !         CONVERGENCE ?
          DEL = US(IJ) - USTOLD
          IF (ABS(DEL) < PCE_GC*US(IJ)) EXIT
          TAUOLD = US(IJ)**2
          USTOLD = US(IJ)
        END DO

        X = TAUWEFF / TAUOLD

        ! protection just in case there is no convergence
        IF (ITER > NITER .and. X >= PMAX) THEN
          CDFG = CDM(U10(IJ))
          US(IJ) = U10(IJ)*SQRT(CDFG)
          Z0MINRST = US(IJ)**2*ALPHA*GM1
          Z0(IJ) = MAX(XNLEV / (EXP(XKUTOP / US(IJ)) - 1.0_JWRB), Z0MINRST)
          Z0B(IJ) = Z0MINRST
        ELSE
          Z0(IJ) = MAX(XNLEV / (EXP(XKUTOP / US(IJ)) - 1.0_JWRB), Z0MIN)
          Z0B(IJ) = Z0(IJ)*SQRT(TAUUNR / TAUOLD)
        END IF

        !       Refine solution

        IF (X < PMAX) THEN

          USNRF = US(IJ)
          Z0NRF = Z0(IJ)
          Z0BNRF = Z0B(IJ)

          USTOLD = US(IJ)
          TAUOLD = MAX(USTOLD**2, TAUWEFF)

          ALPOG = MAX(MIN(Z0B(IJ) / TAUOLD, ALPHAMAX), ALPHAOG)

!$acc loop seq
          DO ITER=1,NITER
            X = MIN(TAUWEFF / TAUOLD, PMAX)
            USTM1 = 1.0_JWRB / MAX(USTOLD, EPSUS)
            Z0VIS = RNUM*USTM1
            HZ0VISO1MX = 0.5_JWRB*Z0VIS / (1.0_JWRB - X)
            Z0B(IJ) = ALPOG*TAUOLD
            Z0(IJ) = HZ0VISO1MX + SQRT(HZ0VISO1MX**2 + Z0B(IJ)**2 / (1.0_JWRB - X))

            XOLOGZ0 = 1.0_JWRB / LOG(XNLEV / Z0(IJ) + 1.0_JWRB)
            F = USTOLD - XKUTOP*XOLOGZ0
            ZZ = 2.0_JWRB*USTM1*(3.0_JWRB*Z0B(IJ)**2 + 0.5_JWRB*Z0VIS*Z0(IJ) - Z0(IJ)**2) / (2.0_JWRB*Z0(IJ)**2*(1.0_JWRB - X) -  &
            & Z0VIS*Z0(IJ))

            DELF = 1.0_JWRB - XKUTOP*XOLOGZ0**2*ZZ

            IF (DELF /= 0.0_JWRB) US(IJ) = USTOLD - F / DELF
            !           CONVERGENCE ?
            TAUNEW = MAX(US(IJ)**2, TAUWEFF)
            US(IJ) = SQRT(TAUNEW)
            DEL = TAUNEW - TAUOLD
            IF (ABS(DEL) < PCE_GC*TAUOLD) EXIT
            TAUOLD = TAUNEW
            USTOLD = US(IJ)

          END DO
          ! protection just in case there is no convergence
          IF (ITER > NITER) THEN
            US(IJ) = USNRF
            Z0(IJ) = Z0NRF
            Z0B(IJ) = Z0BNRF
            USTM1 = 1.0_JWRB / MAX(US(IJ), EPSUS)
            Z0VIS = RNUM*USTM1
            CHRNCK(IJ) = MAX(G*(Z0(IJ) - Z0VIS)*USTM1**2, ALPHAMIN)

          ELSE
            CHRNCK(IJ) = MAX(G*(Z0B(IJ) / SQRT(1.0_JWRB - X)) / MAX(US(IJ), EPSUS)**2, ALPHAMIN)
          END IF

        ELSE
          USTM1 = 1.0_JWRB / MAX(US(IJ), EPSUS)
          Z0VIS = RNUM*USTM1
          CHRNCK(IJ) = MAX(G*(Z0(IJ) - Z0VIS)*USTM1**2, ALPHAMIN)
        END IF



      ELSE

        TAUWEFF = TAUWACT*US2TOTAUW

        IF (LLCAPCHNK) THEN
          CHARNOCK_MIN = CHNKMIN_LOKI(U10(IJ))
          XMIN = 0.15_JWRB*(ALPHA - CHARNOCK_MIN)
          ALPHAOG = CHARNOCK_MIN*GM1
        ELSE
          XMIN = 0.0_JWRB
          ALPHAOG = ALPHA*GM1
        END IF

        XKUTOP = XKAPPA*U10(IJ)

        USTOLD = (1 - IUSFG)*U10(IJ)*SQRT(MIN(ACD + BCD*U10(IJ), CDMAX)) + IUSFG*US(IJ)
        TAUOLD = MAX(USTOLD**2, TAUWEFF)
        US(IJ) = SQRT(TAUOLD)
        USTM1 = 1.0_JWRB / MAX(US(IJ), EPSUS)

!$acc loop seq
        DO ITER=1,NITER
          X = MAX(TAUWACT / TAUOLD, XMIN)
          Z0CH = ALPHAOG*TAUOLD / SQRT(1.0_JWRB - X)
          Z0VIS = RNUM*USTM1
          Z0TOT = Z0CH + Z0VIS

          XOLOGZ0 = 1.0_JWRB / (XLOGXL - LOG(Z0TOT))
          F = US(IJ) - XKUTOP*XOLOGZ0
          ZZ = USTM1*(Z0CH*(2.0_JWRB - TWOXMP1*X) / (1.0_JWRB - X) - Z0VIS) / Z0TOT
          DELF = 1.0_JWRB - XKUTOP*XOLOGZ0**2*ZZ

          IF (DELF /= 0.0_JWRB) US(IJ) = US(IJ) - F / DELF
          TAUNEW = MAX(US(IJ)**2, TAUWEFF)
          US(IJ) = SQRT(TAUNEW)
          IF (TAUNEW == TAUOLD) EXIT
          USTM1 = 1.0_JWRB / MAX(US(IJ), EPSUS)
          TAUOLD = TAUNEW
        END DO

        Z0(IJ) = Z0CH
        Z0B(IJ) = ALPHAOG*TAUOLD
        CHRNCK(IJ) = MAX(G*Z0(IJ)*USTM1**2, ALPHAMIN)


      END IF


      ! =========================================

    ELSE IF (ICODE_WND == 1 .or. ICODE_WND == 2) THEN

      ! [Loki] Removed content of pragma-marked region!
    END IF



!$acc end data
  END SUBROUTINE AIRSEA_LOKI
END MODULE AIRSEA_LOKI_MOD
