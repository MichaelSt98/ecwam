! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE INTPOL_LOKI_MOD
  CONTAINS
  SUBROUTINE INTPOL_LOKI (KIJS, KIJL, FLR, FLA, WAVNUM, UCUR, VCUR, IRA, IJ, DFTH)

    ! ----------------------------------------------------------------------

    !**** *INTPOL* - TRANSFORMATION OF SPECTRA FROM SIGMA TO OMEGA.

    !     S.D.HASSELMANN      MPI            1.1.91
    !     H. GUNTHER          GKSS/ECMWF     1.2.91  MODIFIED FOR CYCLE_4

    !*    PURPOSE.
    !     --------

    !       TRANSFORMATION OF A MOVING COORDINATE SYSTEM TO AN ABSOLUTE
    !       COORDINATE SYSTEM.

    !**   INTERFACE.
    !     ----------

    !       *CALL* *INTPOL (KIJS, KIJL, FLR, FLA, WAVNUM, UCUR, VCUR, IRA)*
    !         *KIJS*  - INDEX OF FIRST GRIDPOINT.
    !         *KIJL*  - INDEX OF LAST GRIDPOINT.
    !         *FLR*   - SPECTRA (INPUT)
    !         *FLA*   - SPECTRA (OUTPUT)
    !         *WAVNUM*- WAVE NUMBER.
    !         *UCUR*  - U-COMPONENT OF THE SURFACE CURRENT
    !         *VCUR*  - V-COMPONENT OF THE SURFACE CURRENT
    !         *IRA*   - = 1 TRANSFORMATION FROM MOVING TO ABSOLUTE COORD.
    !                   =-1 TRANSFORMATION FROM ABSOLUTE TO MOVING COORD.

    !     METHOD.
    !     -------

    !       SCATTERING TO NEIGHBOURING POINT.

    !     EXTERNALS.
    !     ----------

    !       NONE.

    !     REFERENCE.
    !     ----------

    !       NONE.

    ! ----------------------------------------------------------------------

    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWCURR, ONLY: CURRENT_MAX
    USE YOWFRED, ONLY: FR, DFIM, COSTH, SINTH, DELTH, FRATIO, FLOGSPRDM1, FR5
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: G, ZPI, EPSMIN
    USE YOWTEST, ONLY: IU06

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    INTERFACE
      SUBROUTINE ABORT1 ()
      END SUBROUTINE ABORT1
    END INTERFACE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FLR
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: FLA
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: UCUR, VCUR
    INTEGER(KIND=JWIM), INTENT(IN) :: IRA


    INTEGER(KIND=JWIM) :: M, K
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: NFRE_MAX, NEWM, NEWM1, KH
    INTEGER(KIND=JWIM) :: NEWF, NEWFLA, KNEW, zero_current_bin

    REAL(KIND=JWRB) :: FRE0, CDF, ZPI2GM, COEF, FMAX, FREQ, DFREQTH, FR5OFREQ5
    REAL(KIND=JWRB) :: FNEW, GWH
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NFRE) :: DFTH
    REAL(KIND=JWRB) :: OLDFL, WAVN
    REAL(KIND=JWRB) :: FNEF, GWP, GWM

    LOGICAL :: LICE2SEA, zero_current
!$acc routine seq
!$acc data present( FLR, FLA, WAVNUM, UCUR, VCUR, DFTH )
    ! ----------------------------------------------------------------------


    !*    0. INITIAL OUTPUT ARRAY WITH ZERO.
    !        -------------------------------

    FRE0 = FRATIO - 1.0_JWRB
    ZPI2GM = ZPI**2 / G

    COEF = IRA / ZPI

    !!    MAXIMUM EXTEND OF THE HIGH FRQUENCY TAIL THAT CAN BE REMAPPED ONTO THE FREQUENCY DISCRETISATION
    FMAX = FR(NFRE) + (ZPI / G)*FR(NFRE)**2*CURRENT_MAX
    NFRE_MAX = FLOOR(LOG10(FMAX / FR(1))*FLOGSPRDM1) + 1

    !!!??? I believe that in order to be energy conserving, since
    !!! DFIM(1) is only the part above FR(1)
    !!! purpose of trapezoidal rule, one has to take DFTH(1) which
    !!! define on both side of FR(1). Similarly for FR(NFRE)
    !!! Redefine the centered frequency-direction intervals
    CDF = 0.5_JWRB*(FRATIO - 1.0_JWRB / FRATIO)*DELTH
!$acc loop seq
    DO M=1,NFRE
      DFTH(M) = FR(M)*CDF
    END DO

    IF (ABS(IRA) /= 1) THEN
    END IF

    LICE2SEA = .true.
    zero_current = UCUR(IJ) == 0.0_JWRB .AND. VCUR(IJ) == 0.0_JWRB

!$acc loop seq
    DO K=1,NANG
!$acc loop seq
      DO M=1,NFRE
        IF (FLR(IJ, K, M) > EPSMIN) THEN
          LICE2SEA = .false.
          EXIT
        END IF
      END DO
      IF (.not.LICE2SEA) EXIT
    END DO


!$acc loop seq
    DO M=1,NFRE
!$acc loop seq
      DO K=1,NANG
        FLA(IJ, K, M) = 0.0_JWRB
      END DO
    END DO

    !*    1. LOOP OVER FREQUENCIES.
    !        ----------------------

!$acc loop seq
    DO M=1,NFRE_MAX

      IF (M <= NFRE) THEN
        FREQ = FR(M)
        DFREQTH = DFTH(M)
        WAVN = WAVNUM(IJ, M)
      ELSE
        FREQ = FR(NFRE)*FRATIO**(M - NFRE)
        DFREQTH = FREQ*CDF
        WAVN = ZPI2GM*FREQ**2
      END IF

      FR5OFREQ5 = FR5(NFRE) / (FREQ**5)
      IF (zero_current) zero_current_bin = FLOOR(LOG10(FREQ / FR(1))*FLOGSPRDM1) + 1

      !*    1.3 LOOP OVER DIRECTONS.
      !         --------------------

!$acc loop seq
      DO K=1,NANG

        !*    1.3.1 NEW FREQUENCY AND DIRECTION AT ALL GRIDPOINTS.
        !           ----------------------------------------------

        IF (zero_current) THEN
          FNEF = FREQ
          KNEW = K
          NEWF = zero_current_bin
        ELSE
          FNEF = FREQ + COEF*WAVN*(COSTH(K)*VCUR(IJ) + SINTH(K)*UCUR(IJ))
          IF (FNEF > 0.0_JWRB) THEN
            KNEW = K
          ELSE
            KNEW = MOD(K + NANG / 2 - 1, NANG) + 1
            FNEF = -FNEF
          END IF

          !*    1.3.2 NEW FREQUENCY BIN NUMBER AT ALL GRIDPOINTS.
          !           -------------------------------------------

          IF (FNEF <= FR(1) / FRATIO) THEN
            NEWF = -1
          ELSE
            NEWF = FLOOR(LOG10(FNEF / FR(1))*FLOGSPRDM1) + 1
          END IF
        END IF

        !*    1.3.3 INTERPOLATED ENERGY DENSITIES AT ALL GRIDPOINTS.
        !           ------------------------------------------------
        IF (M <= NFRE) THEN
          IF (LICE2SEA) THEN
            OLDFL = 0.0_JWRB
          ELSE
            OLDFL = FLR(IJ, K, M)
          END IF
        ELSE
          IF (LICE2SEA) THEN
            OLDFL = 0.0_JWRB
          ELSE
            !! Extend with a f**-5 tail
            OLDFL = FLR(IJ, K, NFRE)*FR5OFREQ5
          END IF
        END IF

        FNEW = FNEF
        NEWM = NEWF
        IF (NEWM < NFRE .and. NEWM >= 1) THEN
          !!!! DFTH is only defined between 1 and NFRE
          NEWM1 = NEWM + 1
          GWH = DFREQTH / (FR(NEWM1) - FR(NEWM))*OLDFL
          GWM = GWH*(FR(NEWM1) - FNEW) / DFTH(NEWM)
          GWP = GWH*(FNEW - FR(NEWM)) / DFTH(NEWM1)
          NEWFLA = NEWM1
        ELSE IF (NEWM == 0) THEN
          GWH = FRATIO*DFREQTH / (FRE0*FR(1))*OLDFL
          GWP = GWH*(FNEW - FR(1) / FRATIO) / DFTH(1)
          NEWF = -1
          NEWFLA = 1
        ELSE IF (NEWM == NFRE) THEN
          GWH = DFREQTH / (FRE0*FR(NFRE))*OLDFL
          GWM = GWH*(FRATIO*FR(NFRE) - FNEW) / DFTH(NFRE)
          NEWFLA = -1
        ELSE
          NEWF = -1
          NEWFLA = -1
        END IF

        !*    1.3.4 NEW SPECTRUM AT ALL GRIDPOINTS.
        !           -------------------------------

        NEWM = NEWF
        NEWM1 = NEWFLA
        KH = KNEW
        IF (NEWM /= -1) FLA(IJ, KH, NEWM) = FLA(IJ, KH, NEWM) + GWM
        IF (NEWM1 /= -1) FLA(IJ, KH, NEWM1) = FLA(IJ, KH, NEWM1) + GWP

        !*    BRANCH BACK TO 1.3 FOR NEXT DIRECTION.

      END DO

      !*    BRANCH BACK TO 1. FOR NEXT FREQUENCY.

    END DO

!$acc loop seq
    DO M=1,NFRE
!$acc loop seq
      DO K=1,NANG
        FLA(IJ, K, M) = MAX(FLA(IJ, K, M), EPSMIN)
      END DO
    END DO



!$acc end data
  END SUBROUTINE INTPOL_LOKI
END MODULE INTPOL_LOKI_MOD
