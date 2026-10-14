! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE FNDPRT_LOKI_MOD
  CONTAINS
  SUBROUTINE FNDPRT_LOKI (KIJS, KIJL, NPMAX, NPEAK, MIJ, NTHP, NFRP, FLLOW, LLCOSDIFF, FLNOISE, FL1, SWM, SUNASGN, ENE, DIR,  &
  & PER, IJ, KLOC, W2, W1, SPEC, LLW3, PARMEAN_F1D, PARMEAN_TEMP, PARMEAN_SI, PARMEAN_CI, PARMEAN_EM, PARMEAN_FM, PARMEAN_THQ)

    ! ----------------------------------------------------------------------

    !**** *FNDPRT* -  CALCULATE PARTITION MASKS

    !     D.PETTENUZZO     MAY 2011
    !     MODIFIED MAY 2012 TO SUIT ECMWF CODES


    !*    PURPOSE.
    !     --------

    !     FIND ALL THE POINTS IN THE FREQUENCY DIRECTION DOMAIN
    !     WHICH BELONG TO A GIVEN PARTITION IDENTIFIED BY ITS
    !     PEAK

    !**   INTERFACE.
    !     ----------

    !       *CALL* *FNDPRT (KIJS, KIJL, NPMAX,
    !     &                 NPEAK, MIJ, NTHP, NFRP,
    !     &                 FLLOW, LLCOSDIFF, FLNOISE,
    !     &                 FL1, SWM, SUNASGN,
    !     &                 ENE, DIR, PER)

    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *NPMAX*  - MAXIMUM NUMBER OF PARTITIONS
    !          *NPEAK*  - NUMBER OF PEAKS
    !          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
    !          *NTHP*   - DIRECTION INDEX FOR ALL PEAKS
    !          *NFRP*   - FREQUENCY INDEX FOR ALL PEAKS
    !          *FLLOW*  - MINIMUM ENERGY LEVEL
    !          *LLCOSDIFF* LOGICAL (COSDIFF.LT.-0.4)
    !          *FLNOISE* - SPECTRAL NOISE LEVEL
    !          *FL1*    - BLOCK OF SPECTRA
    !          *SWM*    - ORIGINAL SWELL MASK
    !                     IT MIGHT NEED TO BE ADJUSTED
    !                     THIS IS POSSIBLE BECAUSE WE SUPPLY SWELL SPECTRA
    !                     IN WHICH THE WIND SEA HAS BEEN REMOVED. SOME OF
    !                     THE ENERGY LEAKING OUT FROM THE WINDSEA SPECTRUM
    !          *SUNASGN*- SPECTRA OF ALL UNASSIGNED SPECTRAL COMPONENTS
    !          *EME*    - MEAN WAVE ENERGY FOR EACH WAVE SYSTEM
    !          *DIR*    - MEAN WAVE DIRECTION FOR EACH WAVE SYSTEM
    !          *PER*    - MEAN PERIOD BASED ON -1 MOMENT FOR EACH WAVE SYSTEM.


    !     METHOD.
    !     -------

    !       STEEPEST ASCENT

    !     EXTERNALS.
    !     ----------
    !       *PARMEAN* - COMPUTE MEAN PARAMETERS.

    !       NONE.

    ! ----------------------------------------------------------------------

    USE PARMEAN_LOKI_MOD, ONLY: PARMEAN_LOKI
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWCOUT, ONLY: LLPARTITION

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL, NPMAX
    INTEGER(KIND=JWIM), INTENT(INOUT), DIMENSION(KIJL) :: NPEAK
    INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(KIJL) :: MIJ
    INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(KIJL, NPMAX) :: NTHP, NFRP

    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FLLOW
    LOGICAL, INTENT(IN), DIMENSION(KIJL, NANG) :: LLCOSDIFF
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: FLNOISE
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: SWM
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: SUNASGN
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, 0:NPMAX) :: DIR, PER, ENE


    INTEGER(KIND=JWIM) :: ITHC, IFRC
    INTEGER(KIND=JWIM) :: M, K, IP, NITT
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: NANGH, KK, KKMIN, KKMAX
    INTEGER(KIND=JWIM) :: IFRL, ITHL, ITHR
    INTEGER(KIND=JWIM) :: MMIN, MMAX
    INTEGER(KIND=JWIM), INTENT(INOUT), DIMENSION(3*NANG) :: KLOC

    REAL(KIND=JWRB) :: HALF_SECTOR
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NANG, NFRE) :: W2
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: W1
    ! FNDPRT is reached only when partition output is enabled. Loki hoists
    ! this workspace into its output driver even for the inactive branch.
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NANG, NFRE, MERGE(NPMAX, 1, LLPARTITION), KIJL) :: SPEC

    LOGICAL :: LLCHANGE, LLADD
    LOGICAL :: LLADDPART
    LOGICAL, INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: LLW3
!$acc routine seq
    REAL(KIND=JWRB), INTENT(INOUT) :: PARMEAN_F1D(NFRE)
    REAL(KIND=JWRB), INTENT(INOUT) :: PARMEAN_TEMP(NPMAX)
    REAL(KIND=JWRB), INTENT(INOUT) :: PARMEAN_SI(NPMAX)
    REAL(KIND=JWRB), INTENT(INOUT) :: PARMEAN_CI(NPMAX)
    REAL(KIND=JWRB), INTENT(INOUT) :: PARMEAN_EM(NPMAX)
    REAL(KIND=JWRB), INTENT(INOUT) :: PARMEAN_FM(NPMAX)
    REAL(KIND=JWRB), INTENT(INOUT) :: PARMEAN_THQ(NPMAX)
!$acc data present( NPEAK, MIJ, NTHP, NFRP, FLLOW, LLCOSDIFF, FLNOISE, FL1, SWM, SUNASGN, ENE, DIR, PER, KLOC, W2, W1, SPEC,  &
!$acc & LLW3, PARMEAN_F1D, PARMEAN_TEMP, PARMEAN_SI, PARMEAN_CI, PARMEAN_EM, PARMEAN_FM, PARMEAN_THQ )

    ! ----------------------------------------------------------------------


    !     ASSUME THAT ONE HAS ONLY SEARCH A SECTOR IN DIRECTION (not the full 360 degrees)
    !     AROUND THE PEAK
    !!!   THIS IS A STRONG ASSUMPTION ON THE DIRECTIONALITY OF SWELL !!!

    HALF_SECTOR = 75.0_JWRB
    NANGH = NINT((HALF_SECTOR / 360.0_JWRB)*NANG) + 1

    SUNASGN(:, :, :) = FL1(:, :, :)

    !     POINTS BELOW THE MINIMUM LEVEL BELONG TO THE WINDSEA PART
    !     AND CAN BE EXCLUDED FROM THE PARTITIONS OF THE SWELL SPECTRUM
    MMIN = NFRE
    MMAX = 0
!$acc loop seq
    DO M=1,NFRE
!$acc loop seq
      DO K=1,NANG
        IF (FL1(IJ, K, M) <= FLLOW(IJ, K, M)) THEN
          W1(IJ, K, M) = 1._JWRB
          LLW3(IJ, K, M) = .false.
        ELSE
          W1(IJ, K, M) = 0._JWRB
          LLW3(IJ, K, M) = .true.
          MMIN = MIN(M, MMIN)
          MMAX = M
        END IF
      END DO
    END DO


    !     MAIN LOOP ON GRID POINTS
!$acc loop seq
    DO IP=1,NPEAK(IJ)
      ITHC = NTHP(IJ, IP)
      IFRC = NFRP(IJ, IP)

      KKMIN = ITHC - NANGH
      KKMAX = ITHC + NANGH
!$acc loop seq
      DO KK=KKMIN,KKMAX
        KLOC(KK + NANG) = 1 + MOD(NANG + KK - 1, NANG)
      END DO

      !*        1.  SET UP THE W2 MAP
      !             -----------------

!$acc loop seq
      DO M=1,NFRE
!$acc loop seq
        DO K=1,NANG
          W2(K, M) = 0.0_JWRB
        END DO
      END DO

!$acc loop seq
      DO K=ITHC - 1,ITHC + 1
        ITHL = 1 + MOD(NANG + K - 1, NANG)
!$acc loop seq
        DO M=IFRC - 1,IFRC + 1
          IFRL = MAX(1, MIN(NFRE, M))
          IF (W1(IJ, ITHL, IFRL) <= 0.5_JWRB) W2(ITHL, IFRL) = 0.5_JWRB
        END DO
      END DO

      IF (W1(IJ, ITHC, IFRC) < 0.25_JWRB) W2(ITHC, IFRC) = 1.0_JWRB

      !         FIND IF MORE HIGH FREQUENCY BINS HAVE BECOME EXCLUDED
!$acc loop seq
      OUT0: DO M=MMAX,MMIN,-1
!$acc loop seq
        DO KK=KKMIN,KKMAX
          K = KLOC(KK + NANG)
          IF (W1(IJ, K, M) < 1.0_JWRB) THEN
            MMAX = M
            EXIT OUT0
          END IF
        END DO
      END DO OUT0

      !*        2.  ITTERATE SEARCH
      !             ---------------

      NITT = 0

      !     2.a Branch point

200   CONTINUE
      NITT = NITT + 1
      LLCHANGE = .false.

      !         2.b Determine central points

!$acc loop seq
      DO M=MMIN,MIN(MIJ(IJ), MMAX)
        !           by definition bins beyond M=MIJ are never extremas
        !           and bins above MMAX are excluded.
!$acc loop seq
        DO KK=KKMIN,KKMAX
          K = KLOC(KK + NANG)

          IF (LLW3(IJ, K, M)) THEN
            IF (W2(K, M) == 0.5_JWRB .and. W1(IJ, K, M) < 0.5_JWRB) THEN
              LLADD = .true.
!$acc loop seq
              OUT1: DO ITHR=K - 1,K + 1
                ITHL = 1 + MOD(NANG + ITHR - 1, NANG)
!$acc loop seq
                DO IFRL=MAX(1, M - 1),MIN(NFRE, M + 1)
                  IF (W2(ITHL, IFRL) == 0.0_JWRB .and. FL1(IJ, ITHL, IFRL) > FL1(IJ, K, M)) THEN
                    LLADD = .false.
                    EXIT OUT1
                  END IF
                END DO
              END DO OUT1
              IF (LLADD) THEN
                W2(K, M) = 1.0_JWRB
                LLCHANGE = .true.
              END IF
            END IF
          END IF

        END DO
      END DO

      !         2.c Determine peripherical points

!$acc loop seq
      DO M=MMIN,MMAX
!$acc loop seq
        DO KK=KKMIN,KKMAX
          K = KLOC(KK + NANG)

          IF (LLW3(IJ, K, M) .and. W1(IJ, K, M) < 1.0_JWRB) THEN
            IF (W2(K, M) == 0.0_JWRB) THEN
!$acc loop seq
              OUT2: DO ITHR=K - 1,K + 1
                ITHL = 1 + MOD(NANG + ITHR - 1, NANG)
!$acc loop seq
                DO IFRL=MAX(1, M - 1),MIN(NFRE, M + 1)
                  IF (W2(ITHL, IFRL) == 1.0_JWRB) THEN
                    W2(K, M) = 0.5_JWRB
                    LLCHANGE = .true.
                    EXIT OUT2
                  END IF
                END DO
              END DO OUT2
            END IF
          END IF

        END DO
      END DO

      !*        2.d Branch back ?

      IF (LLCHANGE .and. NITT < 25) GO TO 200


      !*        3   UPDATE THE OVERALL MAP

!$acc loop seq
      DO M=1,NFRE
!$acc loop seq
        DO K=1,NANG
          W1(IJ, K, M) = W1(IJ, K, M) + W2(K, M)
        END DO
      END DO

!$acc loop seq
      DO M=1,NFRE
!$acc loop seq
        DO K=1,NANG
          SPEC(K, M, IP, IJ) = FL1(IJ, K, M)*W2(K, M)
        END DO
      END DO

    END DO
    ! loop IP

    !       CHECK THAT UNASSIGNED BINS ARE IN THE WIND SECTOR, OTHERWISE THEY SHOULD BE
    !       ADDED AS A PARTITION IF ABOVE NOISE LEVEL
    LLADDPART = .false.
    IF (NPEAK(IJ) < NPMAX) THEN
!$acc loop seq
      DO K=1,NANG
!$acc loop seq
        DO M=1,NFRE
          IF (LLCOSDIFF(IJ, K) .and. W1(IJ, K, M) <= 0.0_JWRB .and. FL1(IJ, K, M) > FLNOISE(IJ)) THEN
            W2(K, M) = 1.0_JWRB
            W1(IJ, K, M) = 1.0_JWRB
            LLADDPART = .true.
          ELSE
            W2(K, M) = 0.0_JWRB
          END IF
        END DO
      END DO
    END IF

    IF (LLADDPART) THEN
      NPEAK(IJ) = NPEAK(IJ) + 1
      IP = NPEAK(IJ)
!$acc loop seq
      DO M=1,NFRE
!$acc loop seq
        DO K=1,NANG
          SPEC(K, M, IP, IJ) = FL1(IJ, K, M)*W2(K, M)
        END DO
      END DO
    END IF

!$acc loop seq
    DO IP=1,NPEAK(IJ)
!$acc loop seq
      DO M=1,NFRE
!$acc loop seq
        DO K=1,NANG
          IF (SPEC(K, M, IP, IJ) > 0.0_JWRB) SUNASGN(IJ, K, M) = 0.0_JWRB
        END DO
      END DO
    END DO

    ! loop IJ

    CALL PARMEAN_LOKI(KIJS, KIJL, NPMAX, NPEAK, SPEC, ENE, DIR, PER, IJ=IJ, F1D=PARMEAN_F1D, TEMP=PARMEAN_TEMP, SI=PARMEAN_SI,  &
    & CI=PARMEAN_CI, EM=PARMEAN_EM, FM=PARMEAN_FM, THQ=PARMEAN_THQ)


    !     UPDATE SWELL MASK (IF IN THE WIND DIRECTION)
!$acc loop seq
    DO M=1,NFRE
!$acc loop seq
      DO K=1,NANG
        IF (.not.LLCOSDIFF(IJ, K)) THEN
          SWM(IJ, K, M) = SWM(IJ, K, M)*MAX(W1(IJ, K, M), 1.0_JWRB)
        END IF
      END DO
    END DO



!$acc end data
  END SUBROUTINE FNDPRT_LOKI
END MODULE FNDPRT_LOKI_MOD
