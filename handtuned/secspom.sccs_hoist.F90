! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!--------------------------------------------------------------------
!
MODULE SECSPOM_LOKI_MOD
  CONTAINS
  SUBROUTINE SECSPOM_LOKI (F1, F3, KIJS, KIJL, NFRE, NANG, NDEPTH, DEPTHA, DEPTHD, OMSTART, FRAC, MR, OMEGA, DEPTH, AKMEAN, TA,  &
  & TB, TC_QL, TT_4M, TT_4P, IM_P, IM_M, IJ, LL, OMEGA_EXT, F2)
    !
    !--------------------------------------------------------------------
    !
    !*****SECSPOM** COMPUTES SECOND ORDER SPECTRUM IN FREQUENCY SPACE.
    !
    !     P.JANSSEN JULY 2008
    !
    !     PURPOSE
    !     -------
    !             DETERMINES SECOND-ORDER SPECTRUM, BASED ON JANSSEN (2008)
    !             THERE ARE THREE CORRECTIONS:
    !                   1) GENERATION OF SECOND-HARMONICS
    !                   2) QUASI-LINEAR EFFECT
    !                   3) SHIFT OF SPECTRUM BECAUSE OF STOKES FREQUENCY
    !                      CORRECTION.
    !
    !     INTERFACE
    !     ---------
    !             *CALL* *SECSPOM(F1,F3,KIJS,KIJL,NFRE,NANG,NMAX,NDEPTH,
    !                             DEPTHA,DEPTHD,FRAC,MR,DFDTH,OMEGA,
    !                             DEPTH,AKMEAN,TA,TB,TC_QL,
    !                             TT_4M,TT_4P,IM_P,IM_M)*
    !
    !
    !     PARAMETER   TYPE      PURPOSE.
    !     ---------   ----      -------
    !
    !       F1        REAL      2D FREE WAVE SPECTRUM (INPUT)
    !       F3        REAL      BOUND WAVES SPECTRUM (OUTPUT)
    !       KIJS      INTEGER   FIRST GRID POINT
    !       KIJL      INTEGER   LAST GRID POINT
    !       NFRE      INTEGER   NUMBER OF FREQUENCIES
    !       NANG      INTEGER   NUMBER OF DIRECTIONS
    !       NMAX      INTEGER   MAXIMUM INDEX CORRESPONDS TO TWICE THE CUT-OFF
    !                           FREQUENCY
    !       NDEPTH    INTEGER   NUMBER OF ENTRIES IN DEPTH TABLE
    !       OMSTART   REAL      STARTING FREQUENCY OF ORIGINAL FREQUENCY GRID
    !       FRAC      REAL      FRACTIONAL INCREASE IN FREQUENCY SPACE
    !       MR        INTEGER   THINNING FACTOR IN FREQUENCY SPACE
    !       DFDTH     REAL      PRODUCT OF INCREMENT IN FREQUENCY AND DIRECTION
    !       OMEGA     REAL      ANGULAR FREQUENCY ARRAY
    !       DEPTH     REAL      DEPTH ARRAY
    !       AKMEAN    REAL      MEAN WAVENUMBER
    !       TA        REAL      TABLE FOR MINUS INTERACTIONS
    !       TB        REAL      TABLE FOR PLUS INTERACTIONS
    !       TC_QL     REAL      TABLE FOR QUASI-LINEAR INTERACTIONS
    !       TT_4M     REAL      TABLE FOR STOKES FREQUENCY CORRECTION
    !       TT_4P     REAL      TABLE FOR STOKES FREQUENCY CORRECTION
    !       IM_P      INTEGER   TABLE FOR WAVENUMBER M2 PLUS
    !       IM_M      INTEGER   TABLE FOR WAVENUMBER M2 MIN
    !
    !
    !
    !     METHOD
    !     ------
    !             EVALUATE SECOND ORDER SPECTRUM IN FREQUENCY BASED ON
    !             KRASITSKII'S CANONICAL TRANSFORMATION.
    !
    !     EXTERNALS
    !     ---------
    !             NONE
    !
    !     REFERENCES
    !     ----------
    !             V.E. ZAKHAROV, HAMILTONIAN APPROACH (1968)
    !             M.A. SROKOSZ, J.G.R.,91,995-1006 (1986)
    !             P.A.E.M. JANSSEN, JFM (2009)
    !
    !
    !--------------------------------------------------------------------

    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    USE YOWTABL, ONLY: NMAX

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL, NFRE, NANG, NDEPTH, MR
    INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(NFRE, NFRE) :: IM_P, IM_M

    REAL(KIND=JWRB), INTENT(IN) :: DEPTHA, DEPTHD, OMSTART, FRAC
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(NFRE) :: OMEGA
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: DEPTH, AKMEAN
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: F1
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: F3
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(NDEPTH, NANG, NFRE, NFRE) :: TA, TB, TC_QL, TT_4M, TT_4P

    INTEGER(KIND=JWIM) :: M, K, M1, K1, M2_M, M2_P, K2, MP, MM, L, ID
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    INTEGER(KIND=JWIM), INTENT(INOUT), DIMENSION(NANG, NANG) :: LL
    INTEGER(KIND=JWIM) :: JD

    REAL(KIND=JWRB) :: OM0, OM0P, OM0M, OM0H, OM1
    REAL(KIND=JWRB) :: T_4M, T_4P, DELM1, XD, X_MIN, OMRT, XLOGD, OMG5
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NMAX) :: OMEGA_EXT
    REAL(KIND=JWRB) :: XINCR, DF2KP, DF2KM, PSUM
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NMAX) :: F2

    LOGICAL :: LLSAMEDPTH
!$acc routine seq
!$acc data present( F1, F3, OMEGA, DEPTH, AKMEAN, TA, TB, TC_QL, TT_4M, TT_4P, IM_P, IM_M, LL, OMEGA_EXT, F2 )

    ! ----------------------------------------------------------------------


    !
    !***  1. COMPUTATION OF TAIL OF THE SPECTRUM AND INDEX JD
    !     ---------------------------------------------------
    !
    !
    X_MIN = 1.0_JWRB
    XLOGD = LOG(DEPTHD)
    XD = MAX(X_MIN / AKMEAN(IJ), DEPTH(IJ))
    XD = LOG(XD / DEPTHA) / XLOGD + 1.0_JWRB
    ID = NINT(XD)
    ID = MAX(ID, 1)
    JD = MIN(ID, NDEPTH)
    LLSAMEDPTH = .true.
    ID = JD
    LLSAMEDPTH = LLSAMEDPTH .and. JD == JD
    ! [Loki] Removed content of pragma-marked region!

!$acc loop seq
    DO M=1,NFRE
      OMEGA_EXT(M) = OMEGA(M)
!$acc loop seq
      DO K=1,NANG
        F2(IJ, K, M) = F1(IJ, K, M)
      END DO
    END DO

    OMG5 = OMEGA(NFRE)**5
!$acc loop seq
    DO M=NFRE + 1,NMAX
      OM0 = OMSTART*(1.0_JWRB + FRAC)**(MR*M - 1)
      OMEGA_EXT(M) = OM0
      OMRT = OMG5 / OM0**5
!$acc loop seq
      DO K=1,NANG
        F2(IJ, K, M) = OMRT*F1(IJ, K, NFRE)
      END DO
    END DO

!$acc loop seq
    DO K=1,NANG
!$acc loop seq
      DO K1=1,NANG
        L = K - K1
        IF (L > NANG) L = L - NANG
        IF (L < 1) L = L + NANG
        LL(K1, K) = L
      END DO
    END DO

    !
    !***  2. COMPUTATION OF THE 2nd ORDER FREQUENCY SPECTRUM.
    !     ---------------------------------------------------

!$acc loop seq
    DO M=1,NFRE
      OM0 = OMEGA(M)
      OM0H = 0.5_JWRB*OM0
      !!!         MP   = MIN(M+1,NFRE)
      MP = MIN(M + 1, NMAX)
      OM0P = OMEGA_EXT(MP)
      MM = MAX(M - 1, 1)
      OM0M = OMEGA(MM)
      DELM1 = 1. / (OM0P - OM0M)
!$acc loop seq
      DO K=1,NANG
        K2 = K
        DF2KP = F2(IJ, K, MP)*DELM1
        DF2KM = F2(IJ, K, MM)*DELM1
        PSUM = 0.0_JWRB
!$acc loop seq
        DO M1=1,NFRE
          OM1 = OMEGA(M1)
          M2_M = IM_M(M1, M)
          M2_P = IM_P(M1, M)
!$acc loop seq
          DO K1=1,NANG
            L = LL(K1, K)

            IF (LLSAMEDPTH) THEN
              !*                  2.1 OM0-OM1 CASE: SECOND HARMONICS
              !                   ----------------------------------
              !                   OM2 = OM0-OM1
              IF (ABS(OM1) < OM0H) THEN
                !DIR$ PREFERVECTOR
                PSUM = PSUM + TA(ID, L, M1, M)*(F2(IJ, K1, M1)*F2(IJ, K2, M2_M) + F2(IJ, K2, M1)*F2(IJ, K1, M2_M))
              END IF

              !*                  2.2 OM1+OM0 CASE: INFRA-GRAVITY WAVES
              !                   -------------------------------------
              !                   OM2 = OM1+OM0
              !DIR$ PREFERVECTOR
              XINCR = 2.0_JWRB*TB(ID, L, M1, M)*F2(IJ, K2, M2_P)
              !
              !*                  2.3 QUASI-LINEAR EFFECT
              !                   -----------------------
              !DIR$ PREFERVECTOR
              XINCR = XINCR + TC_QL(ID, L, M1, M)*F2(IJ, K, M)
              !
              !*                   2.4 STOKES-FREQUENCY CORRECTION
              !                    -------------------------------
              !DIR$ PREFERVECTOR
              T_4M = TT_4M(ID, L, M1, M)
              T_4P = TT_4P(ID, L, M1, M)
              XINCR = XINCR - (DF2KP*T_4P - DF2KM*T_4M)

            ELSE

              !*                  2.1 OM0-OM1 CASE: SECOND HARMONICS
              !                   ----------------------------------
              !                   OM2 = OM0-OM1
              IF (ABS(OM1) < OM0H) THEN
                !DIR$ PREFERVECTOR
                PSUM = PSUM + TA(JD, L, M1, M)*(F2(IJ, K1, M1)*F2(IJ, K2, M2_M) + F2(IJ, K2, M1)*F2(IJ, K1, M2_M))
              END IF

              !*                  2.2 OM1+OM0 CASE: INFRA-GRAVITY WAVES
              !                   -------------------------------------
              !                   OM2 = OM1+OM0
              !DIR$ PREFERVECTOR
              XINCR = 2._JWRB*TB(JD, L, M1, M)*F2(IJ, K2, M2_P)
              !
              !*                  2.3 QUASI-LINEAR EFFECT
              !                   -----------------------
              !DIR$ PREFERVECTOR
              XINCR = XINCR + TC_QL(JD, L, M1, M)*F2(IJ, K, M)
              !
              !*                   2.4 STOKES-FREQUENCY CORRECTION
              !                    -------------------------------
              !DIR$ PREFERVECTOR
              T_4M = TT_4M(JD, L, M1, M)
              T_4P = TT_4P(JD, L, M1, M)
              XINCR = XINCR - (DF2KP*T_4P - DF2KM*T_4M)
            END IF

            PSUM = PSUM + F2(IJ, K1, M1)*XINCR

          END DO
        END DO
        ! M1
        F3(IJ, K, M) = PSUM
      END DO
      !K
    END DO
    !M

    !--------------------------------------------------------------------



!$acc end data
  END SUBROUTINE SECSPOM_LOKI
END MODULE SECSPOM_LOKI_MOD
