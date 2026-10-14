! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SINPUT_ARD_LOKI_MOD
  CONTAINS
  SUBROUTINE SINPUT_ARD_LOKI (NGST, LLSNEG, KIJS, KIJL, FL1, WAVNUM, CINV, XK2CG, WDWAVE, WSWAVE, UFRIC, Z0M, COSWDIF, SINWDIF2,  &
  & RAORW, WSTAR, RNFAC, FLD, SL, SPOS, XLLWS, IJ, GAM0, DSTAB, COSLP)
    ! ----------------------------------------------------------------------

    !**** *SINPUT_ARD* - COMPUTATION OF INPUT SOURCE FUNCTION.


    !*    PURPOSE.
    !     ---------

    !       COMPUTE THE WIND INPUT SOURCE TRERM BASED ON ARDHUIN ET AL. 2010.

    !       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
    !       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
    !       INPUT SOURCE FUNCTION.
    !
    !       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
    !       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
    !       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
    !       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
    !       FINDS:
    !
    !             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
    !
    !       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
    !       LEVEL.

    !**   INTERFACE.
    !     ----------

    !     *CALL* *SINPUT_ARD (NGST, LLSNEG, KIJS, KIJL, FL1,
    !    &                    WAVNUM, CINV, XK2CG,
    !    &                    WSWAVE, WDWAVE, UFRIC, Z0M,
    !    &                    COSWDIF, SINWDIF2,
    !    &                    RAORW, WSTAR, RNFAC,
    !    &                    FLD, SL, SPOS, XLLWS)
    !         *NGST* - IF = 1 THEN NO GUSTINESS PARAMETERISATION
    !                - IF = 2 THEN GUSTINESS PARAMETERISATION
    !         *LLSNEG- IF TRUE THEN THE NEGATIVE SINPUT (SWELL DAMPING) WILL BE COMPUTED
    !         *KIJS* - INDEX OF FIRST GRIDPOINT.
    !         *KIJL* - INDEX OF LAST GRIDPOINT.
    !          *FL1* - SPECTRUM.
    !       *WAVNUM* - WAVE NUMBER.
    !         *CINV* - INVERSE PHASE VELOCITY.
    !       *XK2CG*  - (WAVNUM)**2 * GROUP SPPED.
    !       *WDWAVE* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
    !                  NOTATION (POINTING ANGLE OF WIND VECTOR,
    !                  CLOCKWISE FROM NORTH).
    !        *UFRIC* - NEW FRICTION VELOCITY IN M/S.
    !        *Z0M* - ROUGHNESS LENGTH IN M.
    !      *COSWDIF* - COS(TH(K)-WDWAVE(IJ))
    !     *SINWDIF2* - SIN(TH(K)-WDWAVE(IJ))**2
    !        *RAORW* - RATIO AIR DENSITY TO WATER DENSITY.
    !        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
    !        *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !          *FLD* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
    !           *SL* - TOTAL SOURCE FUNCTION ARRAY.
    !         *SPOS* - POSITIVE SOURCE FUNCTION ARRAY.
    !       *XLLWS*  - = 1 WHERE SINPUT IS POSITIVE

    !     METHOD.
    !     -------

    !       SEE REFERENCE.


    ! ----------------------------------------------------------------------
    USE YOWWIND, ONLY: WSPMIN
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWCOUP, ONLY: LLCAPCHNK, LLNORMAGAM, LLGCBZ0
    USE YOWFRED, ONLY: FR, TH, DFIM, COSTH, SINTH, ZPIFR, DELTH
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: G, GM1, EPSMIN, EPSUS, ZPI, ACDLIN, BCDLIN
    USE YOWPHYS, ONLY: ZALP, TAUWSHELTER, XKAPPA, BETAMAXOXKAPPA2, RNU, RNUM, SWELLF, SWELLF2, SWELLF3, SWELLF4, SWELLF5,  &
    & SWELLF6, SWELLF7, SWELLF7M1, Z0RAT, Z0TUBMAX, ABMIN, ABMAX, ALPHAMIN, ALPHAMAX
    USE YOWTEST, ONLY: IU06
    USE YOWTABL, ONLY: IAB, SWELLFT

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: NGST
    LOGICAL, INTENT(IN) :: LLSNEG
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM, CINV, XK2CG
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: WDWAVE, WSWAVE, UFRIC, Z0M
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: RAORW, WSTAR, RNFAC
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF, SINWDIF2
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL, SPOS
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: XLLWS


    INTEGER(KIND=JWIM) :: K, M, IND, IGST
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ

    REAL(KIND=JWRB) :: CONSTN
    REAL(KIND=JWRB) :: AVG_GST, ABS_TAUWSHELTER
    REAL(KIND=JWRB) :: CONST1
    REAL(KIND=JWRB) :: ZNZ
    REAL(KIND=JWRB) :: X1, X2, ZLOG, ZLOG1, ZLOG2, ZLOG2X, XV1, XV2, ZBETA1, ZBETA2
    REAL(KIND=JWRB) :: XI, X, DELI1, DELI2
    REAL(KIND=JWRB) :: FU, FUD, NU_AIR, SMOOTH, HFTSWELLF6, Z0TUB
    REAL(KIND=JWRB) :: FAC_NU_AIR, FACM1_NU_AIR
    REAL(KIND=JWRB) :: ARG, DELABM1
    REAL(KIND=JWRB) :: TAUPX, TAUPY
    REAL(KIND=JWRB) :: DSTAB2
    REAL(KIND=JPHOOK) :: SINPUT_ARD_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: CONST, SIG, SIG2, COEF, COEF5, DFIM_SIG2
    REAL(KIND=JWRB) :: CONSTF, CONST11, CONST22
    REAL(KIND=JWRB) :: Z0VIS, Z0NOZ, FWW
    REAL(KIND=JWRB) :: PVISC, PTURB
    REAL(KIND=JWRB) :: ZCN
    REAL(KIND=JWRB) :: SIG_N, UORBT, AORB, TEMP, RE, RE_C, ZORB
    REAL(KIND=JWRB) :: CNSN, SUMF, SUMFSIN2
    REAL(KIND=JWRB) :: CSTRNFAC
    REAL(KIND=JWRB) :: FLP_AVG, SLP_AVG
    REAL(KIND=JWRB) :: ROGOROAIR, AIRD_PVISC, USG2, FLP, SLP
    REAL(KIND=JWRB), DIMENSION(2) :: XSTRESS, YSTRESS
    REAL(KIND=JWRB), DIMENSION(2) :: TAUX, TAUY, USTP, USTPM1, USDIRP, UCN
    REAL(KIND=JWRB), DIMENSION(2) :: UCNZALPD
    REAL(KIND=JWRB) :: XNGAMCONST
    REAL(KIND=JWRB), DIMENSION(2) :: GAMNORMA    ! ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB) :: DSTAB1, TEMP1, TEMP2
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, 2) :: GAM0
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, 2) :: DSTAB
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG) :: COSLP

    LOGICAL :: LTAUWSHELTER
    REAL(KIND=JWRB), PARAMETER :: BG_GUST = 0.0_JWRB    ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
    REAL(KIND=JWRB), PARAMETER :: ONETHIRD = 1.0_JWRB / 3.0_JWRB
    REAL(KIND=JWRB), PARAMETER :: SIG_NMAX = 0.9_JWRB    ! MAX OF RELATIVE STANDARD DEVIATION OF USTAR
    REAL(KIND=JWRB), PARAMETER :: C1 = 1.03E-3_JWRB
    REAL(KIND=JWRB), PARAMETER :: C2 = 0.04E-3_JWRB
    REAL(KIND=JWRB), PARAMETER :: P1 = 1.48_JWRB
    REAL(KIND=JWRB), PARAMETER :: P2 = -0.21_JWRB
    REAL(KIND=JWRB) :: ZCHAR, C_D, DC_DDU, SIG_CONV
    REAL(KIND=JWRB) :: XKAPPAD, U10, C2U10P1, U10P2
    REAL(KIND=JWRB) :: BCD_LOC, U10M1, ZN, WSIGSTAR_Z0VIS
    REAL(KIND=JPHOOK) :: WSIGSTAR_ZHOOK_HANDLE
!$acc routine seq
!$acc data present( FL1, WAVNUM, CINV, XK2CG, WDWAVE, WSWAVE, UFRIC, Z0M, COSWDIF, SINWDIF2, RAORW, WSTAR, RNFAC, FLD, SL, SPOS,  &
!$acc & XLLWS, GAM0, DSTAB, COSLP )
    ! ----------------------------------------------------------------------


    AVG_GST = 1.0_JWRB / NGST
    CONST1 = BETAMAXOXKAPPA2
    CONSTN = DELTH / (XKAPPA*ZPI)

    ABS_TAUWSHELTER = ABS(TAUWSHELTER)
    IF (ABS_TAUWSHELTER == 0.0_JWRB) THEN
      LTAUWSHELTER = .false.
    ELSE
      LTAUWSHELTER = .true.
    END IF


    !     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
    IF (NGST > 1) THEN
      ! [Loki] inlined child subroutine: WSIGSTAR
      ! =========================================


      ! ----------------------------------------------------------------------


      IF (LLGCBZ0 .or. LLNORMAGAM) THEN
        ZN = RNUM

        U10M1 = 1.0_JWRB / MAX(WSWAVE(IJ), WSPMIN)
        ! CHARNOCK:
        WSIGSTAR_Z0VIS = ZN / MAX(UFRIC(IJ), EPSUS)
        ZCHAR = G*(Z0M(IJ) - WSIGSTAR_Z0VIS) / MAX(UFRIC(IJ)**2, EPSUS)
        ZCHAR = MAX(MIN(ZCHAR, ALPHAMAX), ALPHAMIN)

        BCD_LOC = BCDLIN*SQRT(ZCHAR)
        C_D = ACDLIN + BCD_LOC*WSWAVE(IJ)
        DC_DDU = BCD_LOC
        SIG_CONV = 1.0_JWRB + 0.5_JWRB*WSWAVE(IJ) / C_D*DC_DDU
        SIG_N = MIN(SIG_NMAX, SIG_CONV*U10M1*(BG_GUST*UFRIC(IJ)**3 + 0.5_JWRB*XKAPPA*WSTAR(IJ)**3)**ONETHIRD)

      ELSE
        ZN = 0.0_JWRB


        !!! for consistency I have kept the old method, even though the new method above could be used,
        !!! but until LLGCBZ0 or LLNORMAGAM is the default, keep the old scheme whe it is not...
        !
        !       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
        !       BASED ON U*
        !
        XKAPPAD = 1.0_JWRB / XKAPPA
        U10 = UFRIC(IJ)*XKAPPAD*(LOG(10.0_JWRB) - LOG(Z0M(IJ)))
        U10 = MAX(U10, WSPMIN)
        U10M1 = 1.0_JWRB / U10
        C2U10P1 = C2*U10**P1
        U10P2 = U10**P2
        C_D = (C1 + C2U10P1)*U10P2
        DC_DDU = (P2*C1 + (P1 + P2)*C2U10P1)*U10P2*U10M1
        SIG_CONV = 1.0_JWRB + 0.5_JWRB*U10 / C_D*DC_DDU
        SIG_N = MIN(SIG_NMAX, SIG_CONV*U10M1*(BG_GUST*UFRIC(IJ)**3 + 0.5_JWRB*XKAPPA*WSTAR(IJ)**3)**ONETHIRD)

      END IF


      ! =========================================
    END IF


    IF (LLNORMAGAM) THEN
      CSTRNFAC = CONSTN*RNFAC(IJ) / RAORW(IJ)
    END IF
    ! ----------------------------------------------------------------------

    IF (LLSNEG) THEN
      !!!!  only for the negative sinput
      NU_AIR = RNU
      FACM1_NU_AIR = 4.0_JWRB / NU_AIR

      FAC_NU_AIR = RNUM

      FU = ABS(SWELLF3)
      FUD = SWELLF2
      DELABM1 = REAL(IAB) / (ABMAX - ABMIN)


      !       computation of Uorb and Aorb
      UORBT = EPSMIN
      AORB = EPSMIN

!$acc loop seq
      DO M=1,NFRE
        SIG = ZPIFR(M)
        SIG2 = SIG**2
        DFIM_SIG2 = DFIM(M)*SIG2

        K = 1
        TEMP = FL1(IJ, K, M)
!$acc loop seq
        DO K=2,NANG
          TEMP = TEMP + FL1(IJ, K, M)
        END DO

        UORBT = UORBT + DFIM_SIG2*TEMP
        AORB = AORB + DFIM(M)*TEMP
      END DO

      UORBT = 2.0_JWRB*SQRT(UORBT)        ! this is the significant orbital amplitude
      AORB = 2.0_JWRB*SQRT(AORB)        ! this 1/2 Hs
      RE = FACM1_NU_AIR*UORBT*AORB        ! this is the Reynolds number
      Z0VIS = FAC_NU_AIR / MAX(UFRIC(IJ), 0.0001_JWRB)
      Z0TUB = Z0RAT*MIN(Z0TUBMAX, Z0M(IJ))
      Z0NOZ = MAX(Z0VIS, Z0TUB)
      ZORB = AORB / Z0NOZ

      !         compute fww
      XI = (LOG10(MAX(ZORB, 3.0_JWRB)) - ABMIN)*DELABM1
      IND = MIN(IAB - 1, INT(XI))
      DELI1 = MIN(1.0_JWRB, XI - REAL(IND, kind=JWRB))
      DELI2 = 1.0_JWRB - DELI1
      FWW = SWELLFT(IND)*DELI2 + SWELLFT(IND + 1)*DELI1
      TEMP2 = FWW*UORBT

      !       Define the critical Reynolds number
      IF (SWELLF6 == 1.0_JWRB) THEN
        RE_C = SWELLF4
      ELSE
        HFTSWELLF6 = 1.0_JWRB - SWELLF6
        RE_C = SWELLF4*(2.0_JWRB / AORB)**HFTSWELLF6
      END IF

      !       Swell damping weight between viscous and turbulent boundary layer
      IF (SWELLF7 > 0.0_JWRB) THEN
        SMOOTH = 0.5_JWRB*TANH((RE - RE_C)*SWELLF7M1)
        PTURB = 0.5_JWRB + SMOOTH
        PVISC = 0.5_JWRB - SMOOTH
      ELSE
        IF (RE <= RE_C) THEN
          PTURB = 0.0_JWRB
          PVISC = 0.5_JWRB
        ELSE
          PTURB = 0.5_JWRB
          PVISC = 0.0_JWRB
        END IF
      END IF

      AIRD_PVISC = PVISC*RAORW(IJ)

    END IF



    ! Initialisation

    IF (NGST == 1) THEN
      USTP(1) = UFRIC(IJ)
    ELSE IF (NGST == 2) THEN
      USTP(1) = UFRIC(IJ)*(1.0_JWRB + SIG_N)
      USTP(2) = UFRIC(IJ)*(1.0_JWRB - SIG_N)
    END IF

    !... Expressing the IGST loop in this way enables the compiler to
    !... unroll it whilst still retaining correctness for the case
    !... where NGST == 1. This is an important optimisation for GPUs.
!$acc loop seq
    DO IGST=1,2
      IF (IGST <= NGST) THEN
        USTPM1(IGST) = 1.0_JWRB / MAX(USTP(IGST), EPSUS)
      END IF
    END DO

    IF (LTAUWSHELTER) THEN
!$acc loop seq
      DO IGST=1,2
        IF (IGST <= NGST) THEN
          XSTRESS(IGST) = 0.0_JWRB
          YSTRESS(IGST) = 0.0_JWRB
          USG2 = USTP(IGST)**2
          TAUX(IGST) = USG2*SIN(WDWAVE(IJ))
          TAUY(IGST) = USG2*COS(WDWAVE(IJ))
        END IF
      END DO

      ROGOROAIR = G / RAORW(IJ)

    ELSE
!$acc loop seq
      DO K=1,NANG
        COSLP(IJ, K) = COSWDIF(IJ, K)
      END DO
    END IF


    !*    2. MAIN LOOP OVER FREQUENCIES.
    !        ---------------------------

    IF (.not.LLNORMAGAM) THEN
      GAMNORMA(1) = 1.0_JWRB
      GAMNORMA(2) = 1.0_JWRB
    END IF

    IF (.not.LLSNEG) THEN
!$acc loop seq
      DO K=1,NANG
        DSTAB(IJ, K, 1) = 0.0_JWRB
        DSTAB(IJ, K, 2) = 0.0_JWRB
      END DO
    END IF

!$acc loop seq
    DO M=1,NFRE
      SIG = ZPIFR(M)
      SIG2 = SIG**2
      CONST = SIG*CONST1

      IF (LLSNEG) THEN
        COEF = -SWELLF*16._JWRB*SIG2 / G
        COEF5 = -SWELLF5*2._JWRB*SQRT(2._JWRB*NU_AIR*SIG)
      END IF

      IF (LTAUWSHELTER) THEN
!$acc loop seq
        DO IGST=1,2
          IF (IGST <= NGST) THEN
            TAUPX = TAUX(IGST) - ABS_TAUWSHELTER*XSTRESS(IGST)
            TAUPY = TAUY(IGST) - ABS_TAUWSHELTER*YSTRESS(IGST)
            USDIRP(IGST) = ATAN2(TAUPX, TAUPY)
            USTP(IGST) = (TAUPX**2 + TAUPY**2)**0.25_JWRB
            USTPM1(IGST) = 1.0_JWRB / MAX(USTP(IGST), EPSUS)
          END IF
        END DO

        CONSTF = ROGOROAIR*CINV(IJ, M)*DFIM(M)
      END IF


      !*      PRECALCULATE FREQUENCY DEPENDENCE.
      !       ----------------------------------

!$acc loop seq
      DO IGST=1,2
        IF (IGST <= NGST) THEN
          UCN(IGST) = USTP(IGST)*CINV(IJ, M)
          UCNZALPD(IGST) = XKAPPA / (UCN(IGST) + ZALP)
        END IF
      END DO
      ZCN = LOG(WAVNUM(IJ, M)*Z0M(IJ))
      CNSN = CONST*RAORW(IJ)

      !*    2.1 LOOP OVER DIRECTIONS.
      !         ---------------------

!$acc loop seq
      DO K=1,NANG
        XLLWS(IJ, K, M) = 0.0_JWRB
      END DO

      IF (LLNORMAGAM) THEN
        XNGAMCONST = CSTRNFAC*XK2CG(IJ, M)
      END IF

      IF (LLSNEG) THEN
        !       SWELL DAMPING:
        DSTAB1 = COEF5*AIRD_PVISC*WAVNUM(IJ, M)
        TEMP1 = COEF*RAORW(IJ)
      END IF

!$acc loop seq
      DO IGST=1,2
        IF (IGST <= NGST) THEN
!$acc loop seq
          DO K=1,NANG
            IF (LTAUWSHELTER) THEN
              COSLP(IJ, K) = COS(TH(K) - USDIRP(IGST))
            END IF
            GAM0(IJ, K, IGST) = 0.0_JWRB
            IF (COSLP(IJ, K) > 0.01_JWRB) THEN
              X = COSLP(IJ, K)*UCN(IGST)
              ZLOG = ZCN + UCNZALPD(IGST) / COSLP(IJ, K)
              IF (ZLOG < 0.0_JWRB) THEN
                ZLOG2X = ZLOG*ZLOG*X
                GAM0(IJ, K, IGST) = EXP(ZLOG)*ZLOG2X*ZLOG2X*CNSN
                XLLWS(IJ, K, M) = 1.0_JWRB
              END IF
            END IF
          END DO

          IF (LLNORMAGAM) THEN

            SUMF = 0.0_JWRB
            SUMFSIN2 = 0.0_JWRB
!$acc loop seq
            DO K=1,NANG
              SUMF = SUMF + GAM0(IJ, K, IGST)*FL1(IJ, K, M)
              SUMFSIN2 = SUMFSIN2 + GAM0(IJ, K, IGST)*FL1(IJ, K, M)*SINWDIF2(IJ, K)
            END DO

            ZNZ = XNGAMCONST*USTPM1(IGST)
            GAMNORMA(IGST) = (1.0_JWRB + ZNZ*SUMFSIN2) / (1.0_JWRB + ZNZ*SUMF)

          END IF

          IF (LLSNEG) THEN
!$acc loop seq
            DO K=1,NANG
              DSTAB2 = TEMP1*(TEMP2 + (FU + FUD*COSLP(IJ, K))*USTP(IGST))
              DSTAB(IJ, K, IGST) = DSTAB1 + PTURB*DSTAB2
            END DO
          END IF
        END IF
      END DO


      !*    2.2 UPDATE THE SHELTERING STRESS (in any),
      !         AND THEN ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
      !         ---------------------------------------------------------

!$acc loop seq
      DO K=1,NANG

!$acc loop seq
        DO IGST=1,2
          IF (IGST <= NGST) THEN
            ! SLP: only the positive contributions
            SLP = GAM0(IJ, K, IGST)*GAMNORMA(IGST)
            FLP = SLP + DSTAB(IJ, K, IGST)

            SLP = SLP*FL1(IJ, K, M)

            IF (LTAUWSHELTER) THEN
              CONST11 = CONSTF*SINTH(K)
              CONST22 = CONSTF*COSTH(K)
              XSTRESS(IGST) = XSTRESS(IGST) + SLP*CONST11
              YSTRESS(IGST) = YSTRESS(IGST) + SLP*CONST22
            END IF

            IF (IGST == 1) THEN
              SLP_AVG = SLP
              FLP_AVG = FLP
            ELSE
              SLP_AVG = SLP_AVG + SLP
              FLP_AVG = FLP_AVG + FLP
            END IF
          END IF
        END DO

        SPOS(IJ, K, M) = AVG_GST*SLP_AVG

        FLD(IJ, K, M) = AVG_GST*FLP_AVG
        SL(IJ, K, M) = FLD(IJ, K, M)*FL1(IJ, K, M)

      END DO

    END DO
    ! END LOOP OVER FREQUENCIES



!$acc end data
  END SUBROUTINE SINPUT_ARD_LOKI
END MODULE SINPUT_ARD_LOKI_MOD
