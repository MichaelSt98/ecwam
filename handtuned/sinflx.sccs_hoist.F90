! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SINFLX_LOKI_MOD
  CONTAINS
  SUBROUTINE SINFLX_LOKI (ICALL, NCALL, KIJS, KIJL, LUPDTUS, FL1, WAVNUM, CINV, XK2CG, WSWAVE, WDWAVE, AIRD, RAORW, WSTAR,  &
  & CICOVER, COSWDIF, SINWDIF2, FMEAN, HALP, FMEANWS, FLM, UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, PHIWA, FLD, SL, SPOS, MIJ,  &
  & RHOWGDFTH, XLLWS, IJ, WD, FLWD, FD, GAM0, DSTAB, COSLP)

    ! ----------------------------------------------------------------------

    !**** *SINFLX* - UPDATE STRESS AND COMPUTE WIND INPUT SOURCE TERM.

    ! ----------------------------------------------------------------------

    USE NS_GC_LOKI_MOD, ONLY: NS_GC_LOKI
    USE STRESS_GC_LOKI_MOD, ONLY: STRESS_GC_LOKI
    USE CHNKMIN_LOKI_MOD, ONLY: CHNKMIN_LOKI
    USE YOWICE, ONLY: CITHRSH_TAIL
    USE YOWSTAT, ONLY: IPHYS
    USE YOWTABL, ONLY: EPS1, IAB, SWELLFT
    USE YOWTEST, ONLY: IU06
    USE YOWWIND, ONLY: WSPMIN
    USE YOWFRED, ONLY: FR, TH, FR5, DELTH, ZPIFR, DFIM, DFIMOFR, WETAIL, FRTAIL, COSTH, SINTH, FRATIO, FLOGSPRDM1, RHOWG_DFIM,  &
    & FRIC, OMEGA_GC, XK_GC
    USE YOWPCONS, ONLY: G, ZPI, ZPI4GM2, GM1, EPSMIN, EPSUS, ACD, BCD, ACDLIN, BCDLIN, CDMAX, ZPI4GM1
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWCOUP, ONLY: LWCOU, LLCAPCHNK, LLGCBZ0, LLNORMAGAM, X0TAUHF, JTOT_TAUHF, WTAUHF
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPHYS, ONLY: DTHRN_A, DTHRN_U, ALPHAPMAX, XKAPPA, XNLEV, RNU, RNUM, ALPHA, ALPHAMIN, ALPHAMAX, ANG_GC_A, ANG_GC_B,  &
    & ANG_GC_C, ZALP, TAUWSHELTER, BETAMAXOXKAPPA2, SWELLF, SWELLF2, SWELLF3, SWELLF4, SWELLF5, SWELLF6, SWELLF7, SWELLF7M1,  &
    & Z0RAT, Z0TUBMAX, ABMIN, ABMAX, TAILFACTOR, TAILFACTOR_PM, GAMNCONST, SDSBR, ISDSDTH, ISB, IPSAT, SSDSC2, SSDSC4, SSDSC6,  &
    & MICHE, SSDSC3, SSDSBRF1, BRKPBCOEF, SSDSC5, NSDSNTH, INDICESSAT, SATWEIGHTS
    USE YOWWNDG, ONLY: ICODE, ICODE_CPL

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: ICALL    !! CALL NUMBER.
    INTEGER(KIND=JWIM), INTENT(IN) :: NCALL    !! TOTAL NUMBER OF CALLS.
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL    !! GRID POINT INDEXES.

    LOGICAL, INTENT(IN) :: LUPDTUS    !! IF TRUE UFRIC AND Z0M WILL BE UPDATED (CALLING AIRSEA).

    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FL1    !! WAVE SPECTRUM.
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM    !! WAVE NUMBER.
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: CINV    !! INVERSE PHASE VELOCITY.
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: XK2CG    !! (WAVNUM)**2 * GROUP SPPED.
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: WSWAVE    !! WIND SPEED IN M/S.
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: WDWAVE    !! WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC NOTATION.
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: AIRD    !! AIR DENSITY (KG/M**3).
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: RAORW    !! RATIO AIR DENSITY TO WATER DENSITY.
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: WSTAR    !! FREE CONVECTION VELOCITY SCALE (M/S)
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CICOVER    !! SEA ICE COVER.
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF    !! COS(TH(K)-WDWAVE(IJ))
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: SINWDIF2    !! SIN(TH(K)-WDWAVE(IJ))**2
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: FMEAN    !! MEAN FREQUENCY.
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: HALP    !! 1/2 PHILLIPS PARAMETER
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: FMEANWS    !! MEAN FREQUENCY OF THE WINDSEA.
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: FLM    !! SPECTAL DENSITY MINIMUM VALUE
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: UFRIC    !! FRICTION VELOCITY IN M/S.
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: TAUW    !! WAVE STRESS IN (M/S)**2
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: TAUWDIR    !! WAVE STRESS DIRECTION.
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: Z0M    !! ROUGHNESS LENGTH IN M.
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: Z0B    !! BACKGROUND ROUGHNESS LENGTH.
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: CHRNCK    !! CHARNOCK COEFFICIENT.

    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: PHIWA    !! ENERGY FLUX FROM WIND INTO WAVES.
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: FLD    !! DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: SL    !! TOTAL SOURCE FUNCTION ARRAY.
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: SPOS    !!  POSITIVE SINPUT ONLY.

    INTEGER(KIND=JWIM), INTENT(OUT) :: MIJ(KIJL)    !! LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.

    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NFRE) :: RHOWGDFTH    !! WATER DENSITY * G * DF * DTHETA

    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: XLLWS    !! TOTAL WINDSEA MASK FROM INPUT SOURCE TERM.


    INTEGER(KIND=JWIM) :: SINFLX_K
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: IUSFG, ICODE_WND
    INTEGER(KIND=JWIM) :: NGST

    REAL(KIND=JPHOOK) :: SINFLX_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: RNFAC

    LOGICAL :: LLPHIWA, LLSNEG
    INTEGER(KIND=JWIM) :: HALPHAP_K, SINFLX_M
    REAL(KIND=JWRB) :: ZLNFRNFRE
    REAL(KIND=JPHOOK) :: HALPHAP_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: ALPHAP
    REAL(KIND=JWRB) :: XMSS, SINFLX_EM, FM, F1D
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG) :: WD
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FLWD
    INTEGER(KIND=JWIM) :: MEANSQS_LF_M, MEANSQS_LF_K, KFRE
    REAL(KIND=JPHOOK) :: MEANSQS_LF_ZHOOK_HANDLE
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NFRE) :: FD
    REAL(KIND=JWRB) :: TEMP1, SINFLX_TEMP2
    INTEGER(KIND=JWIM) :: FEMEAN_M, FEMEAN_K
    REAL(KIND=JWRB) :: SINFLX_DELT25, DELT2, DEL2
    REAL(KIND=JPHOOK) :: FEMEAN_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: FEMEAN_TEMP2
    INTEGER(KIND=JWIM) :: I, J
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
    REAL(KIND=JPHOOK) :: SINPUT_ZHOOK_HANDLE
    INTEGER(KIND=JWIM) :: SINPUT_K, SINPUT_M, IND, IGST
    REAL(KIND=JWRB) :: CONSTN
    REAL(KIND=JWRB) :: AVG_GST, ABS_TAUWSHELTER
    REAL(KIND=JWRB) :: CONST1
    REAL(KIND=JWRB) :: ZNZ
    REAL(KIND=JWRB) :: X1, X2, ZLOG, ZLOG1, ZLOG2, ZLOG2X, XV1, XV2, ZBETA1, ZBETA2
    REAL(KIND=JWRB) :: SINPUT_XI, SINPUT_X, SINPUT_DELI1, SINPUT_DELI2
    REAL(KIND=JWRB) :: FU, FUD, NU_AIR, SMOOTH, HFTSWELLF6, Z0TUB
    REAL(KIND=JWRB) :: FAC_NU_AIR, FACM1_NU_AIR
    REAL(KIND=JWRB) :: SINPUT_ARG, DELABM1
    REAL(KIND=JWRB) :: TAUPX, TAUPY
    REAL(KIND=JWRB) :: DSTAB2
    REAL(KIND=JPHOOK) :: SINPUT_ARD_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: SINPUT_CONST, SIG, SIG2, COEF, COEF5, DFIM_SIG2
    REAL(KIND=JWRB) :: CONSTF, CONST11, CONST22
    REAL(KIND=JWRB) :: SINPUT_Z0VIS, Z0NOZ, FWW
    REAL(KIND=JWRB) :: PVISC, PTURB
    REAL(KIND=JWRB) :: ZCN
    REAL(KIND=JWRB) :: SIG_N, UORBT, AORB, SINFLX_TEMP, RE, RE_C, ZORB
    REAL(KIND=JWRB) :: CNSN, SUMF, SUMFSIN2
    REAL(KIND=JWRB) :: CSTRNFAC
    REAL(KIND=JWRB) :: FLP_AVG, SLP_AVG
    REAL(KIND=JWRB) :: ROGOROAIR, AIRD_PVISC, USG2, FLP, SLP
    REAL(KIND=JWRB), DIMENSION(2) :: XSTRESS, YSTRESS
    REAL(KIND=JWRB), DIMENSION(2) :: TAUX, TAUY, USTP, USTPM1, USDIRP, UCN
    REAL(KIND=JWRB), DIMENSION(2) :: UCNZALPD
    REAL(KIND=JWRB) :: XNGAMCONST
    REAL(KIND=JWRB), DIMENSION(2) :: GAMNORMA    ! ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB) :: DSTAB1, SINPUT_TEMP1, SINPUT_TEMP2
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
    REAL(KIND=JWRB) :: SINPUT_ZCHAR, C_D, DC_DDU, SIG_CONV
    REAL(KIND=JWRB) :: SINPUT_XKAPPAD, U10, C2U10P1, U10P2
    REAL(KIND=JWRB) :: BCD_LOC, U10M1, ZN, WSIGSTAR_Z0VIS
    REAL(KIND=JPHOOK) :: WSIGSTAR_ZHOOK_HANDLE
    INTEGER(KIND=JWIM) :: FEMEANWS_M, FEMEANWS_K
    REAL(KIND=JWRB) :: FEMEANWS_DELT25, FEMEANWS_DELT2
    REAL(KIND=JPHOOK) :: FEMEANWS_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: FEMEANWS_TEMP2, EM_LOC
    INTEGER(KIND=JWIM) :: FRCUTINDEX_M
    REAL(KIND=JWRB) :: FPMH, FPPM, FM2, FPM, FPM4
    REAL(KIND=JPHOOK) :: FRCUTINDEX_ZHOOK_HANDLE
    INTEGER(KIND=JWIM) :: STRESSO_M, STRESSO_K, STRESSO_I, STRESSO_J, II
    REAL(KIND=JWRB) :: TAUTOUS2
    REAL(KIND=JWRB) :: COSW, FCOSW2
    REAL(KIND=JPHOOK) :: STRESSO_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: STRESSO_XSTRESS, STRESSO_YSTRESS
    REAL(KIND=JWRB) :: TAUHF, PHIHF
    REAL(KIND=JWRB) :: CMRHOWGDFTH
    REAL(KIND=JWRB) :: US2, STRESSO_TAUX, STRESSO_TAUY, STRESSO_TAUPX, STRESSO_TAUPY
    REAL(KIND=JWRB) :: STRESSO_USDIRP, UST
    REAL(KIND=JWRB) :: SUMT, SUMX, SUMY
    LOGICAL :: STRESSO_LTAUWSHELTER
    INTEGER(KIND=JWIM) :: TAU_PHI_HF_J, TAU_PHI_HF_K
    INTEGER(KIND=JWIM) :: NS
    REAL(KIND=JWRB), PARAMETER :: ZSUPMAX = 0.0_JWRB    !  LOG(1.)
    REAL(KIND=JWRB) :: OMEGA, OMEGACC
    REAL(KIND=JWRB) :: X0G
    REAL(KIND=JWRB) :: YC, Y, CM1, ZX, ZARG, STRESSO_ZLOG, ZBETA
    REAL(KIND=JWRB) :: FNC, FNC2
    REAL(KIND=JWRB) :: STRESSO_GAMNORMA    ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB) :: STRESSO_ZNZ, CONFG
    REAL(KIND=JWRB) :: TAU_PHI_HF_COSW, TAU_PHI_HF_FCOSW2
    REAL(KIND=JPHOOK) :: TAU_PHI_HF_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: XKS, OMS
    REAL(KIND=JWRB) :: SQRTZ0OG, ZSUP, ZINF, DELZ
    REAL(KIND=JWRB) :: TAUL, XLOGGZ0, SQRTGZ0
    REAL(KIND=JWRB) :: USTPH
    REAL(KIND=JWRB) :: STRESSO_CONST1, CONST2, CONSTTAU, CONSTPHI
    REAL(KIND=JWRB) :: F1DCOS2, F1DCOS3
    REAL(KIND=JWRB) :: STRESSO_F1D, F1DSIN2
    REAL(KIND=JPHOOK) :: OMEGAGC_ZHOOK_HANDLE
!$acc routine seq
!$acc data present( FL1, WAVNUM, CINV, XK2CG, WSWAVE, WDWAVE, AIRD, RAORW, WSTAR, CICOVER, COSWDIF, SINWDIF2, FMEAN, HALP,  &
!$acc & FMEANWS, FLM, UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, PHIWA, FLD, SL, SPOS, MIJ, RHOWGDFTH, XLLWS, WD, FLWD, FD, GAM0,  &
!$acc & DSTAB, COSLP )

    ! ----------------------------------------------------------------------


    ! UPDATE UFRIC AND Z0M
    IF (ICALL == 1) THEN
      IUSFG = 0
      IF (LWCOU) THEN
        ICODE_WND = ICODE_CPL
      ELSE
        ICODE_WND = ICODE
      END IF
    ELSE
      IUSFG = 1
      ICODE_WND = 3
    END IF

    IF (LLNORMAGAM .and. LLCAPCHNK) THEN
      RNFAC = 1.0_JWRB + DTHRN_A*(1.0_JWRB + TANH(WSWAVE(IJ) - DTHRN_U))
    ELSE
      RNFAC = 1.0_JWRB
    END IF


    IF (LUPDTUS) THEN
      ! increase noise level in the tail
      IF (ICALL == 1) THEN
!$acc loop seq
        DO SINFLX_K=1,NANG
          FL1(IJ, SINFLX_K, NFRE) = MAX(FL1(IJ, SINFLX_K, NFRE), FLM(IJ, SINFLX_K))
        END DO

        IF (LLGCBZ0) THEN
          ! [Loki] inlined child subroutine: HALPHAP
          ! =========================================

          ! ----------------------------------------------------------------------


          ZLNFRNFRE = LOG(FR(NFRE))

          ! Find spectrum in wind direction
!$acc loop seq
          DO HALPHAP_K=1,NANG
            WD(IJ, HALPHAP_K) = 0.5_JWRB + 0.5_JWRB*SIGN(1.0_JWRB, COSWDIF(IJ, HALPHAP_K))
          END DO

!$acc loop seq
          DO SINFLX_M=1,NFRE
!$acc loop seq
            DO HALPHAP_K=1,NANG
              FLWD(IJ, HALPHAP_K, SINFLX_M) = FL1(IJ, HALPHAP_K, SINFLX_M)*WD(IJ, HALPHAP_K)
            END DO
          END DO

          ! [Loki] inlined child subroutine: MEANSQS_LF
          ! =========================================

          ! ----------------------------------------------------------------------

          !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
          !        ------------------------------------------

          KFRE = MIN(NFRE, NFRE)

          XMSS = 0.0_JWRB

          !*    2.2 SHALLOW WATER INTEGRATION.
          !         --------------------------

!$acc loop seq
          DO MEANSQS_LF_M=1,KFRE
            TEMP1 = DFIM(MEANSQS_LF_M)*WAVNUM(IJ, MEANSQS_LF_M)**2
            SINFLX_TEMP2 = 0.0_JWRB
!$acc loop seq
            DO MEANSQS_LF_K=1,NANG
              SINFLX_TEMP2 = SINFLX_TEMP2 + FLWD(IJ, MEANSQS_LF_K, MEANSQS_LF_M)
            END DO
            XMSS = XMSS + TEMP1*SINFLX_TEMP2
          END DO


          ! =========================================

          ! [Loki] inlined child subroutine: FEMEAN
          ! =========================================

          ! ----------------------------------------------------------------------


          !*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
          !        ------------------------------------------------

          SINFLX_EM = 0.0_JWRB
          FM = 0.0_JWRB

          SINFLX_DELT25 = WETAIL*FR(NFRE)*DELTH
          DELT2 = FRTAIL*DELTH

          !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
          !        ------------------------------------------

!$acc loop seq
          DO FEMEAN_M=1,NFRE
            FEMEAN_K = 1
            FEMEAN_TEMP2 = MAX(FLWD(IJ, FEMEAN_K, FEMEAN_M), EPSMIN)
!$acc loop seq
            DO FEMEAN_K=2,NANG
              FEMEAN_TEMP2 = FEMEAN_TEMP2 + MAX(FLWD(IJ, FEMEAN_K, FEMEAN_M), EPSMIN)
            END DO
            SINFLX_EM = SINFLX_EM + FEMEAN_TEMP2*DFIM(FEMEAN_M)
            FM = FM + DFIMOFR(FEMEAN_M)*FEMEAN_TEMP2
          END DO


          !*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
          !*       NORMALIZE WITH TOTAL ENERGY.
          !        ------------------------------------------

          SINFLX_EM = SINFLX_EM + SINFLX_DELT25*FEMEAN_TEMP2
          FM = FM + DELT2*FEMEAN_TEMP2
          FM = SINFLX_EM / FM
          FM = MAX(FM, FR(1))


          ! =========================================

          IF (SINFLX_EM > 0.0_JWRB .and. FM < FR(-2 + NFRE)) THEN
            ALPHAP = XMSS / (ZLNFRNFRE - LOG(FM))
            IF (ALPHAP > ALPHAPMAX) THEN
              ! some odd cases, revert to tail value
              F1D = 0.0_JWRB
!$acc loop seq
              DO HALPHAP_K=1,NANG
                F1D = F1D + FLWD(IJ, HALPHAP_K, NFRE)*DELTH
              END DO
              ALPHAP = ZPI4GM2*FR5(NFRE)*F1D
            END IF
          ELSE
            F1D = 0.0_JWRB
!$acc loop seq
            DO HALPHAP_K=1,NANG
              F1D = F1D + FLWD(IJ, HALPHAP_K, NFRE)*DELTH
            END DO
            ALPHAP = ZPI4GM2*FR5(NFRE)*F1D
          END IF

          !     1/2 ALPHAP:
          HALP(IJ) = 0.5_JWRB*MIN(ALPHAP, ALPHAPMAX)


          ! =========================================
        ELSE
          HALP(IJ) = 0.0_JWRB
        END IF

      END IF

      ! [Loki] inlined child subroutine: AIRSEA
      ! =========================================

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
        COSDIFF = COS(WDWAVE(IJ) - TAUWDIR(IJ))
        TAUWACT = MAX(TAUW(IJ)*COSDIFF, EPSMIN)
        LLCOSDIFF = COSDIFF > 0.9_JWRB

        !  USING THE CG MODEL:
        IF (LLGCBZ0) THEN

          IF (LLCAPCHNK) THEN
            CHARNOCK_MIN = CHNKMIN_LOKI(WSWAVE(IJ))
            ALPHAOG = CHARNOCK_MIN*GM1
          ELSE
            ALPHAOG = 0.0_JWRB
          END IF

          USMAX = MAX(-0.21339_JWRB + 0.093698_JWRB*WSWAVE(IJ) - 0.0020944_JWRB*WSWAVE(IJ)**2 + 5.5091E-5_JWRB*WSWAVE(IJ)**3,  &
          & 0.03_JWRB)
          TAUWEFF = MIN(TAUWACT*US2TOTAUW, USMAX**2)


          IF (IUSFG == 0) THEN
            ALPHAGM1 = ALPHA*GM1
            IF (WSWAVE(IJ) < 1.0_JWRB) THEN
              CDFG = 0.002_JWRB
            ELSE IF (LLCOSDIFF) THEN
              X = MIN(TAUWACT / MAX(UFRIC(IJ), EPSUS)**2, PMAX)
              ZCHAR = MIN(ALPHAGM1*UFRIC(IJ)**2 / SQRT(1.0_JWRB - X), 0.05_JWRB*EXP(-0.05_JWRB*(WSWAVE(IJ) - 35._JWRB)))
              ZCHAR = MIN(ZCHAR, ALPHAMAX)
              CDFG = ACDLIN + BCDLIN*SQRT(ZCHAR)*WSWAVE(IJ)
            ELSE
              CDFG = CDM(WSWAVE(IJ))
            END IF
            UFRIC(IJ) = WSWAVE(IJ)*SQRT(CDFG)
          END IF

          W1 = 0.85_JWRB - 0.05_JWRB*(TANH(10.0_JWRB*(WSWAVE(IJ) - 5.0_JWRB)) + 1.0_JWRB)

          XKUTOP = XKAPPA*WSWAVE(IJ)

          USTOLD = UFRIC(IJ)
          TAUOLD = USTOLD**2

!$acc loop seq
          DO ITER=1,NITER
            !         Z0 IS DERIVED FROM THE NEUTRAL LOG PROFILE: UTOP = (USTAR/XKAPPA)*LOG((XNLEV+Z0)/Z0)
            Z0M(IJ) = MAX(XNLEV / (EXP(MIN(XKUTOP / USTOLD, 50.0_JWRB)) - 1.0_JWRB), Z0MIN)
            ! Viscous kinematic stress nu_air * dU/dz at z=0 of the neutral log profile reduced by factor 25 (0.04)
            TAUV = RNUKAPPAM1*USTOLD / Z0M(IJ)

            ANG_GC = ANG_GC_A + ANG_GC_B*TANH(ANG_GC_C*TAUOLD)

            TAUUNR = STRESS_GC_LOKI(ANG_GC, UFRIC(IJ), Z0M(IJ), Z0MIN, HALP(IJ), RNFAC)

            !         TOTAL kinematic STRESS:
            TAUNEW = TAUWEFF + TAUV + TAUUNR
            USTNEW = SQRT(TAUNEW)
            UFRIC(IJ) = W1*USTOLD + (1.0_JWRB - W1)*USTNEW

            !         CONVERGENCE ?
            DEL = UFRIC(IJ) - USTOLD
            IF (ABS(DEL) < PCE_GC*UFRIC(IJ)) EXIT
            TAUOLD = UFRIC(IJ)**2
            USTOLD = UFRIC(IJ)
          END DO

          X = TAUWEFF / TAUOLD

          ! protection just in case there is no convergence
          IF (ITER > NITER .and. X >= PMAX) THEN
            CDFG = CDM(WSWAVE(IJ))
            UFRIC(IJ) = WSWAVE(IJ)*SQRT(CDFG)
            Z0MINRST = UFRIC(IJ)**2*ALPHA*GM1
            Z0M(IJ) = MAX(XNLEV / (EXP(XKUTOP / UFRIC(IJ)) - 1.0_JWRB), Z0MINRST)
            Z0B(IJ) = Z0MINRST
          ELSE
            Z0M(IJ) = MAX(XNLEV / (EXP(XKUTOP / UFRIC(IJ)) - 1.0_JWRB), Z0MIN)
            Z0B(IJ) = Z0M(IJ)*SQRT(TAUUNR / TAUOLD)
          END IF

          !       Refine solution

          IF (X < PMAX) THEN

            USNRF = UFRIC(IJ)
            Z0NRF = Z0M(IJ)
            Z0BNRF = Z0B(IJ)

            USTOLD = UFRIC(IJ)
            TAUOLD = MAX(USTOLD**2, TAUWEFF)

            ALPOG = MAX(MIN(Z0B(IJ) / TAUOLD, ALPHAMAX), ALPHAOG)

!$acc loop seq
            DO ITER=1,NITER
              X = MIN(TAUWEFF / TAUOLD, PMAX)
              USTM1 = 1.0_JWRB / MAX(USTOLD, EPSUS)
              Z0VIS = RNUM*USTM1
              HZ0VISO1MX = 0.5_JWRB*Z0VIS / (1.0_JWRB - X)
              Z0B(IJ) = ALPOG*TAUOLD
              Z0M(IJ) = HZ0VISO1MX + SQRT(HZ0VISO1MX**2 + Z0B(IJ)**2 / (1.0_JWRB - X))

              XOLOGZ0 = 1.0_JWRB / LOG(XNLEV / Z0M(IJ) + 1.0_JWRB)
              F = USTOLD - XKUTOP*XOLOGZ0
              ZZ = 2.0_JWRB*USTM1*(3.0_JWRB*Z0B(IJ)**2 + 0.5_JWRB*Z0VIS*Z0M(IJ) - Z0M(IJ)**2) / (2.0_JWRB*Z0M(IJ)**2*(1.0_JWRB -  &
              & X) - Z0VIS*Z0M(IJ))

              DELF = 1.0_JWRB - XKUTOP*XOLOGZ0**2*ZZ

              IF (DELF /= 0.0_JWRB) UFRIC(IJ) = USTOLD - F / DELF
              !           CONVERGENCE ?
              TAUNEW = MAX(UFRIC(IJ)**2, TAUWEFF)
              UFRIC(IJ) = SQRT(TAUNEW)
              DEL = TAUNEW - TAUOLD
              IF (ABS(DEL) < PCE_GC*TAUOLD) EXIT
              TAUOLD = TAUNEW
              USTOLD = UFRIC(IJ)

            END DO
            ! protection just in case there is no convergence
            IF (ITER > NITER) THEN
              UFRIC(IJ) = USNRF
              Z0M(IJ) = Z0NRF
              Z0B(IJ) = Z0BNRF
              USTM1 = 1.0_JWRB / MAX(UFRIC(IJ), EPSUS)
              Z0VIS = RNUM*USTM1
              CHRNCK(IJ) = MAX(G*(Z0M(IJ) - Z0VIS)*USTM1**2, ALPHAMIN)

            ELSE
              CHRNCK(IJ) = MAX(G*(Z0B(IJ) / SQRT(1.0_JWRB - X)) / MAX(UFRIC(IJ), EPSUS)**2, ALPHAMIN)
            END IF

          ELSE
            USTM1 = 1.0_JWRB / MAX(UFRIC(IJ), EPSUS)
            Z0VIS = RNUM*USTM1
            CHRNCK(IJ) = MAX(G*(Z0M(IJ) - Z0VIS)*USTM1**2, ALPHAMIN)
          END IF



        ELSE

          TAUWEFF = TAUWACT*US2TOTAUW

          IF (LLCAPCHNK) THEN
            CHARNOCK_MIN = CHNKMIN_LOKI(WSWAVE(IJ))
            XMIN = 0.15_JWRB*(ALPHA - CHARNOCK_MIN)
            ALPHAOG = CHARNOCK_MIN*GM1
          ELSE
            XMIN = 0.0_JWRB
            ALPHAOG = ALPHA*GM1
          END IF

          XKUTOP = XKAPPA*WSWAVE(IJ)

          USTOLD = (1 - IUSFG)*WSWAVE(IJ)*SQRT(MIN(ACD + BCD*WSWAVE(IJ), CDMAX)) + IUSFG*UFRIC(IJ)
          TAUOLD = MAX(USTOLD**2, TAUWEFF)
          UFRIC(IJ) = SQRT(TAUOLD)
          USTM1 = 1.0_JWRB / MAX(UFRIC(IJ), EPSUS)

!$acc loop seq
          DO ITER=1,NITER
            X = MAX(TAUWACT / TAUOLD, XMIN)
            Z0CH = ALPHAOG*TAUOLD / SQRT(1.0_JWRB - X)
            Z0VIS = RNUM*USTM1
            Z0TOT = Z0CH + Z0VIS

            XOLOGZ0 = 1.0_JWRB / (XLOGXL - LOG(Z0TOT))
            F = UFRIC(IJ) - XKUTOP*XOLOGZ0
            ZZ = USTM1*(Z0CH*(2.0_JWRB - TWOXMP1*X) / (1.0_JWRB - X) - Z0VIS) / Z0TOT
            DELF = 1.0_JWRB - XKUTOP*XOLOGZ0**2*ZZ

            IF (DELF /= 0.0_JWRB) UFRIC(IJ) = UFRIC(IJ) - F / DELF
            TAUNEW = MAX(UFRIC(IJ)**2, TAUWEFF)
            UFRIC(IJ) = SQRT(TAUNEW)
            IF (TAUNEW == TAUOLD) EXIT
            USTM1 = 1.0_JWRB / MAX(UFRIC(IJ), EPSUS)
            TAUOLD = TAUNEW
          END DO

          Z0M(IJ) = Z0CH
          Z0B(IJ) = ALPHAOG*TAUOLD
          CHRNCK(IJ) = MAX(G*Z0M(IJ)*USTM1**2, ALPHAMIN)


        END IF


        ! =========================================

      ELSE IF (ICODE_WND == 1 .or. ICODE_WND == 2) THEN

        ! [Loki] Removed content of pragma-marked region!
      END IF


      ! =========================================

    END IF

    ! COMPUTE WIND INPUT
    !!  FLD AND SL ARE INITIALISED IN SINPUT
    IF (ICALL < NCALL) THEN
      NGST = 1
      LLPHIWA = .false.
      LLSNEG = .false.
    ELSE
      NGST = 2
      LLPHIWA = .true.
      LLSNEG = .true.
    END IF

    ! [Loki] inlined child subroutine: SINPUT
    ! =========================================
    ! ----------------------------------------------------------------------


    SELECT CASE (IPHYS)
    CASE (0)
      ! [Loki] Removed content of pragma-marked region!
    CASE (1)
      ! [Loki] inlined child subroutine: SINPUT_ARD
      ! =========================================
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
          SINPUT_ZCHAR = G*(Z0M(IJ) - WSIGSTAR_Z0VIS) / MAX(UFRIC(IJ)**2, EPSUS)
          SINPUT_ZCHAR = MAX(MIN(SINPUT_ZCHAR, ALPHAMAX), ALPHAMIN)

          BCD_LOC = BCDLIN*SQRT(SINPUT_ZCHAR)
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
          SINPUT_XKAPPAD = 1.0_JWRB / XKAPPA
          U10 = UFRIC(IJ)*SINPUT_XKAPPAD*(LOG(10.0_JWRB) - LOG(Z0M(IJ)))
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
        CSTRNFAC = CONSTN*RNFAC / RAORW(IJ)
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
        DO SINPUT_M=1,NFRE
          SIG = ZPIFR(SINPUT_M)
          SIG2 = SIG**2
          DFIM_SIG2 = DFIM(SINPUT_M)*SIG2

          SINPUT_K = 1
          SINFLX_TEMP = FL1(IJ, SINPUT_K, SINPUT_M)
!$acc loop seq
          DO SINPUT_K=2,NANG
            SINFLX_TEMP = SINFLX_TEMP + FL1(IJ, SINPUT_K, SINPUT_M)
          END DO

          UORBT = UORBT + DFIM_SIG2*SINFLX_TEMP
          AORB = AORB + DFIM(SINPUT_M)*SINFLX_TEMP
        END DO

        UORBT = 2.0_JWRB*SQRT(UORBT)          ! this is the significant orbital amplitude
        AORB = 2.0_JWRB*SQRT(AORB)          ! this 1/2 Hs
        RE = FACM1_NU_AIR*UORBT*AORB          ! this is the Reynolds number
        SINPUT_Z0VIS = FAC_NU_AIR / MAX(UFRIC(IJ), 0.0001_JWRB)
        Z0TUB = Z0RAT*MIN(Z0TUBMAX, Z0M(IJ))
        Z0NOZ = MAX(SINPUT_Z0VIS, Z0TUB)
        ZORB = AORB / Z0NOZ

        !         compute fww
        SINPUT_XI = (LOG10(MAX(ZORB, 3.0_JWRB)) - ABMIN)*DELABM1
        IND = MIN(IAB - 1, INT(SINPUT_XI))
        SINPUT_DELI1 = MIN(1.0_JWRB, SINPUT_XI - REAL(IND, kind=JWRB))
        SINPUT_DELI2 = 1.0_JWRB - SINPUT_DELI1
        FWW = SWELLFT(IND)*SINPUT_DELI2 + SWELLFT(IND + 1)*SINPUT_DELI1
        SINPUT_TEMP2 = FWW*UORBT

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
        DO SINPUT_K=1,NANG
          COSLP(IJ, SINPUT_K) = COSWDIF(IJ, SINPUT_K)
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
        DO SINPUT_K=1,NANG
          DSTAB(IJ, SINPUT_K, 1) = 0.0_JWRB
          DSTAB(IJ, SINPUT_K, 2) = 0.0_JWRB
        END DO
      END IF

!$acc loop seq
      DO SINPUT_M=1,NFRE
        SIG = ZPIFR(SINPUT_M)
        SIG2 = SIG**2
        SINPUT_CONST = SIG*CONST1

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

          CONSTF = ROGOROAIR*CINV(IJ, SINPUT_M)*DFIM(SINPUT_M)
        END IF


        !*      PRECALCULATE FREQUENCY DEPENDENCE.
        !       ----------------------------------

!$acc loop seq
        DO IGST=1,2
          IF (IGST <= NGST) THEN
            UCN(IGST) = USTP(IGST)*CINV(IJ, SINPUT_M)
            UCNZALPD(IGST) = XKAPPA / (UCN(IGST) + ZALP)
          END IF
        END DO
        ZCN = LOG(WAVNUM(IJ, SINPUT_M)*Z0M(IJ))
        CNSN = SINPUT_CONST*RAORW(IJ)

        !*    2.1 LOOP OVER DIRECTIONS.
        !         ---------------------

!$acc loop seq
        DO SINPUT_K=1,NANG
          XLLWS(IJ, SINPUT_K, SINPUT_M) = 0.0_JWRB
        END DO

        IF (LLNORMAGAM) THEN
          XNGAMCONST = CSTRNFAC*XK2CG(IJ, SINPUT_M)
        END IF

        IF (LLSNEG) THEN
          !       SWELL DAMPING:
          DSTAB1 = COEF5*AIRD_PVISC*WAVNUM(IJ, SINPUT_M)
          SINPUT_TEMP1 = COEF*RAORW(IJ)
        END IF

!$acc loop seq
        DO IGST=1,2
          IF (IGST <= NGST) THEN
!$acc loop seq
            DO SINPUT_K=1,NANG
              IF (LTAUWSHELTER) THEN
                COSLP(IJ, SINPUT_K) = COS(TH(SINPUT_K) - USDIRP(IGST))
              END IF
              GAM0(IJ, SINPUT_K, IGST) = 0.0_JWRB
              IF (COSLP(IJ, SINPUT_K) > 0.01_JWRB) THEN
                SINPUT_X = COSLP(IJ, SINPUT_K)*UCN(IGST)
                ZLOG = ZCN + UCNZALPD(IGST) / COSLP(IJ, SINPUT_K)
                IF (ZLOG < 0.0_JWRB) THEN
                  ZLOG2X = ZLOG*ZLOG*SINPUT_X
                  GAM0(IJ, SINPUT_K, IGST) = EXP(ZLOG)*ZLOG2X*ZLOG2X*CNSN
                  XLLWS(IJ, SINPUT_K, SINPUT_M) = 1.0_JWRB
                END IF
              END IF
            END DO

            IF (LLNORMAGAM) THEN

              SUMF = 0.0_JWRB
              SUMFSIN2 = 0.0_JWRB
!$acc loop seq
              DO SINPUT_K=1,NANG
                SUMF = SUMF + GAM0(IJ, SINPUT_K, IGST)*FL1(IJ, SINPUT_K, SINPUT_M)
                SUMFSIN2 = SUMFSIN2 + GAM0(IJ, SINPUT_K, IGST)*FL1(IJ, SINPUT_K, SINPUT_M)*SINWDIF2(IJ, SINPUT_K)
              END DO

              ZNZ = XNGAMCONST*USTPM1(IGST)
              GAMNORMA(IGST) = (1.0_JWRB + ZNZ*SUMFSIN2) / (1.0_JWRB + ZNZ*SUMF)

            END IF

            IF (LLSNEG) THEN
!$acc loop seq
              DO SINPUT_K=1,NANG
                DSTAB2 = SINPUT_TEMP1*(SINPUT_TEMP2 + (FU + FUD*COSLP(IJ, SINPUT_K))*USTP(IGST))
                DSTAB(IJ, SINPUT_K, IGST) = DSTAB1 + PTURB*DSTAB2
              END DO
            END IF
          END IF
        END DO


        !*    2.2 UPDATE THE SHELTERING STRESS (in any),
        !         AND THEN ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
        !         ---------------------------------------------------------

!$acc loop seq
        DO SINPUT_K=1,NANG

!$acc loop seq
          DO IGST=1,2
            IF (IGST <= NGST) THEN
              ! SLP: only the positive contributions
              SLP = GAM0(IJ, SINPUT_K, IGST)*GAMNORMA(IGST)
              FLP = SLP + DSTAB(IJ, SINPUT_K, IGST)

              SLP = SLP*FL1(IJ, SINPUT_K, SINPUT_M)

              IF (LTAUWSHELTER) THEN
                CONST11 = CONSTF*SINTH(SINPUT_K)
                CONST22 = CONSTF*COSTH(SINPUT_K)
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

          SPOS(IJ, SINPUT_K, SINPUT_M) = AVG_GST*SLP_AVG

          FLD(IJ, SINPUT_K, SINPUT_M) = AVG_GST*FLP_AVG
          SL(IJ, SINPUT_K, SINPUT_M) = FLD(IJ, SINPUT_K, SINPUT_M)*FL1(IJ, SINPUT_K, SINPUT_M)

        END DO

      END DO
      ! END LOOP OVER FREQUENCIES


      ! =========================================
    END SELECT


    ! =========================================


    ! MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
    ! [Loki] inlined child subroutine: FEMEANWS
    ! =========================================

    ! ----------------------------------------------------------------------


    !*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
    !        ------------------------------------------------

    EM_LOC = EPSMIN
    FMEANWS(IJ) = EPSMIN

    FEMEANWS_DELT25 = WETAIL*FR(NFRE)*DELTH
    FEMEANWS_DELT2 = FRTAIL*DELTH


    !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
    !        ------------------------------------------

!$acc loop seq
    DO FEMEANWS_M=1,NFRE
      FEMEANWS_TEMP2 = 0.0_JWRB
!$acc loop seq
      DO FEMEANWS_K=1,NANG
        FEMEANWS_TEMP2 = FEMEANWS_TEMP2 + XLLWS(IJ, FEMEANWS_K, FEMEANWS_M)*FL1(IJ, FEMEANWS_K, FEMEANWS_M)
      END DO
      EM_LOC = EM_LOC + DFIM(FEMEANWS_M)*FEMEANWS_TEMP2
      FMEANWS(IJ) = FMEANWS(IJ) + DFIMOFR(FEMEANWS_M)*FEMEANWS_TEMP2
    END DO

    !*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
    !*       NORMALIZE WITH TOTAL ENERGY.
    !        ------------------------------------------

    EM_LOC = EM_LOC + FEMEANWS_DELT25*FEMEANWS_TEMP2
    FMEANWS(IJ) = FMEANWS(IJ) + FEMEANWS_DELT2*FEMEANWS_TEMP2
    FMEANWS(IJ) = EM_LOC / FMEANWS(IJ)



    ! =========================================

    ! COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
    ! [Loki] inlined child subroutine: FRCUTINDEX
    ! =========================================

    ! ----------------------------------------------------------------------


    !*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
    !*    FREQUENCIES LE MAX(TAILFACTOR*MAX(FMNWS,FM),TAILFACTOR_PM*FPM),
    !*    WHERE FPM IS THE PIERSON-MOSKOWITZ FREQUENCY BASED ON FRICTION
    !*    VELOCITY. (FPM=G/(FRIC*ZPI*USTAR))
    !     ------------------------------------------------------------

    FPMH = TAILFACTOR / FR(1)
    FPPM = TAILFACTOR_PM*G / (FRIC*ZPIFR(1))

    IF (CICOVER(IJ) <= CITHRSH_TAIL) THEN
      FM2 = MAX(FMEANWS(IJ), FMEAN(IJ))*FPMH
      FPM = FPPM / MAX(UFRIC(IJ), EPSMIN)
      FPM4 = MAX(FM2, FPM)
      MIJ(IJ) = NINT(LOG10(FPM4)*FLOGSPRDM1) + 1
      MIJ(IJ) = MIN(MAX(1, MIJ(IJ)), NFRE)
    ELSE
      MIJ(IJ) = NFRE
    END IF

    !     SET RHOWGDFTH
!$acc loop seq
    DO FRCUTINDEX_M=1,MIJ(IJ)
      RHOWGDFTH(IJ, FRCUTINDEX_M) = RHOWG_DFIM(FRCUTINDEX_M)
    END DO
    IF (MIJ(IJ) /= NFRE) RHOWGDFTH(IJ, MIJ(IJ)) = 0.5_JWRB*RHOWGDFTH(IJ, MIJ(IJ))
!$acc loop seq
    DO FRCUTINDEX_M=MIJ(IJ) + 1,NFRE
      RHOWGDFTH(IJ, FRCUTINDEX_M) = 0.0_JWRB
    END DO


    ! =========================================

    ! UPDATE TAUW
    ! [Loki] inlined child subroutine: STRESSO
    ! =========================================

    ! ----------------------------------------------------------------------


    PHIWA(IJ) = 0.0_JWRB
    STRESSO_XSTRESS = 0.0_JWRB
    STRESSO_YSTRESS = 0.0_JWRB

    !*    CONTRIBUTION TO THE WAVE STRESS FROM THE NEGATIVE PART OF THE WIND INPUT
    !     ------------------------------------------------------------------------

    IF (LLPHIWA) THEN
      !     full energy flux due to negative Sinput (SL-SPOS)
      !     we assume that above NFRE, the contibutions can be neglected
!$acc loop seq
      DO STRESSO_M=1,NFRE
!$acc loop seq
        DO STRESSO_K=1,NANG
          PHIWA(IJ) = PHIWA(IJ) + (SL(IJ, STRESSO_K, STRESSO_M) - SPOS(IJ, STRESSO_K, STRESSO_M))*RHOWG_DFIM(STRESSO_M)
        END DO
      END DO
    END IF

    !*    CALCULATE LOW-FREQUENCY CONTRIBUTION TO STRESS AND ENERGY FLUX (positive sinput).
    !     ---------------------------------------------------------------------------------
!$acc loop seq
    DO STRESSO_M=1,NFRE
      !     THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
      STRESSO_K = 1
      SUMX = SPOS(IJ, STRESSO_K, STRESSO_M)*SINTH(STRESSO_K)
      SUMY = SPOS(IJ, STRESSO_K, STRESSO_M)*COSTH(STRESSO_K)
      SUMT = SPOS(IJ, STRESSO_K, STRESSO_M)
!$acc loop seq
      DO STRESSO_K=2,NANG
        SUMX = SUMX + SPOS(IJ, STRESSO_K, STRESSO_M)*SINTH(STRESSO_K)
        SUMY = SUMY + SPOS(IJ, STRESSO_K, STRESSO_M)*COSTH(STRESSO_K)
        SUMT = SUMT + SPOS(IJ, STRESSO_K, STRESSO_M)
      END DO
      CMRHOWGDFTH = RHOWGDFTH(IJ, STRESSO_M)*CINV(IJ, STRESSO_M)
      STRESSO_XSTRESS = STRESSO_XSTRESS + CMRHOWGDFTH*SUMX
      STRESSO_YSTRESS = STRESSO_YSTRESS + CMRHOWGDFTH*SUMY

      IF (LLPHIWA) THEN
        PHIWA(IJ) = PHIWA(IJ) + RHOWGDFTH(IJ, STRESSO_M)*SUMT
      END IF
    END DO

    !     TAUW is the kinematic wave stress !
    STRESSO_XSTRESS = STRESSO_XSTRESS / MAX(AIRD(IJ), 1.0_JWRB)
    STRESSO_YSTRESS = STRESSO_YSTRESS / MAX(AIRD(IJ), 1.0_JWRB)


    !*    CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS and energy flux (positive sinput).
    !     ----------------------------------------------------------------------------------

    IF (IPHYS == 0 .or. TAUWSHELTER == 0.0_JWRB) THEN
      STRESSO_LTAUWSHELTER = .false.
      STRESSO_USDIRP = WDWAVE(IJ)
      UST = UFRIC(IJ)
    ELSE
      STRESSO_LTAUWSHELTER = .true.
      STRESSO_TAUX = UFRIC(IJ)**2*SIN(WDWAVE(IJ))
      STRESSO_TAUY = UFRIC(IJ)**2*COS(WDWAVE(IJ))
      STRESSO_TAUPX = STRESSO_TAUX - TAUWSHELTER*STRESSO_XSTRESS
      STRESSO_TAUPY = STRESSO_TAUY - TAUWSHELTER*STRESSO_YSTRESS
      STRESSO_USDIRP = ATAN2(STRESSO_TAUPX, STRESSO_TAUPY)
      UST = (STRESSO_TAUPX**2 + STRESSO_TAUPY**2)**0.25_JWRB
    END IF

    ! [Loki] inlined child subroutine: TAU_PHI_HF
    ! =========================================


    ! ----------------------------------------------------------------------


    IF (LLGCBZ0) THEN
      ! [Loki] inlined child subroutine: OMEGAGC
      ! =========================================
      ! ----------------------------------------------------------------------


      NS = NS_GC_LOKI(UFRIC(IJ))
      XKS = XK_GC(NS)
      OMS = OMEGA_GC(NS)


      ! =========================================
    END IF

    !     See INIT_X0TAUHF
    X0G = X0TAUHF*G

    IF (LLPHIWA) USTPH = UST

    !*    COMPUTE THE INTEGRALS
    !     ---------------------

    XLOGGZ0 = LOG(G*Z0M(IJ))
    OMEGACC = MAX(ZPIFR(MIJ(IJ)), X0G / UST)
    SQRTZ0OG = SQRT(Z0M(IJ)*GM1)
    SQRTGZ0 = 1.0_JWRB / SQRTZ0OG
    YC = OMEGACC*SQRTZ0OG
    ZINF = LOG(YC)

    CONSTTAU = ZPI4GM2*FR5(MIJ(IJ))

    TAU_PHI_HF_K = 1
    TAU_PHI_HF_COSW = MAX(COSWDIF(IJ, TAU_PHI_HF_K), 0.0_JWRB)
    TAU_PHI_HF_FCOSW2 = FL1(IJ, TAU_PHI_HF_K, MIJ(IJ))*TAU_PHI_HF_COSW**2
    F1DCOS3 = TAU_PHI_HF_FCOSW2*TAU_PHI_HF_COSW
    F1DCOS2 = TAU_PHI_HF_FCOSW2
    F1DSIN2 = FL1(IJ, TAU_PHI_HF_K, MIJ(IJ))*SINWDIF2(IJ, TAU_PHI_HF_K)
    STRESSO_F1D = FL1(IJ, TAU_PHI_HF_K, MIJ(IJ))
!$acc loop seq
    DO TAU_PHI_HF_K=2,NANG
      TAU_PHI_HF_COSW = MAX(COSWDIF(IJ, TAU_PHI_HF_K), 0.0_JWRB)
      TAU_PHI_HF_FCOSW2 = FL1(IJ, TAU_PHI_HF_K, MIJ(IJ))*TAU_PHI_HF_COSW**2
      F1DCOS3 = F1DCOS3 + TAU_PHI_HF_FCOSW2*TAU_PHI_HF_COSW
      F1DCOS2 = F1DCOS2 + TAU_PHI_HF_FCOSW2
      F1DSIN2 = F1DSIN2 + FL1(IJ, TAU_PHI_HF_K, MIJ(IJ))*SINWDIF2(IJ, TAU_PHI_HF_K)
      STRESSO_F1D = STRESSO_F1D + FL1(IJ, TAU_PHI_HF_K, MIJ(IJ))
    END DO
    F1DCOS3 = DELTH*F1DCOS3
    F1DCOS2 = DELTH*F1DCOS2
    F1DSIN2 = DELTH*F1DSIN2
    STRESSO_F1D = DELTH*STRESSO_F1D

    IF (LLNORMAGAM) THEN
      CONFG = GAMNCONST*FR5(MIJ(IJ))*RNFAC*SQRTGZ0
      STRESSO_CONST1 = CONFG*F1DSIN2
      CONST2 = CONFG*STRESSO_F1D
    ELSE
      STRESSO_CONST1 = 0.0_JWRB
      CONST2 = 0.0_JWRB
    END IF


    !     TAUHF :
    IF (LLGCBZ0) THEN
      ZSUP = MIN(LOG(OMS*SQRTZ0OG), ZSUPMAX)
    ELSE
      ZSUP = ZSUPMAX
    END IF

    TAUL = UST**2
    DELZ = MAX((ZSUP - ZINF) / REAL(JTOT_TAUHF - 1, kind=JWRB), 0.0_JWRB)
    TAUHF = 0.0_JWRB

    ! Intergrals are integrated following a change of variable : Z=LOG(Y)
    IF (STRESSO_LTAUWSHELTER) THEN
!$acc loop seq
      DO TAU_PHI_HF_J=1,JTOT_TAUHF
        Y = EXP(ZINF + REAL(TAU_PHI_HF_J - 1, kind=JWRB)*DELZ)
        OMEGA = Y*SQRTGZ0
        CM1 = OMEGA*GM1
        ZX = UST*CM1 + ZALP
        ZARG = XKAPPA / ZX
        STRESSO_ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
        STRESSO_ZLOG = MIN(STRESSO_ZLOG, 0.0_JWRB)
        ZBETA = STRESSO_ZLOG**4*EXP(STRESSO_ZLOG)
        STRESSO_ZNZ = ZBETA*UST*Y
        STRESSO_GAMNORMA = (1.0_JWRB + STRESSO_CONST1*STRESSO_ZNZ) / (1.0_JWRB + CONST2*STRESSO_ZNZ)
        FNC2 = F1DCOS3*CONSTTAU*ZBETA*TAUL*WTAUHF(TAU_PHI_HF_J)*DELZ*STRESSO_GAMNORMA
        TAUL = MAX(TAUL - TAUWSHELTER*FNC2, 0.0_JWRB)

        UST = SQRT(TAUL)
        TAUHF = TAUHF + FNC2
      END DO
    ELSE
!$acc loop seq
      DO TAU_PHI_HF_J=1,JTOT_TAUHF
        Y = EXP(ZINF + REAL(TAU_PHI_HF_J - 1, kind=JWRB)*DELZ)
        OMEGA = Y*SQRTGZ0
        CM1 = OMEGA*GM1
        ZX = UST*CM1 + ZALP
        ZARG = XKAPPA / ZX
        STRESSO_ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
        STRESSO_ZLOG = MIN(STRESSO_ZLOG, 0.0_JWRB)
        ZBETA = STRESSO_ZLOG**4*EXP(STRESSO_ZLOG)
        FNC2 = ZBETA*WTAUHF(TAU_PHI_HF_J)
        STRESSO_ZNZ = ZBETA*UST*Y
        STRESSO_GAMNORMA = (1.0_JWRB + STRESSO_CONST1*STRESSO_ZNZ) / (1.0_JWRB + CONST2*STRESSO_ZNZ)
        TAUHF = TAUHF + FNC2*STRESSO_GAMNORMA
      END DO
      TAUHF = F1DCOS3*CONSTTAU*TAUL*TAUHF*DELZ
    END IF


    PHIHF = 0.0_JWRB
    IF (LLPHIWA) THEN
      !       PHIHF:
      !       We are neglecting the gravity-capillary contribution
      !       Recompute DELZ over the full interval
      TAUL = USTPH**2
      ZSUP = ZSUPMAX
      DELZ = MAX((ZSUP - ZINF) / REAL(JTOT_TAUHF - 1, kind=JWRB), 0.0_JWRB)

      CONSTPHI = AIRD(IJ)*ZPI4GM1*FR5(MIJ(IJ))

      ! Intergrals are integrated following a change of variable : Z=LOG(Y)
      IF (STRESSO_LTAUWSHELTER) THEN
!$acc loop seq
        DO TAU_PHI_HF_J=1,JTOT_TAUHF
          Y = EXP(ZINF + REAL(TAU_PHI_HF_J - 1, kind=JWRB)*DELZ)
          OMEGA = Y*SQRTGZ0
          CM1 = OMEGA*GM1
          ZX = USTPH*CM1 + ZALP
          ZARG = XKAPPA / ZX
          STRESSO_ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
          STRESSO_ZLOG = MIN(STRESSO_ZLOG, 0.0_JWRB)
          ZBETA = STRESSO_ZLOG**4*EXP(STRESSO_ZLOG)
          STRESSO_ZNZ = ZBETA*UST*Y
          STRESSO_GAMNORMA = (1.0_JWRB + STRESSO_CONST1*STRESSO_ZNZ) / (1.0_JWRB + CONST2*STRESSO_ZNZ)
          FNC2 = ZBETA*TAUL*WTAUHF(TAU_PHI_HF_J)*DELZ*STRESSO_GAMNORMA
          TAUL = MAX(TAUL - TAUWSHELTER*F1DCOS3*CONSTTAU*FNC2, 0.0_JWRB)
          USTPH = SQRT(TAUL)
          PHIHF = PHIHF + FNC2 / Y
        END DO
        PHIHF = F1DCOS2*CONSTPHI*SQRTZ0OG*PHIHF
      ELSE
!$acc loop seq
        DO TAU_PHI_HF_J=1,JTOT_TAUHF
          Y = EXP(ZINF + REAL(TAU_PHI_HF_J - 1, kind=JWRB)*DELZ)
          OMEGA = Y*SQRTGZ0
          CM1 = OMEGA*GM1
          ZX = USTPH*CM1 + ZALP
          ZARG = XKAPPA / ZX
          STRESSO_ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
          STRESSO_ZLOG = MIN(STRESSO_ZLOG, 0.0_JWRB)
          ZBETA = STRESSO_ZLOG**4*EXP(STRESSO_ZLOG)
          STRESSO_ZNZ = ZBETA*UST*Y
          STRESSO_GAMNORMA = (1.0_JWRB + STRESSO_CONST1*STRESSO_ZNZ) / (1.0_JWRB + CONST2*STRESSO_ZNZ)
          FNC2 = ZBETA*WTAUHF(TAU_PHI_HF_J)*STRESSO_GAMNORMA
          PHIHF = PHIHF + FNC2 / Y
        END DO
        PHIHF = F1DCOS2*CONSTPHI*SQRTZ0OG*TAUL*PHIHF*DELZ
      END IF
    END IF


    ! =========================================

    STRESSO_XSTRESS = STRESSO_XSTRESS + TAUHF*SIN(STRESSO_USDIRP)
    STRESSO_YSTRESS = STRESSO_YSTRESS + TAUHF*COS(STRESSO_USDIRP)
    TAUW(IJ) = SQRT(STRESSO_XSTRESS**2 + STRESSO_YSTRESS**2)
    TAUW(IJ) = MAX(TAUW(IJ), 0.0_JWRB)
    TAUWDIR(IJ) = ATAN2(STRESSO_XSTRESS, STRESSO_YSTRESS)

    IF (.not.LLGCBZ0) THEN
      TAUTOUS2 = 1.0_JWRB / (1.0_JWRB + EPS1)
      TAUW(IJ) = MIN(TAUW(IJ), UFRIC(IJ)**2*TAUTOUS2)
    END IF

    IF (LLPHIWA) THEN
      PHIWA(IJ) = PHIWA(IJ) + PHIHF
    END IF


    ! =========================================

    ! ----------------------------------------------------------------------



!$acc end data
  END SUBROUTINE SINFLX_LOKI
END MODULE SINFLX_LOKI_MOD
