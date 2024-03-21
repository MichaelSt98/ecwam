! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SINPUT_ARD_CUF_HOIST_NEW_MOD
  !CONTAINED SUBROUTINES:
  ! - WSIGSTAR
  ! - SINPUT_ARD
  ! - SINPUT_JAN
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE WSIGSTAR_CUF_HOIST_NEW (WSWAVE, UFRIC, Z0M, WSTAR, SIG_N, ACDLIN, ALPHAMAX, ALPHAMIN, BCDLIN,  &
  & EPSUS, G, LLGCBZ0, RNUM, WSPMIN, XKAPPA)
    ! ----------------------------------------------------------------------
    
    !**** *WSIGSTAR* - COMPUTATION OF THE RELATIVE STANDARD DEVIATION OF USTAR.
    
    !*    PURPOSE.
    !     ---------
    
    !     COMPUTES THE STANDARD DEVIATION OF USTAR DUE TO SMALL SCALE GUSTINESS
    !     RELATIVE TO USTAR
    
    !**   INTERFACE.
    !     ----------
    
    !     *CALL* *WSIGSTAR (KIJS, KIJL, WSWAVE, UFRIC, Z0M, WSTAR, SIG_N)
    !             *KIJS*   - INDEX OF FIRST GRIDPOINT.
    !             *KIJL*   - INDEX OF LAST GRIDPOINT.
    !             *WSWAVE* - 10M WIND SPEED (m/s).
    !             *UFRIC*  - NEW FRICTION VELOCITY IN M/S.
    !             *Z0M*    - ROUGHNESS LENGTH IN M.
    !             *WSTAR*  - FREE CONVECTION VELOCITY SCALE (M/S).
    !             *SIG_N*  - ESTINATED RELATIVE STANDARD DEVIATION OF USTAR.
    
    !     METHOD.
    !     -------
    
    !     USE PANOFSKY (1991) TO EXPRESS THE STANDARD DEVIATION OF U10 IN TERMS
    !     USTAR AND  w* THE CONVECTIVE VELOCITY SCALE.
    !     (but with the background gustiness set to 0.)
    !     and USTAR=SQRT(Cd)*U10 to DERIVE THE STANDARD DEVIATION OF USTAR.
    !     WITH CD=A+B*U10 (see below).
    
    !     REFERENCE.
    !     ----------
    
    !     SEE SECTION 3.2.1 OF THE WAM DOCUMENTATION.
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    REAL(KIND=JWRB), INTENT(IN) :: WSWAVE, UFRIC, Z0M, WSTAR
    REAL(KIND=JWRB), INTENT(OUT) :: SIG_N
    
    REAL(KIND=JWRB), PARAMETER :: BG_GUST = 0.0_JWRB    ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
    REAL(KIND=JWRB), PARAMETER :: ONETHIRD = 1.0_JWRB / 3.0_JWRB
    REAL(KIND=JWRB), PARAMETER :: SIG_NMAX = 0.9_JWRB    ! MAX OF RELATIVE STANDARD DEVIATION OF USTAR
    
    REAL(KIND=JWRB), PARAMETER :: LOG10 = LOG(10.0_JWRB)
    REAL(KIND=JWRB), PARAMETER :: C1 = 1.03E-3_JWRB
    REAL(KIND=JWRB), PARAMETER :: C2 = 0.04E-3_JWRB
    REAL(KIND=JWRB), PARAMETER :: P1 = 1.48_JWRB
    REAL(KIND=JWRB), PARAMETER :: P2 = -0.21_JWRB
    
    ! $ loki routine seq
    REAL(KIND=JWRB) :: ZCHAR, C_D, DC_DDU, SIG_CONV
    REAL(KIND=JWRB) :: XKAPPAD, U10, C2U10P1, U10P2
    REAL(KIND=JWRB) :: BCD, U10M1, ZN, Z0VIS
    REAL(KIND=JWRB), INTENT(IN) :: ACDLIN
    REAL(KIND=JWRB), INTENT(IN) :: ALPHAMAX
    REAL(KIND=JWRB), INTENT(IN) :: ALPHAMIN
    REAL(KIND=JWRB), INTENT(IN) :: BCDLIN
    REAL(KIND=JWRB), INTENT(IN) :: EPSUS
    REAL(KIND=JWRB), INTENT(IN) :: G
    LOGICAL, INTENT(IN) :: LLGCBZ0
    REAL(KIND=JWRB), INTENT(IN) :: RNUM
    REAL(KIND=JWRB), INTENT(IN) :: WSPMIN
    REAL(KIND=JWRB), INTENT(IN) :: XKAPPA
!$acc routine seq
    
    
    ! ----------------------------------------------------------------------
    
    
    
    IF (LLGCBZ0) THEN
      ZN = RNUM
      
      U10M1 = 1.0_JWRB / MAX(WSWAVE, WSPMIN)
      ! CHARNOCK:
      Z0VIS = ZN / MAX(UFRIC, EPSUS)
      ZCHAR = G*(Z0M - Z0VIS) / MAX(UFRIC**2, EPSUS)
      ZCHAR = MAX(MIN(ZCHAR, ALPHAMAX), ALPHAMIN)
      
      BCD = BCDLIN*SQRT(ZCHAR)
      C_D = ACDLIN + BCD*WSWAVE
      DC_DDU = BCD
      SIG_CONV = 1.0_JWRB + 0.5_JWRB*WSWAVE / C_D*DC_DDU
      SIG_N = MIN(SIG_NMAX, SIG_CONV*U10M1*(BG_GUST*UFRIC**3 + 0.5_JWRB*XKAPPA*WSTAR**3)**ONETHIRD)
    ELSE
      ZN = 0.0_JWRB
      
      !!! for consistency I have kept the old method, even though the new method above could be used,
      !!! but until LLGCBZ0 is the default, keep the old scheme whe it is not...
      !
      !       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
      !       BASED ON U*
      !
      XKAPPAD = 1.0_JWRB / XKAPPA
      U10 = UFRIC*XKAPPAD*(LOG10 - LOG(Z0M))
      U10 = MAX(U10, WSPMIN)
      U10M1 = 1.0_JWRB / U10
      C2U10P1 = C2*U10**P1
      U10P2 = U10**P2
      C_D = (C1 + C2U10P1)*U10P2
      DC_DDU = (P2*C1 + (P1 + P2)*C2U10P1)*U10P2*U10M1
      SIG_CONV = 1.0_JWRB + 0.5_JWRB*U10 / C_D*DC_DDU
      SIG_N = MIN(SIG_NMAX, SIG_CONV*U10M1*(BG_GUST*UFRIC**3 + 0.5_JWRB*XKAPPA*WSTAR**3)**ONETHIRD)
    END IF
    
    
  END SUBROUTINE WSIGSTAR_CUF_HOIST_NEW
  ATTRIBUTES(DEVICE) SUBROUTINE SINPUT_ARD_CUF_HOIST_NEW (NGST, LLSNEG, KIJS, KIJL, FL1, WAVNUM, CINV, XK2CG, WDWAVE, WSWAVE,  &
  & UFRIC, Z0M, COSWDIF, SINWDIF2, RAORW, WSTAR, RNFAC, FLD, SL, SPOS, XLLWS, ABMAX, ABMIN, ACDLIN, ALPHAMAX, ALPHAMIN, BCDLIN,  &
  & BETAMAXOXKAPPA2, COSTH, DELTH, DFIM, EPSMIN, EPSUS, G, IAB, LLGCBZ0, LLNORMAGAM, NANG, NFRE, RNU, RNUM, SINTH, SWELLF,  &
  & SWELLF2, SWELLF3, SWELLF4, SWELLF5, SWELLF6, SWELLF7, SWELLF7M1, SWELLFT, TAUWSHELTER, TH, WSPMIN, XKAPPA, Z0RAT, Z0TUBMAX,  &
  & ZALP, ZPI, ZPIFR, ICHNK, NCHNK, IJ)
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
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    USE YOWCOUP, ONLY: LLCAPCHNK
    USE YOWFRED, ONLY: FR
    USE YOWPARAM, ONLY: NANG_PARAM
    USE YOWPCONS, ONLY: GM1
    USE YOWPHYS, ONLY: RN1_RN
    USE YOWTEST, ONLY: IU06
    USE YOWSTAT, ONLY: IDAMPING
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
    INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NGST
    LOGICAL, VALUE, INTENT(IN) :: LLSNEG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    REAL(KIND=JWRB), INTENT(IN) :: FL1(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: WAVNUM(KIJL, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: CINV(KIJL, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: XK2CG(KIJL, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: WDWAVE(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: WSWAVE(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: UFRIC(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: Z0M(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: RAORW(KIJL)
    REAL(KIND=JWRB), INTENT(IN) :: RNFAC(KIJL)
    REAL(KIND=JWRB), INTENT(IN) :: WSTAR(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: COSWDIF(KIJL, NANG_loki_param)
    REAL(KIND=JWRB), INTENT(IN) :: SINWDIF2(KIJL, NANG_loki_param)
    REAL(KIND=JWRB), INTENT(OUT) :: FLD(KIJL, NANG_loki_param, NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(OUT) :: SL(KIJL, NANG_loki_param, NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(OUT) :: SPOS(KIJL, NANG_loki_param, NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(OUT) :: XLLWS(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: K
    INTEGER(KIND=JWIM) :: M
    INTEGER(KIND=JWIM) :: IND
    INTEGER(KIND=JWIM) :: IGST
    
    REAL(KIND=JWRB) :: CONSTN
    REAL(KIND=JWRB) :: AVG_GST
    REAL(KIND=JWRB) :: ABS_TAUWSHELTER
    REAL(KIND=JWRB) :: CONST1
    REAL(KIND=JWRB) :: ZNZ
    REAL(KIND=JWRB) :: X1
    REAL(KIND=JWRB) :: X2
    REAL(KIND=JWRB) :: ZLOG
    REAL(KIND=JWRB) :: ZLOG1
    REAL(KIND=JWRB) :: ZLOG2
    REAL(KIND=JWRB) :: ZLOG2X
    REAL(KIND=JWRB) :: XV1
    REAL(KIND=JWRB) :: XV2
    REAL(KIND=JWRB) :: ZBETA1
    REAL(KIND=JWRB) :: ZBETA2
    REAL(KIND=JWRB) :: XI
    REAL(KIND=JWRB) :: X
    REAL(KIND=JWRB) :: DELI1
    REAL(KIND=JWRB) :: DELI2
    REAL(KIND=JWRB) :: FU
    REAL(KIND=JWRB) :: FUD
    REAL(KIND=JWRB) :: NU_AIR
    REAL(KIND=JWRB) :: SMOOTH
    REAL(KIND=JWRB) :: HFTSWELLF6
    REAL(KIND=JWRB) :: Z0TUB
    REAL(KIND=JWRB) :: FAC_NU_AIR
    REAL(KIND=JWRB) :: FACM1_NU_AIR
    REAL(KIND=JWRB) :: ARG
    REAL(KIND=JWRB) :: DELABM1
    REAL(KIND=JWRB) :: TAUPX
    REAL(KIND=JWRB) :: TAUPY
    REAL(KIND=JWRB) :: DSTAB2
    
    REAL(KIND=JWRB) :: SIG2
    REAL(KIND=JWRB) :: COEF
    REAL(KIND=JWRB) :: COEF5
    REAL(KIND=JWRB) :: DFIM_SIG2
    REAL(KIND=JWRB) :: COSLP
    
    REAL(KIND=JWRB) :: XNGAMCONST
    REAL(KIND=JWRB) :: CONSTF
    REAL(KIND=JWRB) :: CONST11
    REAL(KIND=JWRB) :: CONST22
    REAL(KIND=JWRB) :: Z0VIS
    REAL(KIND=JWRB) :: Z0NOZ
    REAL(KIND=JWRB) :: FWW
    REAL(KIND=JWRB) :: PVISC
    REAL(KIND=JWRB) :: PTURB
    REAL(KIND=JWRB) :: ZCN
    REAL(KIND=JWRB) :: SIG_N
    REAL(KIND=JWRB) :: UORBT
    REAL(KIND=JWRB) :: AORB
    REAL(KIND=JWRB) :: TEMP
    REAL(KIND=JWRB) :: RE
    REAL(KIND=JWRB) :: RE_C
    REAL(KIND=JWRB) :: ZORB
    REAL(KIND=JWRB) :: CNSN
    REAL(KIND=JWRB) :: SUMF
    REAL(KIND=JWRB) :: SUMFSIN2
    REAL(KIND=JWRB) :: CSTRNFAC
    REAL(KIND=JWRB) :: FLP_AVG
    REAL(KIND=JWRB) :: SLP_AVG
    REAL(KIND=JWRB) :: ROGOROAIR
    REAL(KIND=JWRB) :: AIRD_PVISC
    REAL(KIND=JWRB) :: DSTAB1
    REAL(KIND=JWRB) :: TEMP1
    REAL(KIND=JWRB) :: TEMP2
    
    REAL(KIND=JWRB) :: XSTRESS(2)
    REAL(KIND=JWRB) :: YSTRESS(2)
    REAL(KIND=JWRB) :: FLP(2)
    REAL(KIND=JWRB) :: SLP(2)
    REAL(KIND=JWRB) :: USG2(2)
    REAL(KIND=JWRB) :: TAUX(2)
    REAL(KIND=JWRB) :: TAUY(2)
    REAL(KIND=JWRB) :: USTP(2)
    REAL(KIND=JWRB) :: USTPM1(2)
    REAL(KIND=JWRB) :: USDIRP(2)
    REAL(KIND=JWRB) :: UCN(2)
    REAL(KIND=JWRB) :: UCNZALPD(2)
    REAL(KIND=JWRB) :: GAMNORMA(2)    ! ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB) :: GAM0(2, NANG_PARAM)
    REAL(KIND=JWRB) :: DSTAB(2, NANG_PARAM)
    
    LOGICAL :: LTAUWSHELTER
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ABMAX
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ABMIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ACDLIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ALPHAMAX
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ALPHAMIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: BCDLIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: BETAMAXOXKAPPA2
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: COSTH(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: DFIM(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPSMIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPSUS
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IAB
    LOGICAL, VALUE, INTENT(IN) :: LLGCBZ0
    LOGICAL, VALUE, INTENT(IN) :: LLNORMAGAM
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: RNU
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: RNUM
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: SINTH(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SWELLF
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SWELLF2
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SWELLF3
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SWELLF4
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SWELLF5
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SWELLF6
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SWELLF7
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SWELLF7M1
    REAL(KIND=JWRB), INTENT(IN) :: SWELLFT(IAB)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: TAUWSHELTER
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: TH(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WSPMIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: XKAPPA
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: Z0RAT
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: Z0TUBMAX
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZALP
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZPI
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: ZPIFR(NFRE_loki_param)
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
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
    
    
    IF (NGST > 1) THEN
      CALL WSIGSTAR_CUF_HOIST_NEW(WSWAVE(IJ, ICHNK), UFRIC(IJ, ICHNK), Z0M(IJ, ICHNK), WSTAR(IJ, ICHNK), SIG_N, ACDLIN,  &
      & ALPHAMAX, ALPHAMIN, BCDLIN, EPSUS, G, LLGCBZ0, RNUM, WSPMIN, XKAPPA)
    END IF
    
    
    IF (LLNORMAGAM) THEN
      CSTRNFAC = CONSTN*RNFAC(IJ) / RAORW(IJ)
    END IF
    
    
    !     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
    
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
      
      DO M=1,NFRE
        SIG2 = ZPIFR(M)**2
        DFIM_SIG2 = DFIM(M)*SIG2
        
        K = 1
        TEMP = FL1(IJ, K, M, ICHNK)
        DO K=2,NANG
          TEMP = TEMP + FL1(IJ, K, M, ICHNK)
        END DO
        
        UORBT = UORBT + DFIM_SIG2*TEMP
        AORB = AORB + DFIM(M)*TEMP
      END DO
      
      UORBT = 2.0_JWRB*SQRT(UORBT)        ! this is the significant orbital amplitude
      AORB = 2.0_JWRB*SQRT(AORB)        ! this 1/2 Hs
      RE = FACM1_NU_AIR*UORBT*AORB        ! this is the Reynolds number
      Z0VIS = FAC_NU_AIR / MAX(UFRIC(IJ, ICHNK), 0.0001_JWRB)
      Z0TUB = Z0RAT*MIN(Z0TUBMAX, Z0M(IJ, ICHNK))
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
      USTP(1) = UFRIC(IJ, ICHNK)
    ELSE
      USTP(1) = UFRIC(IJ, ICHNK)*(1.0_JWRB + SIG_N)
      USTP(2) = UFRIC(IJ, ICHNK)*(1.0_JWRB - SIG_N)
    END IF
    
    DO IGST=1,NGST
      USTPM1(IGST) = 1.0_JWRB / MAX(USTP(IGST), EPSUS)
    END DO
    
    IF (LTAUWSHELTER) THEN
      DO IGST=1,NGST
        XSTRESS(IGST) = 0.0_JWRB
        YSTRESS(IGST) = 0.0_JWRB
        USG2(IGST) = USTP(IGST)**2
        TAUX(IGST) = USG2(IGST)*SIN(WDWAVE(IJ, ICHNK))
        TAUY(IGST) = USG2(IGST)*COS(WDWAVE(IJ, ICHNK))
      END DO
      
      ROGOROAIR = G / RAORW(IJ)
    END IF
    
    
    !*    2. MAIN LOOP OVER FREQUENCIES.
    !        ---------------------------
    
    IF (.not.LLNORMAGAM) THEN
      DO IGST=1,NGST
        GAMNORMA(IGST) = 1.0_JWRB
      END DO
    END IF
    
    IF (.not.LLSNEG) THEN
      DO K=1,NANG
        DO IGST=1,NGST
          DSTAB(IGST, K) = 0.0_JWRB
        END DO
      END DO
    END IF
    
    DO M=1,NFRE
      
      IF (LTAUWSHELTER) THEN
        DO IGST=1,NGST
          TAUPX = TAUX(IGST) - ABS_TAUWSHELTER*XSTRESS(IGST)
          TAUPY = TAUY(IGST) - ABS_TAUWSHELTER*YSTRESS(IGST)
          USDIRP(IGST) = ATAN2(TAUPX, TAUPY)
          USTP(IGST) = (TAUPX**2 + TAUPY**2)**0.25_JWRB
          USTPM1(IGST) = 1.0_JWRB / MAX(USTP(IGST), EPSUS)
        END DO
        
        CONSTF = ROGOROAIR*CINV(IJ, M, ICHNK)*DFIM(M)
      END IF
      
      
      !*      PRECALCULATE FREQUENCY DEPENDENCE.
      !       ----------------------------------
      
      DO IGST=1,NGST
        UCN(IGST) = USTP(IGST)*CINV(IJ, M, ICHNK)
        UCNZALPD(IGST) = XKAPPA / (UCN(IGST) + ZALP)
      END DO
      ZCN = LOG(WAVNUM(IJ, M, ICHNK)*Z0M(IJ, ICHNK))
      CNSN = ZPIFR(M)*CONST1*RAORW(IJ)
      
      !*    2.1 LOOP OVER DIRECTIONS.
      !         ---------------------
      
      DO K=1,NANG
        XLLWS(IJ, K, M, ICHNK) = 0.0_JWRB
      END DO
      
      IF (LLSNEG) THEN
        !       SWELL DAMPING:
        
        SIG2 = ZPIFR(M)**2
        DFIM_SIG2 = DFIM(M)*SIG2
        
        COEF = -SWELLF*16._JWRB*SIG2 / G
        COEF5 = -SWELLF5*2._JWRB*SQRT(2._JWRB*NU_AIR*ZPIFR(M))
        
        DSTAB1 = COEF5*AIRD_PVISC*WAVNUM(IJ, M, ICHNK)
        TEMP1 = COEF*RAORW(IJ)
      END IF
      
      DO K=1,NANG
        DO IGST=1,NGST
          
          SUMF = 0.0_JWRB
          SUMFSIN2 = 0.0_JWRB
          
          IF (LTAUWSHELTER) THEN
            COSLP = COS(TH(K) - USDIRP(IGST))
          ELSE
            COSLP = COSWDIF(IJ, K)
          END IF
          
          GAM0(IGST, K) = 0._JWRB
          IF (COSLP > 0.01_JWRB) THEN
            X = COSLP*UCN(IGST)
            ZLOG = ZCN + UCNZALPD(IGST) / COSLP
            IF (ZLOG < 0.0_JWRB) THEN
              ZLOG2X = ZLOG*ZLOG*X
              GAM0(IGST, K) = EXP(ZLOG)*ZLOG2X*ZLOG2X*CNSN
              XLLWS(IJ, K, M, ICHNK) = 1.0_JWRB
            END IF
          END IF
          
          IF (LLSNEG) THEN
            DSTAB2 = TEMP1*(TEMP2 + (FU + FUD*COSLP)*USTP(IGST))
            DSTAB(IGST, K) = DSTAB1 + PTURB*DSTAB2
          END IF
          
          SUMF = SUMF + GAM0(IGST, K)*FL1(IJ, K, M, ICHNK)
          SUMFSIN2 = SUMFSIN2 + GAM0(IGST, K)*FL1(IJ, K, M, ICHNK)*SINWDIF2(IJ, K)
        END DO
      END DO
      
      IF (LLNORMAGAM) THEN
        
        XNGAMCONST = CSTRNFAC*XK2CG(IJ, M, ICHNK)
        DO IGST=1,NGST
          ZNZ = XNGAMCONST*USTPM1(IGST)
          GAMNORMA(IGST) = (1.0_JWRB + ZNZ*SUMFSIN2) / (1.0_JWRB + ZNZ*SUMF)
        END DO
        
      END IF
      
      
      
      !*    2.2 UPDATE THE SHELTERING STRESS (in any),
      !         AND THEN ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
      !         ---------------------------------------------------------
      
      DO K=1,NANG
        
        DO IGST=1,NGST
          ! SLP: only the positive contributions
          SLP(IGST) = GAM0(IGST, K)*GAMNORMA(IGST)
          FLP(IGST) = SLP(IGST) + DSTAB(IGST, K)
        END DO
        
        DO IGST=1,NGST
          SLP(IGST) = SLP(IGST)*FL1(IJ, K, M, ICHNK)
        END DO
        
        IF (LTAUWSHELTER) THEN
          CONST11 = CONSTF*SINTH(K)
          CONST22 = CONSTF*COSTH(K)
          DO IGST=1,NGST
            XSTRESS(IGST) = XSTRESS(IGST) + SLP(IGST)*CONST11
            YSTRESS(IGST) = YSTRESS(IGST) + SLP(IGST)*CONST22
          END DO
        END IF
        
        IGST = 1
        SLP_AVG = SLP(IGST)
        FLP_AVG = FLP(IGST)
        DO IGST=2,NGST
          SLP_AVG = SLP_AVG + SLP(IGST)
          FLP_AVG = FLP_AVG + FLP(IGST)
        END DO
        
        SPOS(IJ, K, M) = AVG_GST*SLP_AVG
        FLD(IJ, K, M) = AVG_GST*FLP_AVG
        SL(IJ, K, M) = FLD(IJ, K, M)*FL1(IJ, K, M, ICHNK)
        
      END DO
      
    END DO
    
    ! END LOOP OVER FREQUENCIES
    
    
  END SUBROUTINE SINPUT_ARD_CUF_HOIST_NEW
  ATTRIBUTES(DEVICE) SUBROUTINE SINPUT_JAN_CUF_HOIST_NEW (NGST, LLSNEG, KIJS, KIJL, FL1, WAVNUM, CINV, XK2CG, WSWAVE, UFRIC,  &
  & Z0M, COSWDIF, SINWDIF2, RAORW, WSTAR, RNFAC, FLD, SL, SPOS, XLLWS, ACDLIN, ALPHAMAX, ALPHAMIN, BCDLIN, BETAMAXOXKAPPA2,  &
  & DELTH, EPSUS, G, IDAMPING, LLGCBZ0, LLNORMAGAM, NANG, NFRE, RNUM, WSPMIN, XKAPPA, ZALP, ZPI, ZPIFR, ICHNK, NCHNK, IJ)
    ! ----------------------------------------------------------------------
    
    !**** *SINPUT_JAN* - COMPUTATION OF INPUT SOURCE FUNCTION.
    
    !     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
    
    !     OPTIMIZED BY : H. GUENTHER
    
    !     MODIFIED BY :
    !       J-R BIDLOT NOVEMBER 1995
    !       J-R BIDLOT FEBRUARY 1996-97
    !       J-R BIDLOT FEBRUARY 1999 : INTRODUCE ICALL AND NCALL
    !       P.A.E.M. JANSSEN MAY 2000 : INTRODUCE GUSTINESS
    !       J-R BIDLOT FEBRUARY 2001 : MAKE IT FULLY IMPLICIT BY ONLY
    !                                  USING NEW STRESS AND ROUGHNESS.
    !       S. ABDALLA OCTOBER 2001:  INTRODUCTION OF VARIABLE AIR
    !                                 DENSITY AND STABILITY-DEPENDENT
    !                                 WIND GUSTINESS
    !       P.A.E.M. JANSSEN OCTOBER 2008: INTRODUCE DAMPING WHEN WAVES ARE
    !                                      RUNNING FASTER THAN THE WIND.
    !       J-R BIDLOT JANUARY 2013: SHALLOW WATER FORMULATION.
    
    !*    PURPOSE.
    !     ---------
    
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
    
    !     *CALL* *SINPUT_JAN (NGST, LLSNEG, KIJS, KIJL, FL1,
    !    &                    WAVNUM, CINV, XK2CG,
    !    &                    WDWAVE, WSWAVE, UFRIC, Z0M,
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
    !        *UFRIC* - FRICTION VELOCITY IN M/S.
    !        *Z0M*   - ROUGHNESS LENGTH IN M.
    !      *COSWDIF* - COS(TH(K)-WDWAVE(IJ))
    !     *SINWDIF2* - SIN(TH(K)-WDWAVE(IJ))**2
    !        *RAORW* - RATIO AIR DENSITY TO WATER DENSITY
    !        *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
    !          *FLD* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
    !           *SL* - TOTAL SOURCE FUNCTION ARRAY.
    !         *SPOS* - ONLY POSITIVE PART OF INPUT SOURCE FUNCTION ARRAY.
    !        *XLLWS* - 1 WHERE SINPUT IS POSITIVE.
    
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCE.
    
    !     EXTERNALS.
    !     ----------
    
    !       WSIGSTAR.
    
    !     MODIFICATIONS
    !     -------------
    
    !     - REMOVAL OF CALL TO CRAY SPECIFIC FUNCTIONS EXPHF AND ALOGHF
    !       BY THEIR STANDARD FORTRAN EQUIVALENT EXP and ALOGHF
    !     - MODIFIED TO MAKE INTEGRATION SCHEME FULLY IMPLICIT
    !     - INTRODUCTION OF VARIABLE AIR DENSITY
    !     - INTRODUCTION OF WIND GUSTINESS
    
    !     REFERENCE.
    !     ----------
    
    !       P. JANSSEN, J.P.O., 1989.
    !       P. JANSSEN, J.P.O., 1991
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    USE YOWFRED, ONLY: TH
    USE YOWFRED, ONLY: FR, TH
    USE YOWPARAM, ONLY: NANG_PARAM
    USE YOWPCONS, ONLY: GM1
    USE YOWTEST, ONLY: IU06
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
    INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NGST
    LOGICAL, VALUE, INTENT(IN) :: LLSNEG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    REAL(KIND=JWRB), INTENT(IN) :: FL1(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: WAVNUM(KIJL, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: CINV(KIJL, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: XK2CG(KIJL, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: WSWAVE(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: UFRIC(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: Z0M(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: COSWDIF(KIJL, NANG_loki_param)
    REAL(KIND=JWRB), INTENT(IN) :: SINWDIF2(KIJL, NANG_loki_param)
    REAL(KIND=JWRB), INTENT(IN) :: RAORW(KIJL)
    REAL(KIND=JWRB), INTENT(IN) :: RNFAC(KIJL)
    REAL(KIND=JWRB), INTENT(IN) :: WSTAR(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(OUT) :: FLD(KIJL, NANG_loki_param, NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(OUT) :: SL(KIJL, NANG_loki_param, NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(OUT) :: SPOS(KIJL, NANG_loki_param, NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(OUT) :: XLLWS(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: IG
    INTEGER(KIND=JWIM) :: K
    INTEGER(KIND=JWIM) :: M
    INTEGER(KIND=JWIM) :: IGST
    
    REAL(KIND=JWRB) :: CONST1
    REAL(KIND=JWRB) :: CONST3
    REAL(KIND=JWRB) :: XKAPPAD
    REAL(KIND=JWRB) :: CONSTN
    REAL(KIND=JWRB) :: ZNZ
    REAL(KIND=JWRB) :: X
    REAL(KIND=JWRB) :: ZLOG
    REAL(KIND=JWRB) :: ZLOG2X
    REAL(KIND=JWRB) :: ZBETA
    REAL(KIND=JWRB) :: TEMPD
    
    REAL(KIND=JWRB) :: WSIN(2)
    REAL(KIND=JWRB) :: ZTANHKD
    REAL(KIND=JWRB) :: SIG_N
    REAL(KIND=JWRB) :: CNSN
    REAL(KIND=JWRB) :: SUMF
    REAL(KIND=JWRB) :: SUMFSIN2
    REAL(KIND=JWRB) :: CSTRNFAC
    REAL(KIND=JWRB) :: UFAC1
    REAL(KIND=JWRB) :: UFAC2
    REAL(KIND=JWRB) :: GAMNORMA(2)    ! ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB) :: SIGDEV(2)
    REAL(KIND=JWRB) :: US(2)
    REAL(KIND=JWRB) :: Z0(2)
    REAL(KIND=JWRB) :: UCN(2)
    REAL(KIND=JWRB) :: ZCN(2)
    REAL(KIND=JWRB) :: USTPM1(2)
    REAL(KIND=JWRB) :: XVD(2)
    REAL(KIND=JWRB) :: UCND(2)
    REAL(KIND=JWRB) :: CONST3_UCN2(2)
    REAL(KIND=JWRB) :: GAM0(2, NANG_PARAM)
    
    LOGICAL :: LZ
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ACDLIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ALPHAMAX
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ALPHAMIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: BCDLIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: BETAMAXOXKAPPA2
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPSUS
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IDAMPING
    LOGICAL, VALUE, INTENT(IN) :: LLGCBZ0
    LOGICAL, VALUE, INTENT(IN) :: LLNORMAGAM
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: RNUM
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WSPMIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: XKAPPA
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZALP
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZPI
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: ZPIFR(NFRE_loki_param)
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    
    ! ----------------------------------------------------------------------
    
    
    CONST1 = BETAMAXOXKAPPA2
    CONST3 = 2.0_JWRB*XKAPPA / CONST1      ! SEE IDAMPING
    XKAPPAD = 1.E0_JWRB / XKAPPA
    
    CONST3 = IDAMPING*CONST3
    
    CONSTN = DELTH / (XKAPPA*ZPI)
    
    !     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
    
    IF (NGST > 1) THEN
      CALL WSIGSTAR_CUF_HOIST_NEW(WSWAVE(IJ, ICHNK), UFRIC(IJ, ICHNK), Z0M(IJ, ICHNK), WSTAR(IJ, ICHNK), SIG_N, ACDLIN,  &
      & ALPHAMAX, ALPHAMIN, BCDLIN, EPSUS, G, LLGCBZ0, RNUM, WSPMIN, XKAPPA)
    END IF
    
    !     DEFINE WHERE SINPUT WILL BE EVALUATED IN RELATIVE TERM WRT USTAR
    !     DEFINE ALSO THE RELATIVE WEIGHT OF EACH.
    
    IF (NGST == 1) THEN
      WSIN(1) = 1.0_JWRB
      SIGDEV(1) = 1.0_JWRB
    ELSE
      WSIN(1) = 0.5_JWRB
      WSIN(2) = 0.5_JWRB
      SIGDEV(1) = 1.0_JWRB - SIG_N
      SIGDEV(2) = 1.0_JWRB + SIG_N
    END IF
    
    
    IF (NGST == 1) THEN
      US(1) = UFRIC(IJ, ICHNK)
      Z0(1) = Z0M(IJ, ICHNK)
    ELSE
      DO IGST=1,NGST
        US(IGST) = UFRIC(IJ, ICHNK)*SIGDEV(IGST)
        Z0(IGST) = Z0M(IJ, ICHNK)
      END DO
    END IF
    
    DO IGST=1,NGST
      USTPM1(IGST) = 1.0_JWRB / MAX(US(IGST), EPSUS)
    END DO
    
    ! ----------------------------------------------------------------------
    
    !*    2. LOOP OVER FREQUENCIES.
    !        ----------------------
    
    DO M=1,NFRE
      
      !*      PRECALCULATE FREQUENCY DEPENDENCE.
      !       ----------------------------------
      
      ZTANHKD = ZPIFR(M)**2 / (G*WAVNUM(IJ, M, ICHNK))
      CNSN = CONST1*ZPIFR(M)*ZTANHKD*RAORW(IJ)
      
      DO IGST=1,NGST
        UCN(IGST) = US(IGST)*CINV(IJ, M, ICHNK) + ZALP
        CONST3_UCN2(IGST) = CONST3*UCN(IGST)**2
        UCND(IGST) = 1.0_JWRB / UCN(IGST)
        ZCN(IGST) = LOG(WAVNUM(IJ, M, ICHNK)*Z0(IGST))
        XVD(IGST) = 1.0_JWRB / (-US(IGST)*XKAPPAD*ZCN(IGST)*CINV(IJ, M, ICHNK))
      END DO
      
      !*    2.1 LOOP OVER DIRECTIONS.
      !         ---------------------
      
      !       WIND INPUT:
      DO K=1,NANG
        XLLWS(IJ, K, M, ICHNK) = 0.0_JWRB
        
        DO IGST=1,NGST
          
          IF (COSWDIF(IJ, K) > 0.01_JWRB) THEN
            LZ = .true.
            TEMPD = XKAPPA / COSWDIF(IJ, K)
          ELSE
            LZ = .false.
            TEMPD = XKAPPA
          END IF
          
          GAM0(IGST, K) = 0.0_JWRB
          IF (LZ) THEN
            ZLOG = ZCN(IGST) + TEMPD*UCND(IGST)
            IF (ZLOG < 0.0_JWRB) THEN
              X = COSWDIF(IJ, K)*UCN(IGST)
              ZLOG2X = ZLOG*ZLOG*X
              GAM0(IGST, K) = ZLOG2X*ZLOG2X*EXP(ZLOG)*CNSN
              XLLWS(IJ, K, M, ICHNK) = 1.0_JWRB
            END IF
          END IF
        END DO
        
      END DO
      
      
      IF (LLNORMAGAM) THEN
        
        SUMF = 0.0_JWRB
        SUMFSIN2 = 0.0_JWRB
        DO K=1,NANG
          DO IGST=1,NGST
            SUMF = SUMF + GAM0(IGST, K)*FL1(IJ, K, M, ICHNK)
            SUMFSIN2 = SUMFSIN2 + GAM0(IGST, K)*FL1(IJ, K, M, ICHNK)*SINWDIF2(IJ, K)
          END DO
          
          CSTRNFAC = CONSTN*RNFAC(IJ) / RAORW(IJ)
          ZNZ = CSTRNFAC*XK2CG(IJ, M, ICHNK)*USTPM1(IGST)
          GAMNORMA(IGST) = (1.0_JWRB + ZNZ*SUMFSIN2) / (1.0_JWRB + ZNZ*SUMF)
          
        END DO
      ELSE
        DO IGST=1,NGST
          GAMNORMA(IGST) = 1.0_JWRB
        END DO
      END IF
      
      DO K=1,NANG
        UFAC1 = WSIN(1)*GAM0(1, K)*GAMNORMA(1)
        DO IGST=2,NGST
          UFAC1 = UFAC1 + WSIN(IGST)*GAM0(IGST, K)*GAMNORMA(IGST)
        END DO
        
        UFAC2 = 0.0_JWRB
        IF (LLSNEG) THEN
          !         SWELL DAMPING:
          ZBETA = CONST3_UCN2(1)*(COSWDIF(IJ, K) - XVD(1))
          UFAC2 = WSIN(1)*ZBETA
          DO IGST=2,NGST
            ZBETA = CONST3_UCN2(IGST)*(COSWDIF(IJ, K) - XVD(IGST))
            UFAC2 = UFAC2 + WSIN(IGST)*ZBETA
          END DO
        END IF
        
        FLD(IJ, K, M) = UFAC1 + UFAC2*CNSN
        SPOS(IJ, K, M) = UFAC1*FL1(IJ, K, M, ICHNK)
        SL(IJ, K, M) = FLD(IJ, K, M)*FL1(IJ, K, M, ICHNK)
      END DO
      
      !*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
      !         ------------------------------------------------
      
    END DO
    
    
    
  END SUBROUTINE SINPUT_JAN_CUF_HOIST_NEW
END MODULE SINPUT_ARD_CUF_HOIST_NEW_MOD
