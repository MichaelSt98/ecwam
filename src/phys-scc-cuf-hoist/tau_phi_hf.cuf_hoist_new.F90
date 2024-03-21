! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE TAU_PHI_HF_CUF_HOIST_NEW_MOD
  !CONTAINED SUBROUTINES:
  ! - OMEGAGC
  ! - TAU_PHI_HF
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE OMEGAGC_CUF_HOIST_NEW (UST, NS, XKS, OMS, NWAV_GC, OMEGA_GC, SQRTGOSURFT, XKM_GC, XK_GC,  &
  & XLOGKRATIOM1_GC)
    
    !***  DETERMINE THE CUT-OFF ANGULAR FREQUENCY FOR THE GRAV-CAPILLARY WAVES
    !     !!!! rounded to the closest index of XK_GC  !!!!!
    
    !     AUTHOR: PETER JANSSEN
    !     ------
    
    !     REFERENCES:
    !     ----------
    
    !     VIERS PAPER EQ.(29)
    
    !----------------------------------------------------------------------
    
    USE NS_GC_CUF_HOIST_NEW_MOD, ONLY: NS_GC_CUF_HOIST_NEW
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
    REAL(KIND=JWRB), INTENT(IN) :: UST
    INTEGER(KIND=JWIM), INTENT(OUT) :: NS    ! index in array XK_GC corresponding to XKS and OMS
    REAL(KIND=JWRB), INTENT(OUT) :: XKS    ! cut-off wave number
    REAL(KIND=JWRB), INTENT(OUT) :: OMS    ! cut-off angular frequency
    
    INTEGER(KIND=JWIM), INTENT(IN) :: NWAV_GC
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: OMEGA_GC(NWAV_GC)
    REAL(KIND=JWRB), INTENT(IN) :: SQRTGOSURFT
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: XKM_GC(NWAV_GC)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: XK_GC(NWAV_GC)
    REAL(KIND=JWRB), INTENT(IN) :: XLOGKRATIOM1_GC
!$acc routine seq
    
    
    ! ----------------------------------------------------------------------
    
    
    NS = NS_GC_CUF_HOIST_NEW(UST, NWAV_GC, SQRTGOSURFT, XKM_GC, XLOGKRATIOM1_GC)
    XKS = XK_GC(NS)
    OMS = OMEGA_GC(NS)
    
    
  END SUBROUTINE OMEGAGC_CUF_HOIST_NEW
  
  ATTRIBUTES(DEVICE) SUBROUTINE TAU_PHI_HF_CUF_HOIST_NEW (KIJS, KIJL, MIJ, LTAUWSHELTER, UFRIC, Z0M, FL1, AIRD, RNFAC, COSWDIF,  &
  & SINWDIF2, UST, TAUHF, PHIHF, LLPHIHF, DELTH, FR5, G, GAMNCONST, GM1, JTOT_TAUHF, LLGCBZ0, LLNORMAGAM, NANG, NWAV_GC,  &
  & OMEGA_GC, SQRTGOSURFT, TAUWSHELTER, WTAUHF, X0TAUHF, XKAPPA, XKM_GC, XK_GC, XLOGKRATIOM1_GC, ZALP, ZPI4GM1, ZPI4GM2, ZPIFR,  &
  & ICHNK, NCHNK, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *TAU_PHI_HF* - COMPUTATION OF HIGH-FREQUENCY STRESS.
    !                                   HIGH-FREQUENCY ENERGY FLUX.
    
    !     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90
    !     JEAN BIDLOT  ECMWF  JANUARY 2017
    
    !*    PURPOSE.
    !     ---------
    
    !       COMPUTE HIGH-FREQUENCY WAVE STRESS AND ENERGY FLUX
    
    !**   INTERFACE.
    !     ---------
    
    !       *CALL* *TAU_PHI_HF(KIJS, KIJL, MIJ, LTAUWSHELTER, UFRIC, UST, Z0M,
    !                          FL1, AIRD, RNFAC,
    !                          COSWDIF, SINWDIF2,
    !                          UST, TAUHF, PHIHF, LLPHIHF)
    !          *KIJS*         - INDEX OF FIRST GRIDPOINT
    !          *KIJL*         - INDEX OF LAST GRIDPOINT
    !          *MIJ*          - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
    !          *LTAUWSHELTER* - if true then TAUWSHELTER
    !          *FL1*          - WAVE SPECTRUM.
    !          *AIRD*         - AIR DENSITY IN KG/M**3.
    !          *RNFAC*        - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !          *UFRIC*        - FRICTION VELOCITY
    !          *COSWDIF*      - COS(TH(K)-WDWAVE(IJ))
    !          *SINWDIF2*     - SIN(TH(K)-WDWAVE(IJ))**2
    !          *UST*          - REDUCED FRICTION VELOCITY DUE TO SHELTERING
    !          *Z0M*          - ROUGHNESS LENGTH
    !          *TAUHF*        - HIGH-FREQUENCY STRESS
    !          *PHIHF*        - HIGH-FREQUENCY ENERGY FLUX INTO OCEAN
    !          *LLPHIHF*      - TRUE IF PHIHF NEEDS TO COMPUTED
    
    
    !     METHOD.
    !     -------
    
    !       IT NEEDS A CALL TO INIT_X0TAUHF TO INITIALISE
    !       SEE REFERENCE FOR WAVE STRESS CALCULATION.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    USE YOWFRED, ONLY: TH
    USE YOWPARAM, ONLY: NFRE
    USE YOWPCONS, ONLY: ZPI
    USE YOWTEST, ONLY: IU06
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
    INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    INTEGER(KIND=JWIM), INTENT(IN) :: MIJ(KIJL, NCHNK)
    LOGICAL, VALUE, INTENT(IN) :: LTAUWSHELTER
    REAL(KIND=JWRB), INTENT(IN) :: UFRIC(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: Z0M(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: FL1(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: AIRD(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: RNFAC(KIJL)
    REAL(KIND=JWRB), INTENT(IN) :: COSWDIF(KIJL, NANG_loki_param)
    REAL(KIND=JWRB), INTENT(IN) :: SINWDIF2(KIJL, NANG_loki_param)
    REAL(KIND=JWRB), INTENT(INOUT) :: UST(KIJL)
    REAL(KIND=JWRB), INTENT(OUT) :: TAUHF(KIJL)
    REAL(KIND=JWRB), INTENT(OUT) :: PHIHF(KIJL)
    LOGICAL, VALUE, INTENT(IN) :: LLPHIHF
    
    
    INTEGER(KIND=JWIM) :: J
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: K
    INTEGER(KIND=JWIM) :: NS
    
    REAL(KIND=JWRB), PARAMETER :: ZSUPMAX = 0.0_JWRB    !  LOG(1.)
    REAL(KIND=JWRB) :: OMEGA
    REAL(KIND=JWRB) :: OMEGACC
    REAL(KIND=JWRB) :: X0G
    REAL(KIND=JWRB) :: YC
    REAL(KIND=JWRB) :: Y
    REAL(KIND=JWRB) :: CM1
    REAL(KIND=JWRB) :: ZX
    REAL(KIND=JWRB) :: ZARG
    REAL(KIND=JWRB) :: ZLOG
    REAL(KIND=JWRB) :: ZBETA
    REAL(KIND=JWRB) :: FNC
    REAL(KIND=JWRB) :: FNC2
    REAL(KIND=JWRB) :: GAMNORMA    ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB) :: ZNZ
    REAL(KIND=JWRB) :: CONFG
    REAL(KIND=JWRB) :: COSW
    REAL(KIND=JWRB) :: FCOSW2
    
    REAL(KIND=JWRB) :: XKS
    REAL(KIND=JWRB) :: OMS
    REAL(KIND=JWRB) :: SQRTZ0OG
    REAL(KIND=JWRB) :: ZSUP
    REAL(KIND=JWRB) :: ZINF
    REAL(KIND=JWRB) :: DELZ
    REAL(KIND=JWRB) :: TAUL
    REAL(KIND=JWRB) :: XLOGGZ0
    REAL(KIND=JWRB) :: SQRTGZ0
    REAL(KIND=JWRB) :: USTPH
    REAL(KIND=JWRB) :: CONST1
    REAL(KIND=JWRB) :: CONST2
    REAL(KIND=JWRB) :: CONSTTAU
    REAL(KIND=JWRB) :: CONSTPHI
    REAL(KIND=JWRB) :: F1DCOS2
    REAL(KIND=JWRB) :: F1DCOS3
    REAL(KIND=JWRB) :: F1D
    REAL(KIND=JWRB) :: F1DSIN2
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: FR5(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: GAMNCONST
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: GM1
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: JTOT_TAUHF
    LOGICAL, VALUE, INTENT(IN) :: LLGCBZ0
    LOGICAL, VALUE, INTENT(IN) :: LLNORMAGAM
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NWAV_GC
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: OMEGA_GC(NWAV_GC)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SQRTGOSURFT
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: TAUWSHELTER
    REAL(KIND=JWRB), INTENT(IN) :: WTAUHF(JTOT_TAUHF)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: X0TAUHF
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: XKAPPA
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: XKM_GC(NWAV_GC)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: XK_GC(NWAV_GC)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: XLOGKRATIOM1_GC
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZALP
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZPI4GM1
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZPI4GM2
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: ZPIFR(NFRE_loki_param)
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    
    ! ----------------------------------------------------------------------
    
    
    
    IF (LLGCBZ0) THEN
      CALL OMEGAGC_CUF_HOIST_NEW(UFRIC(IJ, ICHNK), NS, XKS, OMS, NWAV_GC, OMEGA_GC(:), SQRTGOSURFT, XKM_GC(:), XK_GC(:),  &
      & XLOGKRATIOM1_GC)
    END IF
    
    !     See INIT_X0TAUHF
    X0G = X0TAUHF*G
    
    IF (LLPHIHF) USTPH = UST(IJ)
    
    !*    COMPUTE THE INTEGRALS
    !     ---------------------
    
    XLOGGZ0 = LOG(G*Z0M(IJ, ICHNK))
    OMEGACC = MAX(ZPIFR(MIJ(IJ, ICHNK)), X0G / UST(IJ))
    SQRTZ0OG = SQRT(Z0M(IJ, ICHNK)*GM1)
    SQRTGZ0 = 1.0_JWRB / SQRTZ0OG
    YC = OMEGACC*SQRTZ0OG
    ZINF = LOG(YC)
    
    CONSTTAU = ZPI4GM2*FR5(MIJ(IJ, ICHNK))
    
    K = 1
    COSW = MAX(COSWDIF(IJ, K), 0.0_JWRB)
    FCOSW2 = FL1(IJ, K, MIJ(IJ, ICHNK), ICHNK)*COSW**2
    F1DCOS3 = FCOSW2*COSW
    F1DCOS2 = FCOSW2
    F1DSIN2 = FL1(IJ, K, MIJ(IJ, ICHNK), ICHNK)*SINWDIF2(IJ, K)
    F1D = FL1(IJ, K, MIJ(IJ, ICHNK), ICHNK)
    DO K=2,NANG
      COSW = MAX(COSWDIF(IJ, K), 0.0_JWRB)
      FCOSW2 = FL1(IJ, K, MIJ(IJ, ICHNK), ICHNK)*COSW**2
      F1DCOS3 = F1DCOS3 + FCOSW2*COSW
      F1DCOS2 = F1DCOS2 + FCOSW2
      F1DSIN2 = F1DSIN2 + FL1(IJ, K, MIJ(IJ, ICHNK), ICHNK)*SINWDIF2(IJ, K)
      F1D = F1D + FL1(IJ, K, MIJ(IJ, ICHNK), ICHNK)
    END DO
    F1DCOS3 = DELTH*F1DCOS3
    F1DCOS2 = DELTH*F1DCOS2
    F1DSIN2 = DELTH*F1DSIN2
    F1D = DELTH*F1D
    
    IF (LLNORMAGAM) THEN
      CONFG = GAMNCONST*FR5(MIJ(IJ, ICHNK))*RNFAC(IJ)*SQRTGZ0
      CONST1 = CONFG*F1DSIN2
      CONST2 = CONFG*F1D
    ELSE
      CONST1 = 0.0_JWRB
      CONST2 = 0.0_JWRB
    END IF
    
    
    !     TAUHF :
    IF (LLGCBZ0) THEN
      ZSUP = MIN(LOG(OMS*SQRTZ0OG), ZSUPMAX)
    ELSE
      ZSUP = ZSUPMAX
    END IF
    
    TAUL = UST(IJ)**2
    DELZ = MAX((ZSUP - ZINF) / REAL(JTOT_TAUHF - 1, kind=JWRB), 0.0_JWRB)
    TAUHF(IJ) = 0.0_JWRB
    
    ! Intergrals are integrated following a change of variable : Z=LOG(Y)
    IF (LTAUWSHELTER) THEN
      DO J=1,JTOT_TAUHF
        Y = EXP(ZINF + REAL(J - 1, kind=JWRB)*DELZ)
        OMEGA = Y*SQRTGZ0
        CM1 = OMEGA*GM1
        ZX = UST(IJ)*CM1 + ZALP
        ZARG = XKAPPA / ZX
        ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
        ZLOG = MIN(ZLOG, 0.0_JWRB)
        ZBETA = ZLOG**4*EXP(ZLOG)
        ZNZ = ZBETA*UST(IJ)*Y
        GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
        FNC2 = F1DCOS3*CONSTTAU*ZBETA*TAUL*WTAUHF(J)*DELZ*GAMNORMA
        TAUL = MAX(TAUL - TAUWSHELTER*FNC2, 0.0_JWRB)
        
        UST(IJ) = SQRT(TAUL)
        TAUHF(IJ) = TAUHF(IJ) + FNC2
      END DO
    ELSE
      DO J=1,JTOT_TAUHF
        Y = EXP(ZINF + REAL(J - 1, kind=JWRB)*DELZ)
        OMEGA = Y*SQRTGZ0
        CM1 = OMEGA*GM1
        ZX = UST(IJ)*CM1 + ZALP
        ZARG = XKAPPA / ZX
        ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
        ZLOG = MIN(ZLOG, 0.0_JWRB)
        ZBETA = ZLOG**4*EXP(ZLOG)
        FNC2 = ZBETA*WTAUHF(J)
        ZNZ = ZBETA*UST(IJ)*Y
        GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
        TAUHF(IJ) = TAUHF(IJ) + FNC2*GAMNORMA
      END DO
      TAUHF(IJ) = F1DCOS3*CONSTTAU*TAUL*TAUHF(IJ)*DELZ
    END IF
    
    
    PHIHF(IJ) = 0.0_JWRB
    IF (LLPHIHF) THEN
      !       PHIHF:
      !       We are neglecting the gravity-capillary contribution
      !       Recompute DELZ over the full interval
      TAUL = USTPH**2
      ZSUP = ZSUPMAX
      DELZ = MAX((ZSUP - ZINF) / REAL(JTOT_TAUHF - 1, kind=JWRB), 0.0_JWRB)
      
      CONSTPHI = AIRD(IJ, ICHNK)*ZPI4GM1*FR5(MIJ(IJ, ICHNK))
      
      ! Intergrals are integrated following a change of variable : Z=LOG(Y)
      IF (LTAUWSHELTER) THEN
        DO J=1,JTOT_TAUHF
          Y = EXP(ZINF + REAL(J - 1, kind=JWRB)*DELZ)
          OMEGA = Y*SQRTGZ0
          CM1 = OMEGA*GM1
          ZX = USTPH*CM1 + ZALP
          ZARG = XKAPPA / ZX
          ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
          ZLOG = MIN(ZLOG, 0.0_JWRB)
          ZBETA = ZLOG**4*EXP(ZLOG)
          ZNZ = ZBETA*UST(IJ)*Y
          GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
          FNC2 = ZBETA*TAUL*WTAUHF(J)*DELZ*GAMNORMA
          TAUL = MAX(TAUL - TAUWSHELTER*F1DCOS3*CONSTTAU*FNC2, 0.0_JWRB)
          USTPH = SQRT(TAUL)
          PHIHF(IJ) = PHIHF(IJ) + FNC2 / Y
        END DO
        PHIHF(IJ) = F1DCOS2*CONSTPHI*SQRTZ0OG*PHIHF(IJ)
      ELSE
        DO J=1,JTOT_TAUHF
          Y = EXP(ZINF + REAL(J - 1, kind=JWRB)*DELZ)
          OMEGA = Y*SQRTGZ0
          CM1 = OMEGA*GM1
          ZX = USTPH*CM1 + ZALP
          ZARG = XKAPPA / ZX
          ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
          ZLOG = MIN(ZLOG, 0.0_JWRB)
          ZBETA = ZLOG**4*EXP(ZLOG)
          ZNZ = ZBETA*UST(IJ)*Y
          GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
          FNC2 = ZBETA*WTAUHF(J)*GAMNORMA
          PHIHF(IJ) = PHIHF(IJ) + FNC2 / Y
        END DO
        PHIHF(IJ) = F1DCOS2*CONSTPHI*SQRTZ0OG*TAUL*PHIHF(IJ)*DELZ
      END IF
    END IF
    
    
    
  END SUBROUTINE TAU_PHI_HF_CUF_HOIST_NEW
END MODULE TAU_PHI_HF_CUF_HOIST_NEW_MOD
MODULE MEANSQS_GC_CUF_HOIST_NEW_MOD
  !CONTAINED SUBROUTINES:
  ! - MEANSQS_GC
  CONTAINS
  SUBROUTINE MEANSQS_GC (XKMSS, KIJS, KIJL, HALP, USTAR, XMSSCG, FRGC)
    
    !***  DETERMINE MSS FOR GRAV-CAP WAVES UP TO WAVE NUMBER XKMSS
    
    !     AUTHOR: PETER JANSSEN
    !     ------
    
    !     REFERENCES:
    !     ----------
    
    !     VIERS PAPER EQ.(29)
    
    !----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: NWAV_GC, XLOGKRATIOM1_GC, XKM_GC, VG_GC, C2OSQRTVG_GC, DELKCC_GC, DELKCC_GC_NS
    USE YOWPCONS, ONLY: G, ZPI, SURFT
    
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE TAU_PHI_HF_MOD, ONLY: OMEGAGC
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    
    REAL(KIND=JWRB), INTENT(IN) :: XKMSS    ! WAVE NUMBER CUT-OFF
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: HALP    ! 1/2 Phillips parameter
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: USTAR    ! friction velocity
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: XMSSCG    ! mean square slope for gravity-capillary waves
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: FRGC    ! Frequency from which the gravity-capillary spectrum is approximated
    
    
    INTEGER(KIND=JWIM) :: IJ, I, NE
    INTEGER(KIND=JWIM), DIMENSION(KIJL) :: NS
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    REAL(KIND=JWRB), DIMENSION(KIJL) :: XKS, OMS, COEF
    
    !     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"
    
    ! ----------------------------------------------------------------------
    
    IF (LHOOK) CALL DR_HOOK('MEANSQS_GC', 0, ZHOOK_HANDLE)
    
    NE = MIN(MAX(NINT(LOG(XKMSS*XKM_GC(1))*XLOGKRATIOM1_GC), 1), NWAV_GC)
    
    DO IJ=KIJS,KIJL
      CALL OMEGAGC(USTAR(IJ), NS(IJ), XKS(IJ), OMS(IJ))
      FRGC(IJ) = OMS(IJ) / ZPI
      IF (XKS(IJ) > XKMSS) THEN
        NS(IJ) = NE
        XMSSCG(IJ) = 0.0_JWRB
      ELSE
        XMSSCG(IJ) = DELKCC_GC_NS(NS(IJ))*XKM_GC(NS(IJ))
      END IF
    END DO
    
    DO IJ=KIJS,KIJL
      DO I=NS(IJ) + 1,NE
        !         ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
        !         BB = COEF(IJ)*SQRT(VG_GC(I))/C_GC(I)**2
        !         mss :  integral of k**2 F(k)  k dk
        XMSSCG(IJ) = XMSSCG(IJ) + DELKCC_GC(I)*XKM_GC(I)
      END DO
      COEF(IJ) = C2OSQRTVG_GC(NS(IJ))*HALP(IJ)
      XMSSCG(IJ) = XMSSCG(IJ)*COEF(IJ)
    END DO
    
    IF (LHOOK) CALL DR_HOOK('MEANSQS_GC', 1, ZHOOK_HANDLE)
    
  END SUBROUTINE MEANSQS_GC
END MODULE MEANSQS_GC_CUF_HOIST_NEW_MOD
