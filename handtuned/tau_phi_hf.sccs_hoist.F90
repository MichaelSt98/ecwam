! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE TAU_PHI_HF_LOKI_MOD
  CONTAINS
  SUBROUTINE TAU_PHI_HF_LOKI (KIJS, KIJL, MIJ, LTAUWSHELTER, UFRIC, Z0M, FL1, AIRD, RNFAC, COSWDIF, SINWDIF2, UST, TAUHF, PHIHF,  &
  & LLPHIHF, IJ)

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
    !          *LLPHIHF*      - TRUE IF PHIHF NEEDS TO BE COMPUTED


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

    USE NS_GC_LOKI_MOD, ONLY: NS_GC_LOKI
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWCOUP, ONLY: X0TAUHF, JTOT_TAUHF, WTAUHF, LLGCBZ0, LLNORMAGAM
    USE YOWFRED, ONLY: ZPIFR, FR5, TH, DELTH, OMEGA_GC, XK_GC
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: G, GM1, ZPI, ZPI4GM1, ZPI4GM2
    USE YOWPHYS, ONLY: ZALP, XKAPPA, TAUWSHELTER, GAMNCONST
    USE YOWTEST, ONLY: IU06

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    INTEGER(KIND=JWIM), INTENT(IN) :: MIJ(KIJL)
    LOGICAL, INTENT(IN) :: LTAUWSHELTER
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: UFRIC, Z0M
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: AIRD, RNFAC
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF, SINWDIF2
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: UST
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: TAUHF, PHIHF
    LOGICAL, INTENT(IN) :: LLPHIHF


    INTEGER(KIND=JWIM) :: TAU_PHI_HF_J, TAU_PHI_HF_K
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: NS

    REAL(KIND=JWRB), PARAMETER :: ZSUPMAX = 0.0_JWRB    !  LOG(1.)
    REAL(KIND=JWRB) :: OMEGA, OMEGACC
    REAL(KIND=JWRB) :: X0G
    REAL(KIND=JWRB) :: YC, Y, CM1, ZX, ZARG, ZLOG, ZBETA
    REAL(KIND=JWRB) :: FNC, FNC2
    REAL(KIND=JWRB) :: GAMNORMA    ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB) :: ZNZ, CONFG
    REAL(KIND=JWRB) :: TAU_PHI_HF_COSW, TAU_PHI_HF_FCOSW2
    REAL(KIND=JPHOOK) :: TAU_PHI_HF_ZHOOK_HANDLE

    REAL(KIND=JWRB) :: XKS, OMS
    REAL(KIND=JWRB) :: SQRTZ0OG, ZSUP, ZINF, DELZ
    REAL(KIND=JWRB) :: TAUL, XLOGGZ0, SQRTGZ0
    REAL(KIND=JWRB) :: USTPH
    REAL(KIND=JWRB) :: CONST1, CONST2, CONSTTAU, CONSTPHI
    REAL(KIND=JWRB) :: F1DCOS2, F1DCOS3
    REAL(KIND=JWRB) :: F1D, F1DSIN2
    REAL(KIND=JPHOOK) :: OMEGAGC_ZHOOK_HANDLE
!$acc routine seq
!$acc data present( MIJ, UFRIC, Z0M, FL1, AIRD, RNFAC, COSWDIF, SINWDIF2, UST, TAUHF, PHIHF )


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

    IF (LLPHIHF) USTPH = UST(IJ)

    !*    COMPUTE THE INTEGRALS
    !     ---------------------

    XLOGGZ0 = LOG(G*Z0M(IJ))
    OMEGACC = MAX(ZPIFR(MIJ(IJ)), X0G / UST(IJ))
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
    F1D = FL1(IJ, TAU_PHI_HF_K, MIJ(IJ))
!$acc loop seq
    DO TAU_PHI_HF_K=2,NANG
      TAU_PHI_HF_COSW = MAX(COSWDIF(IJ, TAU_PHI_HF_K), 0.0_JWRB)
      TAU_PHI_HF_FCOSW2 = FL1(IJ, TAU_PHI_HF_K, MIJ(IJ))*TAU_PHI_HF_COSW**2
      F1DCOS3 = F1DCOS3 + TAU_PHI_HF_FCOSW2*TAU_PHI_HF_COSW
      F1DCOS2 = F1DCOS2 + TAU_PHI_HF_FCOSW2
      F1DSIN2 = F1DSIN2 + FL1(IJ, TAU_PHI_HF_K, MIJ(IJ))*SINWDIF2(IJ, TAU_PHI_HF_K)
      F1D = F1D + FL1(IJ, TAU_PHI_HF_K, MIJ(IJ))
    END DO
    F1DCOS3 = DELTH*F1DCOS3
    F1DCOS2 = DELTH*F1DCOS2
    F1DSIN2 = DELTH*F1DSIN2
    F1D = DELTH*F1D

    IF (LLNORMAGAM) THEN
      CONFG = GAMNCONST*FR5(MIJ(IJ))*RNFAC(IJ)*SQRTGZ0
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
!$acc loop seq
      DO TAU_PHI_HF_J=1,JTOT_TAUHF
        Y = EXP(ZINF + REAL(TAU_PHI_HF_J - 1, kind=JWRB)*DELZ)
        OMEGA = Y*SQRTGZ0
        CM1 = OMEGA*GM1
        ZX = UST(IJ)*CM1 + ZALP
        ZARG = XKAPPA / ZX
        ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
        ZLOG = MIN(ZLOG, 0.0_JWRB)
        ZBETA = ZLOG**4*EXP(ZLOG)
        ZNZ = ZBETA*UST(IJ)*Y
        GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
        FNC2 = F1DCOS3*CONSTTAU*ZBETA*TAUL*WTAUHF(TAU_PHI_HF_J)*DELZ*GAMNORMA
        TAUL = MAX(TAUL - TAUWSHELTER*FNC2, 0.0_JWRB)

        UST(IJ) = SQRT(TAUL)
        TAUHF(IJ) = TAUHF(IJ) + FNC2
      END DO
    ELSE
!$acc loop seq
      DO TAU_PHI_HF_J=1,JTOT_TAUHF
        Y = EXP(ZINF + REAL(TAU_PHI_HF_J - 1, kind=JWRB)*DELZ)
        OMEGA = Y*SQRTGZ0
        CM1 = OMEGA*GM1
        ZX = UST(IJ)*CM1 + ZALP
        ZARG = XKAPPA / ZX
        ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
        ZLOG = MIN(ZLOG, 0.0_JWRB)
        ZBETA = ZLOG**4*EXP(ZLOG)
        FNC2 = ZBETA*WTAUHF(TAU_PHI_HF_J)
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

      CONSTPHI = AIRD(IJ)*ZPI4GM1*FR5(MIJ(IJ))

      ! Intergrals are integrated following a change of variable : Z=LOG(Y)
      IF (LTAUWSHELTER) THEN
!$acc loop seq
        DO TAU_PHI_HF_J=1,JTOT_TAUHF
          Y = EXP(ZINF + REAL(TAU_PHI_HF_J - 1, kind=JWRB)*DELZ)
          OMEGA = Y*SQRTGZ0
          CM1 = OMEGA*GM1
          ZX = USTPH*CM1 + ZALP
          ZARG = XKAPPA / ZX
          ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
          ZLOG = MIN(ZLOG, 0.0_JWRB)
          ZBETA = ZLOG**4*EXP(ZLOG)
          ZNZ = ZBETA*UST(IJ)*Y
          GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
          FNC2 = ZBETA*TAUL*WTAUHF(TAU_PHI_HF_J)*DELZ*GAMNORMA
          TAUL = MAX(TAUL - TAUWSHELTER*F1DCOS3*CONSTTAU*FNC2, 0.0_JWRB)
          USTPH = SQRT(TAUL)
          PHIHF(IJ) = PHIHF(IJ) + FNC2 / Y
        END DO
        PHIHF(IJ) = F1DCOS2*CONSTPHI*SQRTZ0OG*PHIHF(IJ)
      ELSE
!$acc loop seq
        DO TAU_PHI_HF_J=1,JTOT_TAUHF
          Y = EXP(ZINF + REAL(TAU_PHI_HF_J - 1, kind=JWRB)*DELZ)
          OMEGA = Y*SQRTGZ0
          CM1 = OMEGA*GM1
          ZX = USTPH*CM1 + ZALP
          ZARG = XKAPPA / ZX
          ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
          ZLOG = MIN(ZLOG, 0.0_JWRB)
          ZBETA = ZLOG**4*EXP(ZLOG)
          ZNZ = ZBETA*UST(IJ)*Y
          GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
          FNC2 = ZBETA*WTAUHF(TAU_PHI_HF_J)*GAMNORMA
          PHIHF(IJ) = PHIHF(IJ) + FNC2 / Y
        END DO
        PHIHF(IJ) = F1DCOS2*CONSTPHI*SQRTZ0OG*TAUL*PHIHF(IJ)*DELZ
      END IF
    END IF



!$acc end data
  END SUBROUTINE TAU_PHI_HF_LOKI
END MODULE TAU_PHI_HF_LOKI_MOD
