! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDISSIP_ARD_LOKI_MOD
  CONTAINS
  SUBROUTINE SDISSIP_ARD_LOKI (KIJS, KIJL, FL1, FLD, SL, WAVNUM, CGROUP, XK2CG, UFRIC, COSWDIF, RAORW, IJ, FACSAT, FACWTRB,  &
  & TEMP1, BTH0, C_, C_C, DSIP, TRPZ_DSIP, BTH, TEMP2, D, SCUMUL, RENEWALFREQ, WCUMUL)
    ! ----------------------------------------------------------------------

    !**** *SDISSIP_ARD* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

    !     LOTFI AOUF       METEO FRANCE 2013
    !     FABRICE ARDHUIN  IFREMER  2013


    !*    PURPOSE.
    !     --------
    !       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
    !       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
    !       OF DISSIPATION SOURCE FUNCTION.

    !**   INTERFACE.
    !     ----------

    !       *CALL* *SDISSIP_ARD (KIJS, KIJL, FL1, FLD,SL,*
    !                            WAVNUM, CGROUP, XK2CG,
    !                            UFRIC, COSWDIF, RAORW)*
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *FL1*    - SPECTRUM.
    !          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*     - TOTAL SOURCE FUNCTION ARRAY
    !          *WAVNUM* - WAVE NUMBER
    !          *CGROUP* - GROUP SPEED
    !          *XK2CG*  - (WAVE NUMBER)**2 * GROUP SPEED
    !          *UFRIC*  - FRICTION VELOCITY IN M/S.
    !          *RAORW*  - RATIO AIR DENSITY TO WATER DENSITY
    !          *COSWDIF*-  COS(TH(K)-WDWAVE(IJ))


    !     METHOD.
    !     -------

    !       SEE REFERENCES.

    !     EXTERNALS.
    !     ----------

    !       NONE.

    !     REFERENCE.
    !     ----------

    !       ARDHUIN et AL. JPO DOI:10.1175/20110JPO4324.1


    ! ----------------------------------------------------------------------
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWFRED, ONLY: FR, TH, ZPIFR, FRATIO, DELTH
    USE YOWPCONS, ONLY: G, ZPI
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPHYS, ONLY: SDSBR, ISDSDTH, ISB, IPSAT, SSDSC2, SSDSC4, SSDSC6, MICHE, SSDSC3, SSDSBRF1, BRKPBCOEF, SSDSC5, NSDSNTH,  &
    & INDICESSAT, SATWEIGHTS

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM, CGROUP, XK2CG
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: UFRIC, RAORW
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF


    INTEGER(KIND=JWIM) :: K, M, I, J, M2, K2, KK, NANGD
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    !     NDIKCUMUL is the  integer difference in frequency bands
    INTEGER(KIND=JWIM) :: NDIKCUMUL

    REAL(KIND=JWRB) :: TPIINV, TPIINVH, TMP01, TMP02, TMP03
    REAL(KIND=JWRB) :: EPSR, SSDSC6M1, ZCOEF, ZCOEFM1
    REAL(KIND=JWRB) :: XLOGDFRTH, BRLAMBDA

    REAL(KIND=JPHOOK) :: SDISSIP_ARD_ZHOOK_HANDLE


    REAL(KIND=JWRB) :: SSDSC2_SIG
    REAL(KIND=JWRB) :: FACTURB
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: FACSAT
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: FACWTRB
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: TEMP1
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: BTH0    !saturation spectrum
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: C_
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: C_C
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: DSIP
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: TRPZ_DSIP
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: BTH    !saturation spectrum
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: TEMP2
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: D
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: SCUMUL
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: RENEWALFREQ
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, 0:NANG / 2, NFRE) :: WCUMUL
!$acc routine seq
!$acc data present( FL1, FLD, SL, WAVNUM, CGROUP, XK2CG, UFRIC, COSWDIF, RAORW, FACSAT, FACWTRB, TEMP1, BTH0, C_, C_C, DSIP,  &
!$acc & TRPZ_DSIP, BTH, TEMP2, D, SCUMUL, RENEWALFREQ, WCUMUL )

    ! ----------------------------------------------------------------------


    ! INITIALISATION

    EPSR = SQRT(SDSBR)

    TPIINV = 1.0_JWRB / ZPI
    TPIINVH = 0.5_JWRB*TPIINV


    TMP03 = 1.0_JWRB / (SDSBR*MICHE)

    SSDSC6M1 = 1._JWRB - SSDSC6

!$acc loop seq
    DO M=1,NFRE
      FACSAT(IJ, M) = WAVNUM(IJ, M)*TPIINV*XK2CG(IJ, M)
    END DO

    ! COMPUTE SATURATION SPECTRUM
!$acc loop seq
    DO M=1,NFRE
      BTH0(IJ, M) = 0.0_JWRB
!$acc loop seq
      DO K=1,NANG
        BTH(IJ, K, M) = 0.0_JWRB
      END DO
    END DO

!$acc loop seq
    DO M=1,NFRE
!$acc loop seq
      DO K=1,NANG
        ! integrates in directional sector
!$acc loop seq
        DO K2=1,NSDSNTH*2 + 1
          KK = INDICESSAT(K, K2)
          BTH(IJ, K, M) = BTH(IJ, K, M) + SATWEIGHTS(K, K2)*FL1(IJ, KK, M)
        END DO
        BTH(IJ, K, M) = BTH(IJ, K, M)*FACSAT(IJ, M)
        BTH0(IJ, M) = MAX(BTH0(IJ, M), BTH(IJ, K, M))
      END DO
    END DO


    ! SATURATION TERM

!$acc loop seq
    DO M=1,NFRE
      SSDSC2_SIG = SSDSC2*ZPIFR(M)
      ZCOEF = SSDSC2_SIG*SSDSC6
      ZCOEFM1 = SSDSC2_SIG*SSDSC6M1
!$acc loop seq
      DO K=1,NANG
        D(IJ, K, M) =  &
        & ZCOEF*MAX(0._JWRB, BTH0(IJ, M)*TMP03 - SSDSC4)**IPSAT + ZCOEFM1*MAX(0._JWRB, BTH(IJ, K, M)*TMP03 - SSDSC4)**IPSAT
      END DO
    END DO


    ! CUMULATIVE TERM
    IF (SSDSC3 /= 0.0_JWRB) THEN

      NANGD = NANG / 2
      XLOGDFRTH = LOG(FRATIO)*DELTH
      !       l(k,th)=1/(2*pi²)= the breaking crest density
      BRLAMBDA = BRKPBCOEF / (2.0_JWRB*ZPI**2)
      TMP02 = SSDSC3*BRLAMBDA

      !       NDIKCUMUL is the  integer difference in frequency bands
      !       between the "large breakers" and short "wiped-out waves"
      !!! wrong !!???        NDIKCUMUL = NINT(SSDSBRF1/(FRATIO-1.))
      NDIKCUMUL = NINT(-LOG(SSDSBRF1) / LOG(FRATIO))

!$acc loop seq
      DO M=1,NFRE
        C_(IJ, M) = ZPIFR(M) / WAVNUM(IJ, M)
        C_C(IJ, M) = C_(IJ, M)**2
        DSIP(IJ, M) = TMP02*ZPIFR(M)*XLOGDFRTH / CGROUP(IJ, M)          !  coef*dtheta*dk = coef*dtheta*dsigma/cg
      END DO

!$acc loop seq
      DO M2=1,NFRE - NDIKCUMUL
        IF (BTH0(IJ, M2) > SDSBR) THEN
          TEMP1(IJ, M2) = 1.0_JWRB
        ELSE
          TEMP1(IJ, M2) = 0.0_JWRB
        END IF
      END DO
!$acc loop seq
      DO M2=1,NFRE - NDIKCUMUL
!$acc loop seq
        DO K2=1,NANG
          SCUMUL(IJ, K2, M2) = TEMP1(IJ, M2)*MAX(SQRT(BTH(IJ, K2, M2)) - EPSR, 0.0_JWRB)**2
        END DO
      END DO

!$acc loop seq
      DO M=NDIKCUMUL + 1,NFRE
!$acc loop seq
        DO K=1,NANG
          RENEWALFREQ(IJ, K, M) = 0.0_JWRB
        END DO
      END DO


!$acc loop seq
      DO M=NDIKCUMUL + 1,NFRE

        IF (M - NDIKCUMUL >= 3) THEN
          TRPZ_DSIP(IJ, 1) = 0.5_JWRB*DSIP(IJ, 1)
!$acc loop seq
          DO M2=2,M - NDIKCUMUL - 1
            TRPZ_DSIP(IJ, M2) = DSIP(IJ, M2)
          END DO
          TRPZ_DSIP(IJ, M - NDIKCUMUL) = 0.5_JWRB*DSIP(IJ, M - NDIKCUMUL)
        ELSE
!$acc loop seq
          DO M2=1,M - NDIKCUMUL
            TRPZ_DSIP(IJ, M2) = DSIP(IJ, M2)
          END DO
        END IF

!$acc loop seq
        DO M2=1,M - NDIKCUMUL
!$acc loop seq
          DO KK=0,NANGD
            WCUMUL(IJ, KK, M2) =  &
            & SQRT(ABS(C_C(IJ, M) + C_C(IJ, M2) - 2.0_JWRB*C_(IJ, M)*C_(IJ, M2)*COS(KK*DELTH)))*TRPZ_DSIP(IJ, M2)
          END DO
        END DO

!$acc loop seq
        DO K=1,NANG
          ! Correction of saturation level for shallow-water kinematics
          ! Cumulative effect based on lambda   (breaking probability is
          ! the expected rate of sweeping by larger breaking waves)

!$acc loop seq
          DO K2=1,NANG
          END DO

!$acc loop seq
          DO M2=1,M - NDIKCUMUL
!$acc loop seq
            DO K2=1,NANG
              KK = ABS(K2 - K)
              IF (KK > NANGD) KK = KK - NANGD
              ! Integrates over frequencies M2 and directions K2 to
              ! Integration is performed from M2=1 to a frequency lower than M: IK-NDIKCUMUL
              RENEWALFREQ(IJ, K, M) = RENEWALFREQ(IJ, K, M) + WCUMUL(IJ, KK, M2)*SCUMUL(IJ, K2, M2)
            END DO
          END DO
        END DO
      END DO

!$acc loop seq
      DO M=NDIKCUMUL + 1,NFRE
!$acc loop seq
        DO K=1,NANG
          D(IJ, K, M) = D(IJ, K, M) + RENEWALFREQ(IJ, K, M)
        END DO
      END DO
    END IF


    !     WAVE-TURBULENCE INTERACTION TERM
    IF (SSDSC5 /= 0.0_JWRB) THEN
      TMP01 = 2._JWRB*SSDSC5 / G
      FACTURB = TMP01*RAORW(IJ)*UFRIC(IJ)*UFRIC(IJ)
!$acc loop seq
      DO M=1,NFRE
        FACWTRB(IJ, M) = ZPIFR(M)*WAVNUM(IJ, M)*FACTURB
!$acc loop seq
        DO K=1,NANG
          D(IJ, K, M) = D(IJ, K, M) - FACWTRB(IJ, M)*COSWDIF(IJ, K)
        END DO
      END DO
    END IF


    ! ADD ALL CONTRIBUTIONS TO SOURCE TERM
!$acc loop seq
    DO M=1,NFRE
!$acc loop seq
      DO K=1,NANG
        SL(IJ, K, M) = SL(IJ, K, M) + D(IJ, K, M)*FL1(IJ, K, M)
        FLD(IJ, K, M) = FLD(IJ, K, M) + D(IJ, K, M)
      END DO
    END DO



!$acc end data
  END SUBROUTINE SDISSIP_ARD_LOKI
END MODULE SDISSIP_ARD_LOKI_MOD
