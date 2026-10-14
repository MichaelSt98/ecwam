! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDISSIP_LOKI_MOD
  CONTAINS
  SUBROUTINE SDISSIP_LOKI (KIJS, KIJL, FL1, FLD, SL, WAVNUM, CGROUP, XK2CG, UFRIC, COSWDIF, RAORW, IJ, FACSAT, FACWTRB,  &
  & SDISSIP_TEMP1, BTH0, C_, C_C, DSIP, TRPZ_DSIP, BTH, SDISSIP_TEMP2, D, SCUMUL, RENEWALFREQ, WCUMUL)
    ! ----------------------------------------------------------------------

    !**** *SDISSIP* - COMPUTATION OF DEEP WATER DISSIPATION SOURCE FUNCTION.


    !*    PURPOSE.
    !     --------
    !       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
    !       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
    !       OF DISSIPATION SOURCE FUNCTION.

    !**   INTERFACE.
    !     ----------

    !       *CALL* *SDISSIP (KIJS, KIJL, FL1, FLD, SL, *
    !                        WAVNUM, CGROUP, XK2CG,
    !                        EMEAN, F1MEAN, XKMEAN,*
    !                        UFRIC, COSWDIF, RAORW)*
    !         *KIJS* - INDEX OF FIRST GRIDPOINT
    !         *KIJL* - INDEX OF LAST GRIDPOINT
    !         *FL1*  - SPECTRUM.
    !         *FLD*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*  - TOTAL SOURCE FUNCTION ARRAY
    !       *WAVNUM* - WAVE NUMBER
    !       *CGROUP* - GROUP SPEED
    !       *XK2CG*  - (WAVNUM)**2 * GROUP SPEED
    !       *EMEAN* - MEAN ENERGY DENSITY
    !       *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
    !       *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.
    !       *UFRIC*  - FRICTION VELOCITY IN M/S.
    !       *RAORW*  - RATIO AIR DENSITY TO WATER DENSITY
    !       *COSWDIF*-  COS(TH(K)-WDWAVE(IJ))

    ! ----------------------------------------------------------------------
    USE YOWFRED, ONLY: FR, TH, ZPIFR, FRATIO, DELTH
    USE YOWPCONS, ONLY: G, ZPI
    USE YOWPHYS, ONLY: SDSBR, ISDSDTH, ISB, IPSAT, SSDSC2, SSDSC4, SSDSC6, MICHE, SSDSC3, SSDSBRF1, BRKPBCOEF, SSDSC5, NSDSNTH,  &
    & INDICESSAT, SATWEIGHTS
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWSTAT, ONLY: IPHYS

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE
    INTERFACE
      SUBROUTINE SDISSIP_JAN (KIJS, KIJL, FL1, FLD, SL, WAVNUM, EMEAN, F1MEAN, XKMEAN)
        USE parkind_wave, ONLY: jwim, jwrb
        USE YOWPARAM, ONLY: NANG, NFRE
        INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
        REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
        REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL
        REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM
        REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: EMEAN, F1MEAN, XKMEAN
      END SUBROUTINE SDISSIP_JAN
    END INTERFACE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM, CGROUP, XK2CG
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: UFRIC, RAORW
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF

    REAL(KIND=JPHOOK) :: SDISSIP_ZHOOK_HANDLE
    INTEGER(KIND=JWIM) :: SDISSIP_K, SDISSIP_M, SDISSIP_I, SDISSIP_J, M2, K2, KK, NANGD
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: NDIKCUMUL
    REAL(KIND=JWRB) :: TPIINV, TPIINVH, TMP01, TMP02, TMP03
    REAL(KIND=JWRB) :: EPSR, SSDSC6M1, ZCOEF, ZCOEFM1
    REAL(KIND=JWRB) :: XLOGDFRTH, BRLAMBDA
    REAL(KIND=JPHOOK) :: SDISSIP_ARD_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: SSDSC2_SIG
    REAL(KIND=JWRB) :: FACTURB
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: FACSAT
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: FACWTRB
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: SDISSIP_TEMP1
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: BTH0    !saturation spectrum
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: C_
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: C_C
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: DSIP
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: TRPZ_DSIP
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: BTH    !saturation spectrum
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: SDISSIP_TEMP2
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: D
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: SCUMUL
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: RENEWALFREQ
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, 0:NANG / 2, NFRE) :: WCUMUL
!$acc routine seq
!$acc data present( FL1, FLD, SL, WAVNUM, CGROUP, XK2CG, UFRIC, COSWDIF, RAORW, FACSAT, FACWTRB, SDISSIP_TEMP1, BTH0, C_, C_C,  &
!$acc & DSIP, TRPZ_DSIP, BTH, SDISSIP_TEMP2, D, SCUMUL, RENEWALFREQ, WCUMUL )

    ! ----------------------------------------------------------------------


    SELECT CASE (IPHYS)
    CASE (0)
      ! [Loki] Removed content of pragma-marked region!
    CASE (1)
      ! [Loki] inlined child subroutine: SDISSIP_ARD
      ! =========================================

      ! ----------------------------------------------------------------------


      ! INITIALISATION

      EPSR = SQRT(SDSBR)

      TPIINV = 1.0_JWRB / ZPI
      TPIINVH = 0.5_JWRB*TPIINV


      TMP03 = 1.0_JWRB / (SDSBR*MICHE)

      SSDSC6M1 = 1._JWRB - SSDSC6

!$acc loop seq
      DO SDISSIP_M=1,NFRE
        FACSAT(IJ, SDISSIP_M) = WAVNUM(IJ, SDISSIP_M)*TPIINV*XK2CG(IJ, SDISSIP_M)
      END DO

      ! COMPUTE SATURATION SPECTRUM
!$acc loop seq
      DO SDISSIP_M=1,NFRE
        BTH0(IJ, SDISSIP_M) = 0.0_JWRB
!$acc loop seq
        DO SDISSIP_K=1,NANG
          BTH(IJ, SDISSIP_K, SDISSIP_M) = 0.0_JWRB
        END DO
      END DO

!$acc loop seq
      DO SDISSIP_M=1,NFRE
!$acc loop seq
        DO SDISSIP_K=1,NANG
          ! integrates in directional sector
!$acc loop seq
          DO K2=1,NSDSNTH*2 + 1
            KK = INDICESSAT(SDISSIP_K, K2)
            BTH(IJ, SDISSIP_K, SDISSIP_M) = BTH(IJ, SDISSIP_K, SDISSIP_M) + SATWEIGHTS(SDISSIP_K, K2)*FL1(IJ, KK, SDISSIP_M)
          END DO
          BTH(IJ, SDISSIP_K, SDISSIP_M) = BTH(IJ, SDISSIP_K, SDISSIP_M)*FACSAT(IJ, SDISSIP_M)
          BTH0(IJ, SDISSIP_M) = MAX(BTH0(IJ, SDISSIP_M), BTH(IJ, SDISSIP_K, SDISSIP_M))
        END DO
      END DO


      ! SATURATION TERM

!$acc loop seq
      DO SDISSIP_M=1,NFRE
        SSDSC2_SIG = SSDSC2*ZPIFR(SDISSIP_M)
        ZCOEF = SSDSC2_SIG*SSDSC6
        ZCOEFM1 = SSDSC2_SIG*SSDSC6M1
!$acc loop seq
        DO SDISSIP_K=1,NANG
          D(IJ, SDISSIP_K, SDISSIP_M) = ZCOEF*MAX(0._JWRB, BTH0(IJ, SDISSIP_M)*TMP03 - SSDSC4)**IPSAT + ZCOEFM1*MAX(0._JWRB,  &
          & BTH(IJ, SDISSIP_K, SDISSIP_M)*TMP03 - SSDSC4)**IPSAT
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
        DO SDISSIP_M=1,NFRE
          C_(IJ, SDISSIP_M) = ZPIFR(SDISSIP_M) / WAVNUM(IJ, SDISSIP_M)
          C_C(IJ, SDISSIP_M) = C_(IJ, SDISSIP_M)**2
          DSIP(IJ, SDISSIP_M) = TMP02*ZPIFR(SDISSIP_M)*XLOGDFRTH / CGROUP(IJ, SDISSIP_M)            !  coef*dtheta*dk = coef*dtheta*dsigma/cg
        END DO

!$acc loop seq
        DO M2=1,NFRE - NDIKCUMUL
          IF (BTH0(IJ, M2) > SDSBR) THEN
            SDISSIP_TEMP1(IJ, M2) = 1.0_JWRB
          ELSE
            SDISSIP_TEMP1(IJ, M2) = 0.0_JWRB
          END IF
        END DO
!$acc loop seq
        DO M2=1,NFRE - NDIKCUMUL
!$acc loop seq
          DO K2=1,NANG
            SCUMUL(IJ, K2, M2) = SDISSIP_TEMP1(IJ, M2)*MAX(SQRT(BTH(IJ, K2, M2)) - EPSR, 0.0_JWRB)**2
          END DO
        END DO

!$acc loop seq
        DO SDISSIP_M=NDIKCUMUL + 1,NFRE
!$acc loop seq
          DO SDISSIP_K=1,NANG
            RENEWALFREQ(IJ, SDISSIP_K, SDISSIP_M) = 0.0_JWRB
          END DO
        END DO


!$acc loop seq
        DO SDISSIP_M=NDIKCUMUL + 1,NFRE

          IF (SDISSIP_M - NDIKCUMUL >= 3) THEN
            TRPZ_DSIP(IJ, 1) = 0.5_JWRB*DSIP(IJ, 1)
!$acc loop seq
            DO M2=2,SDISSIP_M - NDIKCUMUL - 1
              TRPZ_DSIP(IJ, M2) = DSIP(IJ, M2)
            END DO
            TRPZ_DSIP(IJ, SDISSIP_M - NDIKCUMUL) = 0.5_JWRB*DSIP(IJ, SDISSIP_M - NDIKCUMUL)
          ELSE
!$acc loop seq
            DO M2=1,SDISSIP_M - NDIKCUMUL
              TRPZ_DSIP(IJ, M2) = DSIP(IJ, M2)
            END DO
          END IF

!$acc loop seq
          DO M2=1,SDISSIP_M - NDIKCUMUL
!$acc loop seq
            DO KK=0,NANGD
              WCUMUL(IJ, KK, M2) = SQRT(ABS(C_C(IJ, SDISSIP_M) + C_C(IJ, M2) - 2.0_JWRB*C_(IJ, SDISSIP_M)*C_(IJ, M2)*COS(KK*DELTH &
              & )))*TRPZ_DSIP(IJ, M2)
            END DO
          END DO

!$acc loop seq
          DO SDISSIP_K=1,NANG
            ! Correction of saturation level for shallow-water kinematics
            ! Cumulative effect based on lambda   (breaking probability is
            ! the expected rate of sweeping by larger breaking waves)

!$acc loop seq
            DO K2=1,NANG
            END DO

!$acc loop seq
            DO M2=1,SDISSIP_M - NDIKCUMUL
!$acc loop seq
              DO K2=1,NANG
                KK = ABS(K2 - SDISSIP_K)
                IF (KK > NANGD) KK = KK - NANGD
                ! Integrates over frequencies M2 and directions K2 to
                ! Integration is performed from M2=1 to a frequency lower than M: IK-NDIKCUMUL
                RENEWALFREQ(IJ, SDISSIP_K, SDISSIP_M) =  &
                & RENEWALFREQ(IJ, SDISSIP_K, SDISSIP_M) + WCUMUL(IJ, KK, M2)*SCUMUL(IJ, K2, M2)
              END DO
            END DO
          END DO
        END DO

!$acc loop seq
        DO SDISSIP_M=NDIKCUMUL + 1,NFRE
!$acc loop seq
          DO SDISSIP_K=1,NANG
            D(IJ, SDISSIP_K, SDISSIP_M) = D(IJ, SDISSIP_K, SDISSIP_M) + RENEWALFREQ(IJ, SDISSIP_K, SDISSIP_M)
          END DO
        END DO
      END IF


      !     WAVE-TURBULENCE INTERACTION TERM
      IF (SSDSC5 /= 0.0_JWRB) THEN
        TMP01 = 2._JWRB*SSDSC5 / G
        FACTURB = TMP01*RAORW(IJ)*UFRIC(IJ)*UFRIC(IJ)
!$acc loop seq
        DO SDISSIP_M=1,NFRE
          FACWTRB(IJ, SDISSIP_M) = ZPIFR(SDISSIP_M)*WAVNUM(IJ, SDISSIP_M)*FACTURB
!$acc loop seq
          DO SDISSIP_K=1,NANG
            D(IJ, SDISSIP_K, SDISSIP_M) = D(IJ, SDISSIP_K, SDISSIP_M) - FACWTRB(IJ, SDISSIP_M)*COSWDIF(IJ, SDISSIP_K)
          END DO
        END DO
      END IF


      ! ADD ALL CONTRIBUTIONS TO SOURCE TERM
!$acc loop seq
      DO SDISSIP_M=1,NFRE
!$acc loop seq
        DO SDISSIP_K=1,NANG
          SL(IJ, SDISSIP_K, SDISSIP_M) = SL(IJ, SDISSIP_K, SDISSIP_M) + D(IJ, SDISSIP_K, SDISSIP_M)*FL1(IJ, SDISSIP_K, SDISSIP_M)
          FLD(IJ, SDISSIP_K, SDISSIP_M) = FLD(IJ, SDISSIP_K, SDISSIP_M) + D(IJ, SDISSIP_K, SDISSIP_M)
        END DO
      END DO


      ! =========================================
    END SELECT



!$acc end data
  END SUBROUTINE SDISSIP_LOKI
END MODULE SDISSIP_LOKI_MOD
