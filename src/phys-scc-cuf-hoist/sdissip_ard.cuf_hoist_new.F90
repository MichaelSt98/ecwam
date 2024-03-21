! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDISSIP_ARD_CUF_HOIST_NEW_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE SDISSIP_ARD_CUF_HOIST_NEW (KIJS, KIJL, FL1, FLD, SL, INDEP, WAVNUM, XK2CG, UFRIC, COSWDIF,  &
  & RAORW, CUMULW, G, INDICESSAT, IPSAT, MICHE, NANG, NDEPTH, NDIKCUMUL, NFRE, NSDSNTH, SATWEIGHTS, SDSBR, SSDSC2, SSDSC3,  &
  & SSDSC4, SSDSC5, SSDSC6, ZPI, ZPIFR, ICHNK, NCHNK, IJ)
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
    !                            INDEP, WAVNUM, XK2CG,
    !                            UFRIC, COSWDIF, RAORW)*
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *FL1*    - SPECTRUM.
    !          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*     - TOTAL SOURCE FUNCTION ARRAY
    !          *INDEP*  - DEPTH INDEX
    !          *WAVNUM* - WAVE NUMBER
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
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    USE YOWFRED, ONLY: FR, TH
    USE YOWPARAM, ONLY: NANG_PARAM
    USE YOWPHYS, ONLY: SSDSBRF1, ISDSDTH, ISB, BRKPBCOEF
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
    INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    
    REAL(KIND=JWRB), INTENT(IN) :: FL1(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(INOUT) :: FLD(KIJL, NANG_loki_param, NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(INOUT) :: SL(KIJL, NANG_loki_param, NFRE_loki_param)
    INTEGER(KIND=JWIM), INTENT(IN) :: INDEP(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: WAVNUM(KIJL, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: XK2CG(KIJL, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: UFRIC(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: RAORW(KIJL)
    REAL(KIND=JWRB), INTENT(IN) :: COSWDIF(KIJL, NANG_loki_param)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: K
    INTEGER(KIND=JWIM) :: M
    INTEGER(KIND=JWIM) :: I
    INTEGER(KIND=JWIM) :: J
    INTEGER(KIND=JWIM) :: M2
    INTEGER(KIND=JWIM) :: K2
    INTEGER(KIND=JWIM) :: KK
    
    REAL(KIND=JWRB) :: TPIINV
    REAL(KIND=JWRB) :: TPIINVH
    REAL(KIND=JWRB) :: TMP01
    REAL(KIND=JWRB) :: TMP03
    REAL(KIND=JWRB) :: EPSR
    REAL(KIND=JWRB) :: SSDSC6M1
    REAL(KIND=JWRB) :: ZCOEF
    REAL(KIND=JWRB) :: ZCOEFM1
    
    
    REAL(KIND=JWRB) :: SSDSC2_SIG
    REAL(KIND=JWRB) :: FACTURB
    REAL(KIND=JWRB) :: BTH
    REAL(KIND=JWRB) :: BTH0
    REAL(KIND=JWRB) :: SCUMUL(NANG_PARAM)
    REAL(KIND=JWRB) :: D(NANG_PARAM)
    
    REAL(KIND=JWRB) :: RENEWALFREQ
    INTEGER :: FOO
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: CUMULW(NDEPTH, 0:NANG/2, NFRE_loki_param, NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
    INTEGER(KIND=JWIM), INTENT(IN), DEVICE :: INDICESSAT(NANG_loki_param, NSDSNTH*2+1)
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IPSAT
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: MICHE
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NDEPTH
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NDIKCUMUL
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NSDSNTH
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: SATWEIGHTS(NANG_loki_param, NSDSNTH*2+1)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SDSBR
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SSDSC2
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SSDSC3
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SSDSC4
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SSDSC5
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SSDSC6
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZPI
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: ZPIFR(NFRE_loki_param)
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    ! ----------------------------------------------------------------------
    
    
    ! INITIALISATION
    
    FOO = NDEPTH      ! necessary for Loki ...
    EPSR = SQRT(SDSBR)
    
    TPIINV = 1.0_JWRB / ZPI
    TPIINVH = 0.5_JWRB*TPIINV
    TMP03 = 1.0_JWRB / (SDSBR*MICHE)
    SSDSC6M1 = 1._JWRB - SSDSC6
    
    
    DO M=1,NFRE
      
      ! SATURATION TERM
      SSDSC2_SIG = SSDSC2*ZPIFR(M)
      ZCOEF = SSDSC2_SIG*SSDSC6
      ZCOEFM1 = SSDSC2_SIG*SSDSC6M1
      
      ! COMPUTE SATURATION SPECTRUM
      BTH0 = 0.0_JWRB
      
      DO K=1,NANG
        BTH = 0.0_JWRB
        ! integrates in directional sector
        DO K2=1,NSDSNTH*2 + 1
          KK = INDICESSAT(K, K2)
          BTH = BTH + SATWEIGHTS(K, K2)*FL1(IJ, KK, M, ICHNK)
        END DO
        BTH = BTH*WAVNUM(IJ, M, ICHNK)*TPIINV*XK2CG(IJ, M, ICHNK)
        BTH0 = MAX(BTH0, BTH)
        
        D(K) = ZCOEFM1*MAX(0._JWRB, BTH*TMP03 - SSDSC4)**IPSAT
        
        SCUMUL(K) = MAX(SQRT(ABS(BTH)) - EPSR, 0._JWRB)**2
      END DO
      
      DO K=1,NANG
        ! cumulative term
        D(K) = D(K) + ZCOEF*MAX(0._JWRB, BTH0*TMP03 - SSDSC4)**IPSAT
        IF (BTH0 <= SDSBR) THEN
          SCUMUL(K) = 0._JWRB
        END IF
        
      END DO
      
      IF (M > NDIKCUMUL) THEN
        ! CUMULATIVE TERM
        IF (SSDSC3 /= 0.0_JWRB) THEN
          
          DO K=1,NANG
            ! Correction of saturation level for shallow-water kinematics
            ! Cumulative effect based on lambda   (breaking probability is
            ! the expected rate of sweeping by larger breaking waves)
            
            RENEWALFREQ = 0.0_JWRB
            
            DO M2=1,M - NDIKCUMUL
              DO K2=1,NANG
                KK = ABS(K2 - K)
                IF (KK > NANG / 2) KK = KK - NANG / 2
                ! Integrates over frequencies M2 and directions K2 to
                ! Integration is performed from M2=1 to a frequency lower than M: IK-NDIKCUMUL
                RENEWALFREQ = RENEWALFREQ + CUMULW(INDEP(IJ, ICHNK), KK, M2, M)*SCUMUL(K2)
              END DO
            END DO
            
            D(K) = D(K) + RENEWALFREQ
          END DO
        END IF
      END IF
      
      !       WAVE-TURBULENCE INTERACTION TERM
      IF (SSDSC5 /= 0.0_JWRB) THEN
        TMP01 = 2._JWRB*SSDSC5 / G
        FACTURB = TMP01*RAORW(IJ)*UFRIC(IJ, ICHNK)*UFRIC(IJ, ICHNK)
        DO K=1,NANG
          D(K) = D(K) - ZPIFR(M)*WAVNUM(IJ, M, ICHNK)*FACTURB*COSWDIF(IJ, K)
        END DO
      END IF
      
      
      ! ADD ALL CONTRIBUTIONS TO SOURCE TERM
      DO K=1,NANG
        SL(IJ, K, M) = SL(IJ, K, M) + D(K)*FL1(IJ, K, M, ICHNK)
        FLD(IJ, K, M) = FLD(IJ, K, M) + D(K)
      END DO
    END DO
    
    
    
  END SUBROUTINE SDISSIP_ARD_CUF_HOIST_NEW
END MODULE SDISSIP_ARD_CUF_HOIST_NEW_MOD
