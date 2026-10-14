! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDICE_LOKI_MOD
  CONTAINS
  SUBROUTINE SDICE_LOKI (KIJS, KIJL, FL1, FLD, SL, SLICE, WAVNUM, CGROUP, CICV, CITH, ALPFAC, IJ, ALP, XK2, SDICE3_ALP)
    ! ----------------------------------------------------------------------

    !**** *SDICE* - CALLING OF DIFFERENT SEA ICE ATTENUATION TERMS


    !     JOSH KOUSAL       ECMWF 2023


    !*    PURPOSE.
    !     --------

    !**   INTERFACE.
    !     ----------

    !       *CALL* *SDICE (KIJS, KIJL, FL1, FLD,SL,*
    !                       WAVNUM, CGROUP,
    !                       CICV,CITH, ALPFAC)*
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *FL1*    - SPECTRUM.
    !          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*     - TOTAL SOURCE FUNCTION ARRAY
    !          *SLICE*  - TOTAL SOURCE FUNCTION ARRAY, ICE
    !          *INDEP*  - DEPTH INDEX
    !          *WAVNUM* - WAVE NUMBER
    !          *CGROUP* - GROUP SPEED
    !          *CICV*   - SEA ICE COVER
    !          *CITH*   - SEA ICE THICKNESS
    !          *ALPFAC* - FACTOR TO REDUCE ATTENUATION IN ICE. EQUAL TO 1 IN SOLID ICE

    !     METHOD.
    !     -------

    !       SEE REFERENCES.

    !     EXTERNALS.
    !     ----------

    !       NONE.

    !     REFERENCE.
    !     ----------

    !     Jie Yu * , W. Erick Rogers, David W. Wang, 2022
    !     DOI:10.1016/j.coldregions.2022.103582

    ! ----------------------------------------------------------------------

    USE YOWFRED, ONLY: FR, DFIM
    USE YOWPCONS, ONLY: G, ZPI, EPSMIN
    USE YOWSTAT, ONLY: IDELT, XIMP
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWPARAM, ONLY: NANG, NFRE

    USE YOWICE, ONLY: LCIWA1, LCIWA2, LCIWA3, NICT, NICH, TICMIN, HICMIN, DTIC, DHIC, CIDEAC, ZALPFACB, CDICWA

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: SLICE
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM, CGROUP
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CICV
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CITH
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: ALPFAC

    REAL(KIND=JPHOOK) :: SDICE_ZHOOK_HANDLE
    INTEGER(KIND=JWIM) :: SDICE_K, SDICE_M
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    REAL(KIND=JWRB) :: FLDICE
    REAL(KIND=JWRB) :: SDICE_DELTM, SDICE_DELT5, SDICE_DELT, SDICE_GTEMP1
    INTEGER(KIND=JWIM) :: ICM, SDICE_I, MAXICM
    INTEGER(KIND=JWIM) :: IT, IT1, IH, IH1
    REAL(KIND=JWRB) :: CIDMIN, CIDMAX, CIDMEAN, CIFRGL, CIFRGMT
    REAL(KIND=JWRB) :: SD, SN
    REAL(KIND=JWRB) :: TW, SDICE_X
    REAL(KIND=JWRB) :: A, B, C
    REAL(KIND=JWRB) :: CIDEAC_INT, WT, WT1, WH, WH1
    REAL(KIND=JWRB) :: DINV
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: ALP
    REAL(KIND=JPHOOK) :: SDICE1_ZHOOK_HANDLE
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NFRE) :: XK2
    INTEGER(KIND=JWIM) :: SDICE2_K, SDICE2_M
    REAL(KIND=JWRB) :: EWH
    REAL(KIND=JWRB) :: SDICE2_ALP    !! ALP=SPATIAL ATTENUATION RATE OF ENERGY
    REAL(KIND=JWRB) :: SDICE2_FLDICE
    REAL(KIND=JWRB) :: SDICE2_DELTM, SDICE2_DELT5, SDICE2_DELT, SDICE2_GTEMP1
    REAL(KIND=JPHOOK) :: SDICE2_ZHOOK_HANDLE
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: SDICE3_ALP    !! ALP=SPATIAL ATTENUATION RATE OF ENERGY
    INTEGER(KIND=JWIM) :: IMODEL    !! DAMPING MODEL: 1=FIT TO TEMPELFJORD DATA, 2=Jie Yu 2022
    INTEGER(KIND=JWIM) :: SDICE3_K, SDICE3_M
    REAL(KIND=JWRB) :: SDICE3_FLDICE
    REAL(KIND=JWRB) :: SDICE3_DELTM, SDICE3_DELT5, SDICE3_DELT, SDICE3_GTEMP1
    REAL(KIND=JWRB) :: CDICE
    REAL(KIND=JWRB) :: HICEMAX, HICEMIN
    REAL(KIND=JWRB) :: SDICE_TEMP
    REAL(KIND=JPHOOK) :: SDICE3_ZHOOK_HANDLE
!$acc routine seq
!$acc data present( FL1, FLD, SL, SLICE, WAVNUM, CGROUP, CICV, CITH, ALPFAC, ALP, XK2, SDICE3_ALP )


    ! ----------------------------------------------------------------------


    !     Attenuation of waves in ice (type 1): scattering
    IF (LCIWA1) THEN
      ! [Loki] inlined child subroutine: SDICE1
      ! =========================================


      ! ----------------------------------------------------------------------


      SDICE_DELT = IDELT
      SDICE_DELTM = 1.0_JWRB / SDICE_DELT
      SDICE_DELT5 = XIMP*SDICE_DELT

      !     following  Dumont et al. (2011), eqn (13):
      ! sea ice fragility
      CIFRGL = 0.955_JWRB
      ! sea ice foes minimum size (m)
      CIDMIN = 20.0_JWRB
      ! sea ice fragmentation number (>2)
      CIFRGMT = 2.0_JWRB

      A = 200.0_JWRB
      C = 300.0_JWRB
      !     the maximum number of fragmentation step is set
      !     for the limit value for sea ice = 0 (CIDMAX=A)
      !     to insure that CIDMEAN will be increasing with increasing sea ice cover
      MAXICM = INT(LOG(A / CIDMIN) / LOG(CIFRGMT))
      IF (CITH(IJ) > 0.0_JWRB) THEN
        ! sea ice foes maxmimum size (m)
        !!!       testing making it function of sea ice cover
        CIDMAX = A + C*CICV(IJ)
        ICM = MIN(INT(LOG(CIDMAX / CIDMIN) / LOG(CIFRGMT)), MAXICM)
        SN = 0.0_JWRB
        SD = 0.0_JWRB
!$acc loop seq
        DO SDICE_I=0,ICM
          SDICE_X = (CIFRGMT**2*CIFRGL)**SDICE_I
          SN = SN + SDICE_X*CIDMAX / CIFRGMT**SDICE_I
          SD = SD + SDICE_X
        END DO
        CIDMEAN = SN / SD
        DINV = 1.0_JWRB / CIDMEAN
      ELSE
        DINV = CIDMIN
      END IF


!$acc loop seq
      DO SDICE_M=1,NFRE
        TW = 1.0_JWRB / FR(SDICE_M)
        IT = FLOOR((TW - TICMIN) / DTIC + 1)
        IT = MAX(1, MIN(IT, NICT))
        IT1 = IT + 1
        IT1 = MAX(1, MIN(IT1, NICT))
        WT1 = MAX(MIN(1.0_JWRB, (TW - (TICMIN + (IT - 1)*DTIC)) / DTIC), 0.0_JWRB)
        WT = 1.0_JWRB - WT1
        IF (CITH(IJ) > 0.0_JWRB) THEN
          IH = FLOOR((CITH(IJ) - HICMIN) / DHIC + 1)
          IH = MAX(1, MIN(IH, NICH))
          IH1 = IH + 1
          IH1 = MAX(1, MIN(IH1, NICH))
          WH1 = MAX(MIN(1._JWRB, (CITH(IJ) - (HICMIN + (IH - 1)*DHIC)) / DHIC), 0.0_JWRB)
          WH = 1.0_JWRB - WH1
          CIDEAC_INT = WT*(WH*CIDEAC(IT, IH) + WH1*CIDEAC(IT, IH1)) + WT1*(WH*CIDEAC(IT1, IH) + WH1*CIDEAC(IT1, IH1))
          ALP(IJ, SDICE_M) = EXP(CIDEAC_INT)*DINV*ZALPFACB            ! CICV accounted for below

        ELSE
          ALP(IJ, SDICE_M) = 0.0_JWRB
        END IF
      END DO

!$acc loop seq
      DO SDICE_M=1,NFRE
!$acc loop seq
        DO SDICE_K=1,NANG

          !              apply the source term
          FLDICE = -ALP(IJ, SDICE_M)*CGROUP(IJ, SDICE_M)
          SLICE(IJ, SDICE_K, SDICE_M) = FL1(IJ, SDICE_K, SDICE_M)*FLDICE
          SL(IJ, SDICE_K, SDICE_M) = SL(IJ, SDICE_K, SDICE_M) + CICV(IJ)*SLICE(IJ, SDICE_K, SDICE_M)
          FLD(IJ, SDICE_K, SDICE_M) = FLD(IJ, SDICE_K, SDICE_M) + CICV(IJ)*FLDICE

          !              to be used for wave radiative stress calculation
          SDICE_GTEMP1 = MAX((1.0_JWRB - SDICE_DELT5*FLDICE), 1.0_JWRB)
          SLICE(IJ, SDICE_K, SDICE_M) = SLICE(IJ, SDICE_K, SDICE_M) / SDICE_GTEMP1

        END DO
      END DO


      ! =========================================
    END IF

    !     Attenuation of waves in ice (type 2): bottom friction
    IF (LCIWA2) THEN
      ! [Loki] inlined child subroutine: SDICE2
      ! =========================================


      ! ----------------------------------------------------------------------


      SDICE2_DELT = IDELT
      SDICE2_DELTM = 1.0_JWRB / SDICE2_DELT
      SDICE2_DELT5 = XIMP*SDICE2_DELT

!$acc loop seq
      DO SDICE2_M=1,NFRE
!$acc loop seq
        DO SDICE2_K=1,NANG

          EWH = 4.0_JWRB*SQRT(MAX(EPSMIN, FL1(IJ, SDICE2_K, SDICE2_M)*DFIM(SDICE2_M)))
          XK2(SDICE2_M) = WAVNUM(IJ, SDICE2_M)**2
          SDICE2_ALP = CDICWA*XK2(SDICE2_M)*EWH*ZALPFACB

          !              apply the source term
          SDICE2_FLDICE = -SDICE2_ALP*CGROUP(IJ, SDICE2_M)
          SLICE(IJ, SDICE2_K, SDICE2_M) = FL1(IJ, SDICE2_K, SDICE2_M)*SDICE2_FLDICE
          SL(IJ, SDICE2_K, SDICE2_M) = SL(IJ, SDICE2_K, SDICE2_M) + CICV(IJ)*SLICE(IJ, SDICE2_K, SDICE2_M)
          FLD(IJ, SDICE2_K, SDICE2_M) = FLD(IJ, SDICE2_K, SDICE2_M) + CICV(IJ)*SDICE2_FLDICE

          !              to be used for wave radiative stress calculation
          SDICE2_GTEMP1 = MAX((1.0_JWRB - SDICE2_DELT5*SDICE2_FLDICE), 1.0_JWRB)
          SLICE(IJ, SDICE2_K, SDICE2_M) = SLICE(IJ, SDICE2_K, SDICE2_M) / SDICE2_GTEMP1

        END DO
      END DO


      ! =========================================
    END IF

    !     Attenuation of waves in ice (type 3): viscous friction
    IF (LCIWA3) THEN
      ! [Loki] inlined child subroutine: SDICE3
      ! =========================================


      ! ----------------------------------------------------------------------


      SDICE3_DELT = IDELT
      SDICE3_DELTM = 1.0_JWRB / SDICE3_DELT
      SDICE3_DELT5 = XIMP*SDICE3_DELT

      IMODEL = 2
      HICEMAX = 4.0_JWRB
      HICEMIN = 0.1_JWRB

      SELECT CASE (IMODEL)
      CASE (1)
        !          Best fit w Tempelfjorde obs from Lotfi Aouf
        CDICE = 0.0656_JWRB

!$acc loop seq
        DO SDICE3_M=1,NFRE
          SDICE3_ALP(IJ, SDICE3_M) = (CDICE*CITH(IJ)*WAVNUM(IJ, SDICE3_M)**2)*ALPFAC(IJ)
        END DO

      CASE (2)
        !          Jie Yu, W. Erik Rogers, David W. Wang 2022
        CDICE = 0.1274_JWRB*(ZPI / SQRT(G))**4.5_JWRB

!$acc loop seq
        DO SDICE3_M=1,NFRE
          SDICE3_ALP(IJ, SDICE3_M) = (2._JWRB*CDICE*(CITH(IJ)**1.25_JWRB)*(FR(SDICE3_M)**4.5_JWRB))*ALPFAC(IJ)*ZALPFACB
        END DO

      END SELECT

!$acc loop seq
      DO SDICE3_M=1,NFRE
!$acc loop seq
        DO SDICE3_K=1,NANG

          !              calculate source term increments
          SDICE3_FLDICE = -SDICE3_ALP(IJ, SDICE3_M)*CGROUP(IJ, SDICE3_M)
          SLICE(IJ, SDICE3_K, SDICE3_M) = FL1(IJ, SDICE3_K, SDICE3_M)*SDICE3_FLDICE

          !              apply the source term (old way to ensure bit-identicality)
          SDICE_TEMP = -CICV(IJ)*SDICE3_ALP(IJ, SDICE3_M)*CGROUP(IJ, SDICE3_M)
          SL(IJ, SDICE3_K, SDICE3_M) = SL(IJ, SDICE3_K, SDICE3_M) + FL1(IJ, SDICE3_K, SDICE3_M)*SDICE_TEMP
          FLD(IJ, SDICE3_K, SDICE3_M) = FLD(IJ, SDICE3_K, SDICE3_M) + SDICE_TEMP

          !              apply the source term (new way)
          !                SL(IJ,K,M)     =  SL(IJ,K,M)  + CICV(IJ)*SLICE(IJ,K,M)
          !                FLD(IJ,K,M)    =  FLD(IJ,K,M) + CICV(IJ)*FLDICE
          !
          !              to be used for wave radiative stress calculation
          SDICE3_GTEMP1 = MAX((1.0_JWRB - SDICE3_DELT5*SDICE3_FLDICE), 1.0_JWRB)
          SLICE(IJ, SDICE3_K, SDICE3_M) = SLICE(IJ, SDICE3_K, SDICE3_M) / SDICE3_GTEMP1

        END DO
      END DO


      ! =========================================
    END IF



!$acc end data
  END SUBROUTINE SDICE_LOKI
END MODULE SDICE_LOKI_MOD
