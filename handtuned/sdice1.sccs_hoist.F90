! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDICE1_LOKI_MOD
  CONTAINS
  SUBROUTINE SDICE1_LOKI (KIJS, KIJL, FL1, FLD, SL, SLICE, CGROUP, CICV, CITH, IJ, ALP)
    ! ----------------------------------------------------------------------

    !**** *SDICE1* - COMPUTATION OF SEA ICE ATTENUATION DUE TO SCATTERING (CAME FROM CIWAF)


    !     JEAN BIDLOT       ECMWF ~2012
    !     JOSH KOUSAL       ECMWF 2023

    !*    PURPOSE.
    !     --------

    !**   INTERFACE.
    !     ----------

    !       *CALL* *SDICE1 (KIJS, KIJL, FL1, FLD,SL,
    !                       WAVNUM, CGROUP,
    !                       CICV,CITH)*
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *FL1*    - SPECTRUM.
    !          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*     - TOTAL SOURCE FUNCTION ARRAY
    !          *SLICE*  - TOTAL SOURCE FUNCTION ARRAY, ICE
    !          *WAVNUM* - WAVE NUMBER
    !          *CGROUP* - GROUP SPEED
    !          *CICV*   - SEA ICE COVER
    !          *CITH*   - SEA ICE THICKNESS

    !     METHOD.
    !     -------

    !       SEE REFERENCES.

    !     EXTERNALS.
    !     ----------

    !       NONE.

    !     REFERENCE.
    !     ----------

    !     KOHOUT AND MEYLAN, 2008: JGR, 113, doi:10.1029/2007JC004424
    !     M.J. Doble, J.-R. Bidlot / Ocean Modelling 70 (2013), 166-173

    ! ----------------------------------------------------------------------

    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWFRED, ONLY: FR, DFIM
    USE YOWICE, ONLY: NICT, NICH, TICMIN, HICMIN, DTIC, DHIC, CIDEAC, ZALPFACB
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: G, ZPI, EPSMIN
    USE YOWSTAT, ONLY: IDELT, XIMP

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: SLICE
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: CGROUP
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CICV
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CITH

    INTEGER(KIND=JWIM) :: K, M
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ

    REAL(KIND=JWRB) :: FLDICE
    REAL(KIND=JWRB) :: DELTM, DELT5, DELT, GTEMP1

    INTEGER(KIND=JWIM) :: ICM, I, MAXICM
    INTEGER(KIND=JWIM) :: IT, IT1, IH, IH1

    REAL(KIND=JWRB) :: CIDMIN, CIDMAX, CIDMEAN, CIFRGL, CIFRGMT
    REAL(KIND=JWRB) :: SD, SN
    REAL(KIND=JWRB) :: TW, X
    REAL(KIND=JWRB) :: A, B, C
    REAL(KIND=JWRB) :: CIDEAC_INT, WT, WT1, WH, WH1
    REAL(KIND=JWRB) :: DINV
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: ALP

    REAL(KIND=JPHOOK) :: SDICE1_ZHOOK_HANDLE
!$acc routine seq
!$acc data present( FL1, FLD, SL, SLICE, CGROUP, CICV, CITH, ALP )


    ! ----------------------------------------------------------------------


    DELT = IDELT
    DELTM = 1.0_JWRB / DELT
    DELT5 = XIMP*DELT

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
      DO I=0,ICM
        X = (CIFRGMT**2*CIFRGL)**I
        SN = SN + X*CIDMAX / CIFRGMT**I
        SD = SD + X
      END DO
      CIDMEAN = SN / SD
      DINV = 1.0_JWRB / CIDMEAN
    ELSE
      DINV = CIDMIN
    END IF


!$acc loop seq
    DO M=1,NFRE
      TW = 1.0_JWRB / FR(M)
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
        ALP(IJ, M) = EXP(CIDEAC_INT)*DINV*ZALPFACB          ! CICV accounted for below

      ELSE
        ALP(IJ, M) = 0.0_JWRB
      END IF
    END DO

!$acc loop seq
    DO M=1,NFRE
!$acc loop seq
      DO K=1,NANG

        !              apply the source term
        FLDICE = -ALP(IJ, M)*CGROUP(IJ, M)
        SLICE(IJ, K, M) = FL1(IJ, K, M)*FLDICE
        SL(IJ, K, M) = SL(IJ, K, M) + CICV(IJ)*SLICE(IJ, K, M)
        FLD(IJ, K, M) = FLD(IJ, K, M) + CICV(IJ)*FLDICE

        !              to be used for wave radiative stress calculation
        GTEMP1 = MAX((1.0_JWRB - DELT5*FLDICE), 1.0_JWRB)
        SLICE(IJ, K, M) = SLICE(IJ, K, M) / GTEMP1

      END DO
    END DO



!$acc end data
  END SUBROUTINE SDICE1_LOKI
END MODULE SDICE1_LOKI_MOD
