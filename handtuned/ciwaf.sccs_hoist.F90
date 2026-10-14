! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE CIWAF_LOKI_MOD
  CONTAINS
  SUBROUTINE CIWAF_LOKI (KIJS, KIJL, CGROUP, CICOVER, CITHICK, CIWA, IJ)

    ! ----------------------------------------------------------------------

    !**** *CIWAF* - COMPUTE SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
    !               SCATTERING.

    !*    PURPOSE.
    !     --------

    !       CIWAF COMPUTES SEA ICE WAVE ATTENUATION FACTORS

    !**   INTERFACE.
    !     ----------

    !       *CALL* *CIWAF (KIJS, KIJL, CGROUP, CICOVER, CITHICK, CIWA)

    !          *KIJS*     - INDEX OF FIRST POINT.
    !          *KIJL*     - INDEX OF LAST POINT.
    !          *CGROUP*   - GROUP SPEED.
    !          *CICOVER*  - SEA ICE COVER.
    !          *CITHICK*  - SEA ICE THICKNESS.
    !          *CIWA*     - SEA ICE WAVE ATTENUATION FACTOR.

    !     METHOD.
    !     -------

    !     EXTERNALS.
    !     ----------

    !     REFERENCES.
    !     -----------

    !     KOHOUT AND MEYLAN, 2008: JGR, 113, doi:10.1029/2007JC004424


    ! ----------------------------------------------------------------------

    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWFRED, ONLY: FR
    USE YOWICE, ONLY: NICT, NICH, TICMIN, HICMIN, DTIC, DHIC, CIDEAC, CIBLOCK
    USE YOWPARAM, ONLY: NFRE
    USE YOWPCONS, ONLY: EPSMIN
    USE YOWSTAT, ONLY: IDELT

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: CGROUP
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CICOVER
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CITHICK
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NFRE) :: CIWA


    INTEGER(KIND=JWIM) :: ICM, I, MAXICM
    INTEGER(KIND=JWIM) :: M
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: IT, IT1, IH, IH1

    REAL(KIND=JWRB) :: CIDMIN, CIDMAX, CIDMEAN, CIFRGL, CIFRGMT
    REAL(KIND=JWRB) :: SD, SN
    REAL(KIND=JWRB) :: TW, X
    REAL(KIND=JWRB) :: A, B, C
    REAL(KIND=JWRB) :: CIDEAC_INT, WT, WT1, WH, WH1
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    REAL(KIND=JWRB) :: DINV
    REAL(KIND=JWRB) :: ALP
!$acc routine seq
!$acc data present( CGROUP, CICOVER, CITHICK, CIWA )

    ! ----------------------------------------------------------------------


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
    IF (CITHICK(IJ) > 0.0_JWRB) THEN
      ! sea ice foes maxmimum size (m)
      !!!       testing making it function of sea ice cover
      CIDMAX = A + C*CICOVER(IJ)
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
      IF (CITHICK(IJ) > 0.0_JWRB) THEN
        IH = FLOOR((CITHICK(IJ) - HICMIN) / DHIC + 1)
        IH = MAX(1, MIN(IH, NICH))
        IH1 = IH + 1
        IH1 = MAX(1, MIN(IH1, NICH))
        WH1 = MAX(MIN(1._JWRB, (CITHICK(IJ) - (HICMIN + (IH - 1)*DHIC)) / DHIC), 0.0_JWRB)
        WH = 1.0_JWRB - WH1
        CIDEAC_INT = WT*(WH*CIDEAC(IT, IH) + WH1*CIDEAC(IT, IH1)) + WT1*(WH*CIDEAC(IT1, IH) + WH1*CIDEAC(IT1, IH1))
        !!!            ALP(IJ,M)=CICOVER(IJ)*CIDEAC_INT*DINV(IJ)
        ALP = CICOVER(IJ)*EXP(CIDEAC_INT)*DINV
      ELSE
        ALP = 0.0_JWRB
      END IF

      X = ALP*CGROUP(IJ, M)*IDELT
      IF (X < EPSMIN) THEN
        CIWA(IJ, M) = 1.0_JWRB
      ELSE IF (CICOVER(IJ) > CIBLOCK) THEN
        CIWA(IJ, M) = 0.0_JWRB
      ELSE
        CIWA(IJ, M) = EXP(-MIN(X, 50.0_JWRB))
      END IF
    END DO



!$acc end data
  END SUBROUTINE CIWAF_LOKI
END MODULE CIWAF_LOKI_MOD
