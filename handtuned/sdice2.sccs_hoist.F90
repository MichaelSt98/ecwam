! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDICE2_LOKI_MOD
  CONTAINS
  SUBROUTINE SDICE2_LOKI (KIJS, KIJL, FL1, FLD, SL, SLICE, WAVNUM, CGROUP, CICV, IJ, XK2)
    ! ----------------------------------------------------------------------

    !**** *SDICE2* - COMPUTE SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
    !                BOTTOM FRICTION (CAME FROM CIWABR)


    !     JEAN BIDLOT       ECMWF ~2012
    !     JOSH KOUSAL       ECMWF 2023


    !*    PURPOSE.
    !     --------

    !**   INTERFACE.
    !     ----------

    !       *CALL* *SDICE2 (KIJS, KIJL, FL1, FLD,SL,*
    !                       WAVNUM, CGROUP,
    !                       CICV)*
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *FL1*    - SPECTRUM.
    !          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*     - TOTAL SOURCE FUNCTION ARRAY
    !          *SLICE*  - TOTAL SOURCE FUNCTION ARRAY, ICE
    !          *WAVNUM* - WAVE NUMBER
    !          *CGROUP* - GROUP SPEED
    !          *CICV*   - SEA ICE COVER

    !     METHOD.
    !     -------

    !       SEE REFERENCES.

    !     EXTERNALS.
    !     ----------

    !       NONE.

    !     REFERENCE.
    !     ----------

    !     KOHOUT A., M. MEYLAN, D PLEW, 2011: ANNALS OF GLACIOLOGY, 2011.
    !     M.J. Doble, J.-R. Bidlot / Ocean Modelling 70 (2013), 166-173

    ! ----------------------------------------------------------------------

    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWFRED, ONLY: DFIM
    USE YOWICE, ONLY: CDICWA, ZALPFACB
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: EPSMIN
    USE YOWSTAT, ONLY: IDELT, XIMP

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: SLICE
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM, CGROUP
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CICV

    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NFRE) :: XK2

    INTEGER(KIND=JWIM) :: SDICE2_K, SDICE2_M
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    REAL(KIND=JWRB) :: EWH
    REAL(KIND=JWRB) :: SDICE2_ALP    !! ALP=SPATIAL ATTENUATION RATE OF ENERGY
    REAL(KIND=JWRB) :: SDICE2_FLDICE
    REAL(KIND=JWRB) :: SDICE2_DELTM, SDICE2_DELT5, SDICE2_DELT, SDICE2_GTEMP1

    REAL(KIND=JPHOOK) :: SDICE2_ZHOOK_HANDLE
!$acc routine seq
!$acc data present( FL1, FLD, SL, SLICE, WAVNUM, CGROUP, CICV, XK2 )


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



!$acc end data
  END SUBROUTINE SDICE2_LOKI
END MODULE SDICE2_LOKI_MOD
