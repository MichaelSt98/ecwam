! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDICE3_LOKI_MOD
  CONTAINS
  SUBROUTINE SDICE3_LOKI (KIJS, KIJL, FL1, FLD, SL, SLICE, WAVNUM, CGROUP, CICV, CITH, ALPFAC, IJ, SDICE3_ALP)
    ! ----------------------------------------------------------------------

    !**** *SDICE3* - COMPUTATION OF SEA ICE ATTENUATION DUE TO VISCOUS FRICTION


    !     LOTFI AOUF       METEO FRANCE 2023
    !     JOSH KOUSAL      ECMWF 2023


    !*    PURPOSE.
    !     --------

    !**   INTERFACE.
    !     ----------

    !       *CALL* *SDICE3 (KIJS, KIJL, FL1, FLD,SL,
    !                       WAVNUM, CGROUP,
    !                       CICV,CITH, ALPFAC)*
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

    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWFRED, ONLY: FR
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: G, ZPI
    USE YOWICE, ONLY: ZALPFACB
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
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CITH
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: ALPFAC

    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: SDICE3_ALP    !! ALP=SPATIAL ATTENUATION RATE OF ENERGY

    INTEGER(KIND=JWIM) :: IMODEL    !! DAMPING MODEL: 1=FIT TO TEMPELFJORD DATA, 2=Jie Yu 2022
    INTEGER(KIND=JWIM) :: SDICE3_K, SDICE3_M
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    REAL(KIND=JWRB) :: SDICE3_FLDICE

    REAL(KIND=JWRB) :: SDICE3_DELTM, SDICE3_DELT5, SDICE3_DELT, SDICE3_GTEMP1

    REAL(KIND=JWRB) :: CDICE
    REAL(KIND=JWRB) :: HICEMAX, HICEMIN

    REAL(KIND=JWRB) :: TEMP

    REAL(KIND=JPHOOK) :: SDICE3_ZHOOK_HANDLE
!$acc routine seq
!$acc data present( FL1, FLD, SL, SLICE, WAVNUM, CGROUP, CICV, CITH, ALPFAC, SDICE3_ALP )


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
        TEMP = -CICV(IJ)*SDICE3_ALP(IJ, SDICE3_M)*CGROUP(IJ, SDICE3_M)
        SL(IJ, SDICE3_K, SDICE3_M) = SL(IJ, SDICE3_K, SDICE3_M) + FL1(IJ, SDICE3_K, SDICE3_M)*TEMP
        FLD(IJ, SDICE3_K, SDICE3_M) = FLD(IJ, SDICE3_K, SDICE3_M) + TEMP

        !              apply the source term (new way)
        !                SL(IJ,K,M)     =  SL(IJ,K,M)  + CICV(IJ)*SLICE(IJ,K,M)
        !                FLD(IJ,K,M)    =  FLD(IJ,K,M) + CICV(IJ)*FLDICE
        !
        !              to be used for wave radiative stress calculation
        SDICE3_GTEMP1 = MAX((1.0_JWRB - SDICE3_DELT5*SDICE3_FLDICE), 1.0_JWRB)
        SLICE(IJ, SDICE3_K, SDICE3_M) = SLICE(IJ, SDICE3_K, SDICE3_M) / SDICE3_GTEMP1

      END DO
    END DO



!$acc end data
  END SUBROUTINE SDICE3_LOKI
END MODULE SDICE3_LOKI_MOD
