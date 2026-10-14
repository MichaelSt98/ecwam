! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE STOKESTRN_LOKI_MOD
  CONTAINS
  SUBROUTINE STOKESTRN_LOKI (KIJS, KIJL, FL1, STOKFAC, WSWAVE, WDWAVE, CICOVER, USTOKES, VSTOKES, STRNMS, NEMOUSTOKES,  &
  & NEMOVSTOKES, NEMOSTRN, IJ)

    ! ----------------------------------------------------------------------

    !**** *STOKESTRN* - WRAPPER TO CALL STOKESDRIFT and CIMSSTRN

    !*    PURPOSE.
    !     --------

    !**   INTERFACE.
    !     ----------

    !       *CALL* *STOKESTRN (KIJS, KIJL, FL1, WAVNUM, STOKFAC, DEPTH, FF_NOW, INTFLDS, WAM2NEMO)

    !          *KIJS*    - INDEX OF FIRST GRIDPOINT.
    !          *KIJL*    - INDEX OF LAST GRIDPOINT.
    !          *FL1*     - SPECTRUM(INPUT).
    !          *WAVNUM*  - WAVE NUMBER.
    !          *STOKFAC* - STOKES DRIFT FACTOR.
    !          *DEPTH*   - WATER DEPTH.
    !          *FF_NOW*  - FORCING FIELDS AT CURRENT TIME.
    !          *INTFLDS* - INTEGRATED/DERIVED PARAMETERS
    !          *WAM2NEMO*- WAVE FIELDS PASSED TO NEMO

    ! ----------------------------------------------------------------------

    USE YOWPCONS, ONLY: G, ZPI
    USE YOWFRED, ONLY: FR, DFIM, DELTH, TH, DFIM_SIM, FRATIO, COSTH, SINTH
    USE YOWICE, ONLY: LICERUN, LWAMRSETCI, CITHRSH
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU, JWRO
    USE YOWDRVTYPE, ONLY: FORCING_FIELDS, INTGT_PARAM_FIELDS, WAVE2OCEAN

    USE YOWCOUP, ONLY: LWCOU, LWNEMOCOU, LWNEMOCOUSEND, LWNEMOCOUSTK, LWNEMOCOUSTRN
    USE YOWPARAM, ONLY: NANG, NFRE, NFRE_ODD

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE
    INTERFACE
      SUBROUTINE CIMSSTRN (KIJS, KIJL, FL1, WAVNUM, DEPTH, CITHICK, STRN)
        USE parkind_wave, ONLY: jwim, jwrb
        USE YOWPARAM, ONLY: NANG, NFRE
        INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
        REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
        REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM
        REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: DEPTH
        REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CITHICK
        REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: STRN
      END SUBROUTINE CIMSSTRN
    END INTERFACE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: STOKFAC
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: WSWAVE, WDWAVE, CICOVER
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: USTOKES, VSTOKES, STRNMS
    REAL(KIND=JWRO), INTENT(INOUT), DIMENSION(KIJL) :: NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN


    REAL(KIND=JPHOOK) :: STOKESTRN_ZHOOK_HANDLE
    INTEGER(KIND=JWIM) :: STOKESTRN_M, STOKESTRN_K
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    REAL(KIND=JWRB), PARAMETER :: STMAX = 1.5_JWRB    ! maximum magnitude (this is for safety when coupled)
    REAL(KIND=JWRB) :: STOKESTRN_CONST, FAC, FAC1, FAC2, FAC3
    REAL(KIND=JPHOOK) :: STOKESDRIFT_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: STFAC
!$acc routine seq
!$acc data present( FL1, STOKFAC, WSWAVE, WDWAVE, CICOVER, USTOKES, VSTOKES, STRNMS, NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN )

    ! ----------------------------------------------------------------------


    ! [Loki] inlined child subroutine: STOKESDRIFT
    ! =========================================

    ! ----------------------------------------------------------------------



    !***  1. DETERMINE STOKE DRIFT VECTOR.
    !     --------------------------------

    STOKESTRN_CONST = 2.0_JWRB*DELTH*ZPI**3 / G*FR(NFRE_ODD)**4

    !***  1.1 PERFORM INTEGRATION.
    !     ------------------------

    USTOKES(IJ) = 0.0_JWRB
    VSTOKES(IJ) = 0.0_JWRB

!$acc loop seq
    DO STOKESTRN_M=1,NFRE_ODD
      STFAC = STOKFAC(IJ, STOKESTRN_M)*DFIM_SIM(STOKESTRN_M)
!$acc loop seq
      DO STOKESTRN_K=1,NANG
        FAC3 = STFAC*FL1(IJ, STOKESTRN_K, STOKESTRN_M)
        USTOKES(IJ) = USTOKES(IJ) + FAC3*SINTH(STOKESTRN_K)
        VSTOKES(IJ) = VSTOKES(IJ) + FAC3*COSTH(STOKESTRN_K)
      END DO
    END DO

    !***  1.2 ADD CONTRIBUTION OF UNRESOLVED WAVES.
    !     -----------------------------------------

!$acc loop seq
    DO STOKESTRN_K=1,NANG
      FAC1 = STOKESTRN_CONST*SINTH(STOKESTRN_K)
      FAC2 = STOKESTRN_CONST*COSTH(STOKESTRN_K)
      USTOKES(IJ) = USTOKES(IJ) + FAC1*FL1(IJ, STOKESTRN_K, NFRE_ODD)
      VSTOKES(IJ) = VSTOKES(IJ) + FAC2*FL1(IJ, STOKESTRN_K, NFRE_ODD)
    END DO


    !***  1.3 Sea Ice exception
    !     ---------------------
    IF (LICERUN .and. LWAMRSETCI) THEN
      IF (CICOVER(IJ) > CITHRSH) THEN
        USTOKES(IJ) = 0.016_JWRB*WSWAVE(IJ)*SIN(WDWAVE(IJ))*(1.0_JWRB - CICOVER(IJ))
        VSTOKES(IJ) = 0.016_JWRB*WSWAVE(IJ)*COS(WDWAVE(IJ))*(1.0_JWRB - CICOVER(IJ))
      END IF
    END IF

    !***  1.4 Protection
    !     --------------

    USTOKES(IJ) = MIN(MAX(USTOKES(IJ), -STMAX), STMAX)
    VSTOKES(IJ) = MIN(MAX(VSTOKES(IJ), -STMAX), STMAX)


    ! =========================================

    ! [Loki] Removed content of pragma-marked region!


    IF (LWNEMOCOU .and. (LWNEMOCOUSEND .and. LWCOU .or. .not.LWCOU)) THEN
      IF (LWNEMOCOUSTK) THEN
        NEMOUSTOKES(IJ) = USTOKES(IJ)
        NEMOVSTOKES(IJ) = VSTOKES(IJ)
      ELSE
        NEMOUSTOKES(IJ) = 0.0_JWRO
        NEMOVSTOKES(IJ) = 0.0_JWRO
      END IF

      IF (LWNEMOCOUSTRN) NEMOSTRN(IJ) = STRNMS(IJ)
    END IF


    ! ----------------------------------------------------------------------


!$acc end data
  END SUBROUTINE STOKESTRN_LOKI
END MODULE STOKESTRN_LOKI_MOD
