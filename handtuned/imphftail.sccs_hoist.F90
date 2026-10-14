! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE IMPHFTAIL_LOKI_MOD
  CONTAINS
  SUBROUTINE IMPHFTAIL_LOKI (KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1, IJ)
    ! ----------------------------------------------------------------------

    !**** *IMPHFTAIL* - IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM


    !*    PURPOSE.
    !     --------

    !     IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM ABOVE FREQUENCY INDEX MIJ


    !**   INTERFACE.
    !     ----------

    !       *CALL* *IMPHFTAIL (KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1)
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT
    !          *KIJL*    - INDEX OF LAST GRIDPOINT
    !          *MIJ*     - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
    !          *FLM*     - SPECTAL DENSITY MINIMUM VALUE
    !          *WAVNUM*  - WAVENUMBER
    !          *XK2CG*   - (WAVNUM)**2 * GROUP SPEED
    !          *FL1*     - SPECTRUM (INPUT AND OUTPUT).

    !     METHOD.
    !     -------

    !     EXTERNALS.
    !     ---------

    !     REFERENCE.
    !     ----------

    ! ----------------------------------------------------------------------
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOMHOOK, ONLY: LHOOK, JPHOOK
    USE YOWPARAM, ONLY: NANG, NFRE

    ! ----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(KIJL) :: MIJ
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: FLM
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM, XK2CG
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FL1


    INTEGER(KIND=JWIM) :: IMPHFTAIL_K, IMPHFTAIL_M
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ

    REAL(KIND=JWRB) :: AKM1, TFAC
    REAL(KIND=JPHOOK) :: IMPHFTAIL_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: IMPHFTAIL_TEMP1, IMPHFTAIL_TEMP2
!$acc routine seq
!$acc data present( MIJ, FLM, WAVNUM, XK2CG, FL1 )

    ! ----------------------------------------------------------------------


    !*    DIAGNOSTIC TAIL.
    !     ----------------

    IMPHFTAIL_TEMP1 = 1.0_JWRB / XK2CG(IJ, MIJ(IJ)) / WAVNUM(IJ, MIJ(IJ))

!$acc loop seq
    DO IMPHFTAIL_M=MIJ(IJ) + 1,NFRE
      IMPHFTAIL_TEMP2 = 1.0_JWRB / XK2CG(IJ, IMPHFTAIL_M) / WAVNUM(IJ, IMPHFTAIL_M)
      IMPHFTAIL_TEMP2 = IMPHFTAIL_TEMP2 / IMPHFTAIL_TEMP1

      !*    MERGE TAIL INTO SPECTRA.
      !     ------------------------
!$acc loop seq
      DO IMPHFTAIL_K=1,NANG
        TFAC = FL1(IJ, IMPHFTAIL_K, MIJ(IJ))
        FL1(IJ, IMPHFTAIL_K, IMPHFTAIL_M) = MAX(IMPHFTAIL_TEMP2*TFAC, FLM(IJ, IMPHFTAIL_K))
      END DO
    END DO

    ! ----------------------------------------------------------------------



!$acc end data
  END SUBROUTINE IMPHFTAIL_LOKI
END MODULE IMPHFTAIL_LOKI_MOD
