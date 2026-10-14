! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDEPTHLIM_LOKI_MOD
  CONTAINS
  SUBROUTINE SDEPTHLIM_LOKI (KIJS, KIJL, EMAXDPT, FL1, IJ)
    ! ----------------------------------------------------------------------
    !     J. BIDLOT    ECMWF   NOVEMBER 2017

    !*    PURPOSE.
    !     --------
    !     LIMITS THE SPECTRAL VARIANCE SUCH THAT THE TOTAL VARIANCE
    !     DOES NOT EXCEED THE MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH

    !**   INTERFACE.
    !     ----------
    !     *CALL* *SDEPTHLIM((KIJS, KIJL, EMAXDPT, FL1)
    !          *KIJS*   - LOCAL INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - LOCAL  INDEX OF LAST GRIDPOIN
    !          *EMAXDPT - MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
    !          *FL1*    - SPECTRUM.


    !     METHOD.
    !     -------

    !     EXTERNALS.
    !     ----------

    !     REFERENCE.
    !     ----------
    !     NONE

    ! ----------------------------------------------------------------------
    USE YOWFRED, ONLY: FR, DFIM, DELTH, WETAIL
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: EPSMIN

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: EMAXDPT
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FL1

    INTEGER(KIND=JWIM) :: SDEPTHLIM_K, SDEPTHLIM_M
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    REAL(KIND=JPHOOK) :: SDEPTHLIM_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: EM
    LOGICAL :: LLEPSMIN
    INTEGER(KIND=JWIM) :: SEMEAN_M, SEMEAN_K
    REAL(KIND=JWRB) :: DELT25
    REAL(KIND=JPHOOK) :: SEMEAN_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: SDEPTHLIM_TEMP
!$acc routine seq
!$acc data present( EMAXDPT, FL1 )

    ! ----------------------------------------------------------------------


    LLEPSMIN = .true.
    ! [Loki] inlined child subroutine: SEMEAN
    ! =========================================

    ! ----------------------------------------------------------------------


    !*    1. INITIALISE ENERGY ARRAY.
    !        ------------------------

    IF (LLEPSMIN) THEN
      EM = EPSMIN
    ELSE
      EM = 0.0_JWRB
    END IF

    ! ----------------------------------------------------------------------

    !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
    !        -----------------------------------------

!$acc loop seq
    DO SEMEAN_M=1,NFRE
      SEMEAN_K = 1
      SDEPTHLIM_TEMP = FL1(IJ, SEMEAN_K, SEMEAN_M)
!$acc loop seq
      DO SEMEAN_K=2,NANG
        SDEPTHLIM_TEMP = SDEPTHLIM_TEMP + FL1(IJ, SEMEAN_K, SEMEAN_M)
      END DO
      EM = EM + DFIM(SEMEAN_M)*SDEPTHLIM_TEMP
    END DO

    ! ----------------------------------------------------------------------

    !*    3. ADD TAIL ENERGY.
    !        ----------------

    DELT25 = WETAIL*FR(NFRE)*DELTH
    EM = EM + DELT25*SDEPTHLIM_TEMP


    ! =========================================

    EM = MIN(EMAXDPT(IJ) / EM, 1.0_JWRB)

!$acc loop seq
    DO SDEPTHLIM_M=1,NFRE
!$acc loop seq
      DO SDEPTHLIM_K=1,NANG
        FL1(IJ, SDEPTHLIM_K, SDEPTHLIM_M) = MAX(FL1(IJ, SDEPTHLIM_K, SDEPTHLIM_M)*EM, EPSMIN)
      END DO
    END DO



!$acc end data
  END SUBROUTINE SDEPTHLIM_LOKI
END MODULE SDEPTHLIM_LOKI_MOD
