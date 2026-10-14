! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE WDIRSPREAD_LOKI_MOD
  CONTAINS
  SUBROUTINE WDIRSPREAD_LOKI (KIJS, KIJL, F, EMEAN, LLPEAKF, WDIRSPRD, IJ, F1D)

    ! ----------------------------------------------------------------------

    !**** *WDIRSPREAD* - COMPUTATION OF THE DIRECTIONAL SPREAD

    !     JEAN BIDLOT    ECMWF           MARCH 2000

    !*    PURPOSE.
    !     --------

    !     TO COMPUTE THE DIRECTIONAL SPREAD EITHER FOR THE WHOLE FREQUENCY
    !     RANGE OR BASED ON THE PEAK FREQUENCY.

    !**   INTERFACE.
    !     ----------

    !       *CALL* *WDIRSPREAD (KIJS, KIJL, F, LLPEAKF, WDIRSPRD)*
    !          *KIJS*     - INDEX OF FIRST GRIDPOINT.
    !          *KIJL*     - INDEX OF LAST GRIDPOINT.
    !          *F*        - SPECTRUM.
    !          *EMEAN*    - MEAN WAVE ENERGY FOR THE WAVE SYSTEM (INPUT).
    !          *LLPEAKF*  - TRUE IF THE DIRECTIONAL SPREAD IS BASED ON
    !                       THE PEAK FREQUENCY.
    !          *WDIRSPRD* - DIRECTIONAL SPREAD.

    !     METHOD.
    !     -------

    !       NONE.

    !     EXTERNALS.
    !     ----------

    !      *SCOSFL*
    !      *PEAKFRI*

    !     REFERENCE.
    !     ----------

    !       ORIGINALLY REQUESTED BY LUIGI CAVALERI. WE WILL TRY TO PRODUCE
    !       SOME FORM OF REFERENCE.

    ! ----------------------------------------------------------------------

    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWFRED, ONLY: FR, DFIM, DELTH, WETAIL, TH, COSTH, SINTH
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: EPSMIN

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    USE, INTRINSIC :: IEEE_EXCEPTIONS

    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: F
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: EMEAN
    LOGICAL, INTENT(IN) :: LLPEAKF
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: WDIRSPRD


    INTEGER(KIND=JWIM) :: M
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: IFRINDEX
    REAL(KIND=JWRB) :: COEF_FR, ONE
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    REAL(KIND=JWRB) :: TEMP
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NFRE) :: F1D

    LOGICAL :: LL_HALT_INVALID
    INTEGER(KIND=JWIM) :: K, PEAKFRI_M
    REAL(KIND=JPHOOK) :: PEAKFRI_ZHOOK_HANDLE
    INTEGER(KIND=JWIM) :: SCOSFL_K
    REAL(KIND=JPHOOK) :: SCOSFL_ZHOOK_HANDLE
    REAL(KIND=JWRB) :: MEANDIR, SI, CI
!$acc routine seq
!$acc data present( F, EMEAN, WDIRSPRD, F1D )

    ! ----------------------------------------------------------------------


    ! Turn off Floating-Point-Exceptions in this scope to avoid FPE_INVALID in optimized code
    !   with branch prediction. It is safe to do so as DIV_BY_ZERO is protected.





    !*    1. INITIALIZE ARRAYS
    !        -----------------

    ONE = 1.0_JWRB
    COEF_FR = WETAIL*FR(NFRE)

    WDIRSPRD(IJ) = 0.0_JWRB

    IF (LLPEAKF) THEN
      !     COMPUTATION IS BASED ON THE PEAK FREQUENCY
      ! [Loki] inlined child subroutine: PEAKFRI
      ! =========================================

      ! ----------------------------------------------------------------------


      !*    1. INITIALIZE ARRAYS
      !        -----------------

      TEMP = 0.0_JWRB
      IFRINDEX = NFRE

      !*    2. LOOP OVER FREQUENCIES
      !        ---------------------

!$acc loop seq
      DO PEAKFRI_M=1,NFRE

        !*    2.1 COMPUTE 1-D SPECTRUM
        !        ---------------------

        F1D(IJ, PEAKFRI_M) = 0.0_JWRB
!$acc loop seq
        DO K=1,NANG
          F1D(IJ, PEAKFRI_M) = F1D(IJ, PEAKFRI_M) + F(IJ, K, PEAKFRI_M)*DELTH
        END DO

        !*    2.2 DEFINE PEAK FREQUENCY
        !         ---------------------

        IF (TEMP < F1D(IJ, PEAKFRI_M)) THEN
          TEMP = F1D(IJ, PEAKFRI_M)
          IFRINDEX = PEAKFRI_M
        END IF

      END DO


      ! =========================================
      ! [Loki] inlined child subroutine: SCOSFL
      ! =========================================

      ! ----------------------------------------------------------------------


      !*    1. INITIALISE ARRAYS
      !        -----------------

      SI = 0.0_JWRB
      CI = 0.0_JWRB
      WDIRSPRD(IJ) = 0.0_JWRB

      ! ----------------------------------------------------------------------

      !*    2. FIND MEAN DIRECTION
      !        -------------------

!$acc loop seq
      DO SCOSFL_K=1,NANG
        SI = SI + SINTH(SCOSFL_K)*F(IJ, SCOSFL_K, IFRINDEX)
        CI = CI + COSTH(SCOSFL_K)*F(IJ, SCOSFL_K, IFRINDEX)
      END DO

      IF (CI == 0.0_JWRB .and. SI == 0.0_JWRB) THEN
        MEANDIR = 0.0_JWRB
      ELSE
        MEANDIR = ATAN2(SI, CI)
      END IF

      !*    3. COMPUTE FOURRIER COEFFICIENT
      !        -----------------------------

!$acc loop seq
      DO SCOSFL_K=1,NANG
        WDIRSPRD(IJ) = WDIRSPRD(IJ) + COS(TH(SCOSFL_K) - MEANDIR)*F(IJ, SCOSFL_K, IFRINDEX)
      END DO
      WDIRSPRD(IJ) = DELTH*WDIRSPRD(IJ)


      ! =========================================
      IF (TEMP > 0.0_JWRB) THEN
        WDIRSPRD(IJ) = MIN(WDIRSPRD(IJ) / TEMP, ONE)
      ELSE
        WDIRSPRD(IJ) = ONE
      END IF

    ELSE
      !     COMPUTATION IS BASED ON THE WHOLE FREQUENCY RANGE
!$acc loop seq
      DO M=1,NFRE
        IFRINDEX = M
        ! [Loki] inlined child subroutine: SCOSFL
        ! =========================================

        ! ----------------------------------------------------------------------


        !*    1. INITIALISE ARRAYS
        !        -----------------

        SI = 0.0_JWRB
        CI = 0.0_JWRB
        TEMP = 0.0_JWRB

        ! ----------------------------------------------------------------------

        !*    2. FIND MEAN DIRECTION
        !        -------------------

!$acc loop seq
        DO SCOSFL_K=1,NANG
          SI = SI + SINTH(SCOSFL_K)*F(IJ, SCOSFL_K, IFRINDEX)
          CI = CI + COSTH(SCOSFL_K)*F(IJ, SCOSFL_K, IFRINDEX)
        END DO

        IF (CI == 0.0_JWRB .and. SI == 0.0_JWRB) THEN
          MEANDIR = 0.0_JWRB
        ELSE
          MEANDIR = ATAN2(SI, CI)
        END IF

        !*    3. COMPUTE FOURRIER COEFFICIENT
        !        -----------------------------

!$acc loop seq
        DO SCOSFL_K=1,NANG
          TEMP = TEMP + COS(TH(SCOSFL_K) - MEANDIR)*F(IJ, SCOSFL_K, IFRINDEX)
        END DO
        TEMP = DELTH*TEMP


        ! =========================================
        WDIRSPRD(IJ) = WDIRSPRD(IJ) + TEMP*DFIM(M)
      END DO

      !       ADD TAIL CONTRIBUTION
      !       note that it uses TEMP because the last call to SCOSFL is made
      !         the division by delth is needed since dfim contains delth
      WDIRSPRD(IJ) = WDIRSPRD(IJ) / DELTH + TEMP*COEF_FR

      !       NORMALISATION
      IF (EMEAN(IJ) > EPSMIN) THEN
        WDIRSPRD(IJ) = MIN(WDIRSPRD(IJ) / EMEAN(IJ), ONE)
      ELSE
        WDIRSPRD(IJ) = ONE
      END IF

    END IF

    !     COMPUTE THE ACTUAL SPREAD

    WDIRSPRD(IJ) = SQRT(2.0_JWRB*(ONE - WDIRSPRD(IJ)))







!$acc end data
  END SUBROUTINE WDIRSPREAD_LOKI
END MODULE WDIRSPREAD_LOKI_MOD
