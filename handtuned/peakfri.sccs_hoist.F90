! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE PEAKFRI_LOKI_MOD
  CONTAINS
  SUBROUTINE PEAKFRI_LOKI (KIJS, KIJL, F, IPEAKF, EPEAKF, F1D, IJ)

    ! ----------------------------------------------------------------------

    !**** *PEAKFRI* - COMPUTATION OF INDEX OF THE PEAK FREQUENCY AND
    !                 THE CORRESPONDING VALUE OF THE 1-D SPECTRAM
    !                 AT EACH GRID POINT.

    !     JEAN BIDLOT    ECMWF           MARCH 2000

    !*    PURPOSE.
    !     --------

    !       COMPUTE PEAK FREQUENCY INDEX AT EACH GRID POINT.
    !       IF NOT PEAK FOUND THAN EPEAKF WILL BE SET TO 0. AND
    !       IPEAKF TO NFRE

    !**   INTERFACE.
    !     ----------

    !       *CALL* *PEAKFRI (KIJS, KIJL, F, IPEAKF, EPEAKF, F1D)*
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT
    !          *KIJL*    - INDEX OF LAST GRIDPOINT
    !          *F*       - SPECTRUM.
    !          *IPEAKF   - INDEX OF PEAK FREQUENCY
    !          *EPEAKF*  - 1-D SPECTRUM AT PEAK FREQUENCY
    !          *F1D*     - 1D FREQUENCY SPECTRUM

    !     METHOD.
    !     -------

    !       NONE.

    !     EXTERNALS.
    !     ----------

    !       NONE.

    !     REFERENCE.
    !     ----------

    !       NONE.

    ! ----------------------------------------------------------------------

    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWFRED, ONLY: DELTH, FR
    USE YOWPARAM, ONLY: NANG, NFRE

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN) :: F(KIJL, NANG, NFRE)
    INTEGER(KIND=JWIM), INTENT(OUT) :: IPEAKF(KIJL)
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: EPEAKF
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NFRE) :: F1D

    INTEGER(KIND=JWIM) :: K, PEAKFRI_M
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    REAL(KIND=JPHOOK) :: PEAKFRI_ZHOOK_HANDLE
!$acc routine seq
!$acc data present( F, IPEAKF, EPEAKF, F1D )

    ! ----------------------------------------------------------------------


    !*    1. INITIALIZE ARRAYS
    !        -----------------

    EPEAKF(IJ) = 0.0_JWRB
    IPEAKF(IJ) = NFRE

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

      IF (EPEAKF(IJ) < F1D(IJ, PEAKFRI_M)) THEN
        EPEAKF(IJ) = F1D(IJ, PEAKFRI_M)
        IPEAKF(IJ) = PEAKFRI_M
      END IF

    END DO



!$acc end data
  END SUBROUTINE PEAKFRI_LOKI
END MODULE PEAKFRI_LOKI_MOD
