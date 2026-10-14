! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE MEANSQS_LF_LOKI_MOD
  CONTAINS
  SUBROUTINE MEANSQS_LF_LOKI (NFRE_EFF, KIJS, KIJL, F, WAVNUM, XMSS, IJ, FD)

    ! ----------------------------------------------------------------------

    !**** *MEANSQS_LF* - COMPUTATION OF RESOLVED MEAN SQUARE SLOPE UP TO WAVE FREQUENCY FR(NFRE_EFF)

    !*    PURPOSE.
    !     --------

    !       COMPUTE MEAN SQUARE SLOPE AT EACH GRID POINT FOR THE RESOLVED SPECTRUM UP TO FREQUENCY FR(NFRE_EFF).
    !       NFRE_EFF <= NFRE

    !**   INTERFACE.
    !     ----------

    !       *CALL* *MEANSQS_LF (NFRE_EFF, KIJS, KIJL, F, WAVNUM, XMSS)*
    !              *XKMSS*   - WAVE NUMBER CUT OFF
    !              *KIJS*    - INDEX OF FIRST GRIDPOINT
    !              *KIJL*    - INDEX OF LAST GRIDPOINT
    !              *F*       - SPECTRUM.
    !              *WAVNUM*  - WAVE NUMBER.
    !              *XMSS*    - MEAN SQUARE SLOPE (OUTPUT).

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

    USE YOWPCONS, ONLY: G, GM1, ZPI
    USE YOWFRED, ONLY: FR, ZPIFR, DFIM
    USE YOWPARAM, ONLY: NANG, NFRE

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    ! ----------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: NFRE_EFF
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: F
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: XMSS


    INTEGER(KIND=JWIM) :: MEANSQS_LF_M, MEANSQS_LF_K, KFRE
    INTEGER(KIND=JWIM), INTENT(IN) :: IJ

    REAL(KIND=JPHOOK) :: MEANSQS_LF_ZHOOK_HANDLE
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NFRE) :: FD
    REAL(KIND=JWRB) :: TEMP1, TEMP2
!$acc routine seq
!$acc data present( F, WAVNUM, XMSS, FD )

    ! ----------------------------------------------------------------------

    !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
    !        ------------------------------------------

    KFRE = MIN(NFRE_EFF, NFRE)

    XMSS(IJ) = 0.0_JWRB

    !*    2.2 SHALLOW WATER INTEGRATION.
    !         --------------------------

!$acc loop seq
    DO MEANSQS_LF_M=1,KFRE
      TEMP1 = DFIM(MEANSQS_LF_M)*WAVNUM(IJ, MEANSQS_LF_M)**2
      TEMP2 = 0.0_JWRB
!$acc loop seq
      DO MEANSQS_LF_K=1,NANG
        TEMP2 = TEMP2 + F(IJ, MEANSQS_LF_K, MEANSQS_LF_M)
      END DO
      XMSS(IJ) = XMSS(IJ) + TEMP1*TEMP2
    END DO



!$acc end data
  END SUBROUTINE MEANSQS_LF_LOKI
END MODULE MEANSQS_LF_LOKI_MOD
