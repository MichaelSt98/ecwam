! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
ATTRIBUTES(DEVICE) SUBROUTINE CIWABR_FC (KIJS, KIJL, CICOVER, FL1, WAVNUM, CGROUP, CIWAB, CDICWA, DFIM, EPSMIN, IDELT, LICERUN,  &
& LMASKICE, NANG, NFRE, ICHNK, NCHNK, IJ)
  
  ! ----------------------------------------------------------------------
  
  !**** *CIWABR* - COMPUTE SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
  !                BOTTOM FRICTION.
  
  !*    PURPOSE.
  !     --------
  
  !       CIWABR COMPUTES SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
  !              BOTTOM FRICTION.
  
  !**   INTERFACE.
  !     ----------
  
  !       *CALL* *CIWABR (KIJS,KIJL,CICOVER,FL1,CIWAB)
  
  !          *KIJS*     - INDEX OF FIRST POINT.
  !          *KIJL*     - INDEX OF LAST POINT.
  !          *CICOVER*  -SEA ICE COVER.
  !          *FL1*      -ENERGY SPECTRUM.
  !          *CIWAB*    -SEA ICE WAVE ATTENUATION FACTOR DUE TO ICE FLOE BOTTOM FRICTION
  
  !     METHOD.
  !     -------
  
  !     EXTERNALS.
  !     ----------
  
  !     REFERENCES.
  !     -----------
  
  !     KOHOUT A., M. MEYLAN, D PLEW, 2011: ANNALS OF GLACIOLOGY, 2011.
  
  
  ! ----------------------------------------------------------------------
  
  USE PARKIND_WAVE, ONLY: JWRU, JWIM, JWRB
  
  USE YOWFRED, ONLY: FR, DELTH
  USE YOWPCONS, ONLY: ZPI4GM2, G, ZPI
  
  
  ! ----------------------------------------------------------------------
  IMPLICIT NONE
  
  INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
  INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: CICOVER(:, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: FL1(:, :, :, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: WAVNUM(:, :, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: CGROUP(:, :, :)
  REAL(KIND=JWRB), TARGET, INTENT(OUT) :: CIWAB(:, :, :)
  
  
  INTEGER(KIND=JWIM) :: K
  INTEGER(KIND=JWIM) :: M
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
  REAL(KIND=JWRB) :: EWH
  REAL(KIND=JWRB) :: X
  REAL(KIND=JWRB) :: ALP
  REAL(KIND=JWRB) :: XK2
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: CDICWA
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIM(:)
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPSMIN
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IDELT
  LOGICAL, VALUE, INTENT(IN) :: LICERUN
  LOGICAL, VALUE, INTENT(IN) :: LMASKICE
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
  INTEGER, VALUE, INTENT(IN) :: NCHNK
  
  ! ----------------------------------------------------------------------
  
  
  
  IF (.not.LICERUN .or. LMASKICE) THEN
    
    DO M=1,NFRE
      DO K=1,NANG
        CIWAB(IJ, K, M) = 1.0_JWRB
      END DO
    END DO
    
  ELSE
    
    DO M=1,NFRE
      DO K=1,NANG
        EWH = 4.0_JWRB*SQRT(MAX(EPSMIN, FL1(IJ, K, M, ICHNK)*DFIM(M)))
        XK2 = WAVNUM(IJ, M, ICHNK)**2
        ALP = CDICWA*XK2*EWH
        X = ALP*CGROUP(IJ, M, ICHNK)*IDELT
        CIWAB(IJ, K, M) = 1.0_JWRB - CICOVER(IJ, ICHNK)*(1.0_JWRB - EXP(-MIN(X, 50.0_JWRB)))
      END DO
    END DO
    
  END IF
  
  
  
END SUBROUTINE CIWABR_FC
