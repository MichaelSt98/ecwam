! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
ATTRIBUTES(DEVICE) SUBROUTINE HALPHAP_FC (KIJS, KIJL, WAVNUM, COSWDIF, FL1, HALP, ALPHAPMAX, DELTH, DFIM, DFIMOFR, EPSMIN, FR,  &
& FR5, FRTAIL, NANG, NFRE, WETAIL, ZPI4GM2, ICHNK, NCHNK, IJ)
  
  ! ----------------------------------------------------------------------
  
  !**** *HALPHAP* - COMPUTATION OF 1/2 PHILLIPS PARAMETER
  
  
  !**   INTERFACE.
  !     ----------
  
  !       *CALL* *HALPHAP(KIJS, KIJL, WAVNUM, UDIR, FL1, HALP)
  !          *KIJS*   - INDEX OF FIRST GRIDPOINT
  !          *KIJL*   - INDEX OF LAST GRIDPOINT
  !          *WAVNUM* - WAVE NUMBER
  !          *COSWDIF*- COSINE ( WIND SPEED DIRECTION - WAVE DIRECTIONS)
  !          *FL1*    - SPECTRA
  !          *HALP*   - 1/2 PHILLIPS PARAMETER
  
  !     METHOD.
  !     -------
  
  ! ----------------------------------------------------------------------
  
  USE PARKIND_WAVE, ONLY: JWRU, JWIM, JWRB
  
  USE YOWFRED, ONLY: TH
  USE YOWPARAM, ONLY: NANG_PARAM
  USE YOWPCONS, ONLY: G, ZPI
  
  
  ! ----------------------------------------------------------------------
  
  IMPLICIT NONE
  
  INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
  INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: WAVNUM(:, :, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: COSWDIF(:, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: FL1(:, :, :, :)
  REAL(KIND=JWRB), TARGET, INTENT(OUT) :: HALP(:)
  
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
  INTEGER(KIND=JWIM) :: K
  INTEGER(KIND=JWIM) :: M
  
  REAL(KIND=JWRB) :: ZLNFRNFRE
  REAL(KIND=JWRB) :: DELT25
  REAL(KIND=JWRB) :: DELT2
  REAL(KIND=JWRB) :: DEL2
  REAL(KIND=JWRB) :: TEMP1
  REAL(KIND=JWRB) :: TEMP2
  REAL(KIND=JWRB) :: ALPHAP
  REAL(KIND=JWRB) :: XMSS
  REAL(KIND=JWRB) :: EM
  REAL(KIND=JWRB) :: FM
  REAL(KIND=JWRB) :: F1D
  REAL(KIND=JWRB) :: FLWD(NANG_PARAM)
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: ALPHAPMAX
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIM(:)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIMOFR(:)
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPSMIN
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: FR(:)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: FR5(:)
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: FRTAIL
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: WETAIL
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZPI4GM2
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
  INTEGER, VALUE, INTENT(IN) :: NCHNK
  
  ! ----------------------------------------------------------------------
  
  
  ZLNFRNFRE = LOG(FR(NFRE))
  
  DELT25 = WETAIL*FR(NFRE)*DELTH
  DELT2 = FRTAIL*DELTH
  
  ! Find spectrum in wind direction
  
  DO M=1,NFRE
    DO K=1,NANG
      FLWD(K) = FL1(IJ, K, M, ICHNK)*0.5_JWRB + 0.5_JWRB*SIGN(1.0_JWRB, COSWDIF(IJ, K))
    END DO
    
    XMSS = 0._JWRB
    TEMP1 = DFIM(M)*WAVNUM(IJ, M, ICHNK)**2
    TEMP2 = 0.0_JWRB
    DO K=1,NANG
      TEMP2 = TEMP2 + FLWD(K)
    END DO
    XMSS = XMSS + TEMP1*TEMP2
    
    K = 1
    EM = 0._JWRB
    FM = 0._JWRB
    TEMP2 = MAX(FLWD(K), EPSMIN)
    DO K=2,NANG
      TEMP2 = TEMP2 + MAX(FLWD(K), EPSMIN)
    END DO
    EM = EM + TEMP2*DFIM(M)
    FM = FM + DFIMOFR(M)*TEMP2
  END DO
  
  DO K=1,NANG
    FLWD(K) = FL1(IJ, K, NFRE, ICHNK)*0.5_JWRB + 0.5_JWRB*SIGN(1.0_JWRB, COSWDIF(IJ, K))
  END DO
  
  EM = EM + DELT25*TEMP2
  FM = FM + DELT2*TEMP2
  FM = EM / FM
  FM = MAX(FM, FR(1))
  
  IF (EM > 0.0_JWRB .and. FM < FR(-2 + NFRE)) THEN
    ALPHAP = XMSS / (ZLNFRNFRE - LOG(FM))
    IF (ALPHAP > ALPHAPMAX) THEN
      ! some odd cases, revert to tail value
      F1D = 0.0_JWRB
      DO K=1,NANG
        F1D = F1D + FLWD(K)*DELTH
      END DO
      ALPHAP = ZPI4GM2*FR5(NFRE)*F1D
    END IF
  ELSE
    F1D = 0.0_JWRB
    DO K=1,NANG
      F1D = F1D + FLWD(K)*DELTH
    END DO
    ALPHAP = ZPI4GM2*FR5(NFRE)*F1D
  END IF
  
  !     1/2 ALPHAP:
  HALP(IJ) = 0.5_JWRB*MIN(ALPHAP, ALPHAPMAX)
  
  
  
END SUBROUTINE HALPHAP_FC
