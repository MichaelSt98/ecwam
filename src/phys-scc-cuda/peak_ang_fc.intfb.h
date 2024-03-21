INTERFACE
  SUBROUTINE PEAK_ANG_FC (KIJS, KIJL, FL1, XNU, SIG_TH, COSTH, DELTH, DFIM, DFIMFR, DFIMFR2, FR, FRATIO, NANG, NFRE, SINTH, TH,  &
  & WETAIL, WP1TAIL, WP2TAIL, ICHNK, NCHNK, IJ)
    USE PARKIND_WAVE, ONLY: JWRU, JWIM, JWRB
    
    USE YOWFRED, ONLY: DFIMOFR
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FL1(:, :, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: XNU(:)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: SIG_TH(:)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: COSTH(:)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIM(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIMFR(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIMFR2(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FR(:)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: FRATIO
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: SINTH(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: TH(:)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WETAIL
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WP1TAIL
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WP2TAIL
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
  END SUBROUTINE PEAK_ANG_FC
END INTERFACE