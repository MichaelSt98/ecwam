INTERFACE
  SUBROUTINE FKMEAN_FC (KIJS, KIJL, FL1, WAVNUM, EM, FM1, F1, AK, XK, DELTH, DFIM, DFIMFR, DFIMOFR, EPSMIN, FR, FRTAIL, G, NANG,  &
  & NFRE, WETAIL, WP1TAIL, ZPI, ICHNK, NCHNK, IJ)
    USE PARKIND_WAVE, ONLY: JWRU, JWIM, JWRB
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FL1(:, :, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: WAVNUM(:, :, :)
    
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: EM(:)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: FM1(:)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: F1(:)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: AK(:)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: XK(:)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIM(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIMFR(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIMOFR(:)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPSMIN
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FR(:)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: FRTAIL
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WETAIL
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WP1TAIL
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZPI
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
  END SUBROUTINE FKMEAN_FC
END INTERFACE