INTERFACE
  SUBROUTINE Z0WAVE_FC (KIJS, KIJL, US, TAUW, UTOP, Z0, Z0B, CHRNCK, ALPHA, ALPHAMIN, CHNKMIN_U, EPS1, G, GM1, LLCAPCHNK, ICHNK,  &
  & NCHNK, IJ)
    USE PARKIND_WAVE, ONLY: JWRU, JWIM, JWRB
    
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTERFACE
      FUNCTION CHNKMIN_FC (U10)
        USE parkind_wave, ONLY: jwrb
        REAL(KIND=JWRB) :: CHNKMIN
        REAL(KIND=JWRB), INTENT(IN) :: U10
      END FUNCTION CHNKMIN_FC
    END INTERFACE
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: US(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: TAUW(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: UTOP(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: Z0(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: Z0B(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: CHRNCK(:, :)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ALPHA
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ALPHAMIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: CHNKMIN_U
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPS1
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: GM1
    LOGICAL, VALUE, INTENT(IN) :: LLCAPCHNK
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
  END SUBROUTINE Z0WAVE_FC
END INTERFACE