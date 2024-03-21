INTERFACE
  SUBROUTINE SETICE_FC (KIJS, KIJL, FL1, CICOVER, COSWDIF, CITHRSH, EPSMIN, FLMIN, NANG, NFRE, ICHNK, NCHNK, IJ)
    USE PARKIND_WAVE, ONLY: JWRU, JWIM, JWRB
    
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: FL1(:, :, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: CICOVER(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: COSWDIF(:, :)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: CITHRSH
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPSMIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: FLMIN
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
  END SUBROUTINE SETICE_FC
END INTERFACE