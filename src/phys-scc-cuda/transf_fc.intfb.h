INTERFACE
  SUBROUTINE TRANSF_FC (XK, D, DKMAX, G)
    USE PARKIND_WAVE, ONLY: JWRU, JWIM, JWRB
    
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    REAL(KIND=JWRB), INTENT(IN) :: XK, D
    REAL(KIND=JWRB), INTENT(IN) :: DKMAX
    REAL(KIND=JWRB), INTENT(IN) :: G
!$acc routine seq
  END SUBROUTINE TRANSF_FC
END INTERFACE