INTERFACE
  SUBROUTINE TRANSF_SNL_FC (XK0, D, XNU, SIG_TH, BATHYMAX, DKMAX, G, XKDMIN)
    USE PARKIND_WAVE, ONLY: JWRU, JWIM, JWRB
    
    
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    REAL(KIND=JWRB), INTENT(IN) :: XK0, D, XNU, SIG_TH
    
    
    REAL(KIND=JWRB), INTENT(IN) :: BATHYMAX
    REAL(KIND=JWRB), INTENT(IN) :: DKMAX
    REAL(KIND=JWRB), INTENT(IN) :: G
    REAL(KIND=JWRB), INTENT(IN) :: XKDMIN
!$acc routine seq
  END SUBROUTINE TRANSF_SNL_FC
END INTERFACE