MODULE TRANSF_SNL_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE TRANSF_SNL_fc (XK0, D, XNU, SIG_TH, BATHYMAX, DKMAX, G, XKDMIN)
    USE PARKIND_WAVE, ONLY: JWRB
    
    
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: XK0
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: D
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: XNU
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SIG_TH
    
    
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: BATHYMAX
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DKMAX
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: XKDMIN
!$acc routine seq
    INTERFACE
      SUBROUTINE TRANSF_SNL_iso_c (XK0, D, XNU, SIG_TH, BATHYMAX, DKMAX, G, XKDMIN) BIND(c, name="transf_snl_c_launch")
        implicit none
        REAL, VALUE :: XK0
        REAL, VALUE :: D
        REAL, VALUE :: XNU
        REAL, VALUE :: SIG_TH
        REAL, VALUE :: BATHYMAX
        REAL, VALUE :: DKMAX
        REAL, VALUE :: G
        REAL, VALUE :: XKDMIN
      END SUBROUTINE TRANSF_SNL_iso_c
    END INTERFACE
!$acc host_data use_device
    CALL TRANSF_SNL_iso_c(XK0, D, XNU, SIG_TH, BATHYMAX, DKMAX, G, XKDMIN)
!$acc end host_data
  END SUBROUTINE TRANSF_SNL_fc
END MODULE TRANSF_SNL_FC_MOD
