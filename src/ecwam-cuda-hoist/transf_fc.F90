MODULE TRANSF_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE TRANSF_fc (XK, D, DKMAX, G)
    USE PARKIND_WAVE, ONLY: JWRB
    
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: XK
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: D
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DKMAX
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
!$acc routine seq
    INTERFACE
      SUBROUTINE TRANSF_iso_c (XK, D, DKMAX, G) BIND(c, name="transf_c_launch")
        implicit none
        REAL, VALUE :: XK
        REAL, VALUE :: D
        REAL, VALUE :: DKMAX
        REAL, VALUE :: G
      END SUBROUTINE TRANSF_iso_c
    END INTERFACE
!$acc host_data use_device
    CALL TRANSF_iso_c(XK, D, DKMAX, G)
!$acc end host_data
  END SUBROUTINE TRANSF_fc
END MODULE TRANSF_FC_MOD
