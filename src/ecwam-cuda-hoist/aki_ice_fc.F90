MODULE AKI_ICE_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE AKI_ICE_fc (G, XK, DEPTH, RHOW, CITH)
    USE PARKIND_WAVE, ONLY: JWRB
    
    IMPLICIT NONE
    
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: XK
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DEPTH
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: RHOW
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: CITH
    
    !     ICE PROPERTIES (assumed fixed for now)
    
    !     RELATIVE ERROR LIMIT OF NEWTONS METHOD.
    !     MAXIMUM WAVE NUMBER
    
!$acc routine seq
    INTERFACE
      SUBROUTINE AKI_ICE_iso_c (G, XK, DEPTH, RHOW, CITH) BIND(c, name="aki_ice_c_launch")
        implicit none
        REAL, VALUE :: G
        REAL, VALUE :: XK
        REAL, VALUE :: DEPTH
        REAL, VALUE :: RHOW
        REAL, VALUE :: CITH
      END SUBROUTINE AKI_ICE_iso_c
    END INTERFACE
!$acc host_data use_device
    CALL AKI_ICE_iso_c(G, XK, DEPTH, RHOW, CITH)
!$acc end host_data
  END SUBROUTINE AKI_ICE_fc
END MODULE AKI_ICE_FC_MOD
