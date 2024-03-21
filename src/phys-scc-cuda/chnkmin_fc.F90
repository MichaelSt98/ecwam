MODULE CHNKMIN_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE CHNKMIN_fc (U10, ALPHA, ALPHAMIN, CHNKMIN_U)
    USE PARKIND_WAVE, ONLY: JWRB
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    REAL(KIND=JWRB), INTENT(IN) :: U10
    REAL(KIND=JWRB), INTENT(IN) :: ALPHA
    REAL(KIND=JWRB), INTENT(IN) :: ALPHAMIN
    REAL(KIND=JWRB), INTENT(IN) :: CHNKMIN_U
!$acc routine seq
    INTERFACE
      SUBROUTINE CHNKMIN_iso_c (U10, ALPHA, ALPHAMIN, CHNKMIN_U) BIND(c, name="chnkmin_c_launch")
        implicit none
        REAL, VALUE :: U10
        REAL, VALUE :: ALPHA
        REAL, VALUE :: ALPHAMIN
        REAL, VALUE :: CHNKMIN_U
      END SUBROUTINE CHNKMIN_iso_c
    END INTERFACE
!$acc host_data use_device
    CALL CHNKMIN_iso_c(U10, ALPHA, ALPHAMIN, CHNKMIN_U)
!$acc end host_data
  END SUBROUTINE CHNKMIN_fc
END MODULE CHNKMIN_FC_MOD
