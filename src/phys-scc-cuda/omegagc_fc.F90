MODULE OMEGAGC_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE OMEGAGC_fc (UST, NS, XKS, OMS, NWAV_GC, OMEGA_GC, SQRTGOSURFT, XKM_GC, XK_GC, XLOGKRATIOM1_GC)
    USE iso_c_binding, ONLY: c_loc
    USE PARKIND_WAVE, ONLY: JWIM, JWRB
    
    
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTERFACE
      FUNCTION NS_GC (USTAR)
        USE parkind_wave, ONLY: jwrb
        INTEGER :: NS_GC
        REAL(KIND=JWRB), INTENT(IN) :: USTAR
      END FUNCTION NS_GC
    END INTERFACE
    REAL(KIND=JWRB), INTENT(IN) :: UST
    INTEGER(KIND=JWIM), INTENT(OUT) :: NS    ! index in array XK_GC corresponding to XKS and OMS
    REAL(KIND=JWRB), INTENT(OUT) :: XKS    ! cut-off wave number
    REAL(KIND=JWRB), INTENT(OUT) :: OMS    ! cut-off angular frequency
    
    INTEGER(KIND=JWIM), INTENT(IN) :: NWAV_GC
    REAL(KIND=JWRB), INTENT(IN) :: SQRTGOSURFT
    REAL(KIND=JWRB), INTENT(IN) :: XLOGKRATIOM1_GC
!$acc routine seq
    INTERFACE
      SUBROUTINE OMEGAGC_iso_c (UST, NS, XKS, OMS, NWAV_GC, OMEGA_GC, SQRTGOSURFT, XKM_GC, XK_GC, XLOGKRATIOM1_GC) &
      &  BIND(c, name="omegagc_c_launch")
        USE iso_c_binding, ONLY: c_int, c_ptr
        implicit none
        REAL, VALUE :: UST
        INTEGER(KIND=c_int) :: NS
        REAL :: XKS
        REAL :: OMS
        INTEGER(KIND=c_int), VALUE :: NWAV_GC
        TYPE(c_ptr), VALUE :: OMEGA_GC
        REAL, VALUE :: SQRTGOSURFT
        TYPE(c_ptr), VALUE :: XKM_GC
        TYPE(c_ptr), VALUE :: XK_GC
        REAL, VALUE :: XLOGKRATIOM1_GC
      END SUBROUTINE OMEGAGC_iso_c
    END INTERFACE
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: OMEGA_GC(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: XKM_GC(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: XK_GC(:)
!$acc host_data use_device( OMEGA_GC, XKM_GC, XK_GC )
    CALL OMEGAGC_iso_c(UST, NS, XKS, OMS, NWAV_GC, c_loc(OMEGA_GC), SQRTGOSURFT, c_loc(XKM_GC), c_loc(XK_GC), XLOGKRATIOM1_GC)
!$acc end host_data
  END SUBROUTINE OMEGAGC_fc
END MODULE OMEGAGC_FC_MOD
