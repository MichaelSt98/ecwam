MODULE NS_GC_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE NS_GC_fc (USTAR, NWAV_GC, SQRTGOSURFT, XKM_GC, XLOGKRATIOM1_GC)
    USE iso_c_binding, ONLY: c_loc
    USE PARKIND_WAVE, ONLY: JWIM, JWRB
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    REAL(KIND=JWRB), INTENT(IN) :: USTAR
    
    INTEGER(KIND=JWIM), INTENT(IN) :: NWAV_GC
    REAL(KIND=JWRB), INTENT(IN) :: SQRTGOSURFT
    REAL(KIND=JWRB), INTENT(IN) :: XLOGKRATIOM1_GC
!$acc routine seq
    INTERFACE
      SUBROUTINE NS_GC_iso_c (USTAR, NWAV_GC, SQRTGOSURFT, XKM_GC, XLOGKRATIOM1_GC) BIND(c, name="ns_gc_c_launch")
        USE iso_c_binding, ONLY: c_int, c_ptr
        implicit none
        REAL, VALUE :: USTAR
        INTEGER(KIND=c_int), VALUE :: NWAV_GC
        REAL, VALUE :: SQRTGOSURFT
        TYPE(c_ptr), VALUE :: XKM_GC
        REAL, VALUE :: XLOGKRATIOM1_GC
      END SUBROUTINE NS_GC_iso_c
    END INTERFACE
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: XKM_GC(:)
!$acc host_data use_device( XKM_GC )
    CALL NS_GC_iso_c(USTAR, NWAV_GC, SQRTGOSURFT, c_loc(XKM_GC), XLOGKRATIOM1_GC)
!$acc end host_data
  END SUBROUTINE NS_GC_fc
END MODULE NS_GC_FC_MOD
