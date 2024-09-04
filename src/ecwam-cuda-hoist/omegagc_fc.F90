MODULE OMEGAGC_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE OMEGAGC_fc (KIJS, KIJL, UST, NS, XKS, OMS, NWAV_GC, OMEGA_GC, SQRTGOSURFT, XKM_GC, XK_GC, XLOGKRATIOM1_GC, ICHNK,  &
  & NCHNK, IJ)
    USE iso_c_binding, ONLY: c_loc
    USE PARKIND_WAVE, ONLY: JWRB, JWIM
    
    
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTERFACE
      FUNCTION NS_GC (USTAR)
        USE parkind_wave, ONLY: jwrb
        INTEGER :: NS_GC
        REAL(KIND=JWRB), INTENT(IN) :: USTAR
      END FUNCTION NS_GC
    END INTERFACE
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NWAV_GC
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SQRTGOSURFT
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: XLOGKRATIOM1_GC
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NCHNK
    INTERFACE
      SUBROUTINE OMEGAGC_iso_c (KIJS, KIJL, UST, NS, XKS, OMS, NWAV_GC, OMEGA_GC, SQRTGOSURFT, XKM_GC, XK_GC, XLOGKRATIOM1_GC,  &
      & ICHNK, NCHNK, IJ) BIND(c, name="omegagc_c_launch")
        USE iso_c_binding, ONLY: c_int, c_ptr
        implicit none
        INTEGER(KIND=c_int), VALUE :: KIJS
        INTEGER(KIND=c_int), VALUE :: KIJL
        TYPE(c_ptr), VALUE :: UST
        TYPE(c_ptr), VALUE :: NS
        TYPE(c_ptr), VALUE :: XKS
        TYPE(c_ptr), VALUE :: OMS
        INTEGER(KIND=c_int), VALUE :: NWAV_GC
        TYPE(c_ptr), VALUE :: OMEGA_GC
        REAL, VALUE :: SQRTGOSURFT
        TYPE(c_ptr), VALUE :: XKM_GC
        TYPE(c_ptr), VALUE :: XK_GC
        REAL, VALUE :: XLOGKRATIOM1_GC
        INTEGER(KIND=c_int), VALUE :: ICHNK
        INTEGER(KIND=c_int), VALUE :: NCHNK
        INTEGER(KIND=c_int), VALUE :: IJ
      END SUBROUTINE OMEGAGC_iso_c
    END INTERFACE
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: UST(:, :)
    INTEGER(KIND=JWIM), TARGET, INTENT(OUT) :: NS(:)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: XKS(:)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: OMS(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: OMEGA_GC(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: XKM_GC(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: XK_GC(:)
!$acc host_data use_device( UST, NS, XKS, OMS, OMEGA_GC, XKM_GC, XK_GC )
    CALL OMEGAGC_iso_c(KIJS, KIJL, c_loc(UST), c_loc(NS), c_loc(XKS), c_loc(OMS), NWAV_GC, c_loc(OMEGA_GC), SQRTGOSURFT,  &
    & c_loc(XKM_GC), c_loc(XK_GC), XLOGKRATIOM1_GC, ICHNK, NCHNK, IJ)
!$acc end host_data
  END SUBROUTINE OMEGAGC_fc
END MODULE OMEGAGC_FC_MOD
