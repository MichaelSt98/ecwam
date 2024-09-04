MODULE FEMEAN_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE FEMEAN_fc (KIJS, KIJL, F, EM, FM, DELTH, DFIM, DFIMOFR, EPSMIN, FR, FRTAIL, NANG, NFRE, WETAIL, ICHNK, NCHNK, IJ,  &
  & TEMP2)
    USE iso_c_binding, ONLY: c_loc
    USE PARKIND_WAVE, ONLY: JWRB, JWIM
    
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPSMIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: FRTAIL
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WETAIL
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NCHNK
    INTERFACE
      SUBROUTINE FEMEAN_iso_c (KIJS, KIJL, F, EM, FM, DELTH, DFIM, DFIMOFR, EPSMIN, FR, FRTAIL, NANG, NFRE, WETAIL, ICHNK,  &
      & NCHNK, IJ, TEMP2) BIND(c, name="femean_c_launch")
        USE iso_c_binding, ONLY: c_int, c_ptr
        implicit none
        INTEGER(KIND=c_int), VALUE :: KIJS
        INTEGER(KIND=c_int), VALUE :: KIJL
        TYPE(c_ptr), VALUE :: F
        TYPE(c_ptr), VALUE :: EM
        TYPE(c_ptr), VALUE :: FM
        REAL, VALUE :: DELTH
        TYPE(c_ptr), VALUE :: DFIM
        TYPE(c_ptr), VALUE :: DFIMOFR
        REAL, VALUE :: EPSMIN
        TYPE(c_ptr), VALUE :: FR
        REAL, VALUE :: FRTAIL
        INTEGER(KIND=c_int), VALUE :: NANG
        INTEGER(KIND=c_int), VALUE :: NFRE
        REAL, VALUE :: WETAIL
        INTEGER(KIND=c_int), VALUE :: ICHNK
        INTEGER(KIND=c_int), VALUE :: NCHNK
        INTEGER(KIND=c_int), VALUE :: IJ
        TYPE(c_ptr), VALUE :: TEMP2
      END SUBROUTINE FEMEAN_iso_c
    END INTERFACE
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: F(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: EM(:)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: FM(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIM(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIMOFR(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FR(:)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: TEMP2(:, :)
!$acc host_data use_device( F, EM, FM, DFIM, DFIMOFR, FR, TEMP2 )
    CALL FEMEAN_iso_c(KIJS, KIJL, c_loc(F), c_loc(EM), c_loc(FM), DELTH, c_loc(DFIM), c_loc(DFIMOFR), EPSMIN, c_loc(FR), FRTAIL,  &
    & NANG, NFRE, WETAIL, ICHNK, NCHNK, IJ, c_loc(TEMP2))
!$acc end host_data
  END SUBROUTINE FEMEAN_fc
END MODULE FEMEAN_FC_MOD
