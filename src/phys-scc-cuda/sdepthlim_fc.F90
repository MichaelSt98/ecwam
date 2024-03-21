MODULE SDEPTHLIM_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE SDEPTHLIM_fc (KIJS, KIJL, EMAXDPT, FL1, DELTH, DFIM, EPSMIN, FR, NANG, NFRE, WETAIL, ICHNK, NCHNK, IJ)
    USE iso_c_binding, ONLY: c_loc
    USE PARKIND_WAVE, ONLY: JWIM, JWRB
    
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPSMIN
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WETAIL
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    INTERFACE
      SUBROUTINE SDEPTHLIM_iso_c (KIJS, KIJL, EMAXDPT, FL1, DELTH, DFIM, EPSMIN, FR, NANG, NFRE, WETAIL, ICHNK, NCHNK, IJ) &
      &  BIND(c, name="sdepthlim_c_launch")
        USE iso_c_binding, ONLY: c_int, c_ptr
        implicit none
        INTEGER(KIND=c_int), VALUE :: KIJS
        INTEGER(KIND=c_int), VALUE :: KIJL
        TYPE(c_ptr), VALUE :: EMAXDPT
        TYPE(c_ptr), VALUE :: FL1
        REAL, VALUE :: DELTH
        TYPE(c_ptr), VALUE :: DFIM
        REAL, VALUE :: EPSMIN
        TYPE(c_ptr), VALUE :: FR
        INTEGER(KIND=c_int), VALUE :: NANG
        INTEGER(KIND=c_int), VALUE :: NFRE
        REAL, VALUE :: WETAIL
        INTEGER(KIND=c_int), VALUE :: ICHNK
        INTEGER(KIND=c_int), VALUE :: NCHNK
        INTEGER(KIND=c_int), VALUE :: IJ
      END SUBROUTINE SDEPTHLIM_iso_c
    END INTERFACE
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: EMAXDPT(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: FL1(:, :, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIM(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FR(:)
!$acc host_data use_device( EMAXDPT, FL1, DFIM, FR )
    CALL SDEPTHLIM_iso_c(KIJS, KIJL, c_loc(EMAXDPT), c_loc(FL1), DELTH, c_loc(DFIM), EPSMIN, c_loc(FR), NANG, NFRE, WETAIL,  &
    & ICHNK, NCHNK, IJ)
!$acc end host_data
  END SUBROUTINE SDEPTHLIM_fc
END MODULE SDEPTHLIM_FC_MOD
