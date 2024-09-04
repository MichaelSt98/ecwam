MODULE SEMEAN_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE SEMEAN_fc (FL1, KIJS, KIJL, EM, LLEPSMIN, DELTH, DFIM, EPSMIN, FR, NANG, NFRE, WETAIL, ICHNK, NCHNK, IJ, TEMP)
    USE iso_c_binding, ONLY: c_loc
    USE PARKIND_WAVE, ONLY: JWRB, JWIM
    
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    LOGICAL, VALUE, INTENT(IN) :: LLEPSMIN
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPSMIN
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WETAIL
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NCHNK
    INTERFACE
      SUBROUTINE SEMEAN_iso_c (FL1, KIJS, KIJL, EM, LLEPSMIN, DELTH, DFIM, EPSMIN, FR, NANG, NFRE, WETAIL, ICHNK, NCHNK, IJ, TEMP &
      & ) BIND(c, name="semean_c_launch")
        USE iso_c_binding, ONLY: c_int, c_ptr
        implicit none
        TYPE(c_ptr), VALUE :: FL1
        INTEGER(KIND=c_int), VALUE :: KIJS
        INTEGER(KIND=c_int), VALUE :: KIJL
        TYPE(c_ptr), VALUE :: EM
        LOGICAL, VALUE :: LLEPSMIN
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
        TYPE(c_ptr), VALUE :: TEMP
      END SUBROUTINE SEMEAN_iso_c
    END INTERFACE
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FL1(:, :, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: EM(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIM(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FR(:)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: TEMP(:, :)
!$acc host_data use_device( FL1, EM, DFIM, FR, TEMP )
    CALL SEMEAN_iso_c(c_loc(FL1), KIJS, KIJL, c_loc(EM), LLEPSMIN, DELTH, c_loc(DFIM), EPSMIN, c_loc(FR), NANG, NFRE, WETAIL,  &
    & ICHNK, NCHNK, IJ, c_loc(TEMP))
!$acc end host_data
  END SUBROUTINE SEMEAN_fc
END MODULE SEMEAN_FC_MOD
