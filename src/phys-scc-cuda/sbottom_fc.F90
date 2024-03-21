MODULE SBOTTOM_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE SBOTTOM_fc (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, BATHYMAX, GM1, NANG, NFRE_RED, ICHNK, NCHNK, IJ)
    USE iso_c_binding, ONLY: c_loc
    USE PARKIND_WAVE, ONLY: JWIM, JWRB
    
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: BATHYMAX
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: GM1
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE_RED
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    INTERFACE
      SUBROUTINE SBOTTOM_iso_c (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, BATHYMAX, GM1, NANG, NFRE_RED, ICHNK, NCHNK, IJ) &
      &  BIND(c, name="sbottom_c_launch")
        USE iso_c_binding, ONLY: c_int, c_ptr
        implicit none
        INTEGER(KIND=c_int), VALUE :: KIJS
        INTEGER(KIND=c_int), VALUE :: KIJL
        TYPE(c_ptr), VALUE :: FL1
        TYPE(c_ptr), VALUE :: FLD
        TYPE(c_ptr), VALUE :: SL
        TYPE(c_ptr), VALUE :: WAVNUM
        TYPE(c_ptr), VALUE :: DEPTH
        REAL, VALUE :: BATHYMAX
        REAL, VALUE :: GM1
        INTEGER(KIND=c_int), VALUE :: NANG
        INTEGER(KIND=c_int), VALUE :: NFRE_RED
        INTEGER(KIND=c_int), VALUE :: ICHNK
        INTEGER(KIND=c_int), VALUE :: NCHNK
        INTEGER(KIND=c_int), VALUE :: IJ
      END SUBROUTINE SBOTTOM_iso_c
    END INTERFACE
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FL1(:, :, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: FLD(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: SL(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: WAVNUM(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DEPTH(:, :)
!$acc host_data use_device( FL1, FLD, SL, WAVNUM, DEPTH )
    CALL SBOTTOM_iso_c(KIJS, KIJL, c_loc(FL1), c_loc(FLD), c_loc(SL), c_loc(WAVNUM), c_loc(DEPTH), BATHYMAX, GM1, NANG,  &
    & NFRE_RED, ICHNK, NCHNK, IJ)
!$acc end host_data
  END SUBROUTINE SBOTTOM_fc
END MODULE SBOTTOM_FC_MOD
