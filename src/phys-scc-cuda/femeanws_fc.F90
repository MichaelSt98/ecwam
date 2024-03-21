MODULE FEMEANWS_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE FEMEANWS_fc (KIJS, KIJL, FL1, XLLWS, FM, EM, DELTH, DFIM, DFIMOFR, EPSMIN, FR, FRTAIL, NANG, NFRE, WETAIL, ICHNK,  &
  & NCHNK, IJ)
    USE iso_c_binding, ONLY: c_loc
    USE PARKIND_WAVE, ONLY: JWIM, JWRB
    
    
    
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
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    INTERFACE
      SUBROUTINE FEMEANWS_iso_c (KIJS, KIJL, FL1, XLLWS, FM, EM, DELTH, DFIM, DFIMOFR, EPSMIN, FR, FRTAIL, NANG, NFRE, WETAIL,  &
      & ICHNK, NCHNK, IJ) BIND(c, name="femeanws_c_launch")
        USE iso_c_binding, ONLY: c_int, c_ptr
        implicit none
        INTEGER(KIND=c_int), VALUE :: KIJS
        INTEGER(KIND=c_int), VALUE :: KIJL
        TYPE(c_ptr), VALUE :: FL1
        TYPE(c_ptr), VALUE :: XLLWS
        TYPE(c_ptr), VALUE :: FM
        TYPE(c_ptr), VALUE :: EM
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
      END SUBROUTINE FEMEANWS_iso_c
    END INTERFACE
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FL1(:, :, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: XLLWS(:, :, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: FM(:)
    REAL(KIND=JWRB), OPTIONAL, TARGET, INTENT(OUT) :: EM(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIM(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIMOFR(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FR(:)
!$acc host_data use_device( FL1, XLLWS, FM, EM, DFIM, DFIMOFR, FR )
    CALL FEMEANWS_iso_c(KIJS, KIJL, c_loc(FL1), c_loc(XLLWS), c_loc(FM), c_loc(EM), DELTH, c_loc(DFIM), c_loc(DFIMOFR), EPSMIN,  &
    & c_loc(FR), FRTAIL, NANG, NFRE, WETAIL, ICHNK, NCHNK, IJ)
!$acc end host_data
  END SUBROUTINE FEMEANWS_fc
END MODULE FEMEANWS_FC_MOD
