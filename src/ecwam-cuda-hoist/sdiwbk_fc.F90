MODULE SDIWBK_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE SDIWBK_fc (KIJS, KIJL, FL1, FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN, LBIWBK, NANG, NFRE, NFRE_RED, ICHNK, NCHNK, IJ,  &
  & SDS)
    USE iso_c_binding, ONLY: c_loc
    USE PARKIND_WAVE, ONLY: JWRB, JWIM
    
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    LOGICAL, VALUE, INTENT(IN) :: LBIWBK
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE_RED
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NCHNK
    INTERFACE
      SUBROUTINE SDIWBK_iso_c (KIJS, KIJL, FL1, FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN, LBIWBK, NANG, NFRE, NFRE_RED, ICHNK,  &
      & NCHNK, IJ, SDS) BIND(c, name="sdiwbk_c_launch")
        USE iso_c_binding, ONLY: c_int, c_ptr
        implicit none
        INTEGER(KIND=c_int), VALUE :: KIJS
        INTEGER(KIND=c_int), VALUE :: KIJL
        TYPE(c_ptr), VALUE :: FL1
        TYPE(c_ptr), VALUE :: FLD
        TYPE(c_ptr), VALUE :: SL
        TYPE(c_ptr), VALUE :: DEPTH
        TYPE(c_ptr), VALUE :: EMAXDPT
        TYPE(c_ptr), VALUE :: EMEAN
        TYPE(c_ptr), VALUE :: F1MEAN
        LOGICAL, VALUE :: LBIWBK
        INTEGER(KIND=c_int), VALUE :: NANG
        INTEGER(KIND=c_int), VALUE :: NFRE
        INTEGER(KIND=c_int), VALUE :: NFRE_RED
        INTEGER(KIND=c_int), VALUE :: ICHNK
        INTEGER(KIND=c_int), VALUE :: NCHNK
        INTEGER(KIND=c_int), VALUE :: IJ
        TYPE(c_ptr), VALUE :: SDS
      END SUBROUTINE SDIWBK_iso_c
    END INTERFACE
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FL1(:, :, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: FLD(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: SL(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DEPTH(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: EMAXDPT(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: EMEAN(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: F1MEAN(:)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: SDS(:, :)
!$acc host_data use_device( FL1, FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN, SDS )
    CALL SDIWBK_iso_c(KIJS, KIJL, c_loc(FL1), c_loc(FLD), c_loc(SL), c_loc(DEPTH), c_loc(EMAXDPT), c_loc(EMEAN), c_loc(F1MEAN),  &
    & LBIWBK, NANG, NFRE, NFRE_RED, ICHNK, NCHNK, IJ, c_loc(SDS))
!$acc end host_data
  END SUBROUTINE SDIWBK_fc
END MODULE SDIWBK_FC_MOD
