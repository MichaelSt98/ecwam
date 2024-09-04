MODULE SDISSIP_JAN_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE SDISSIP_JAN_fc (KIJS, KIJL, FL1, FLD, SL, WAVNUM, EMEAN, F1MEAN, XKMEAN, CDIS, CDISVIS, DELTA_SDIS, NANG, NFRE,  &
  & RNU, ZPI, ICHNK, NCHNK, IJ, TEMP1, SDS, X, XK2)
    USE iso_c_binding, ONLY: c_loc
    USE PARKIND_WAVE, ONLY: JWRB, JWIM
    
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: CDIS
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: CDISVIS
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTA_SDIS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: RNU
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZPI
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NCHNK
    INTERFACE
      SUBROUTINE SDISSIP_JAN_iso_c (KIJS, KIJL, FL1, FLD, SL, WAVNUM, EMEAN, F1MEAN, XKMEAN, CDIS, CDISVIS, DELTA_SDIS, NANG,  &
      & NFRE, RNU, ZPI, ICHNK, NCHNK, IJ, TEMP1, SDS, X, XK2) BIND(c, name="sdissip_jan_c_launch")
        USE iso_c_binding, ONLY: c_int, c_ptr
        implicit none
        INTEGER(KIND=c_int), VALUE :: KIJS
        INTEGER(KIND=c_int), VALUE :: KIJL
        TYPE(c_ptr), VALUE :: FL1
        TYPE(c_ptr), VALUE :: FLD
        TYPE(c_ptr), VALUE :: SL
        TYPE(c_ptr), VALUE :: WAVNUM
        TYPE(c_ptr), VALUE :: EMEAN
        TYPE(c_ptr), VALUE :: F1MEAN
        TYPE(c_ptr), VALUE :: XKMEAN
        REAL, VALUE :: CDIS
        REAL, VALUE :: CDISVIS
        REAL, VALUE :: DELTA_SDIS
        INTEGER(KIND=c_int), VALUE :: NANG
        INTEGER(KIND=c_int), VALUE :: NFRE
        REAL, VALUE :: RNU
        REAL, VALUE :: ZPI
        INTEGER(KIND=c_int), VALUE :: ICHNK
        INTEGER(KIND=c_int), VALUE :: NCHNK
        INTEGER(KIND=c_int), VALUE :: IJ
        TYPE(c_ptr), VALUE :: TEMP1
        TYPE(c_ptr), VALUE :: SDS
        TYPE(c_ptr), VALUE :: X
        TYPE(c_ptr), VALUE :: XK2
      END SUBROUTINE SDISSIP_JAN_iso_c
    END INTERFACE
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FL1(:, :, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: FLD(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: SL(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: WAVNUM(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: EMEAN(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: F1MEAN(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: XKMEAN(:)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: TEMP1(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: SDS(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: X(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: XK2(:, :)
!$acc host_data use_device( FL1, FLD, SL, WAVNUM, EMEAN, F1MEAN, XKMEAN, TEMP1, SDS, X, XK2 )
    CALL SDISSIP_JAN_iso_c(KIJS, KIJL, c_loc(FL1), c_loc(FLD), c_loc(SL), c_loc(WAVNUM), c_loc(EMEAN), c_loc(F1MEAN),  &
    & c_loc(XKMEAN), CDIS, CDISVIS, DELTA_SDIS, NANG, NFRE, RNU, ZPI, ICHNK, NCHNK, IJ, c_loc(TEMP1), c_loc(SDS), c_loc(X),  &
    & c_loc(XK2))
!$acc end host_data
  END SUBROUTINE SDISSIP_JAN_fc
END MODULE SDISSIP_JAN_FC_MOD
