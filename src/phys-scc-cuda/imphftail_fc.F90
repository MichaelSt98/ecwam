MODULE IMPHFTAIL_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE IMPHFTAIL_fc (KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1, NANG, NFRE, ICHNK, NCHNK, IJ)
    USE iso_c_binding, ONLY: c_loc
    USE PARKIND_WAVE, ONLY: JWIM, JWRB
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    INTERFACE
      SUBROUTINE IMPHFTAIL_iso_c (KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1, NANG, NFRE, ICHNK, NCHNK, IJ) &
      &  BIND(c, name="imphftail_c_launch")
        USE iso_c_binding, ONLY: c_int, c_ptr
        implicit none
        INTEGER(KIND=c_int), VALUE :: KIJS
        INTEGER(KIND=c_int), VALUE :: KIJL
        TYPE(c_ptr), VALUE :: MIJ
        TYPE(c_ptr), VALUE :: FLM
        TYPE(c_ptr), VALUE :: WAVNUM
        TYPE(c_ptr), VALUE :: XK2CG
        TYPE(c_ptr), VALUE :: FL1
        INTEGER(KIND=c_int), VALUE :: NANG
        INTEGER(KIND=c_int), VALUE :: NFRE
        INTEGER(KIND=c_int), VALUE :: ICHNK
        INTEGER(KIND=c_int), VALUE :: NCHNK
        INTEGER(KIND=c_int), VALUE :: IJ
      END SUBROUTINE IMPHFTAIL_iso_c
    END INTERFACE
    INTEGER(KIND=JWIM), TARGET, INTENT(IN) :: MIJ(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FLM(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: WAVNUM(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: XK2CG(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: FL1(:, :, :, :)
!$acc host_data use_device( MIJ, FLM, WAVNUM, XK2CG, FL1 )
    CALL IMPHFTAIL_iso_c(KIJS, KIJL, c_loc(MIJ), c_loc(FLM), c_loc(WAVNUM), c_loc(XK2CG), c_loc(FL1), NANG, NFRE, ICHNK, NCHNK,  &
    & IJ)
!$acc end host_data
  END SUBROUTINE IMPHFTAIL_fc
END MODULE IMPHFTAIL_FC_MOD
