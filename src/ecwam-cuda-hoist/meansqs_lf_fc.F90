MODULE MEANSQS_LF_FC_MOD
  USE iso_c_binding
  CONTAINS
  SUBROUTINE MEANSQS_LF_fc (NFRE_EFF, KIJS, KIJL, F, WAVNUM, XMSS, DFIM, NANG, NFRE, ICHNK, NCHNK, IJ, FD, TEMP1, TEMP2)
    USE iso_c_binding, ONLY: c_loc
    USE PARKIND_WAVE, ONLY: JWRB, JWIM
    
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE_EFF
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NCHNK
    INTERFACE
      SUBROUTINE MEANSQS_LF_iso_c (NFRE_EFF, KIJS, KIJL, F, WAVNUM, XMSS, DFIM, NANG, NFRE, ICHNK, NCHNK, IJ, FD, TEMP1, TEMP2) &
      &  BIND(c, name="meansqs_lf_c_launch")
        USE iso_c_binding, ONLY: c_int, c_ptr
        implicit none
        INTEGER(KIND=c_int), VALUE :: NFRE_EFF
        INTEGER(KIND=c_int), VALUE :: KIJS
        INTEGER(KIND=c_int), VALUE :: KIJL
        TYPE(c_ptr), VALUE :: F
        TYPE(c_ptr), VALUE :: WAVNUM
        TYPE(c_ptr), VALUE :: XMSS
        TYPE(c_ptr), VALUE :: DFIM
        INTEGER(KIND=c_int), VALUE :: NANG
        INTEGER(KIND=c_int), VALUE :: NFRE
        INTEGER(KIND=c_int), VALUE :: ICHNK
        INTEGER(KIND=c_int), VALUE :: NCHNK
        INTEGER(KIND=c_int), VALUE :: IJ
        TYPE(c_ptr), VALUE :: FD
        TYPE(c_ptr), VALUE :: TEMP1
        TYPE(c_ptr), VALUE :: TEMP2
      END SUBROUTINE MEANSQS_LF_iso_c
    END INTERFACE
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: F(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: WAVNUM(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(OUT) :: XMSS(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIM(:)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: FD(:)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: TEMP1(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: TEMP2(:, :)
!$acc host_data use_device( F, WAVNUM, XMSS, DFIM, FD, TEMP1, TEMP2 )
    CALL MEANSQS_LF_iso_c(NFRE_EFF, KIJS, KIJL, c_loc(F), c_loc(WAVNUM), c_loc(XMSS), c_loc(DFIM), NANG, NFRE, ICHNK, NCHNK, IJ,  &
    & c_loc(FD), c_loc(TEMP1), c_loc(TEMP2))
!$acc end host_data
  END SUBROUTINE MEANSQS_LF_fc
END MODULE MEANSQS_LF_FC_MOD
