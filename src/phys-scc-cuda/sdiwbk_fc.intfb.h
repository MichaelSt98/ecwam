INTERFACE
  SUBROUTINE SDIWBK_FC (KIJS, KIJL, FL1, FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN, LBIWBK, NANG, NFRE_RED, ICHNK, NCHNK, IJ)
    USE PARKIND_WAVE, ONLY: JWRU, JWIM, JWRB
    
    USE YOWPARAM, ONLY: NFRE
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: FL1(:, :, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: FLD(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: SL(:, :, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DEPTH(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: EMAXDPT(:, :)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: EMEAN(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: F1MEAN(:)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    LOGICAL, VALUE, INTENT(IN) :: LBIWBK
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE_RED
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
  END SUBROUTINE SDIWBK_FC
END INTERFACE