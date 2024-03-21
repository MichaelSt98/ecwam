INTERFACE
  SUBROUTINE STRESS_GC_FC (ANG_GC, USTAR, Z0, Z0MIN, HALP, RNFAC, BETAMAXOXKAPPA2, BMAXOKAP, C2OSQRTVG_GC, CM_GC, DELKCC_GC_NS,  &
  & DELKCC_OMXKM3_GC, EPSUS, LLNORMAGAM, NWAV_GC, OM3GMKM_GC, OMXKM3_GC, RN1_RN, SQRTGOSURFT, XKAPPA, XKMSQRTVGOC2_GC, XKM_GC,  &
  & XK_GC, XLOGKRATIOM1_GC, ZALP)
    USE PARKIND_WAVE, ONLY: JWRU, JWIM, JWRB
    
    USE YOWFRED, ONLY: OMEGA_GC
    USE YOWPCONS, ONLY: SURFT, G
    
    
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTERFACE
      FUNCTION NS_GC_FC (USTAR)
        USE parkind_wave, ONLY: jwrb
        INTEGER :: NS_GC
        REAL(KIND=JWRB), INTENT(IN) :: USTAR
      END FUNCTION NS_GC_FC
    END INTERFACE
    
    REAL(KIND=JWRB), INTENT(IN) :: ANG_GC    ! factor to account for angular spreading of the input.
    REAL(KIND=JWRB), INTENT(IN) :: USTAR    ! friction velocity
    REAL(KIND=JWRB), INTENT(IN) :: Z0    !  surface roughness
    REAL(KIND=JWRB), INTENT(IN) :: Z0MIN    ! minimum surface roughness
    REAL(KIND=JWRB), INTENT(IN) :: HALP    ! 1/2 Phillips parameter
    REAL(KIND=JWRB), INTENT(IN) :: RNFAC    ! wind dependent factor used in the growth renormalisation
    
    
    
    ! XLAMBDA = 1.0_JWRB + XLAMA * TANH(XLAMB * USTAR**NLAM)
    
    REAL(KIND=JWRB), INTENT(IN) :: BETAMAXOXKAPPA2
    REAL(KIND=JWRB), INTENT(IN) :: BMAXOKAP
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: C2OSQRTVG_GC(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: CM_GC(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DELKCC_GC_NS(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: DELKCC_OMXKM3_GC(:)
    REAL(KIND=JWRB), INTENT(IN) :: EPSUS
    LOGICAL, INTENT(IN) :: LLNORMAGAM
    INTEGER(KIND=JWIM), INTENT(IN) :: NWAV_GC
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: OM3GMKM_GC(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: OMXKM3_GC(:)
    REAL(KIND=JWRB), INTENT(IN) :: RN1_RN
    REAL(KIND=JWRB), INTENT(IN) :: SQRTGOSURFT
    REAL(KIND=JWRB), INTENT(IN) :: XKAPPA
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: XKMSQRTVGOC2_GC(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: XKM_GC(:)
    REAL(KIND=JWRB), TARGET, INTENT(IN) :: XK_GC(:)
    REAL(KIND=JWRB), INTENT(IN) :: XLOGKRATIOM1_GC
    REAL(KIND=JWRB), INTENT(IN) :: ZALP
!$acc routine seq
  END SUBROUTINE STRESS_GC_FC
END INTERFACE