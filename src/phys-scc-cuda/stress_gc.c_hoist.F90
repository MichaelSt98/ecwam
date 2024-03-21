! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
ATTRIBUTES(DEVICE) FUNCTION STRESS_GC_FC (ANG_GC, USTAR, Z0, Z0MIN, HALP, RNFAC, BETAMAXOXKAPPA2, BMAXOKAP, C2OSQRTVG_GC, CM_GC,  &
& DELKCC_GC_NS, DELKCC_OMXKM3_GC, EPSUS, LLNORMAGAM, NWAV_GC, OM3GMKM_GC, OMXKM3_GC, RN1_RN, SQRTGOSURFT, XKAPPA,  &
& XKMSQRTVGOC2_GC, XKM_GC, XK_GC, XLOGKRATIOM1_GC, ZALP)
  
  !***  DETERMINE WAVE INDUCED STRESS FOR GRAV-CAP WAVES
  
  !     AUTHOR: PETER JANSSEN
  !     ------
  
  !     REFERENCES:
  !     ----------
  
  !     VIERS PAPER EQ.(29)
  !     FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
  
  !----------------------------------------------------------------------
  
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
  
  REAL(KIND=JWRB) :: STRESS_GC_FC
  REAL(KIND=JWRB), INTENT(IN) :: ANG_GC  ! factor to account for angular spreading of the input.
  REAL(KIND=JWRB), INTENT(IN) :: USTAR  ! friction velocity
  REAL(KIND=JWRB), INTENT(IN) :: Z0  !  surface roughness
  REAL(KIND=JWRB), INTENT(IN) :: Z0MIN  ! minimum surface roughness
  REAL(KIND=JWRB), INTENT(IN) :: HALP  ! 1/2 Phillips parameter
  REAL(KIND=JWRB), INTENT(IN) :: RNFAC  ! wind dependent factor used in the growth renormalisation
  
  
  INTEGER(KIND=JWIM) :: NS
  INTEGER(KIND=JWIM) :: I
  
  REAL(KIND=JWRB) :: XLAMBDA  ! Correction factor in the wave growth for gravity-capillary waves
  ! XLAMBDA = 1.0_JWRB + XLAMA * TANH(XLAMB * USTAR**NLAM)
  REAL(KIND=JWRB), PARAMETER :: XLAMA = 0.25_JWRB
  REAL(KIND=JWRB), PARAMETER :: XLAMB = 4.0_JWRB
  INTEGER(KIND=JWIM), PARAMETER :: NLAM = 4
  
  REAL(KIND=JWRB) :: TAUWCG_MIN
  REAL(KIND=JWRB) :: TAUWCG
  REAL(KIND=JWRB) :: ZABHRC
  REAL(KIND=JWRB) :: X, XLOG, ZLOG, ZLOG2X
  REAL(KIND=JWRB) :: CONST, ZN
  REAL(KIND=JWRB) :: GAMNORMA  ! RENORMALISATION FACTOR OF THE GROWTH RATE
  REAL(KIND=JWRB) :: GAM_W
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
  
  !     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
  
  ! ----------------------------------------------------------------------
  
  
  !*    1.0  DETERMINE GRAV_CAP SPECTRUM, TAUWHF.
  !          ------------------------------------
  
  !     FIND NS:
  NS = NS_GC_FC(USTAR, NWAV_GC, SQRTGOSURFT, XKM_GC, XLOGKRATIOM1_GC)
  
  TAUWCG_MIN = (USTAR*(Z0MIN / Z0))**2
  
  XLAMBDA = 1.0_JWRB + XLAMA*TANH(XLAMB*USTAR**NLAM)
  
  ZABHRC = ANG_GC*BETAMAXOXKAPPA2*HALP*C2OSQRTVG_GC(NS)
  IF (LLNORMAGAM) THEN
    CONST = RNFAC*BMAXOKAP*HALP*C2OSQRTVG_GC(NS) / MAX(USTAR, EPSUS)
  ELSE
    CONST = 0.0_JWRB
  END IF
  
  DO I=NS,NWAV_GC
    !       GROWTHRATE BY WIND WITHOUT the multiplicative factor representing the ratio of air density to water density (eps)
    !       and BETAMAXOXKAPPA2
    X = USTAR*CM_GC(I)
    XLOG = LOG(XK_GC(I)*Z0) + XKAPPA / (X + ZALP)
    ZLOG = XLOG - LOG(XLAMBDA)
    ZLOG = MIN(ZLOG, 0.0_JWRB)
    ZLOG2X = ZLOG*ZLOG*X
  END DO
  
  GAM_W = ZLOG2X*ZLOG2X*EXP(XLOG)*OM3GMKM_GC(NS)
  ZN = CONST*XKMSQRTVGOC2_GC(NS)*GAM_W
  GAMNORMA = (1.0_JWRB + RN1_RN*ZN) / (1.0_JWRB + ZN)
  TAUWCG = GAM_W*DELKCC_GC_NS(NS)*OMXKM3_GC(NS)*GAMNORMA
  DO I=NS + 1,NWAV_GC
    !       ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
    !       BB = HALP * C2OSQRTVG_GC(NS)*SQRT(VG_GC(I))/C_GC(I)**2
    !       Tauwcg : (rhow * g /rhoa) * integral of (1/c) * gammma * F(k)  k dk
    !       with omega=g*k and omega=k*c,  then
    !       Tauwcg : (rhow /rhoa) * integral of omega * gammma * F(k)  k dk
    !       but gamma is computed wihtout the rhoa/rhow factor so
    !       Tauwcg : integral of omega * gammma_wam * F(k)  k dk
    !       It should be done in vector form with actual directional spreading information
    !       It simplified here by using the ANG_GC factor.
    GAM_W = ZLOG2X*ZLOG2X*EXP(XLOG)*OM3GMKM_GC(I)
    ZN = CONST*XKMSQRTVGOC2_GC(I)*GAM_W
    GAMNORMA = (1.0_JWRB + RN1_RN*ZN) / (1.0_JWRB + ZN)
    TAUWCG = TAUWCG + GAM_W*DELKCC_OMXKM3_GC(I)*GAMNORMA
  END DO
  STRESS_GC_FC = MAX(ZABHRC*TAUWCG, TAUWCG_MIN)
  
  
END FUNCTION STRESS_GC_FC
