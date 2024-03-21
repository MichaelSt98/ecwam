! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE FRCUTINDEX_CUF_HOIST_NEW_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE FRCUTINDEX_CUF_HOIST_NEW (KIJS, KIJL, FM, FMWS, UFRIC, CICOVER, MIJ, RHOWGDFTH, CITHRSH_TAIL,  &
  & EPSMIN, FLOGSPRDM1, FR, FRIC, G, NFRE, RHOWG_DFIM, TAILFACTOR, TAILFACTOR_PM, ZPIFR, ICHNK, NCHNK, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *FRCUTINDEX* - RETURNS THE LAST FREQUENCY INDEX OF
    !                    PROGNOSTIC PART OF SPECTRUM.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *FRCUTINDEX (KIJS, KIJL, FM, FMWS, CICOVER, MIJ, RHOWGDFTH)
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *FM*     - MEAN FREQUENCY
    !          *FMWS*   - MEAN FREQUENCY OF WINDSEA
    !          *UFRIC*  - FRICTION VELOCITY IN M/S
    !          *CICOVER*- CICOVER
    !          *MIJ*    - LAST FREQUENCY INDEX for imposing high frequency tail
    !          *RHOWGDFTH - WATER DENSITY * G * DF * DTHETA
    !                       FOR TRAPEZOIDAL INTEGRATION BETWEEN FR(1) and FR(MIJ)
    !                       !!!!!!!!  RHOWGDFTH=0 FOR FR > FR(MIJ)
    
    
    !     METHOD.
    !     -------
    
    !*    COMPUTES LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
    !*    FREQUENCIES LE 2.5*MAX(FMWS,FM).
    
    
    !!! be aware that if this is NOT used, for iphys=1, the cumulative dissipation has to be
    !!! re-activated (see module yowphys) !!!
    
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    USE YOWFRED, ONLY: DELTH, DFIM, FRATIO
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    INTEGER(KIND=JWIM), INTENT(OUT) :: MIJ(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: FM(KIJL)
    REAL(KIND=JWRB), INTENT(IN) :: FMWS(KIJL)
    REAL(KIND=JWRB), INTENT(IN) :: UFRIC(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: CICOVER(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(OUT) :: RHOWGDFTH(KIJL, NFRE_loki_param)
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: M
    
    REAL(KIND=JWRB) :: FPMH
    REAL(KIND=JWRB) :: FPPM
    REAL(KIND=JWRB) :: FM2
    REAL(KIND=JWRB) :: FPM
    REAL(KIND=JWRB) :: FPM4
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: CITHRSH_TAIL
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPSMIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: FLOGSPRDM1
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: FR(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: FRIC
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: RHOWG_DFIM(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: TAILFACTOR
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: TAILFACTOR_PM
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: ZPIFR(NFRE_loki_param)
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    
    ! ----------------------------------------------------------------------
    
    
    !*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
    !*    FREQUENCIES LE MAX(TAILFACTOR*MAX(FMNWS,FM),TAILFACTOR_PM*FPM),
    !*    WHERE FPM IS THE PIERSON-MOSKOWITZ FREQUENCY BASED ON FRICTION
    !*    VELOCITY. (FPM=G/(FRIC*ZPI*USTAR))
    !     ------------------------------------------------------------
    
    FPMH = TAILFACTOR / FR(1)
    FPPM = TAILFACTOR_PM*G / (FRIC*ZPIFR(1))
    
    
    IF (CICOVER(IJ, ICHNK) <= CITHRSH_TAIL) THEN
      FM2 = MAX(FMWS(IJ), FM(IJ))*FPMH
      FPM = FPPM / MAX(UFRIC(IJ, ICHNK), EPSMIN)
      FPM4 = MAX(FM2, FPM)
      MIJ(IJ, ICHNK) = NINT(LOG10(FPM4)*FLOGSPRDM1) + 1
      MIJ(IJ, ICHNK) = MIN(MAX(1, MIJ(IJ, ICHNK)), NFRE)
    ELSE
      MIJ(IJ, ICHNK) = NFRE
    END IF
    
    !     SET RHOWGDFTH
    DO M=1,MIJ(IJ, ICHNK)
      RHOWGDFTH(IJ, M) = RHOWG_DFIM(M)
    END DO
    IF (MIJ(IJ, ICHNK) /= NFRE) RHOWGDFTH(IJ, MIJ(IJ, ICHNK)) = 0.5_JWRB*RHOWGDFTH(IJ, MIJ(IJ, ICHNK))
    DO M=MIJ(IJ, ICHNK) + 1,NFRE
      RHOWGDFTH(IJ, M) = 0.0_JWRB
    END DO
    
    
    
  END SUBROUTINE FRCUTINDEX_CUF_HOIST_NEW
END MODULE FRCUTINDEX_CUF_HOIST_NEW_MOD
