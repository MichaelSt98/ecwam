! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE STRESSO_CUF_HOIST_NEW_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE STRESSO_CUF_HOIST_NEW (KIJS, KIJL, MIJ, RHOWGDFTH, FL1, SL, SPOS, CINV, WDWAVE, UFRIC, Z0M,  &
  & AIRD, RNFAC, COSWDIF, SINWDIF2, TAUW, TAUWDIR, PHIWA, LLPHIWA, COSTH, DELTH, EPS1, FR5, G, GAMNCONST, GM1, IPHYS,  &
  & JTOT_TAUHF, LLGCBZ0, LLNORMAGAM, NANG, NFRE, NWAV_GC, OMEGA_GC, RHOWG_DFIM, SINTH, SQRTGOSURFT, TAUWSHELTER, WTAUHF,  &
  & X0TAUHF, XKAPPA, XKM_GC, XK_GC, XLOGKRATIOM1_GC, ZALP, ZPI4GM1, ZPI4GM2, ZPIFR, ICHNK, NCHNK, IJ, TAUHF, PHIHF, UST)
    
    ! ----------------------------------------------------------------------
    
    !**** *STRESSO* - COMPUTATION OF WAVE STRESS.
    
    !     H. GUNTHER      GKSS/ECMWF NOVEMBER   1989 CODE MOVED FROM SINPUT.
    !     P.A.E.M. JANSSEN     KNMI  AUGUST     1990
    !     J. BIDLOT            ECMWF FEBRUARY   1996-97
    !     S. ABDALLA           ECMWF OCTOBER    2001 INTRODUCTION OF VARIABLE
    !                                                AIR DENSITY
    !     P.A.E.M. JANSSEN     ECMWF            2011  ADD FLUX CALULATIONS
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION
    
    !**   INTERFACE.
    !     ----------
    
    !        *CALL* *STRESSO (KIJS, KIJL, MIJ, RHOWGDFTH,
    !                         FL1, SL, SPOS,
    !    &                    CINV,
    !    &                    WDWAVE, UFRIC, Z0M, AIRD, RNFAC,
    !    &                    COSWDIF, SINWDIF2,
    !    &                    TAUW, TAUWDIR, PHIWA)*
    !         *KIJS*        - INDEX OF FIRST GRIDPOINT.
    !         *KIJL*        - INDEX OF LAST GRIDPOINT.
    !         *MIJ*         - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
    !         *RHOWGDFTH    - WATER DENSITY * G * DF * DTHETA
    !         *FL1*         - WAVE SPECTRUM.
    !         *SL*          - WIND INPUT SOURCE FUNCTION ARRAY (positive and negative contributions).
    !         *SPOS*        - POSITIVE WIND INPUT SOURCE FUNCTION ARRAY.
    !         *CINV*        - INVERSE PHASE VELOCITY.
    !         *WDWAVE*      - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
    !                         NOTATION (POINTING ANGLE OF WIND VECTOR,
    !                         CLOCKWISE FROM NORTH).
    !         *UFRIC*       - FRICTION VELOCITY IN M/S.
    !         *Z0M*         - ROUGHNESS LENGTH IN M.
    !         *AIRD*        - AIR DENSITY IN KG/M**3.
    !         *RNFAC*       - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !         *COSWDIF*     - COS(TH(K)-WDWAVE(IJ))
    !         *SINWDIF2*    - SIN(TH(K)-WDWAVE(IJ))**2
    !         *TAUW*        - KINEMATIC WAVE STRESS IN (M/S)**2
    !         *TAUWDIR*     - KINEMATIC WAVE STRESS DIRECTION
    !         *PHIWA*       - ENERGY FLUX FROM WIND INTO WAVES INTEGRATED
    !                         OVER THE FULL FREQUENCY RANGE.
    !         *LLPHIWA*     - TRUE IF PHIWA NEEDS TO BE COMPUTED
    
    !     METHOD.
    !     -------
    
    !       THE INPUT SOURCE FUNCTION IS INTEGRATED OVER FREQUENCY
    !       AND DIRECTIONS.
    !       BECAUSE ARRAY *SPOS* IS USED, ONLY THE INPUT SOURCE
    !       HAS TO BE STORED IN *SPOS* (CALL FIRST SINPUT, THEN
    !       STRESSO, AND THEN THE REST OF THE SOURCE FUNCTIONS)
    
    !     REFERENCE.
    !     ----------
    !       P. JANSSEN,
    
    ! ----------------------------------------------------------------------
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    USE YOWFRED, ONLY: FR, TH
    
    USE TAU_PHI_HF_CUF_HOIST_NEW_MOD, ONLY: TAU_PHI_HF_CUF_HOIST_NEW
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
    INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    INTEGER(KIND=JWIM), INTENT(IN) :: MIJ(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: RHOWGDFTH(KIJL, NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(IN) :: FL1(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: SL(KIJL, NANG_loki_param, NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(IN) :: SPOS(KIJL, NANG_loki_param, NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(IN) :: CINV(KIJL, NFRE_loki_param, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: WDWAVE(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: UFRIC(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: Z0M(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: AIRD(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(IN) :: RNFAC(KIJL)
    REAL(KIND=JWRB), INTENT(IN) :: COSWDIF(KIJL, NANG_loki_param)
    REAL(KIND=JWRB), INTENT(IN) :: SINWDIF2(KIJL, NANG_loki_param)
    REAL(KIND=JWRB), INTENT(OUT) :: TAUW(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(OUT) :: TAUWDIR(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(OUT) :: PHIWA(KIJL)
    LOGICAL, VALUE, INTENT(IN) :: LLPHIWA
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: M
    INTEGER(KIND=JWIM) :: K
    INTEGER(KIND=JWIM) :: I
    INTEGER(KIND=JWIM) :: J
    INTEGER(KIND=JWIM) :: II
    
    REAL(KIND=JWRB) :: TAUTOUS2
    REAL(KIND=JWRB) :: COSW
    REAL(KIND=JWRB) :: FCOSW2
    REAL(KIND=JWRB) :: XSTRESS
    REAL(KIND=JWRB) :: YSTRESS
    REAL(KIND=JWRB), INTENT(INOUT) :: TAUHF(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(INOUT) :: PHIHF(KIJL, NCHNK)
    REAL(KIND=JWRB) :: USDIRP
    REAL(KIND=JWRB), INTENT(INOUT) :: UST(KIJL, NCHNK)
    
    REAL(KIND=JWRB) :: CMRHOWGDFTH
    REAL(KIND=JWRB) :: TAUX
    REAL(KIND=JWRB) :: TAUY
    REAL(KIND=JWRB) :: TAUPX
    REAL(KIND=JWRB) :: TAUPY
    REAL(KIND=JWRB) :: SUMT
    REAL(KIND=JWRB) :: SUMX
    REAL(KIND=JWRB) :: SUMY
    
    LOGICAL :: LTAUWSHELTER
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: COSTH(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPS1
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: FR5(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: GAMNCONST
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: GM1
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IPHYS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: JTOT_TAUHF
    LOGICAL, VALUE, INTENT(IN) :: LLGCBZ0
    LOGICAL, VALUE, INTENT(IN) :: LLNORMAGAM
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NWAV_GC
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: OMEGA_GC(NWAV_GC)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: RHOWG_DFIM(NFRE_loki_param)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: SINTH(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: SQRTGOSURFT
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: TAUWSHELTER
    REAL(KIND=JWRB), INTENT(IN) :: WTAUHF(JTOT_TAUHF)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: X0TAUHF
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: XKAPPA
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: XKM_GC(NWAV_GC)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: XK_GC(NWAV_GC)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: XLOGKRATIOM1_GC
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZALP
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZPI4GM1
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZPI4GM2
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: ZPIFR(NFRE_loki_param)
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    
    ! ----------------------------------------------------------------------
    
    
    
    PHIWA(IJ) = 0.0_JWRB
    XSTRESS = 0.0_JWRB
    YSTRESS = 0.0_JWRB
    
    !*    CONTRIBUTION TO THE WAVE STRESS FROM THE NEGATIVE PART OF THE WIND INPUT
    !     ------------------------------------------------------------------------
    
    IF (LLPHIWA) THEN
      !     full energy flux due to negative Sinput (SL-SPOS)
      !     we assume that above NFRE, the contibutions can be neglected
      DO M=1,NFRE
        DO K=1,NANG
          PHIWA(IJ) = PHIWA(IJ) + (SL(IJ, K, M) - SPOS(IJ, K, M))*RHOWG_DFIM(M)
        END DO
      END DO
    END IF
    
    !*    CALCULATE LOW-FREQUENCY CONTRIBUTION TO STRESS AND ENERGY FLUX (positive sinput).
    !     ---------------------------------------------------------------------------------
    DO M=1,NFRE
      !     THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
      K = 1
      SUMX = SPOS(IJ, K, M)*SINTH(K)
      SUMY = SPOS(IJ, K, M)*COSTH(K)
      DO K=2,NANG
        SUMX = SUMX + SPOS(IJ, K, M)*SINTH(K)
        SUMY = SUMY + SPOS(IJ, K, M)*COSTH(K)
      END DO
      CMRHOWGDFTH = RHOWGDFTH(IJ, M)*CINV(IJ, M, ICHNK)
      XSTRESS = XSTRESS + CMRHOWGDFTH*SUMX
      YSTRESS = YSTRESS + CMRHOWGDFTH*SUMY
    END DO
    
    !     TAUW is the kinematic wave stress !
    XSTRESS = XSTRESS / MAX(AIRD(IJ, ICHNK), 1.0_JWRB)
    YSTRESS = YSTRESS / MAX(AIRD(IJ, ICHNK), 1.0_JWRB)
    
    IF (LLPHIWA) THEN
      DO M=1,NFRE
        !       THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
        K = 1
        SUMT = SPOS(IJ, K, M)
        DO K=2,NANG
          SUMT = SUMT + SPOS(IJ, K, M)
        END DO
        PHIWA(IJ) = PHIWA(IJ) + RHOWGDFTH(IJ, M)*SUMT
      END DO
    END IF
    
    !*    CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS and energy flux (positive sinput).
    !     ----------------------------------------------------------------------------------
    
    IF (IPHYS == 0 .or. TAUWSHELTER == 0.0_JWRB) THEN
      LTAUWSHELTER = .false.
      USDIRP = WDWAVE(IJ, ICHNK)
      UST(IJ, ICHNK) = UFRIC(IJ, ICHNK)
    ELSE
      LTAUWSHELTER = .true.
      TAUX = UFRIC(IJ, ICHNK)**2*SIN(WDWAVE(IJ, ICHNK))
      TAUY = UFRIC(IJ, ICHNK)**2*COS(WDWAVE(IJ, ICHNK))
      TAUPX = TAUX - TAUWSHELTER*XSTRESS
      TAUPY = TAUY - TAUWSHELTER*YSTRESS
      USDIRP = ATAN2(TAUPX, TAUPY)
      UST(IJ, ICHNK) = (TAUPX**2 + TAUPY**2)**0.25_JWRB
    END IF
    
    
    CALL TAU_PHI_HF_CUF_HOIST_NEW(KIJS, KIJL, MIJ(:, :), LTAUWSHELTER, UFRIC(:, :), Z0M(:, :), FL1(:, :, :, :), AIRD(:, :),  &
    & RNFAC(:), COSWDIF(:, :), SINWDIF2(:, :), UST(:, ICHNK), TAUHF(:, ICHNK), PHIHF(:, ICHNK), LLPHIWA, DELTH, FR5(:), G,  &
    & GAMNCONST, GM1, JTOT_TAUHF, LLGCBZ0, LLNORMAGAM, NANG, NWAV_GC, OMEGA_GC(:), SQRTGOSURFT, TAUWSHELTER, WTAUHF(:), X0TAUHF,  &
    & XKAPPA, XKM_GC(:), XK_GC(:), XLOGKRATIOM1_GC, ZALP, ZPI4GM1, ZPI4GM2, ZPIFR(:), ICHNK, NCHNK, IJ)
    
    
    XSTRESS = XSTRESS + TAUHF(IJ, ICHNK)*SIN(USDIRP)
    YSTRESS = YSTRESS + TAUHF(IJ, ICHNK)*COS(USDIRP)
    TAUW(IJ, ICHNK) = SQRT(XSTRESS**2 + YSTRESS**2)
    TAUW(IJ, ICHNK) = MAX(TAUW(IJ, ICHNK), 0.0_JWRB)
    TAUWDIR(IJ, ICHNK) = ATAN2(XSTRESS, YSTRESS)
    
    IF (.not.LLGCBZ0) THEN
      TAUTOUS2 = 1.0_JWRB / (1.0_JWRB + EPS1)
      TAUW(IJ, ICHNK) = MIN(TAUW(IJ, ICHNK), UFRIC(IJ, ICHNK)**2*TAUTOUS2)
    END IF
    
    IF (LLPHIWA) THEN
      PHIWA(IJ) = PHIWA(IJ) + PHIHF(IJ, ICHNK)
    END IF
    
    
    
  END SUBROUTINE STRESSO_CUF_HOIST_NEW
END MODULE STRESSO_CUF_HOIST_NEW_MOD
