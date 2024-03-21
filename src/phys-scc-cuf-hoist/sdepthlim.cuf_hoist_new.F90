! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDEPTHLIM_CUF_HOIST_NEW_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE SDEPTHLIM_CUF_HOIST_NEW (KIJS, KIJL, EMAXDPT, FL1, DELTH, DFIM, EPSMIN, FR, NANG, NFRE, WETAIL,  &
  & ICHNK, NCHNK, IJ)
    ! ----------------------------------------------------------------------
    !     J. BIDLOT    ECMWF   NOVEMBER 2017
    
    !*    PURPOSE.
    !     --------
    !     LIMITS THE SPECTRAL VARIANCE SUCH THAT THE TOTAL VARIANCE
    !     DOES NOT EXCEED THE MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
    
    !**   INTERFACE.
    !     ----------
    !     *CALL* *SDEPTHLIM((KIJS, KIJL, EMAXDPT, FL1)
    !          *KIJS*   - LOCAL INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - LOCAL  INDEX OF LAST GRIDPOIN
    !          *EMAXDPT - MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
    !          *FL1*    - SPECTRUM.
    
    
    !     METHOD.
    !     -------
    
    !     EXTERNALS.
    !     ----------
    
    !     REFERENCE.
    !     ----------
    !     NONE
    
    ! ----------------------------------------------------------------------
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
    INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
    REAL(KIND=JWRB), INTENT(IN) :: EMAXDPT(KIJL, NCHNK)
    REAL(KIND=JWRB), INTENT(INOUT) :: FL1(KIJL, NANG_loki_param, NFRE_loki_param, NCHNK)
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: K
    INTEGER(KIND=JWIM) :: M
    REAL(KIND=JWRB) :: DELT25
    REAL(KIND=JWRB) :: EM
    REAL(KIND=JWRB) :: TEMP
    LOGICAL :: LLEPSMIN
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: DFIM(NFRE_loki_param)
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: EPSMIN
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: FR(NFRE_loki_param)
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE
    REAL(KIND=JWRB), VALUE, INTENT(IN) :: WETAIL
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
    INTEGER, VALUE, INTENT(IN) :: NCHNK
    
    ! ----------------------------------------------------------------------
    
    
    
    EM = EPSMIN
    DO M=1,NFRE
      K = 1
      TEMP = FL1(IJ, K, M, ICHNK)
      DO K=2,NANG
        TEMP = TEMP + FL1(IJ, K, M, ICHNK)
      END DO
      EM = EM + DFIM(M)*TEMP
    END DO
    ! ----------------------------------------------------------------------
    
    !*    3. ADD TAIL ENERGY.
    !        ----------------
    
    DELT25 = WETAIL*FR(NFRE)*DELTH
    EM = EM + DELT25*TEMP
    
    EM = MIN(EMAXDPT(IJ, ICHNK) / EM, 1.0_JWRB)
    
    DO M=1,NFRE
      DO K=1,NANG
        FL1(IJ, K, M, ICHNK) = MAX(FL1(IJ, K, M, ICHNK)*EM, EPSMIN)
      END DO
    END DO
    
    
    
  END SUBROUTINE SDEPTHLIM_CUF_HOIST_NEW
END MODULE SDEPTHLIM_CUF_HOIST_NEW_MOD
