! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE NS_GC_CUF_HOIST_NEW_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) FUNCTION NS_GC_CUF_HOIST_NEW (USTAR, NWAV_GC, SQRTGOSURFT, XKM_GC, XLOGKRATIOM1_GC)
    
    ! ----------------------------------------------------------------------
    
    !**** *NS_GC* - FUNCTION TO DETERMINE THE CUT-OFF ANGULAR FREQUENCY INDEX
    !               FOR THE GRAVITY-CAPILLARY MODEL
    !               !!!! rounded to the closest index of XK_GC  !!!!!
    
    !**   INTERFACE.
    !     ----------
    
    !       *FUNCTION* *NS_GC (USTAR)*
    
    !       *USTAR*  - FRICTION VELOCITY.
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRU, JWRB
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER :: NS_GC_CUF_HOIST_NEW
    REAL(KIND=JWRB), INTENT(IN) :: USTAR
    
    REAL(KIND=JWRB) :: Y, XKS
    INTEGER(KIND=JWIM), INTENT(IN) :: NWAV_GC
    REAL(KIND=JWRB), INTENT(IN) :: SQRTGOSURFT
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: XKM_GC(NWAV_GC)
    REAL(KIND=JWRB), INTENT(IN) :: XLOGKRATIOM1_GC
!$acc routine seq
    
    ! ----------------------------------------------------------------------
    
    
    !!!Y = 1.0_JWRB/(1.48_JWRB+2.05_JWRB*UST)
    !!!Y = (1.0_JWRB + UST**2)/(1.0_JWRB+10.0_JWRB*UST**2)
    
    XKS = SQRTGOSURFT / (1.48_JWRB + 2.05_JWRB*USTAR)
    
    NS_GC_CUF_HOIST_NEW = MIN(INT(LOG(XKS*XKM_GC(1))*XLOGKRATIOM1_GC) + 1, NWAV_GC - 1)
    
    
  END FUNCTION NS_GC_CUF_HOIST_NEW
END MODULE NS_GC_CUF_HOIST_NEW_MOD
