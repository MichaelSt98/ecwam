! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
ATTRIBUTES(DEVICE) FUNCTION CHNKMIN_FC (U10, ALPHA, ALPHAMIN, CHNKMIN_U) RESULT(CHNKMIN)
  
  ! ----------------------------------------------------------------------
  
  !**** *CHNKMIN* - FUNCTION TO COMPUTE THE MINMUM CHARNOCK
  
  !*    PURPOSE.
  !     -------
  
  
  !**   INTERFACE.
  !     ----------
  
  !       *FUNCTION* *CHNKMIN (U10)*
  
  !     METHOD.
  !     -------
  
  !     CHNKMIN = ALPHAMIN + (ALPHA-ALPHAMIN)*0.5_JWRB*(1.0_JWRB-TANH(U10-A))
  
  !     EXTERNALS.
  !     ----------
  
  !       NONE.
  
  !     REFERENCE.
  !     ----------
  
  !       NONE.
  
  ! ----------------------------------------------------------------------
  
  USE PARKIND_WAVE, ONLY: JWRB, JWIM, JWRU
  
  USE YOMHOOK, ONLY: LHOOK, JPHOOK
  ! ----------------------------------------------------------------------
  
  IMPLICIT NONE
  
  REAL(KIND=JWRB) :: CHNKMIN
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: U10
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: ALPHA
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: ALPHAMIN
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: CHNKMIN_U
!$acc routine seq
  
  ! ----------------------------------------------------------------------
  
  
  CHNKMIN = ALPHAMIN + (ALPHA - ALPHAMIN)*0.5_JWRB*(1.0_JWRB - TANH(U10 - CHNKMIN_U))
  
  
END FUNCTION CHNKMIN_FC
