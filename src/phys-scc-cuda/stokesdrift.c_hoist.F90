! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
ATTRIBUTES(DEVICE) SUBROUTINE STOKESDRIFT_FC (KIJS, KIJL, FL1, STOKFAC, WSWAVE, WDWAVE, CICOVER, USTOKES, VSTOKES, CITHRSH,  &
& COSTH, DELTH, DFIM_SIM, FR, G, LICERUN, LWAMRSETCI, NANG, NFRE_ODD, SINTH, ZPI, ICHNK, NCHNK, IJ)
  
  !
  !***  *STOKESDRIFT*   DETERMINES THE STOKES DRIFT
  !
  !     PETER JANSSEN MARCH 2009
  !
  !     PURPOSE.
  !     --------
  !
  !              DETERMINATION OF STOKES DRIFT VECTOR
  !
  !     INTERFACE.
  !     ----------
  !              *CALL*  *STOKESDRIFT(KIJS, KIJL, FL1, STOKFAC, WSWAVE,WDWAVE,CICOVER,USTOKES,VSTOKES)*
  !
  !                       INPUT:
  !                            *KIJS*   - FIRST GRIDPOINT
  !                            *KIJL*   - LAST GRIDPOINT
  !                            *FL1*    - 2-D SPECTRUM
  !                            *STOKFAC*- FACTOR TO COMPUTE THE STOKES DRIFT
  !                            Auxilliary fields to specify Stokes when model sea ice cover the blocking threshold
  !                            as 0.016*WSWAVE, aligned in the wind direction
  !                            *WSWAVE* - WIND SPEED IN M/S.
  !                            *WDWAVE* - WIND DIRECTION IN RADIANS.
  !                            *CICOVER*- SEA ICE COVER.
  !
  !                       OUTPUT:
  !                            *USTOKES*   - U-COMPONENT STOKES DRIFT
  !                            *VSTOKES*   - V-COMPONENT STOKES DRIFT
  !
  !     METHOD.
  !     -------
  !              DETERMINE U- AND V-COMPONENT OF STOKES DRIFT FOLLOWING
  !              K.E. KENYON, J.G.R., 74, 6991-6994
  !
  !     EXTERNALS.
  !     ----------
  !              NONE
  !
  !
  !-----------------------------------------------------------------------
  
  USE PARKIND_WAVE, ONLY: JWRU, JWIM, JWRB
  
  USE YOWFRED, ONLY: TH, DFIM, FRATIO
  USE YOWPARAM, ONLY: NFRE
  
  
  ! ----------------------------------------------------------------------
  IMPLICIT NONE
  
  INTEGER(KIND=JWIM), PARAMETER :: NANG_loki_param = 24
  INTEGER(KIND=JWIM), PARAMETER :: NFRE_loki_param = 36
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJL
  
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: FL1(:, :, :, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: STOKFAC(:, :, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: WSWAVE(:, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: WDWAVE(:, :)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: CICOVER(:, :)
  REAL(KIND=JWRB), TARGET, INTENT(OUT) :: USTOKES(:, :)
  REAL(KIND=JWRB), TARGET, INTENT(OUT) :: VSTOKES(:, :)
  
  
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
  INTEGER(KIND=JWIM) :: M
  INTEGER(KIND=JWIM) :: K
  
  REAL(KIND=JWRB), PARAMETER :: STMAX = 1.5_JWRB  ! maximum magnitude (this is for safety when coupled)
  REAL(KIND=JWRB) :: CONST
  REAL(KIND=JWRB) :: FAC
  REAL(KIND=JWRB) :: FAC1
  REAL(KIND=JWRB) :: FAC2
  REAL(KIND=JWRB) :: FAC3
  REAL(KIND=JWRB) :: STFAC
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: CITHRSH
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: COSTH(:)
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: DELTH
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: DFIM_SIM(:)
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: FR(:)
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: G
  LOGICAL, VALUE, INTENT(IN) :: LICERUN
  LOGICAL, VALUE, INTENT(IN) :: LWAMRSETCI
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NANG
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NFRE_ODD
  REAL(KIND=JWRB), TARGET, INTENT(IN) :: SINTH(:)
  REAL(KIND=JWRB), VALUE, INTENT(IN) :: ZPI
  INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: ICHNK
  INTEGER, VALUE, INTENT(IN) :: NCHNK
  
  ! ----------------------------------------------------------------------
  
  
  
  !***  1. DETERMINE STOKE DRIFT VECTOR.
  !     --------------------------------
  
  CONST = 2.0_JWRB*DELTH*ZPI**3 / G*FR(NFRE_ODD)**4
  
  !***  1.1 PERFORM INTEGRATION.
  !     ------------------------
  
  
  USTOKES(IJ, ICHNK) = 0.0_JWRB
  VSTOKES(IJ, ICHNK) = 0.0_JWRB
  
  DO M=1,NFRE_ODD
    STFAC = STOKFAC(IJ, M, ICHNK)*DFIM_SIM(M)
    DO K=1,NANG
      FAC3 = STFAC*FL1(IJ, K, M, ICHNK)
      USTOKES(IJ, ICHNK) = USTOKES(IJ, ICHNK) + FAC3*SINTH(K)
      VSTOKES(IJ, ICHNK) = VSTOKES(IJ, ICHNK) + FAC3*COSTH(K)
    END DO
  END DO
  
  !***  1.2 ADD CONTRIBUTION OF UNRESOLVED WAVES.
  !     -----------------------------------------
  
  DO K=1,NANG
    FAC1 = CONST*SINTH(K)
    FAC2 = CONST*COSTH(K)
    USTOKES(IJ, ICHNK) = USTOKES(IJ, ICHNK) + FAC1*FL1(IJ, K, NFRE_ODD, ICHNK)
    VSTOKES(IJ, ICHNK) = VSTOKES(IJ, ICHNK) + FAC2*FL1(IJ, K, NFRE_ODD, ICHNK)
  END DO
  
  
  !***  1.3 Sea Ice exception
  !     ---------------------
  IF (LICERUN .and. LWAMRSETCI) THEN
    IF (CICOVER(IJ, ICHNK) > CITHRSH) THEN
      USTOKES(IJ, ICHNK) = 0.016_JWRB*WSWAVE(IJ, ICHNK)*SIN(WDWAVE(IJ, ICHNK))*(1.0_JWRB - CICOVER(IJ, ICHNK))
      VSTOKES(IJ, ICHNK) = 0.016_JWRB*WSWAVE(IJ, ICHNK)*COS(WDWAVE(IJ, ICHNK))*(1.0_JWRB - CICOVER(IJ, ICHNK))
    END IF
  END IF
  
  !***  1.4 Protection
  !     --------------
  
  USTOKES(IJ, ICHNK) = MIN(MAX(USTOKES(IJ, ICHNK), -STMAX), STMAX)
  VSTOKES(IJ, ICHNK) = MIN(MAX(VSTOKES(IJ, ICHNK), -STMAX), STMAX)
  
  
  
END SUBROUTINE STOKESDRIFT_FC
