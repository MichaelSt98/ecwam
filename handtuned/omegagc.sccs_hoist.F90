! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE OMEGAGC_LOKI_MOD
  CONTAINS
  SUBROUTINE OMEGAGC_LOKI (KIJS, KIJL, UST, NS, XKS, OMS, IJ)

    !***  DETERMINE THE CUT-OFF ANGULAR FREQUENCY FOR THE GRAV-CAPILLARY WAVES
    !     !!!! rounded to the closest index of XK_GC  !!!!!

    !     AUTHOR: PETER JANSSEN
    !     ------

    !     REFERENCES:
    !     ----------

    !     VIERS PAPER EQ.(29)

    !----------------------------------------------------------------------

    USE NS_GC_LOKI_MOD, ONLY: NS_GC_LOKI
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU

    USE YOWFRED, ONLY: OMEGA_GC, XK_GC

    USE YOMHOOK, ONLY: LHOOK, JPHOOK

    !----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: UST
    INTEGER(KIND=JWIM), INTENT(OUT), DIMENSION(KIJL) :: NS    ! index in array XK_GC corresponding to XKS and OMS
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: XKS    ! cut-off wave number
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: OMS    ! cut-off angular frequency


    INTEGER(KIND=JWIM), INTENT(IN) :: IJ
    REAL(KIND=JPHOOK) :: OMEGAGC_ZHOOK_HANDLE
!$acc routine seq
!$acc data present( UST, NS, XKS, OMS )
    ! ----------------------------------------------------------------------


    NS(IJ) = NS_GC_LOKI(UST(IJ))
    XKS(IJ) = XK_GC(NS(IJ))
    OMS(IJ) = OMEGA_GC(NS(IJ))



!$acc end data
  END SUBROUTINE OMEGAGC_LOKI
END MODULE OMEGAGC_LOKI_MOD
