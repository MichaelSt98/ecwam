! Shared ordered CPU/GPU refraction plane updates. P2_OUT selects output storage.
!$loki loop-fusion
DO IJ = PROPAGS2_IJS, PROPAGS2_IJL
  P2_OUT(IJ,K,M) = (1.0_JWRB-SUMWN(IJ,K,M))*F1(IJ,K,M)
ENDDO
!$loki loop-unroll
DO IC=1,2
  IF (LLWLONN(K,M,IC)) THEN
    !$loki loop-fusion
    DO IJ = PROPAGS2_IJS, PROPAGS2_IJL
      P2_OUT(IJ,K,M) = P2_OUT(IJ,K,M) + WLONN(IJ,K,M,IC)*F1(KLON(IJ,IC),K,M)
    ENDDO
  ENDIF
ENDDO
#if defined(ECWAM_PROPAG_FACTORED_PAIRS) && !defined(_OPENACC)
IF (P2_FACTORED) THEN
#include "propags2_pairs_plane.h"
ELSE
#endif
!$loki loop-unroll
DO ICL=1,2
  !$loki loop-unroll
  DO IC=1,2
    IF (LLWLATN(K,M,IC,ICL)) THEN
      !$loki loop-fusion
      DO IJ = PROPAGS2_IJS, PROPAGS2_IJL
        P2_OUT(IJ,K,M) = P2_OUT(IJ,K,M) + WLATN(IJ,K,M,IC,ICL)*F1(KLAT(IJ,IC,ICL),K,M)
      ENDDO
    ENDIF
  ENDDO
  !$loki loop-unroll
  DO ICR=1,4
    IF (LLWCORN(K,M,ICR,ICL)) THEN
      !$loki loop-fusion
      DO IJ = PROPAGS2_IJS, PROPAGS2_IJL
        P2_OUT(IJ,K,M) = P2_OUT(IJ,K,M) + WCORN(IJ,K,M,ICR,ICL)*F1(KCOR(IJ,KCR(K,ICR),ICL),K,M)
      ENDDO
    ENDIF
  ENDDO
ENDDO
#if defined(ECWAM_PROPAG_FACTORED_PAIRS) && !defined(_OPENACC)
ENDIF
#endif
LSPECTRAL = LLWKPMN(K,M,-1) .AND. LLWMPMN(K,M,-1) .AND. &
            LLWKPMN(K,M,1) .AND. LLWMPMN(K,M,1)
IF (LSPECTRAL) THEN
  !$loki loop-fusion
!DIR$ VECTOR UNALIGNED
!DIR$ UNROLL(2)
  DO IJ = PROPAGS2_IJS, PROPAGS2_IJL
    ZUPDATE = P2_OUT(IJ,K,M)
    ZUPDATE = ZUPDATE + WKPMN(IJ,K,M,-1)*F1(IJ,KPM(K,-1),M)
    ZUPDATE = ZUPDATE + WMPMN(IJ,K,M,-1)*F1(IJ,K,MPM(M,-1))
    ZUPDATE = ZUPDATE + WKPMN(IJ,K,M,1)*F1(IJ,KPM(K,1),M)
    ZUPDATE = ZUPDATE + WMPMN(IJ,K,M,1)*F1(IJ,K,MPM(M,1))
    P2_OUT(IJ,K,M) = ZUPDATE
  ENDDO
ELSE
  !$loki loop-unroll
  DO IC=-1,1,2
    IF (LLWKPMN(K,M,IC)) THEN
      !$loki loop-fusion
      DO IJ = PROPAGS2_IJS, PROPAGS2_IJL
        P2_OUT(IJ,K,M) = P2_OUT(IJ,K,M) + WKPMN(IJ,K,M,IC)*F1(IJ,KPM(K,IC),M)
      ENDDO
    ENDIF
    IF (LLWMPMN(K,M,IC)) THEN
      !$loki loop-fusion
      DO IJ = PROPAGS2_IJS, PROPAGS2_IJL
        P2_OUT(IJ,K,M) = P2_OUT(IJ,K,M) + WMPMN(IJ,K,M,IC)*F1(IJ,K,MPM(M,IC))
      ENDDO
    ENDIF
  ENDDO
ENDIF
