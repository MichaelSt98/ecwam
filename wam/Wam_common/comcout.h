C*    *COMMON* *COUT*  OUTPUT POINTS INDICES AND FLAGS.
C
      CHARACTER*10 COUTT(MOUTT), COUTLST
      PARAMETER (JPPFLAG=21)
      LOGICAL FFLAG(JPPFLAG),FFLAG20,FFLAG21,FFLAG25,FFLAG26,FFLAG27,
     1        PFLAG(JPPFLAG),PFLAG20,PFLAG21,PFLAG25,PFLAG26,PFLAG27,
     2        GFLAG(JPPFLAG),GFLAG20,GFLAG21,GFLAG25,GFLAG26,GFLAG27,
     3        LFDB
      COMMON /COUT/ NGOUT, IGAR(MOUTP), IJAR(MOUTP),
     1              NOUTT, COUTT, COUTLST,
     2              FFLAG, FFLAG20, FFLAG21, FFLAG25, FFLAG26, FFLAG27,
     3              PFLAG, PFLAG20, PFLAG21, PFLAG25, PFLAG26, PFLAG27,
     4              GFLAG, GFLAG20, GFLAG21, GFLAG25, GFLAG26, GFLAG27,
     5              LFDB,IPFGTBL(JPPFLAG+1)
