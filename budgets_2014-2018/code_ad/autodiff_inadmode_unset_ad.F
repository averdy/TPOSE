#include "AUTODIFF_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif

      SUBROUTINE ADAUTODIFF_INADMODE_UNSET( myThid )
C     *==========================================================*
C     | SUBROUTINE ADAUTODIFF_INADMODE_UNSET
C     *==========================================================*

      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "AUTODIFF_PARAMS.h"
#ifdef ALLOW_SEAICE
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
#endif
#ifdef ALLOW_CTRL
#include "ctrl.h"
#endif

C     == Routine arguments ==
C     myThid :: my Thread Id number
      INTEGER myThid
      CHARACTER*(MAX_LEN_MBUF) msgBuf

CAV(
      viscAhd=5.E2
      viscAhz=5.E2
      diffKhT=1.E2
      diffKhS=1.E2
CAV)

      inAdMode  = inAdFalse

      useKPP    = useKPPinFwdMode
      useGMRedi = useGMRediInFwdMode
      useSEAICE = useSEAICEinFwdMode
      useGGL90    = useGGL90inFwdMode
      useSALT_PLUME    = useSALT_PLUMEinFwdMode
#if (defined (ALLOW_CTRL) && defined (ECCO_CTRL_DEPRECATED))
      ctrlSmoothCorrel2D    = useSmoothCorrel2DinFwdMode
#endif
#ifdef ALLOW_SEAICE
      IF ( SEAICEuseFREEDRIFTswitchInAd ) THEN
        SEAICEuseFREEDRIFT = SEAICEuseFREEDRIFTinFwdMode
        SEAICEuseLSR       = .NOT.SEAICEuseFREEDRIFT
      ENDIF
      IF ( SEAICEuseDYNAMICSswitchInAd ) THEN
        SEAICEuseDYNAMICS  = SEAICEuseDYNAMICSinFwdMode
      ENDIF
      SEAICEadjMODE=0
#endif /* ALLOW_SEAICE */
      viscFacAdj=1. _d 0

      IF ( debugLevel.GE.debLevC ) THEN
       WRITE(msgBuf,'(A,L5)')
     &      'S/R ADAUTODIFF_INADMODE_UNSET: useKPP    =', useKPP
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       WRITE(msgBuf,'(A,L5)')
     &      'S/R ADAUTODIFF_INADMODE_UNSET: useGMRedi =', useGMRedi
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       WRITE(msgBuf,'(A,L5)')
     &      'S/R ADAUTODIFF_INADMODE_UNSET: useSEAICE =', useSEAICE
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       WRITE(msgBuf,'(A,L5)')
     &      'S/R ADAUTODIFF_INADMODE_UNSET: useGGL90    =', useGGL90
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       WRITE(msgBuf,'(A,L5)')
     &      'S/R ADAUTODIFF_INADMODE_UNSET: useSALT_PLUME    =',
     &      useSALT_PLUME
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )

       WRITE(msgBuf,'(A,1PE21.14)')
     &      'S/R ADAUTODIFF_INADMODE_UNSET: viscFacAdj =',viscFacAdj
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )

#if (defined (ALLOW_CTRL) && defined (ECCO_CTRL_DEPRECATED))
       WRITE(msgBuf,'(A,L5)')
     &      'S/R ADAUTODIFF_INADMODE_UNSET: ctrlSmoothCorrel2D    =',
     &       ctrlSmoothCorrel2D
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
#endif

#ifdef ALLOW_SEAICE
       IF ( SEAICEuseFREEDRIFTswitchInAd ) THEN
         WRITE(msgBuf,'(2A,L5)') 'S/R ADAUTODIFF_INADMODE_UNSET: ',
     &        'SEAICEuseFREEDRIFT =', SEAICEuseFREEDRIFT
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
         WRITE(msgBuf,'(2A,L5)') 'S/R ADAUTODIFF_INADMODE_UNSET: ',
     &        'SEAICEuseLSR       =', SEAICEuseLSR
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
       ENDIF
       IF ( SEAICEuseDYNAMICSswitchInAd ) THEN
         WRITE(msgBuf,'(2A,L5)') 'S/R ADAUTODIFF_INADMODE_UNSET: ',
     &        'SEAICEuseDYNAMICS =', SEAICEuseDYNAMICS
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
       ENDIF
       IF ( SEAICEapproxLevInAd.NE.0 ) THEN
         WRITE(msgBuf,'(2A,I2)') 'S/R ADAUTODIFF_INADMODE_UNSET: ',
     &        'SEAICEadjMODE =', SEAICEadjMODE
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
       ENDIF
#endif /* ALLOW_SEAICE */

      ENDIF

      RETURN
      END
