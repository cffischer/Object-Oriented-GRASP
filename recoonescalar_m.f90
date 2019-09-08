!*******************************************************************
!                                                                  *
  MODULE recoonescalar_m
!                                                                  *
!   Third-level module for evaluating the reduced matrix elements
!   of a one particle operator for configurations in jj-coupling.  *
!                                                                  *
!                                                                  *
!                                                                  *
!   Structure by A. Senchuk                      September 2019    *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def,   ONLY:  NNNW
      USE CONS_C
      USE m_C
      USE orb_C,   ONLY: NW, NAK
      USE dumx_C,  ONLY: JLIS, JC1S, JC2S
      USE trk_C


    PRIVATE   ! By default  
!-------------------------------------------------------------
!   E x p l i c i t   P u b l i c  D e c l a r a t i o n s     
!-------------------------------------------------------------
    PUBLIC :: RECOONESCALAR

  CONTAINS

!*******************************************************************
!                                                                  *
      SUBROUTINE RECOONESCALAR(NS,JA1,JA2,JA3,JA4,KA,IAT)
!                                                                  *
!     -------------  SECTION REC    SUBPROGRAM 05  --------------  *
!                                                                  *
!     NO SUBROUTINE CALLED                                         *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: NS, JA1, JA2, JA3, JA4, KA
      INTEGER, INTENT(OUT) :: IAT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IA1, IA2, IJ, IJ1, IJ2, J, NPEELGG
!-----------------------------------------------
      IAT=1
      IF(NPEEL == 1 .AND. NS == -1)RETURN
      IF(NS == -1) THEN
         NPEELGG = NPEEL
      ELSE
         NPEELGG = NS
      END IF
      IJ1=JLIST(JA1)
      IJ2=JLIST(JA2)
      IF(JA1 == 1.AND.JA2 == 2) GO TO 1
      IF(KA /= 0)GO TO 5
!
!  CASES WHEN :          KA = 0
!                  OR    JA1 = JA2
!                  OR    JA1 = 1    JA2 = 2
!
    1 DO I=1,NPEELGG
        IJ=JLIST(I)
        IF(I < NPEELGG-1) THEN
         IF(JJC1(I) /= JJC2(I))IAT=0
        END IF
        IF(KA /= 0) THEN
          IF(I == JA1) CYCLE
          IF(I == JA2) CYCLE
        END IF
        DO J=1,3
          IF(JJQ1(J,IJ) /= JJQ2(J,IJ))IAT=0
        END DO
      END DO
      RETURN
!
!  OTHER CASES
!
    5 CONTINUE
      DO I=1,NPEELGG
        IJ=JLIST(I)
        IF(I < NPEELGG-1) THEN
          IA1=JA1-1
          IA2=JA2-1
          IF(JA1 == 1)IA1=JA1
          IF(I >= IA1.AND.I < IA2)GO TO 7
          IF(JJC1(I) /= JJC2(I))IAT=0
        END IF
    7   IF(I == JA1) CYCLE
        IF(I == JA2) CYCLE
        IF((KA == 2).AND.(I == JA3)) CYCLE
        IF((KA == 3).AND.(I == JA3)) CYCLE
        IF((KA == 3).AND.(I == JA4)) CYCLE
        DO J=1,3
          IF(JJQ1(J,IJ) /= JJQ2(J,IJ))IAT=0
        END DO
      END DO
      RETURN
      END SUBROUTINE RECOONESCALAR

    END MODULE recoonescalar_m
