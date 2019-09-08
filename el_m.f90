!*******************************************************************
    MODULE el
!                                                                  *
!      Matrix elements for two particle operater                   *
!          < N1 N2 | Op | N'1 N'2>                                 *
!                                                                  *
!      SUBROUTINE CALLED: COULOM,GG1122,ITREXG,IXJTIK,PERKO2,      *
!                         RECO,RECO2,SIXJ,SPEAK,WW1                *
!      Structure by C. Froese Fischer                May     2019  *
!      Modification by A. Senchuk                  September 2019  *    
!                                                                  *
!*******************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, HALF, EPS
      USE m_C,             ONLY: JLIST, NPEEL
      USE orb_C,           ONLY: NAK
      USE trk_C
      USE el3n
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE reco_I
      USE reco2_I
      USE perko2_I
      USE itrig_I
      USE itrexg_I
      USE ixjtik_I
      USE snrc_I
      USE speak_I
      USE coulom_I
      USE ww1_I
      USE sixj_I
      USE cxk_I
      USE talk_I
      USE gg1122_I

   PRIVATE  ! By default

!-------------------------------------------------------------
!   E x p l i c i t   P u b l i c  D e c l a r a t i o n s     
!-------------------------------------------------------------
   PUBLIC :: EL1
   PUBLIC :: EL2
   PUBLIC :: EL3
   PUBLIC :: EL4
   PUBLIC :: EL5
   
   CONTAINS 

!*******************************************************************
!                                                                  *
      SUBROUTINE EL1(JJA,JJB,JA,JB,IIRE,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 03  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                  N'2 = N2        *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   Restructure by C. Froese Fischer                May      2019  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,IIRE,ICOLBREI 
!      DIMENSION CONE(7,20),S(12),IS(4),KAPS(4),KS(4)
!      DIMENSION PMGG(30),RAGG(30),J(2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: II,IA,IB,IAT,IP1,IP2,IP3,IG1,IG2,IG3,IKK,I1,I2,I3,I4,& 
                 IFAZ,J12,IBRD,IBRE,KRA,KRA1,L1,L2,MU,N,NU,ND1,ND2,   &
                 NE1,NE2,NUP1
      INTEGER, DIMENSION(2) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(DOUBLE)          :: QM1,QM2,QM3,QM4,AA,AB,A1,BB,SI,RECC,RAG
      REAL(DOUBLE), DIMENSION(12)   :: S
      REAL(DOUBLE), DIMENSION(30)   :: PMGG,RAGG
      REAL(DOUBLE), DIMENSION(7,20) :: CONE 
!-----------------------------------------------
      IF(JA /= JB)GO TO 9
!
      IF(JJA /= JJB) RETURN
!
!     THE CASE 1111   + + - -
!
      IF(IIRE /= 0) THEN
        CALL RECO(JA,JA,JA,JA,0,IAT)
        IF(IAT /= 0)RETURN
      END IF
      CALL PERKO2(JA,JA,JA,JA,1)
      QM1=HALF
      QM2=HALF
      QM3=-HALF
      QM4=-HALF
      IA=JLIST(JA)
      J(1)=ID1(3)
      IP2=ITREXG(J(1),J(1),J(1),J(1),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      L1=(J(1)+1)/2
      IP1=IP2
      IG1=IG2
      IF (ICOLBREI == 2) THEN
        IS(1)=IA
        IS(2)=IA
        IS(3)=IA
        IS(4)=IA
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0)RETURN 
      END IF
      DO I2=IP1,IG1,2
        KRA=(I2-1)/2
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L1,L1,L1,L1,ID1(5),ID1(5),ID1(5),ID1(5),KRA,A1)
          IF(DABS(A1) < EPS)CYCLE
          A1=-A1*HALF
        END IF
        AB=ZERO
        DO I3=IP2,IG2,2
          J12=(I3-1)/2
          IF(IXJTIK(J(1),J(1),KRA*2,J(1),J(1),J12*2) == 0)CYCLE
          CALL WW1(IK1,BK1,ID1,BD1,J12,QM1,QM2,QM3,QM4,AA)
          IF(DABS(AA) < EPS)CYCLE
          CALL SIXJ(J(1),J(1),KRA*2,J(1),J(1),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=IK1(3)+J12+KRA
          IF((IFAZ/2)*2 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
!
!     RECOUPLING COEFFICIENTS
!
        IF (ICOLBREI == 1) THEN
          BB=AB*A1
          BB=BB/DSQRT(DBLE(IK1(6)+1))
          IF(DABS(BB) > EPS)CALL SPEAK(JJA,JJB,IA,IA,IA,IA,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          N=(KRA-ND1)/2+1
          IF(((KRA-ND1)/2)*2 == (KRA-ND1)) THEN
            CALL CXK(S,IS,KAPS,KRA,KRA,3,1)
            IF(DABS(S(1)) > EPS) THEN
              BB =-HALF*S(1)*AB/DSQRT(DBLE(IK1(6)+1))
              IF(DABS(BB) > EPS)CALL TALK(JJA,JJB,KRA,IA,IA,IA,IA,4,BB)
            END IF
          END IF
        END IF
      END DO
      RETURN
!  ............................................................
    9 IF(NPEEL <= 1)RETURN
      IF(IIRE /= 0) THEN
        CALL RECO(JA,JB,JB,JB,1,IAT)
        IF(IAT == 0)RETURN
      END IF
      IA=JLIST(JA)
      IB=JLIST(JB)
      QM1=HALF
      QM2=-HALF
      QM3=HALF
      QM4=-HALF
      CALL PERKO2(JA,JB,JA,JA,2)
      J(1)=ID1(3)
      J(2)=ID2(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      IP1=ITREXG(J(1),J(1),J(2),J(2),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      IP3=IP1
      IG3=IG1
      DO I4=IP1,IG1,2
        KRA=(I4-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        RAGG(KRA1)=ZERO
        PMGG(KRA1)=ZERO
        CALL RECO2(JA,JB,KRA*2,0,IAT,RECC)
        IF(IAT == 0) CYCLE
        CALL GG1122(KRA,KRA,QM1,QM2,QM3,QM4,RAG)
        IF(DABS(RAG) < EPS) CYCLE
        RAGG(KRA1)=RAG
        CALL RECO2(JA,JB,KRA*2,1,IAT,RECC)
        PMGG(KRA1)=RECC
      END DO
! * * *                      * * *                      * * *
!     CASES 1212   + + - -        TRANSFORM TO  1122   + - + -
!           2121                                1122
!
      IF (ICOLBREI == 2) THEN
        IS(1)=IA
        IS(2)=IB
        IS(3)=IA
        IS(4)=IB
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        DO II=1,20
           CONE(1,II)=ZERO
           CONE(2,II)=ZERO
           CONE(3,II)=ZERO
           CONE(4,II)=ZERO
           CONE(5,II)=ZERO
           CONE(6,II)=ZERO
           CONE(7,II)=ZERO
        END DO
        IF(IBRD == 0 .AND. IBRE == 0)RETURN 
      END IF
      DO I1=IP1,IG1,2
        KRA=(I1-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L1,L2,L1,L2,ID1(5),ID2(5),ID1(5),ID2(5),KRA,AA)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA/DSQRT(DBLE(I1))
          IF(DABS(AA) > EPS) CALL SPEAK(JJA,JJB,IA,IB,IA,IB,KRA,AA)
        ELSE IF (ICOLBREI == 2) THEN
          N=(KRA-ND1)/2+1
          IF(((KRA-ND1)/2)*2 == (KRA-ND1)) THEN
            CALL CXK(S,IS,KAPS,KRA,KRA,3,1)
            IF(DABS(S(1)) > EPS) THEN
              BB=S(1)*PMGG(KRA1)*RAGG(KRA1)/DSQRT(DBLE(I1))
              IF(DABS(BB) > EPS)CALL TALK(JJA,JJB,KRA,IA,IA,IB,IB,4,BB)
            END IF
          END IF
        END IF
      END DO
! * * *                      * * *                      * * *
!     CASES 1221   + + - -        TRANSFORM TO  1122   + - + -
!           2112                                1122
!
      IP2=ITREXG(J(1),J(2),J(1),J(2),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
        IF(KRA > 30)GO TO 10
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L1,L2,L2,L1,ID1(5),ID2(5),ID2(5),ID1(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
        END IF
        AB=ZERO
        DO I3=IP3,IG3,2
          J12=(I3-1)/2
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(2),KRA*2,J(2),J(1),J12*2) == 0)CYCLE
          CALL SIXJ(J(1),J(2),KRA*2,J(2),J(1),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          AB=AB+AA
        END DO
        IF (ICOLBREI == 1) THEN
          BB=A1*AB
          IF(DABS(BB) > EPS)CALL SPEAK(JJA,JJB,IA,IB,IB,IA,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA 
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+1) /= 0)) THEN
              IF(NU > 0) THEN
                N=(NU-NE1)/2+1
                CALL CXK(S,IS,KAPS,NU,KRA,4,2)
                DO MU = 1,3
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU-1) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-NE1)/2+1
                IF(N <= NE2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,4,2)
                  DO MU = 1,3
                    CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-NE1)/2+1
                IF(N < NE2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,4,2)
                  DO MU = 1,7
                    CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,NE2
           NU=NE1+2*(N-1)
           CALL TALK(JJA,JJB,NU,IB,IA,IB,IA,5,CONE(1,N))
           CALL TALK(JJA,JJB,NU,IA,IB,IB,IA,5,CONE(2,N))
           CALL TALK(JJA,JJB,NU,IA,IB,IA,IB,5,CONE(3,N))
           IF(N == NE2) CYCLE
           NUP1=NU+1
           CALL TALK(JJA,JJB,NUP1,IA,IB,IA,IB,6,CONE(4,N))
           CALL TALK(JJA,JJB,NUP1,IB,IA,IB,IA,6,CONE(5,N))
           CALL TALK(JJA,JJB,NUP1,IA,IB,IB,IA,6,CONE(6,N))
           CALL TALK(JJA,JJB,NUP1,IB,IA,IA,IB,6,CONE(7,N))
        END DO
      END IF
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL1  PMGG RAGG')
      STOP
      END SUBROUTINE EL1
!*******************************************************************
!                                                                  *
      SUBROUTINE EL2(JJA,JJB,JA,JB,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 04  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 2        *
!                                              N'2 = N2 + 2        *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,GG1122,ITREXG,IXJTIK,PERKO2,       *
!                        RECO,RECO2,SIXJ,SPEAK                     *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   Restructure by C. Froese Fischer                May      2019  *
!                                                                  *
!*******************************************************************

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,ICOLBREI
!      DIMENSION J(2)
!      DIMENSION COND(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IAT,IA,IB,IBRD,IBRE,IP1,IP2,IG1,IG2,IKK,II,I2,I3, &
                 IFAZ,IFAZ1,IFAZFRCS,J12,JAA,JBB, &
                 KRA,L1,L2,N,ND1,ND2,NE1,NE2,NUP1,NU,MU
      INTEGER, DIMENSION(2) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(DOUBLE)          :: AA,A1,AB,BB,QM1,QM2,QM3,QM4,RECC,SI
      REAL(DOUBLE), DIMENSION(12)    :: S
      REAL(DOUBLE), DIMENSION(12,20) :: COND
!-----------------------------------------------
      IF(NPEEL <= 1)RETURN
      IF(JA > JB) THEN
        JAA=JB
        JBB=JA
      ELSE
        JAA=JA
        JBB=JB
      END IF
      CALL RECO(JAA,JBB,JBB,JBB,1,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      QM1=HALF
      QM2=HALF
      QM3=-HALF
      QM4=-HALF
      CALL PERKO2(JA,JB,JA,JA,2)
      J(1)=ID1(3)
      J(2)=ID2(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      IP1=ITREXG(J(1),J(1),J(2),J(2),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
! * * *                      * * *                      * * *
!     THE CASE 1122   + + - -
!
      IP2=ITREXG(J(1),J(2),J(1),J(2),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      IF (ICOLBREI == 2) THEN
        IS(1)=IA
        IS(2)=IA
        IS(3)=IB
        IS(4)=IB
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0)RETURN
        DO II=1,20
          COND(1,II) =ZERO
          COND(2,II) =ZERO
          COND(3,II) =ZERO
          COND(4,II) =ZERO
          COND(5,II) =ZERO
          COND(6,II) =ZERO
          COND(7,II) =ZERO
          COND(8,II) =ZERO
          COND(9,II) =ZERO
          COND(10,II)=ZERO
          COND(11,II)=ZERO
          COND(12,II)=ZERO
        END DO
      END IF
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L1,L1,L2,L2,ID1(5),ID1(5),ID2(5),ID2(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
          A1=-HALF*A1
        END IF
        AB=ZERO
          DO I3=IP1,IG1,2
            J12=(I3-1)/2
            CALL RECO2(JAA,JBB,J12*2,0,IAT,RECC)
            IF(IAT /= 0) THEN
              IF(IXJTIK(J(1),J(2),KRA*2,J(2),J(1),J12*2) /= 0) THEN
                CALL GG1122(J12,J12,QM1,QM2,QM3,QM4,AA)
                IF(DABS(AA) > EPS) THEN
                  CALL RECO2(JAA,JBB,J12*2,1,IAT,RECC)
                  AA=AA*RECC
                  CALL SIXJ(J(1),J(2),KRA*2,J(2),J(1),J12*2,0,SI)
                  AA=AA*SI*DSQRT(DBLE(I3))
                  IFAZ=IK1(3)+IK2(3)+KRA*2+J12*2
                  IF((IFAZ/4)*4 /= IFAZ)AA=-AA
                  AB=AB+AA
                END IF
              END IF
            END IF
          END DO
!
!       TRANSFORM FANO & RACAH PHASE CONVENTION
!       TO CONDON & SHORTLEY PHASE CONVENTION
!
        IFAZFRCS=1
        IFAZ1=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
        IF((IFAZ1/4)*4 /= IFAZ1)IFAZFRCS=-IFAZFRCS
!
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(DABS(BB) > EPS)CALL SPEAK(JJA,JJB,IA,IA,IB,IB,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND.   &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU-1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(N <= ND2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+3) /= 0) .AND.   &
               (ITRIG(KS(2),KS(4),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-ND1)/2+1
                IF(N < ND2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                  DO MU = 1,12
                    COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI .EQ. 2) THEN
        DO N = 1,ND2
          NU=ND1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
          IF(N == ND2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
        END DO
      END IF
      RETURN
      END SUBROUTINE EL2
!*******************************************************************
!                                                                  *
      SUBROUTINE EL3(JJA,JJB,JA,JB,JC,JD,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 05  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!                                                                  *
!     SUBROUTINE CALLED: EL31,EL32,EL33                            *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   Restructured by C. Froese Fischer               May      2019  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,ICOLBREI
!-----------------------------------------------
      IF(NPEEL <= 1)RETURN
      IF(JB == JD) THEN
        IF(JA == JB.OR.JC == JB) THEN
          IF(JA == JC)GO TO 10
          IF(JC /= JB) THEN
            CALL EL32(JJA,JJB,JC,JA,JA,JB,JC,JD,ICOLBREI)
          ELSE
            CALL EL31(JJA,JJB,JC,JA,JA,JB,JC,JD,ICOLBREI)
          END IF
        ELSE
          CALL EL33(JJA,JJB,JC,JA,JB,1,JA,JB,JC,JD,ICOLBREI)
        END IF
        RETURN
      ELSE IF(JA == JC) THEN
        IF(JB == JA.OR.JD == JA) THEN
          IF(JB == JD)GO TO 10
          IF(JD /= JA) THEN
            CALL EL32(JJA,JJB,JD,JB,JA,JB,JC,JD,ICOLBREI)
          ELSE
            CALL EL31(JJA,JJB,JD,JB,JA,JB,JC,JD,ICOLBREI)
          END IF
        ELSE
          CALL EL33(JJA,JJB,JD,JB,JA,1,JA,JB,JC,JD,ICOLBREI)
        END IF
        RETURN
      ELSE IF(JA == JD) THEN
        IF(JB == JA.OR.JC == JA) THEN
          IF(JB == JC)GO TO 10
          IF(JC /= JD) THEN
             CALL EL32(JJA,JJB,JC,JB,JA,JB,JC,JD,ICOLBREI)
          ELSE
             CALL EL31(JJA,JJB,JC,JB,JA,JB,JC,JD,ICOLBREI)
          END IF
        ELSE
          CALL EL33(JJA,JJB,JC,JB,JA,2,JA,JB,JD,JC,ICOLBREI)
        END IF
        RETURN
      ELSE IF(JB == JC) THEN
        IF(JA == JB.OR.JD == JB) THEN
          IF(JA == JD)GO TO 10
          IF(JD /= JB) THEN
            CALL EL32(JJA,JJB,JD,JA,JA,JB,JC,JD,ICOLBREI)
          ELSE
            CALL EL31(JJA,JJB,JD,JA,JA,JB,JC,JD,ICOLBREI)
          END IF
        ELSE
          CALL EL33(JJA,JJB,JD,JA,JB,2,JA,JB,JD,JC,ICOLBREI)
        END IF
        RETURN
      END IF
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL3  PMGG RAGG')
      STOP
   END SUBROUTINE EL3
!*******************************************************************
!                                                                  *
      SUBROUTINE EL4(JJA,JJB,JA,JB,JC,JD,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 09  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :      N'1 = N1 +- 1        *
!                                             N'2 = N2 +- 1        *
!                                             N'3 = N3 -+ 2        *
!                                                                  *
!     SUBROUTINE CALLED: EL41                                      *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE m_C,            ONLY: NPEEL
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE el41_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,ICOLBREI
!-----------------------------------------------
      IF(NPEEL <= 2)RETURN
      IF(JA == JB) THEN
        CALL EL41(JJA,JJB,JC,JD,JA,1,JA,JB,JC,JD,ICOLBREI)
      ELSE IF(JC == JD) THEN
        CALL EL41(JJA,JJB,JA,JB,JC,2,JA,JB,JC,JD,ICOLBREI)
      ELSE
        WRITE(99,100)
        STOP
      END IF
      RETURN
  100 FORMAT(5X,'ERRO IN EL4 ')
      END SUBROUTINE EL4
!*******************************************************************
!                                                                  *
      SUBROUTINE EL5(JJA,JJB,JA,JB,JC,JD,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 11  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :    N'1 = N1 (+-) 1        *
!                                           N'2 = N2 (+-) 1        *
!                                           N'3 = N3 (+-) 1        *
!                                           N'4 = N4 (+-) 1        *
!                                                                  *
!      SUBROUTINE CALLED: EL51,EL52,EL53                           *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE m_C,          ONLY: NPEEL
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE el51_I
      USE el52_I
      USE el53_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,ICOLBREI
!-----------------------------------------------
      IF(NPEEL <= 3)RETURN
      IF(JB < JC) THEN
        CALL EL51(JJA,JJB,JA,JB,JC,JD,1,ICOLBREI)
      ELSE IF(JA > JD.AND.JB > JD) THEN
        CALL EL51(JJA,JJB,JC,JD,JA,JB,2,ICOLBREI)
      ELSE IF(JB > JC.AND.JB < JD.AND.JA < JC) THEN
        CALL EL52(JJA,JJB,JA,JC,JB,JD,1,ICOLBREI)
      ELSE IF(JB > JC.AND.JB > JD.AND.JA > JC) THEN
        CALL EL52(JJA,JJB,JC,JA,JD,JB,2,ICOLBREI)
      ELSE IF(JB > JC.AND.JB > JD.AND.JA < JC) THEN
        CALL EL53(JJA,JJB,JA,JC,JD,JB,1,ICOLBREI)
      ELSE IF(JB > JC.AND.JB < JD.AND.JA > JC) THEN
        CALL EL53(JJA,JJB,JC,JA,JB,JD,2,ICOLBREI)
      ELSE
        WRITE(99,100)
        STOP
      END IF
      RETURN
  100 FORMAT(5X,'ERRO IN EL5 ')
      END SUBROUTINE EL5

 END MODULE el
      
