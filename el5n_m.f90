!*******************************************************************
 MODULE el5n
!                                                                  *
!      Module comtaining subroutines called by subroutine "el5"    *
!      found in the "el" module.                                   *
!                                                                  *
!      SUBROUTINE CALLED: COULOM,GG1122,ITREXG,IXJTIK,PERKO2,      *
!                         RECO,RECO2,SIXJ,SPEAK,WW1                *
!      Structure by A. Senchuk                     September 2019  *   
!                                                                  *
!*******************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, HALF, EPS
      USE m_C,             ONLY: NQ1, JLIST, NPEEL
      USE orb_C,           ONLY: NAK
      USE trk_C

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE reco_I
      USE reco4_I
      USE perko2_I
      USE snrc_I
      USE gg1234_I
      USE itrexg_I
      USE ixjtik_I
      USE coulom_I
      USE sixj_I
      USE speak_I
      USE itrig_I
      USE cxk_I
      USE talk_I

 PRIVATE  ! By default

!-------------------------------------------------------------
!   E x p l i c i t   P u b l i c  D e c l a r a t i o n s     
!-------------------------------------------------------------
      PUBLIC :: EL51
      PUBLIC :: EL52
      PUBLIC :: EL53

      
 CONTAINS

!*******************************************************************
!                                                                  *
      SUBROUTINE EL51(JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 12  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 1234, 2134, 1243, 2134    *
!                                                   ( IREZ = 1),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 - 1   *
!                                                   N'3 = N3 + 1   *
!                                                   N'4 = N4 + 1   *
!     AND    3412, 4321, 3421, 4312                 ( IREZ = 2),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 + 1   *
!                                                   N'2 = N2 + 1   *
!                                                   N'3 = N3 - 1   *
!                                                   N'4 = N4 - 1   *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,GG1234,ITREXG,IXJTIK,PERKO2,       *
!                        RECO,REC4,SIXJ,SPEAK                      *
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
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI
!      DIMENSION PMGG(30),J(4)
!      DIMENSION COND(12,20),CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA,IB,IC,ID,IBRD,IBRE,II,IP1,IP2,IG1,IG2,IKK,I2,I3,  &
                 I4,IFAZ,IFAZP,IFAZFRCS,INN,IAT,KRA,KRA1,L1,L2,L3,L4, &
                 J12,JB1,JD1,ND1,ND2,NE1,NE2,N,NN,NU,NUP1,MU
      INTEGER, DIMENSION(4) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(DOUBLE)          :: AA,AB,A1,BB,QM1,QM2,QM3,QM4,RAG,RECC,SI
      REAL(DOUBLE), DIMENSION(12) :: S
      REAL(DOUBLE), DIMENSION(30) :: PMGG
      REAL(DOUBLE), DIMENSION(12,20) :: COND,CONE
!-----------------------------------------------
      IF(NPEEL <= 3)RETURN
      CALL RECO(JA,JD,JC,JB,3,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      ID=JLIST(JD)
      IF(IREZ == 2) THEN
        QM1=-HALF
        QM2=-HALF
        QM3=HALF
        QM4=HALF
      ELSE
        QM1=HALF
        QM2=HALF
        QM3=-HALF
        QM4=-HALF
      END IF
      CALL PERKO2(JA,JB,JC,JD,4)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      J(4)=ID4(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      L4=(J(4)+1)/2
      IF (ICOLBREI == 2) THEN
        IF(IREZ == 1) THEN
          IS(1)=IA
          IS(2)=IB
          IS(3)=IC
          IS(4)=ID
        ELSE IF(IREZ == 2) THEN
          IS(1)=IC
          IS(2)=ID
          IS(3)=IA
          IS(4)=IB
        END IF
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0 .AND. IBRE <= 0)RETURN
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
          CONE(1,II) =ZERO
          CONE(2,II) =ZERO
          CONE(3,II) =ZERO
          CONE(4,II) =ZERO
          CONE(5,II) =ZERO
          CONE(6,II) =ZERO
          CONE(7,II) =ZERO
          CONE(8,II) =ZERO
          CONE(9,II) =ZERO
          CONE(10,II)=ZERO
          CONE(11,II)=ZERO
          CONE(12,II)=ZERO
        END DO
      END IF
      CALL GG1234(IK1,IK2,IK3,IK4,BK1,BK2,BK3,BK4,ID1,ID2,   &
                  ID3,ID4,BD1,BD2,BD3,BD4,QM1,QM2,QM3,QM4,RAG)
      IF(DABS(RAG) < EPS) RETURN
      IP1=ITREXG(J(1),J(2),J(3),J(4),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      DO I4=IP1,IG1,2
        KRA=(I4-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        PMGG(KRA1)=ZERO
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,0,IAT,RECC)
        IF(IAT == 0) CYCLE
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,1,IAT,RECC)
        PMGG(KRA1)=RECC
      END DO
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      IFAZFRCS=1
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)  &
      +IK3(5)*IK3(4)-ID3(5)*ID3(4)+IK4(5)*IK4(4)-ID4(5)*ID4(4)
      IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
      NN=0
      JB1=JB-1
      IFAZP=1
      DO II=JA,JB1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
      NN=0
      JD1=JD-1
      DO II=JC,JD1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
! * * *                      * * *                      * * *
!     CASES 1234   + + - -
!           2134                  TRANSFORM TO  1234   + + - -
!                                                    (IREZ = 1)
!     OR
!     CASES 3412   + + - -        TRANSFORM TO  1234   - - + +
!           3421                                1234
!                                                    (IREZ = 2)
      IP2=ITREXG(J(1),J(3),J(2),J(4),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
!
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L3,L4,L1,L2,ID3(5),ID4(5),ID1(5),ID2(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L2,L3,L4,ID1(5),ID2(5),ID3(5),ID4(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
!
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(1)+J(4)-2*J12
          IF(IREZ == 2)IFAZ=J(2)+J(3)-2*J12
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(3),KRA*2,J(4),J(2),J12*2) == 0) CYCLE
          AA=AA*RAG
          IF(DABS(AA) < EPS) CYCLE
          CALL SIXJ(J(1),J(3),KRA*2,J(4),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=J(2)+J(3)+2*J12+2*KRA
          IF(IREZ == 2)IFAZ=J(1)+J(4)-2*J12+2*KRA
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=-AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,IB,IC,ID,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IC,ID,IA,IB,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
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
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-ND1)/2+1
                IF(N < ND2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                  DO MU = 1,12
                    COND(MU,N)=COND(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
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
! * * *                      * * *                      * * *
!     CASES 1243   + + - -        TRANSFORM TO  1234   + + - -
!           2134                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 3421   + + - -        TRANSFORM TO  1234   - - + +
!           4321                                1234
!                                                    (IREZ = 2)
      IP2=ITREXG(J(1),J(4),J(2),J(3),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
!
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L3,L4,L2,L1,ID3(5),ID4(5),ID2(5),ID1(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L2,L4,L3,ID1(5),ID2(5),ID4(5),ID3(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
!
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(1)-J(4)
          IF(IREZ == 2)IFAZ=J(3)-J(2)
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAG
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(4),KRA*2,J(3),J(2),J12*2) == 0) CYCLE
          CALL SIXJ(J(1),J(4),KRA*2,J(3),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=J(2)+J(3)+2*J(4)+2*KRA
          IF(IREZ == 2)IFAZ=J(1)+2*J(2)+J(4)+4*J12+2*KRA
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,IB,ID,IC,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IC,ID,IB,IA,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+1) /= 0) .AND.   &
               (ITRIG(KS(2),KS(3),NU+NU+1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU-1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(N <= NE2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
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
                  CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                  DO MU = 1,12
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
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(2),IS(3),1,CONE(1,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(3),IS(2),1,CONE(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(3),IS(2),1,CONE(3,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(2),IS(3),1,CONE(4,N))
          IF(N == NE2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(2),IS(3),2,CONE(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(1),IS(4),2,CONE(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(3),IS(2),2,CONE(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(4),IS(1),2,CONE(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(3),IS(2),2,CONE(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(1),IS(4),2,CONE(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(2),IS(3),2,CONE(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(4),IS(1),2,CONE(12,N))
        END DO
      END IF
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL51  PMGG RAGG')
      STOP
      END SUBROUTINE EL51
   
!*******************************************************************
!                                                                  *
      SUBROUTINE EL52(JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 13  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 1324, 3142, 1342, 3124    *
!                                                   ( IREZ = 1),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                   N'3 = N3 - 1   *
!                                                   N'4 = N4 + 1   *
!     AND    2413, 4231, 2431, 4213                 ( IREZ = 2),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 + 1   *
!                                                   N'2 = N2 - 1   *
!                                                   N'3 = N3 + 1   *
!                                                   N'4 = N4 - 1   *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,GG1234,ITREXG,IXJTIK,PERKO2,       *
!                        RECO,REC4,SIXJ,SPEAK                      *
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
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI
!      DIMENSION PMGG(30),J(4)
!     DIMENSION COND(12,20),CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA,IB,IC,ID,IBRD,IBRE,II,IP1,IP2,IG1,IG2,IKK,I2,I3,  &
                 I4,IFAZ,IFAZP,IFAZFRCS,INN,IAT,KRA,KRA1,L1,L2,L3,L4, &
                 J12,JB1,JD1,ND1,ND2,NE1,NE2,N,NN,NU,NUP1,MU
      INTEGER, DIMENSION(4) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(DOUBLE)          :: AA,AB,A1,BB,QM1,QM2,QM3,QM4,RAG,RECC,SI
      REAL(DOUBLE), DIMENSION(12) :: S
      REAL(DOUBLE), DIMENSION(30) :: PMGG
      REAL(DOUBLE), DIMENSION(12,20) :: COND,CONE
!-----------------------------------------------
      IF(NPEEL <= 3)RETURN
      CALL RECO(JA,JD,JC,JB,3,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      ID=JLIST(JD)
      IF(IREZ == 2) THEN
        QM1=-HALF
        QM2=HALF
        QM3=-HALF
        QM4=HALF
      ELSE
        QM1=HALF
        QM2=-HALF
        QM3=HALF
        QM4=-HALF
      END IF
      CALL PERKO2(JA,JB,JC,JD,4)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      J(4)=ID4(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      L4=(J(4)+1)/2
      IF (ICOLBREI == 2) THEN
        IF(IREZ == 1) THEN
          IS(1)=IA
          IS(2)=IC
          IS(3)=IB
          IS(4)=ID
        ELSE IF(IREZ == 2) THEN
          IS(1)=IB
          IS(2)=ID
          IS(3)=IA
          IS(4)=IC
        END IF
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0 .AND. IBRE <= 0)RETURN
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
          CONE(1,II) =ZERO
          CONE(2,II) =ZERO
          CONE(3,II) =ZERO
          CONE(4,II) =ZERO
          CONE(5,II) =ZERO
          CONE(6,II) =ZERO
          CONE(7,II) =ZERO
          CONE(8,II) =ZERO
          CONE(9,II) =ZERO
          CONE(10,II)=ZERO
          CONE(11,II)=ZERO
          CONE(12,II)=ZERO
        END DO
      END IF
      CALL GG1234(IK1,IK2,IK3,IK4,BK1,BK2,BK3,BK4,ID1,ID2,  &
      ID3,ID4,BD1,BD2,BD3,BD4,QM1,QM2,QM3,QM4,RAG)
      IF(DABS(RAG) < EPS) RETURN
      IP1=ITREXG(J(1),J(2),J(3),J(4),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      DO I4=IP1,IG1,2
        KRA=(I4-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        PMGG(KRA1)=ZERO
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,0,IAT,RECC)
        IF(IAT == 0) CYCLE
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,1,IAT,RECC)
        PMGG(KRA1)=RECC
      END DO
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      IFAZFRCS=1
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)  &
      +IK3(5)*IK3(4)-ID3(5)*ID3(4)+IK4(5)*IK4(4)-ID4(5)*ID4(4)
      IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
      NN=0
      JB1=JB-1
      IFAZP=1
      DO II=JA,JB1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
      NN=0
      JD1=JD-1
      DO II=JC,JD1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
! * * *                      * * *                      * * *
!     CASES 1324   + + - -        TRANSFORM TO  1234   + - + -
!           1342                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2413   + + - -        TRANSFORM TO  1234   - + - +
!           4231                                1234
!                                                    (IREZ = 2)
      DO I3=IP1,IG1,2
        KRA=(I3-1)/2
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L2,L4,L1,L3,ID2(5),ID4(5),ID1(5),ID3(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L3,L2,L4,ID1(5),ID3(5),ID2(5),ID4(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        AA=PMGG(KRA1)
        IF(DABS(AA) < EPS) CYCLE
        AA=AA*RAG
        IF(DABS(AA) < EPS) CYCLE
        AA=AA/DSQRT(DBLE(I3))
        AB=AA*DBLE(IFAZP)
        IF(IREZ == 2) THEN
          IFAZ=J(1)+J(2)+J(4)+J(3)+4*KRA
          IF((IFAZ/4)*4 /= IFAZ)AB=-AB
        END IF
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,IC,IB,ID,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IB,ID,IA,IC,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU >  0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
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
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
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
                    COND(MU,N)=COND(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
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
! * * *                      * * *                      * * *
!     CASES 1342   + + - -        TRANSFORM TO  1234   + - + -
!           3124                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2431   + + - -        TRANSFORM TO  1234   - + - +
!           4213                                1234
!                                                    (IREZ = 2)
      IP2=ITREXG(J(1),J(4),J(2),J(3),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
!
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L2,L4,L3,L1,ID2(5),ID4(5),ID3(5),ID1(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L3,L4,L2,ID1(5),ID3(5),ID4(5),ID2(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
!
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(1)-J(3)+2*J12
          IF(IREZ == 2)IFAZ=J(4)-J(2)-2*J12
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAG
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(4),KRA*2,J(3),J(2),J12*2) == 0) CYCLE
          CALL SIXJ(J(1),J(4),KRA*2,J(3),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=2*J(3)-4*KRA-4*J12
          IF(IREZ == 2)IFAZ=J(1)+J(2)+J(3)-J(4)-4*KRA
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=-AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,IC,ID,IB,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IB,ID,IC,IA,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU-1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(N <= NE2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
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
                  CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                  DO MU = 1,12
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
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(2),IS(3),1,CONE(1,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(3),IS(2),1,CONE(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(3),IS(2),1,CONE(3,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(2),IS(3),1,CONE(4,N))
          IF(N == NE2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(2),IS(3),2,CONE(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(1),IS(4),2,CONE(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(3),IS(2),2,CONE(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(4),IS(1),2,CONE(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(3),IS(2),2,CONE(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(1),IS(4),2,CONE(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(2),IS(3),2,CONE(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(4),IS(1),2,CONE(12,N))
        END DO
      END IF
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL52  PMGG RAGG')
      STOP
      END SUBROUTINE EL52

   
!*******************************************************************
!                                                                  *
      SUBROUTINE EL53(JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 14  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 1423, 4132, 1432, 4123    *
!                                                   ( IREZ = 1),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                   N'3 = N3 + 1   *
!                                                   N'4 = N4 - 1   *
!     AND    2314, 3241, 2341, 3214                 ( IREZ = 2),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 + 1   *
!                                                   N'2 = N2 - 1   *
!                                                   N'3 = N3 - 1   *
!                                                   N'4 = N4 + 1   *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,GG1234,ITREXG,IXJTIK,PERKO2,       *
!                      RECO,REC4,SIXJ,SPEAK                        *
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
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI
!      DIMENSION PMGG(30),J(4)
!      DIMENSION COND(12,20),CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA,IB,IC,ID,IBRD,IBRE,II,IP1,IP2,IG1,IG2,IKK,I2,I3,  &
                 I4,IFAZ,IFAZP,IFAZFRCS,INN,IAT,KRA,KRA1,L1,L2,L3,L4, &
                 J12,JB1,JD1,ND1,ND2,NE1,NE2,N,NN,NU,NUP1,MU
      INTEGER, DIMENSION(4) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(DOUBLE)          :: AA,AB,A1,BB,QM1,QM2,QM3,QM4,RAG,RECC,SI
      REAL(DOUBLE), DIMENSION(12) :: S
      REAL(DOUBLE), DIMENSION(30) :: PMGG
      REAL(DOUBLE), DIMENSION(12,20) :: COND,CONE
!-----------------------------------------------
      IF(NPEEL <= 3)RETURN
      CALL RECO(JA,JD,JC,JB,3,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      ID=JLIST(JD)
      IF(IREZ == 2) THEN
        QM1=-HALF
        QM2=HALF
        QM3=HALF
        QM4=-HALF
      ELSE
        QM1=HALF
        QM2=-HALF
        QM3=-HALF
        QM4=HALF
      END IF
      CALL PERKO2(JA,JB,JC,JD,4)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      J(4)=ID4(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      L4=(J(4)+1)/2
      IF (ICOLBREI == 2) THEN
        IF(IREZ == 1) THEN
          IS(1)=IA
          IS(2)=ID
          IS(3)=IB
          IS(4)=IC
        ELSE IF(IREZ == 2) THEN
          IS(1)=IB
          IS(2)=IC
          IS(3)=IA
          IS(4)=ID
        END IF
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0 .AND. IBRE <= 0)RETURN
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
          CONE(1,II) =ZERO
          CONE(2,II) =ZERO
          CONE(3,II) =ZERO
          CONE(4,II) =ZERO
          CONE(5,II) =ZERO
          CONE(6,II) =ZERO
          CONE(7,II) =ZERO
          CONE(8,II) =ZERO
          CONE(9,II) =ZERO
          CONE(10,II)=ZERO
          CONE(11,II)=ZERO
          CONE(12,II)=ZERO
        END DO
      END IF
      CALL GG1234(IK1,IK2,IK3,IK4,BK1,BK2,BK3,BK4,ID1,ID2,  &
      ID3,ID4,BD1,BD2,BD3,BD4,QM1,QM2,QM3,QM4,RAG)
      IF(DABS(RAG) < EPS) RETURN
      IP1=ITREXG(J(1),J(2),J(3),J(4),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      DO I4=IP1,IG1,2
        KRA=(I4-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        PMGG(KRA1)=ZERO
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,0,IAT,RECC)
        IF(IAT == 0) CYCLE
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,1,IAT,RECC)
        PMGG(KRA1)=RECC
      END DO
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      IFAZFRCS=1
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4) &
      +IK3(5)*IK3(4)-ID3(5)*ID3(4)+IK4(5)*IK4(4)-ID4(5)*ID4(4)
      IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
      NN=0
      JB1=JB-1
      IFAZP=1
      DO II=JA,JB1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
      NN=0
      JD1=JD-1
      DO II=JC,JD1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
! * * *                      * * *                      * * *
!     CASES 1423   + + - -        TRANSFORM TO  1234   + - - +
!           4132                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2314   + + - -        TRANSFORM TO  1234   - + + -
!           3241                                1234
!                                                    (IREZ = 2)
      DO I3=IP1,IG1,2
        KRA=(I3-1)/2
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L2,L3,L1,L4,ID2(5),ID3(5),ID1(5),ID4(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L4,L2,L3,ID1(5),ID4(5),ID2(5),ID3(5),KRA,A1)
          END IF
        IF(DABS(A1) < EPS) CYCLE
        END IF
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        AA=PMGG(KRA1)
        IF(DABS(AA) < EPS) CYCLE
        AA=AA*RAG
        IF(DABS(AA) < EPS) CYCLE
        AA=AA/DSQRT(DBLE(I3))
        IFAZ=J(4)+J(3)-2*KRA+2
        IF(IREZ == 2)IFAZ=J(1)+J(2)-2*KRA+2
        IF((IFAZ/4)*4 /= IFAZ)AA=-AA
        AB=AA*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,ID,IB,IC,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IB,IC,IA,ID,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND. &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
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
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-ND1)/2+1
                IF(N < ND2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                  DO MU = 1,12
                    COND(MU,N)=COND(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
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
! * * *                      * * *                      * * *
!     CASES 1432   + + - -        TRANSFORM TO  1234   + - - +
!           4132                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2341   + + - -        TRANSFORM TO  1234   - + + -
!           3214                                1234
!                                                    (IREZ = 2)
      IP2=ITREXG(J(1),J(3),J(4),J(2),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L2,L3,L4,L1,ID2(5),ID3(5),ID4(5),ID1(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L4,L3,L2,ID1(5),ID4(5),ID3(5),ID2(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(1)+J(4)+2*J12
          IF(IREZ == 2)IFAZ=J(2)+J(3)+2*J12
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAG
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(3),KRA*2,J(4),J(2),J12*2) == 0) CYCLE
          CALL SIXJ(J(1),J(3),KRA*2,J(4),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=J(3)-J(4)-4*KRA+2*J12
          IF(IREZ == 2)IFAZ=J(1)-J(2)-4*KRA-2*J12
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,ID,IC,IB,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IB,IC,ID,IA,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU-1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(N <= NE2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
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
                  CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                  DO MU = 1,12
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
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(2),IS(3),1,CONE(1,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(3),IS(2),1,CONE(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(3),IS(2),1,CONE(3,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(2),IS(3),1,CONE(4,N))
          IF(N == NE2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(2),IS(3),2,CONE(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(1),IS(4),2,CONE(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(3),IS(2),2,CONE(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(4),IS(1),2,CONE(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(3),IS(2),2,CONE(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(1),IS(4),2,CONE(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(2),IS(3),2,CONE(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(4),IS(1),2,CONE(12,N))
        END DO
      END IF
      RETURN
10    WRITE(99,100)KRA1
100   FORMAT(5X,'ERRO IN EL53  PMGG RAGG KRA1=',I100)
      STOP
      END SUBROUTINE EL53

 END MODULE el5n
