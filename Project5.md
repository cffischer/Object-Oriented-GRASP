This project restructures the librangular90 files (rangular90.f90) into a more object oriented arrangement as indicated in librang-module.tgz.

## How to rearrange subroutines into modules

An initial attempt has been as follows. Nearly all routines in the library are subroutines. The “call graph” is than closely related by the “calls” that are being made. Calls made from applications are “top level” routines and turned out to be few in number. Three are listed here along with their calls.

1) oneparticlejj.f90 with calls to
     CALL ONEPARTICLEJJ1(NS,KA,JA,JB,JA1,JA2,TCOEFF)
     CALL ONEPARTICLEJJ2(NS,KA,JA1,JA2,TCOEFF)
     CALL RECOP00(NS,JA,JA,KA,IAT)
     CALL RECOP1(NS,JA,KA,0,IAT,REC)
     CALL RECOP2(NS,IA,IB,KS1,KS2,KA,1,IAT,REC)
     CALL PERKO2(JA,JA,JA,JA,1)
     CALL WJ1(IK1,BK1,ID1,BD1,KA,QM1,QM2,WJ)
     CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A2)
     CALL RMEAJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S1)

2) onescalar.f90  with calls to
     CALL ONESCALAR1(NS,JA,JB,JA1,JA2,TCOEFF)
     CALL ONESCALAR2(JA,JB,JA1,JA2,TCOEFF)
     CALL RECOONESCALAR(NS,JA,JA,JA,JA,0,IAT)
     CALL PERKO2(JA,JA,JA,JA,1)
     CALL WJ1(IK1,BK1,ID1,BD1,0,QM1,QM2,WJ)
     CALL RECO2(JAA,JBB,ID2(3),0,IAT,REC)
     CALL GG12(IK1,IK2,BK1,BK2,ID1,ID2,BD1,BD2,QM1,QM2,WW)

3) el.f90  with calls to
     CALL EL1(JA,JB,JA1,JB1,1,ICOLBREI)
     CALL EL2(JA,JB,JA1,JA2,ICOLBREI)
     CALL EL3(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
     CALL EL4(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
     CALL EL5(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
     CALL CORD (JA,JB,JA1,IPCA,JB1

So we notice that some routines are only call by a single top-level routine and these are then routines that should be "contained" in a single module for top levels.  Others are called several so I suppose they need to be in a separate level 2  module. But if you follow through and look at the calls of the second level routines that are called only by one second level routine, you see that some need be be "contained". For example EL3 call  EL31,EL32,EL33. So we have

       EL CONTAINS EL3 which CONTAINS  EL31,EL32,EL33

The analysis is not complete but is essential for organizing the library.

##  Testing the new modules.  

Note that, rearranging the library into modules should not affect the execution of the application although maybe some changes may be desirable and with regard to USE statements.  

## Documentaion

An important goal as well is to include documentation, as best we can.