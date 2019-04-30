************************************************************************
*                                                                      *
      PROGRAM RNUCLEUS
*                                                                      *
*   Generates the isotope data file for the GRASP92 suite of codes.    *
*                                                                      *
*   Call(s) to: GETCPR, GETYN, LENGTH, OPENFL.                         *
*                                                                      *
*   Written by Farid A. Parpia.           Last revision: 16 Oct 1994   *
*                                                                      *
*   Update: 2014-02-13 - Jon Grumer, Lund University, Sweden           *
*           Grid parameters are written to isodata so other routines   *
*           may read this file instead of asking the user (non-default *
*           grid parameters).                                          *
*   Update: 2016-01-01 - Jörgen Ekman, Malmö University, Sweden        *
*           rms radius from I. Angeli and K.P. Marinova,               * 
*           Atomic Data and Nuclear Data Tables 99 (2013) 69–95,       *
*           is read in by default. Parametrization is only used in     *
*           cases when tabulated rms radius is missing.                *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      CHARACTER*256 FILNAM
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
      LOGICAL GETYN,YES
*     
      ! JG (Lund, 2013) ---------------------------------------------------------
      INCLUDE 'parameters.def'       ! Grasp2k parameter file (src/lib/def/)
      CHARACTER*1 GRIDANSW,GRIDANSW2 ! Answers to the default grid questions
      DOUBLE PRECISION RNT, H, HP    ! Grid parameters to be included in isodata
      INTEGER N                      ! ...
      ! -------------------------------------------------------------------------

      ! JE Jan 2016
      DOUBLE PRECISION RRMSVEC(120,300)
*     
      EXTERNAL CONSTS
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
      COMMON/iounit/istdi,istdo,istde

cb alpha constant from lib/lib92/setcon.f
cb AUMAMU from lib/lib92/setcon.f
cb
cb    DATA EMEAMU /5.48579903D-04/
cb   :     ALFAI  /137.0359895D 00/
      COMMON/DEF2/C
     :      /DEF11/FMTOAU,AUMAMU
cb
      ALFAI = 1.0D 00/C
cb EMEAMU: Electron mass in amu
      AUMAMU = EMEAMU
cb end alpha constant

! JE Jan 2016
      RRMSVEC(:,:) = 0.0D 00
      RRMSVEC(1,1) = 0.8783D 00
      RRMSVEC(1,2) = 2.1421D 00
      RRMSVEC(1,3) = 1.7591D 00
      RRMSVEC(2,3) = 1.9661D 00
      RRMSVEC(2,4) = 1.6755D 00
      RRMSVEC(2,6) = 2.0660D 00
      RRMSVEC(2,8) = 1.9239D 00
      RRMSVEC(3,6) = 2.5890D 00
      RRMSVEC(3,7) = 2.4440D 00
      RRMSVEC(3,8) = 2.3390D 00
      RRMSVEC(3,9) = 2.2450D 00
      RRMSVEC(3,11) = 2.4820D 00
      RRMSVEC(4,7) = 2.6460D 00
      RRMSVEC(4,9) = 2.5190D 00
      RRMSVEC(4,10) = 2.3550D 00
      RRMSVEC(4,11) = 2.4630D 00
      RRMSVEC(5,10) = 2.4277D 00
      RRMSVEC(5,11) = 2.4060D 00
      RRMSVEC(6,12) = 2.4702D 00
      RRMSVEC(6,13) = 2.4614D 00
      RRMSVEC(6,14) = 2.5025D 00
      RRMSVEC(7,14) = 2.5582D 00
      RRMSVEC(7,15) = 2.6058D 00
      RRMSVEC(8,16) = 2.6991D 00
      RRMSVEC(8,17) = 2.6932D 00
      RRMSVEC(8,18) = 2.7726D 00
      RRMSVEC(9,19) = 2.8976D 00
      RRMSVEC(10,17) = 3.0413D 00
      RRMSVEC(10,18) = 2.9714D 00
      RRMSVEC(10,19) = 3.0082D 00
      RRMSVEC(10,20) = 3.0055D 00
      RRMSVEC(10,21) = 2.9695D 00
      RRMSVEC(10,22) = 2.9525D 00
      RRMSVEC(10,23) = 2.9104D 00
      RRMSVEC(10,24) = 2.9007D 00
      RRMSVEC(10,25) = 2.9316D 00
      RRMSVEC(10,26) = 2.9251D 00
      RRMSVEC(10,28) = 2.9642D 00
      RRMSVEC(11,20) = 2.9718D 00
      RRMSVEC(11,21) = 3.0136D 00
      RRMSVEC(11,22) = 2.9852D 00
      RRMSVEC(11,23) = 2.9936D 00
      RRMSVEC(11,24) = 2.9735D 00
      RRMSVEC(11,25) = 2.9769D 00
      RRMSVEC(11,26) = 2.9928D 00
      RRMSVEC(11,27) = 3.0136D 00
      RRMSVEC(11,28) = 3.0400D 00
      RRMSVEC(11,29) = 3.0922D 00
      RRMSVEC(11,30) = 3.1180D 00
      RRMSVEC(11,31) = 3.1704D 00
      RRMSVEC(12,24) = 3.0570D 00
      RRMSVEC(12,25) = 3.0284D 00
      RRMSVEC(12,26) = 3.0337D 00
      RRMSVEC(13,27) = 3.0610D 00
      RRMSVEC(14,28) = 3.1224D 00
      RRMSVEC(14,29) = 3.1176D 00
      RRMSVEC(14,30) = 3.1336D 00
      RRMSVEC(15,31) = 3.1889D 00
      RRMSVEC(16,32) = 3.2611D 00
      RRMSVEC(16,34) = 3.2847D 00
      RRMSVEC(16,36) = 3.2985D 00
      RRMSVEC(17,35) = 3.3654D 00
      RRMSVEC(17,37) = 3.3840D 00
      RRMSVEC(18,32) = 3.3468D 00
      RRMSVEC(18,33) = 3.3438D 00
      RRMSVEC(18,34) = 3.3654D 00
      RRMSVEC(18,35) = 3.3636D 00
      RRMSVEC(18,36) = 3.3905D 00
      RRMSVEC(18,37) = 3.3908D 00
      RRMSVEC(18,38) = 3.4028D 00
      RRMSVEC(18,39) = 3.4093D 00
      RRMSVEC(18,40) = 3.4274D 00
      RRMSVEC(18,41) = 3.4251D 00
      RRMSVEC(18,43) = 3.4414D 00
      RRMSVEC(18,42) = 3.4354D 00
      RRMSVEC(18,44) = 3.4454D 00
      RRMSVEC(18,46) = 3.4377D 00
      RRMSVEC(19,38) = 3.4264D 00
      RRMSVEC(19,39) = 3.4349D 00
      RRMSVEC(19,40) = 3.4381D 00
      RRMSVEC(19,41) = 3.4518D 00
      RRMSVEC(19,42) = 3.4517D 00
      RRMSVEC(19,43) = 3.4556D 00
      RRMSVEC(19,44) = 3.4563D 00
      RRMSVEC(19,45) = 3.4605D 00
      RRMSVEC(19,46) = 3.4558D 00
      RRMSVEC(19,47) = 3.4534D 00
      RRMSVEC(20,39) = 3.4595D 00
      RRMSVEC(20,40) = 3.4776D 00
      RRMSVEC(20,41) = 3.4780D 00
      RRMSVEC(20,42) = 3.5081D 00
      RRMSVEC(20,43) = 3.4954D 00
      RRMSVEC(20,44) = 3.5179D 00
      RRMSVEC(20,45) = 3.4944D 00
      RRMSVEC(20,46) = 3.4953D 00
      RRMSVEC(20,47) = 3.4783D 00
      RRMSVEC(20,48) = 3.4771D 00
      RRMSVEC(20,50) = 3.5168D 00
      RRMSVEC(21,42) = 3.5702D 00
      RRMSVEC(21,43) = 3.5575D 00
      RRMSVEC(21,44) = 3.5432D 00
      RRMSVEC(21,45) = 3.5459D 00
      RRMSVEC(21,46) = 3.5243D 00
      RRMSVEC(22,44) = 3.6115D 00
      RRMSVEC(22,45) = 3.5939D 00
      RRMSVEC(22,46) = 3.6070D 00
      RRMSVEC(22,47) = 3.5962D 00
      RRMSVEC(22,48) = 3.5921D 00
      RRMSVEC(22,49) = 3.5733D 00
      RRMSVEC(22,50) = 3.5704D 00
      RRMSVEC(23,51) = 3.6002D 00
      RRMSVEC(24,50) = 3.6588D 00
      RRMSVEC(24,52) = 3.6452D 00
      RRMSVEC(24,53) = 3.6511D 00
      RRMSVEC(24,54) = 3.6885D 00
      RRMSVEC(25,50) = 3.7120D 00
      RRMSVEC(25,51) = 3.7026D 00
      RRMSVEC(25,52) = 3.6706D 00
      RRMSVEC(25,53) = 3.6662D 00
      RRMSVEC(25,54) = 3.6834D 00
      RRMSVEC(25,55) = 3.7057D 00
      RRMSVEC(25,56) = 3.7146D 00
      RRMSVEC(26,54) = 3.6933D 00
      RRMSVEC(26,56) = 3.7377D 00
      RRMSVEC(26,57) = 3.7532D 00
      RRMSVEC(26,58) = 3.7745D 00
      RRMSVEC(27,59) = 3.7875D 00
      RRMSVEC(28,58) = 3.7757D 00
      RRMSVEC(28,60) = 3.8118D 00
      RRMSVEC(28,61) = 3.8225D 00
      RRMSVEC(28,62) = 3.8399D 00
      RRMSVEC(28,64) = 3.8572D 00
      RRMSVEC(29,63) = 3.8823D 00
      RRMSVEC(29,65) = 3.9022D 00
      RRMSVEC(30,64) = 3.9283D 00
      RRMSVEC(30,66) = 3.9491D 00
      RRMSVEC(30,67) = 3.9530D 00
      RRMSVEC(30,68) = 3.9658D 00
      RRMSVEC(30,70) = 3.9845D 00
      RRMSVEC(31,69) = 3.9973D 00
      RRMSVEC(31,71) = 4.0118D 00
      RRMSVEC(32,70) = 4.0414D 00
      RRMSVEC(32,72) = 4.0576D 00
      RRMSVEC(32,73) = 4.0632D 00
      RRMSVEC(32,74) = 4.0742D 00
      RRMSVEC(32,76) = 4.0811D 00
      RRMSVEC(33,75) = 4.0968D 00
      RRMSVEC(34,74) = 4.0700D 00
      RRMSVEC(34,76) = 4.1395D 00
      RRMSVEC(34,77) = 4.1395D 00
      RRMSVEC(34,78) = 4.1406D 00
      RRMSVEC(34,80) = 4.1400D 00
      RRMSVEC(34,82) = 4.1400D 00
      RRMSVEC(35,79) = 4.1629D 00
      RRMSVEC(35,81) = 4.1599D 00
      RRMSVEC(36,72) = 4.1635D 00
      RRMSVEC(36,74) = 4.1870D 00
      RRMSVEC(36,75) = 4.2097D 00
      RRMSVEC(36,76) = 4.2020D 00
      RRMSVEC(36,77) = 4.2082D 00
      RRMSVEC(36,78) = 4.2038D 00
      RRMSVEC(36,79) = 4.2034D 00
      RRMSVEC(36,80) = 4.1970D 00
      RRMSVEC(36,81) = 4.1952D 00
      RRMSVEC(36,82) = 4.1919D 00
      RRMSVEC(36,83) = 4.1871D 00
      RRMSVEC(36,84) = 4.1884D 00
      RRMSVEC(36,85) = 4.1846D 00
      RRMSVEC(36,86) = 4.1835D 00
      RRMSVEC(36,87) = 4.1984D 00
      RRMSVEC(36,88) = 4.2171D 00
      RRMSVEC(36,89) = 4.2286D 00
      RRMSVEC(36,90) = 4.2423D 00
      RRMSVEC(36,91) = 4.2543D 00
      RRMSVEC(36,92) = 4.2724D 00
      RRMSVEC(36,93) = 4.2794D 00
      RRMSVEC(36,94) = 4.3002D 00
      RRMSVEC(36,95) = 4.3067D 00
      RRMSVEC(36,96) = 4.3267D 00
      RRMSVEC(37,76) = 4.2273D 00
      RRMSVEC(37,77) = 4.2356D 00
      RRMSVEC(37,78) = 4.2385D 00
      RRMSVEC(37,79) = 4.2284D 00
      RRMSVEC(37,80) = 4.2271D 00
      RRMSVEC(37,81) = 4.2213D 00
      RRMSVEC(37,82) = 4.2160D 00
      RRMSVEC(37,83) = 4.2058D 00
      RRMSVEC(37,84) = 4.1999D 00
      RRMSVEC(37,85) = 4.2036D 00
      RRMSVEC(37,86) = 4.2025D 00
      RRMSVEC(37,87) = 4.1989D 00
      RRMSVEC(37,88) = 4.2170D 00
      RRMSVEC(37,89) = 4.2391D 00
      RRMSVEC(37,90) = 4.2554D 00
      RRMSVEC(37,91) = 4.2723D 00
      RRMSVEC(37,92) = 4.2903D 00
      RRMSVEC(37,93) = 4.3048D 00
      RRMSVEC(37,94) = 4.3184D 00
      RRMSVEC(37,95) = 4.3391D 00
      RRMSVEC(37,96) = 4.3501D 00
      RRMSVEC(37,97) = 4.4231D 00
      RRMSVEC(37,98) = 4.4336D 00
      RRMSVEC(38,77) = 4.2569D 00
      RRMSVEC(38,78) = 4.2561D 00
      RRMSVEC(38,79) = 4.2586D 00
      RRMSVEC(38,80) = 4.2562D 00
      RRMSVEC(38,81) = 4.2547D 00
      RRMSVEC(38,82) = 4.2478D 00
      RRMSVEC(38,83) = 4.2455D 00
      RRMSVEC(38,84) = 4.2394D 00
      RRMSVEC(38,85) = 4.2304D 00
      RRMSVEC(38,86) = 4.2307D 00
      RRMSVEC(38,87) = 4.2249D 00
      RRMSVEC(38,88) = 4.2240D 00
      RRMSVEC(38,89) = 4.2407D 00
      RRMSVEC(38,90) = 4.2611D 00
      RRMSVEC(38,91) = 4.2740D 00
      RRMSVEC(38,92) = 4.2924D 00
      RRMSVEC(38,93) = 4.3026D 00
      RRMSVEC(38,94) = 4.3191D 00
      RRMSVEC(38,95) = 4.3305D 00
      RRMSVEC(38,96) = 4.3522D 00
      RRMSVEC(38,97) = 4.3625D 00
      RRMSVEC(38,98) = 4.4377D 00
      RRMSVEC(38,99) = 4.4495D 00
      RRMSVEC(38,100) = 4.4640D 00
      RRMSVEC(39,86) = 4.2513D 00
      RRMSVEC(39,87) = 4.2498D 00
      RRMSVEC(39,88) = 4.2441D 00
      RRMSVEC(39,89) = 4.2430D 00
      RRMSVEC(39,90) = 4.2573D 00
      RRMSVEC(39,92) = 4.2887D 00
      RRMSVEC(39,93) = 4.3052D 00
      RRMSVEC(39,94) = 4.3142D 00
      RRMSVEC(39,95) = 4.3284D 00
      RRMSVEC(39,96) = 4.3402D 00
      RRMSVEC(39,97) = 4.3580D 00
      RRMSVEC(39,98) = 4.3711D 00
      RRMSVEC(39,99) = 4.4658D 00
      RRMSVEC(39,100) = 4.4705D 00
      RRMSVEC(39,101) = 4.4863D 00
      RRMSVEC(39,102) = 4.4911D 00
      RRMSVEC(40,87) = 4.2789D 00
      RRMSVEC(40,88) = 4.2787D 00
      RRMSVEC(40,89) = 4.2706D 00
      RRMSVEC(40,90) = 4.2694D 00
      RRMSVEC(40,91) = 4.2845D 00
      RRMSVEC(40,92) = 4.3057D 00
      RRMSVEC(40,94) = 4.3320D 00
      RRMSVEC(40,96) = 4.3512D 00
      RRMSVEC(40,97) = 4.3792D 00
      RRMSVEC(40,98) = 4.4012D 00
      RRMSVEC(40,99) = 4.4156D 00
      RRMSVEC(40,100) = 4.4891D 00
      RRMSVEC(40,101) = 4.5119D 00
      RRMSVEC(40,102) = 4.5292D 00
      RRMSVEC(41,90) = 4.2891D 00
      RRMSVEC(41,91) = 4.2878D 00
      RRMSVEC(41,92) = 4.3026D 00
      RRMSVEC(41,93) = 4.3240D 00
      RRMSVEC(41,99) = 4.4062D 00
      RRMSVEC(41,101) = 4.4861D 00
      RRMSVEC(41,103) = 4.5097D 00
      RRMSVEC(42,90) = 4.3265D 00
      RRMSVEC(42,91) = 4.3182D 00
      RRMSVEC(42,92) = 4.3151D 00
      RRMSVEC(42,94) = 4.3529D 00
      RRMSVEC(42,95) = 4.3628D 00
      RRMSVEC(42,96) = 4.3847D 00
      RRMSVEC(42,97) = 4.3880D 00
      RRMSVEC(42,98) = 4.4091D 00
      RRMSVEC(42,100) = 4.4468D 00
      RRMSVEC(42,102) = 4.4914D 00
      RRMSVEC(42,103) = 4.5145D 00
      RRMSVEC(42,104) = 4.5249D 00
      RRMSVEC(42,105) = 4.5389D 00
      RRMSVEC(42,106) = 4.5490D 00
      RRMSVEC(42,108) = 4.5602D 00
      RRMSVEC(44,96) = 4.3908D 00
      RRMSVEC(44,98) = 4.4229D 00
      RRMSVEC(44,99) = 4.4338D 00
      RRMSVEC(44,100) = 4.4531D 00
      RRMSVEC(44,101) = 4.4606D 00
      RRMSVEC(44,102) = 4.4809D 00
      RRMSVEC(44,104) = 4.5098D 00
      RRMSVEC(45,103) = 4.4945D 00
      RRMSVEC(46,102) = 4.4827D 00
      RRMSVEC(46,104) = 4.5078D 00
      RRMSVEC(46,105) = 4.5150D 00
      RRMSVEC(46,106) = 4.5318D 00
      RRMSVEC(46,108) = 4.5563D 00
      RRMSVEC(46,110) = 4.5782D 00
      RRMSVEC(47,101) = 4.4799D 00
      RRMSVEC(47,103) = 4.5036D 00
      RRMSVEC(47,104) = 4.5119D 00
      RRMSVEC(47,105) = 4.5269D 00
      RRMSVEC(47,107) = 4.5454D 00
      RRMSVEC(47,109) = 4.5638D 00
      RRMSVEC(48,102) = 4.4810D 00
      RRMSVEC(48,103) = 4.4951D 00
      RRMSVEC(48,104) = 4.5122D 00
      RRMSVEC(48,105) = 4.5216D 00
      RRMSVEC(48,106) = 4.5383D 00
      RRMSVEC(48,107) = 4.5466D 00
      RRMSVEC(48,108) = 4.5577D 00
      RRMSVEC(48,109) = 4.5601D 00
      RRMSVEC(48,110) = 4.5765D 00
      RRMSVEC(48,111) = 4.5845D 00
      RRMSVEC(48,112) = 4.5944D 00
      RRMSVEC(48,113) = 4.6012D 00
      RRMSVEC(48,114) = 4.6087D 00
      RRMSVEC(48,115) = 4.6114D 00
      RRMSVEC(48,116) = 4.6203D 00
      RRMSVEC(48,117) = 4.6136D 00
      RRMSVEC(48,118) = 4.6246D 00
      RRMSVEC(48,120) = 4.6300D 00
      RRMSVEC(49,104) = 4.5184D 00
      RRMSVEC(49,105) = 4.5311D 00
      RRMSVEC(49,106) = 4.5375D 00
      RRMSVEC(49,107) = 4.5494D 00
      RRMSVEC(49,108) = 4.5571D 00
      RRMSVEC(49,109) = 4.5685D 00
      RRMSVEC(49,110) = 4.5742D 00
      RRMSVEC(49,111) = 4.5856D 00
      RRMSVEC(49,112) = 4.5907D 00
      RRMSVEC(49,113) = 4.6010D 00
      RRMSVEC(49,114) = 4.6056D 00
      RRMSVEC(49,115) = 4.6156D 00
      RRMSVEC(49,116) = 4.6211D 00
      RRMSVEC(49,117) = 4.6292D 00
      RRMSVEC(49,118) = 4.6335D 00
      RRMSVEC(49,119) = 4.6407D 00
      RRMSVEC(49,120) = 4.6443D 00
      RRMSVEC(49,121) = 4.6505D 00
      RRMSVEC(49,122) = 4.6534D 00
      RRMSVEC(49,123) = 4.6594D 00
      RRMSVEC(49,124) = 4.6625D 00
      RRMSVEC(49,125) = 4.6670D 00
      RRMSVEC(49,126) = 4.6702D 00
      RRMSVEC(49,127) = 4.6733D 00
      RRMSVEC(50,108) = 4.5605D 00
      RRMSVEC(50,109) = 4.5679D 00
      RRMSVEC(50,110) = 4.5785D 00
      RRMSVEC(50,111) = 4.5836D 00
      RRMSVEC(50,112) = 4.5948D 00
      RRMSVEC(50,113) = 4.6015D 00
      RRMSVEC(50,114) = 4.6099D 00
      RRMSVEC(50,115) = 4.6148D 00
      RRMSVEC(50,116) = 4.6250D 00
      RRMSVEC(50,117) = 4.6302D 00
      RRMSVEC(50,118) = 4.6393D 00
      RRMSVEC(50,119) = 4.6438D 00
      RRMSVEC(50,120) = 4.6519D 00
      RRMSVEC(50,121) = 4.6566D 00
      RRMSVEC(50,122) = 4.6634D 00
      RRMSVEC(50,123) = 4.6665D 00
      RRMSVEC(50,124) = 4.6735D 00
      RRMSVEC(50,125) = 4.6765D 00
      RRMSVEC(50,126) = 4.6833D 00
      RRMSVEC(50,127) = 4.6867D 00
      RRMSVEC(50,128) = 4.6921D 00
      RRMSVEC(50,129) = 4.6934D 00
      RRMSVEC(50,130) = 4.7019D 00
      RRMSVEC(50,131) = 4.7078D 00
      RRMSVEC(50,132) = 4.7093D 00
      RRMSVEC(51,121) = 4.6802D 00
      RRMSVEC(51,123) = 4.6879D 00
      RRMSVEC(52,116) = 4.6847D 00
      RRMSVEC(52,118) = 4.6956D 00
      RRMSVEC(52,120) = 4.7038D 00
      RRMSVEC(52,122) = 4.7095D 00
      RRMSVEC(52,123) = 4.7117D 00
      RRMSVEC(52,124) = 4.7183D 00
      RRMSVEC(52,125) = 4.7204D 00
      RRMSVEC(52,126) = 4.7266D 00
      RRMSVEC(52,128) = 4.7346D 00
      RRMSVEC(52,130) = 4.7423D 00
      RRMSVEC(52,132) = 4.7500D 00
      RRMSVEC(52,134) = 4.7569D 00
      RRMSVEC(52,136) = 4.7815D 00
      RRMSVEC(53,127) = 4.7500D 00
      RRMSVEC(54,116) = 4.7211D 00
      RRMSVEC(54,118) = 4.7387D 00
      RRMSVEC(54,120) = 4.7509D 00
      RRMSVEC(54,122) = 4.7590D 00
      RRMSVEC(54,124) = 4.7661D 00
      RRMSVEC(54,126) = 4.7722D 00
      RRMSVEC(54,127) = 4.7747D 00
      RRMSVEC(54,128) = 4.7774D 00
      RRMSVEC(54,129) = 4.7775D 00
      RRMSVEC(54,130) = 4.7818D 00
      RRMSVEC(54,131) = 4.7808D 00
      RRMSVEC(54,132) = 4.7859D 00
      RRMSVEC(54,133) = 4.7831D 00
      RRMSVEC(54,134) = 4.7899D 00
      RRMSVEC(54,136) = 4.7964D 00
      RRMSVEC(54,137) = 4.8094D 00
      RRMSVEC(54,138) = 4.8279D 00
      RRMSVEC(54,139) = 4.8409D 00
      RRMSVEC(54,140) = 4.8566D 00
      RRMSVEC(54,141) = 4.8694D 00
      RRMSVEC(54,142) = 4.8841D 00
      RRMSVEC(54,143) = 4.8942D 00
      RRMSVEC(54,144) = 4.9082D 00
      RRMSVEC(54,146) = 4.9315D 00
      RRMSVEC(55,118) = 4.7832D 00
      RRMSVEC(55,119) = 4.7896D 00
      RRMSVEC(55,120) = 4.7915D 00
      RRMSVEC(55,121) = 4.7769D 00
      RRMSVEC(55,122) = 4.7773D 00
      RRMSVEC(55,123) = 4.7820D 00
      RRMSVEC(55,124) = 4.7828D 00
      RRMSVEC(55,125) = 4.7880D 00
      RRMSVEC(55,126) = 4.7872D 00
      RRMSVEC(55,127) = 4.7936D 00
      RRMSVEC(55,128) = 4.7921D 00
      RRMSVEC(55,129) = 4.7981D 00
      RRMSVEC(55,130) = 4.7992D 00
      RRMSVEC(55,131) = 4.8026D 00
      RRMSVEC(55,132) = 4.8002D 00
      RRMSVEC(55,133) = 4.8041D 00
      RRMSVEC(55,134) = 4.8031D 00
      RRMSVEC(55,135) = 4.8067D 00
      RRMSVEC(55,136) = 4.8059D 00
      RRMSVEC(55,137) = 4.8128D 00
      RRMSVEC(55,138) = 4.8255D 00
      RRMSVEC(55,139) = 4.8422D 00
      RRMSVEC(55,140) = 4.8554D 00
      RRMSVEC(55,141) = 4.8689D 00
      RRMSVEC(55,142) = 4.8825D 00
      RRMSVEC(55,143) = 4.8965D 00
      RRMSVEC(55,144) = 4.9055D 00
      RRMSVEC(55,145) = 4.9188D 00
      RRMSVEC(55,146) = 4.9281D 00
      RRMSVEC(56,120) = 4.8092D 00
      RRMSVEC(56,121) = 4.8176D 00
      RRMSVEC(56,122) = 4.8153D 00
      RRMSVEC(56,123) = 4.8135D 00
      RRMSVEC(56,124) = 4.8185D 00
      RRMSVEC(56,125) = 4.8177D 00
      RRMSVEC(56,126) = 4.8221D 00
      RRMSVEC(56,127) = 4.8204D 00
      RRMSVEC(56,128) = 4.8255D 00
      RRMSVEC(56,129) = 4.8248D 00
      RRMSVEC(56,130) = 4.8283D 00
      RRMSVEC(56,131) = 4.8276D 00
      RRMSVEC(56,132) = 4.8303D 00
      RRMSVEC(56,133) = 4.8286D 00
      RRMSVEC(56,134) = 4.8322D 00
      RRMSVEC(56,135) = 4.8294D 00
      RRMSVEC(56,136) = 4.8334D 00
      RRMSVEC(56,137) = 4.8314D 00
      RRMSVEC(56,138) = 4.8378D 00
      RRMSVEC(56,139) = 4.8513D 00
      RRMSVEC(56,140) = 4.8684D 00
      RRMSVEC(56,141) = 4.8807D 00
      RRMSVEC(56,142) = 4.8953D 00
      RRMSVEC(56,143) = 4.9087D 00
      RRMSVEC(56,144) = 4.9236D 00
      RRMSVEC(56,145) = 4.9345D 00
      RRMSVEC(56,146) = 4.9479D 00
      RRMSVEC(56,148) = 4.9731D 00
      RRMSVEC(57,135) = 4.8488D 00
      RRMSVEC(57,137) = 4.8496D 00
      RRMSVEC(57,138) = 4.8473D 00
      RRMSVEC(57,139) = 4.8550D 00
      RRMSVEC(58,136) = 4.8739D 00
      RRMSVEC(58,138) = 4.8737D 00
      RRMSVEC(58,140) = 4.8771D 00
      RRMSVEC(58,142) = 4.9063D 00
      RRMSVEC(58,144) = 4.9303D 00
      RRMSVEC(58,146) = 4.9590D 00
      RRMSVEC(58,148) = 4.9893D 00
      RRMSVEC(59,141) = 4.8919D 00
      RRMSVEC(60,132) = 4.9174D 00
      RRMSVEC(60,134) = 4.9128D 00
      RRMSVEC(60,135) = 4.9086D 00
      RRMSVEC(60,136) = 4.9111D 00
      RRMSVEC(60,137) = 4.9080D 00
      RRMSVEC(60,138) = 4.9123D 00
      RRMSVEC(60,139) = 4.9076D 00
      RRMSVEC(60,140) = 4.9101D 00
      RRMSVEC(60,141) = 4.9057D 00
      RRMSVEC(60,142) = 4.9123D 00
      RRMSVEC(60,143) = 4.9254D 00
      RRMSVEC(60,144) = 4.9421D 00
      RRMSVEC(60,145) = 4.9535D 00
      RRMSVEC(60,146) = 4.9696D 00
      RRMSVEC(60,148) = 4.9999D 00
      RRMSVEC(60,150) = 5.0400D 00
      RRMSVEC(62,138) = 4.9599D 00
      RRMSVEC(62,139) = 4.9556D 00
      RRMSVEC(62,140) = 4.9565D 00
      RRMSVEC(62,141) = 4.9517D 00
      RRMSVEC(62,142) = 4.9518D 00
      RRMSVEC(62,143) = 4.9479D 00
      RRMSVEC(62,144) = 4.9524D 00
      RRMSVEC(62,145) = 4.9651D 00
      RRMSVEC(62,146) = 4.9808D 00
      RRMSVEC(62,147) = 4.9892D 00
      RRMSVEC(62,148) = 5.0042D 00
      RRMSVEC(62,149) = 5.0134D 00
      RRMSVEC(62,150) = 5.0387D 00
      RRMSVEC(62,151) = 5.0550D 00
      RRMSVEC(62,152) = 5.0819D 00
      RRMSVEC(62,153) = 5.0925D 00
      RRMSVEC(62,154) = 5.1053D 00
      RRMSVEC(63,137) = 4.9762D 00
      RRMSVEC(63,138) = 4.9779D 00
      RRMSVEC(63,139) = 4.9760D 00
      RRMSVEC(63,140) = 4.9695D 00
      RRMSVEC(63,141) = 4.9697D 00
      RRMSVEC(63,142) = 4.9607D 00
      RRMSVEC(63,143) = 4.9636D 00
      RRMSVEC(63,144) = 4.9612D 00
      RRMSVEC(63,145) = 4.9663D 00
      RRMSVEC(63,146) = 4.9789D 00
      RRMSVEC(63,147) = 4.9938D 00
      RRMSVEC(63,148) = 5.0045D 00
      RRMSVEC(63,149) = 5.0202D 00
      RRMSVEC(63,150) = 5.0296D 00
      RRMSVEC(63,151) = 5.0522D 00
      RRMSVEC(63,152) = 5.1064D 00
      RRMSVEC(63,153) = 5.1115D 00
      RRMSVEC(63,154) = 5.1239D 00
      RRMSVEC(63,155) = 5.1221D 00
      RRMSVEC(63,156) = 5.1264D 00
      RRMSVEC(63,157) = 5.1351D 00
      RRMSVEC(63,158) = 5.1413D 00
      RRMSVEC(63,159) = 5.1498D 00
      RRMSVEC(64,145) = 4.9786D 00
      RRMSVEC(64,146) = 4.9801D 00
      RRMSVEC(64,148) = 5.0080D 00
      RRMSVEC(64,150) = 5.0342D 00
      RRMSVEC(64,152) = 5.0774D 00
      RRMSVEC(64,154) = 5.1223D 00
      RRMSVEC(64,155) = 5.1319D 00
      RRMSVEC(64,156) = 5.1420D 00
      RRMSVEC(64,157) = 5.1449D 00
      RRMSVEC(64,158) = 5.1569D 00
      RRMSVEC(64,160) = 5.1734D 00
      RRMSVEC(65,147) = 4.9201D 00
      RRMSVEC(65,148) = 4.9291D 00
      RRMSVEC(65,149) = 4.9427D 00
      RRMSVEC(65,150) = 4.9499D 00
      RRMSVEC(65,151) = 4.9630D 00
      RRMSVEC(65,152) = 4.9689D 00
      RRMSVEC(65,153) = 4.9950D 00
      RRMSVEC(65,154) = 5.0333D 00
      RRMSVEC(65,155) = 5.0391D 00
      RRMSVEC(65,157) = 5.0489D 00
      RRMSVEC(65,159) = 5.0600D 00
      RRMSVEC(66,146) = 5.0438D 00
      RRMSVEC(66,148) = 5.0455D 00
      RRMSVEC(66,149) = 5.0567D 00
      RRMSVEC(66,150) = 5.0706D 00
      RRMSVEC(66,151) = 5.0801D 00
      RRMSVEC(66,152) = 5.0950D 00
      RRMSVEC(66,153) = 5.1035D 00
      RRMSVEC(66,154) = 5.1241D 00
      RRMSVEC(66,155) = 5.1457D 00
      RRMSVEC(66,156) = 5.1622D 00
      RRMSVEC(66,157) = 5.1709D 00
      RRMSVEC(66,158) = 5.1815D 00
      RRMSVEC(66,159) = 5.1825D 00
      RRMSVEC(66,160) = 5.1951D 00
      RRMSVEC(66,161) = 5.1962D 00
      RRMSVEC(66,162) = 5.2074D 00
      RRMSVEC(66,163) = 5.2099D 00
      RRMSVEC(66,164) = 5.2218D 00
      RRMSVEC(67,151) = 5.0398D 00
      RRMSVEC(67,152) = 5.0614D 00
      RRMSVEC(67,153) = 5.0760D 00
      RRMSVEC(67,154) = 5.0856D 00
      RRMSVEC(67,155) = 5.1076D 00
      RRMSVEC(67,156) = 5.1156D 00
      RRMSVEC(67,157) = 5.1535D 00
      RRMSVEC(67,158) = 5.1571D 00
      RRMSVEC(67,159) = 5.1675D 00
      RRMSVEC(67,160) = 5.1662D 00
      RRMSVEC(67,161) = 5.1785D 00
      RRMSVEC(67,162) = 5.1817D 00
      RRMSVEC(67,163) = 5.1907D 00
      RRMSVEC(67,165) = 5.2022D 00
      RRMSVEC(68,150) = 5.0548D 00
      RRMSVEC(68,152) = 5.0843D 00
      RRMSVEC(68,154) = 5.1129D 00
      RRMSVEC(68,156) = 5.1429D 00
      RRMSVEC(68,158) = 5.1761D 00
      RRMSVEC(68,160) = 5.2045D 00
      RRMSVEC(68,162) = 5.2246D 00
      RRMSVEC(68,164) = 5.2389D 00
      RRMSVEC(68,166) = 5.2516D 00
      RRMSVEC(68,167) = 5.2560D 00
      RRMSVEC(68,168) = 5.2644D 00
      RRMSVEC(68,170) = 5.2789D 00
      RRMSVEC(69,153) = 5.0643D 00
      RRMSVEC(69,154) = 5.0755D 00
      RRMSVEC(69,156) = 5.0976D 00
      RRMSVEC(69,157) = 5.1140D 00
      RRMSVEC(69,158) = 5.1235D 00
      RRMSVEC(69,159) = 5.1392D 00
      RRMSVEC(69,160) = 5.1504D 00
      RRMSVEC(69,161) = 5.1616D 00
      RRMSVEC(69,162) = 5.1713D 00
      RRMSVEC(69,163) = 5.1849D 00
      RRMSVEC(69,164) = 5.1906D 00
      RRMSVEC(69,165) = 5.2004D 00
      RRMSVEC(69,166) = 5.2046D 00
      RRMSVEC(69,167) = 5.2129D 00
      RRMSVEC(69,168) = 5.2170D 00
      RRMSVEC(69,169) = 5.2256D 00
      RRMSVEC(69,170) = 5.2303D 00
      RRMSVEC(69,171) = 5.2388D 00
      RRMSVEC(69,172) = 5.2411D 00
      RRMSVEC(70,152) = 5.0423D 00
      RRMSVEC(70,154) = 5.0875D 00
      RRMSVEC(70,155) = 5.1040D 00
      RRMSVEC(70,156) = 5.1219D 00
      RRMSVEC(70,157) = 5.1324D 00
      RRMSVEC(70,158) = 5.1498D 00
      RRMSVEC(70,159) = 5.1629D 00
      RRMSVEC(70,160) = 5.1781D 00
      RRMSVEC(70,161) = 5.1889D 00
      RRMSVEC(70,162) = 5.2054D 00
      RRMSVEC(70,163) = 5.2157D 00
      RRMSVEC(70,164) = 5.2307D 00
      RRMSVEC(70,165) = 5.2399D 00
      RRMSVEC(70,166) = 5.2525D 00
      RRMSVEC(70,167) = 5.2621D 00
      RRMSVEC(70,168) = 5.2702D 00
      RRMSVEC(70,169) = 5.2771D 00
      RRMSVEC(70,170) = 5.2853D 00
      RRMSVEC(70,171) = 5.2906D 00
      RRMSVEC(70,172) = 5.2995D 00
      RRMSVEC(70,173) = 5.3046D 00
      RRMSVEC(70,174) = 5.3108D 00
      RRMSVEC(70,175) = 5.3135D 00
      RRMSVEC(70,176) = 5.3215D 00
      RRMSVEC(71,161) = 5.2293D 00
      RRMSVEC(71,162) = 5.2398D 00
      RRMSVEC(71,163) = 5.2567D 00
      RRMSVEC(71,164) = 5.2677D 00
      RRMSVEC(71,165) = 5.2830D 00
      RRMSVEC(71,166) = 5.2972D 00
      RRMSVEC(71,167) = 5.3108D 00
      RRMSVEC(71,168) = 5.3227D 00
      RRMSVEC(71,169) = 5.3290D 00
      RRMSVEC(71,170) = 5.3364D 00
      RRMSVEC(71,171) = 5.3436D 00
      RRMSVEC(71,172) = 5.3486D 00
      RRMSVEC(71,173) = 5.3577D 00
      RRMSVEC(71,174) = 5.3634D 00
      RRMSVEC(71,175) = 5.3700D 00
      RRMSVEC(71,176) = 5.3739D 00
      RRMSVEC(71,177) = 5.3815D 00
      RRMSVEC(71,178) = 5.3857D 00
      RRMSVEC(71,179) = 5.3917D 00
      RRMSVEC(72,170) = 5.2898D 00
      RRMSVEC(72,171) = 5.3041D 00
      RRMSVEC(72,172) = 5.3065D 00
      RRMSVEC(72,173) = 5.3140D 00
      RRMSVEC(72,174) = 5.3201D 00
      RRMSVEC(72,175) = 5.3191D 00
      RRMSVEC(72,176) = 5.3286D 00
      RRMSVEC(72,177) = 5.3309D 00
      RRMSVEC(72,178) = 5.3371D 00
      RRMSVEC(72,179) = 5.3408D 00
      RRMSVEC(72,180) = 5.3470D 00
      RRMSVEC(72,182) = 5.3516D 00
      RRMSVEC(73,181) = 5.3507D 00
      RRMSVEC(74,180) = 5.3491D 00
      RRMSVEC(74,182) = 5.3559D 00
      RRMSVEC(74,183) = 5.3611D 00
      RRMSVEC(74,184) = 5.3658D 00
      RRMSVEC(74,186) = 5.3743D 00
      RRMSVEC(75,185) = 5.3596D 00
      RRMSVEC(75,187) = 5.3698D 00
      RRMSVEC(76,184) = 5.3823D 00
      RRMSVEC(76,186) = 5.3909D 00
      RRMSVEC(76,187) = 5.3933D 00
      RRMSVEC(76,188) = 5.3993D 00
      RRMSVEC(76,189) = 5.4016D 00
      RRMSVEC(76,190) = 5.4062D 00
      RRMSVEC(76,192) = 5.4126D 00
      RRMSVEC(77,182) = 5.3705D 00
      RRMSVEC(77,183) = 5.3780D 00
      RRMSVEC(77,184) = 5.3805D 00
      RRMSVEC(77,185) = 5.3854D 00
      RRMSVEC(77,186) = 5.3900D 00
      RRMSVEC(77,187) = 5.3812D 00
      RRMSVEC(77,188) = 5.3838D 00
      RRMSVEC(77,189) = 5.3898D 00
      RRMSVEC(77,191) = 5.3968D 00
      RRMSVEC(77,193) = 5.4032D 00
      RRMSVEC(78,178) = 5.3728D 00
      RRMSVEC(78,179) = 5.3915D 00
      RRMSVEC(78,180) = 5.3891D 00
      RRMSVEC(78,181) = 5.3996D 00
      RRMSVEC(78,182) = 5.3969D 00
      RRMSVEC(78,183) = 5.4038D 00
      RRMSVEC(78,184) = 5.4015D 00
      RRMSVEC(78,185) = 5.4148D 00
      RRMSVEC(78,186) = 5.4037D 00
      RRMSVEC(78,187) = 5.4063D 00
      RRMSVEC(78,188) = 5.4053D 00
      RRMSVEC(78,189) = 5.4060D 00
      RRMSVEC(78,190) = 5.4108D 00
      RRMSVEC(78,191) = 5.4102D 00
      RRMSVEC(78,192) = 5.4169D 00
      RRMSVEC(78,193) = 5.4191D 00
      RRMSVEC(78,194) = 5.4236D 00
      RRMSVEC(78,195) = 5.4270D 00
      RRMSVEC(78,196) = 5.4307D 00
      RRMSVEC(78,198) = 5.4383D 00
      RRMSVEC(79,183) = 5.4247D 00
      RRMSVEC(79,184) = 5.4306D 00
      RRMSVEC(79,185) = 5.4296D 00
      RRMSVEC(79,186) = 5.4354D 00
      RRMSVEC(79,187) = 5.4018D 00
      RRMSVEC(79,188) = 5.4049D 00
      RRMSVEC(79,189) = 5.4084D 00
      RRMSVEC(79,190) = 5.4109D 00
      RRMSVEC(79,191) = 5.4147D 00
      RRMSVEC(79,192) = 5.4179D 00
      RRMSVEC(79,193) = 5.4221D 00
      RRMSVEC(79,194) = 5.4252D 00
      RRMSVEC(79,195) = 5.4298D 00
      RRMSVEC(79,196) = 5.4332D 00
      RRMSVEC(79,197) = 5.4371D 00
      RRMSVEC(79,198) = 5.4400D 00
      RRMSVEC(79,199) = 5.4454D 00
      RRMSVEC(80,181) = 5.4364D 00
      RRMSVEC(80,182) = 5.3833D 00
      RRMSVEC(80,183) = 5.4405D 00
      RRMSVEC(80,184) = 5.3949D 00
      RRMSVEC(80,185) = 5.4397D 00
      RRMSVEC(80,186) = 5.4017D 00
      RRMSVEC(80,187) = 5.4046D 00
      RRMSVEC(80,188) = 5.4085D 00
      RRMSVEC(80,189) = 5.4100D 00
      RRMSVEC(80,190) = 5.4158D 00
      RRMSVEC(80,191) = 5.4171D 00
      RRMSVEC(80,192) = 5.4232D 00
      RRMSVEC(80,193) = 5.4238D 00
      RRMSVEC(80,194) = 5.4309D 00
      RRMSVEC(80,195) = 5.4345D 00
      RRMSVEC(80,196) = 5.4385D 00
      RRMSVEC(80,197) = 5.4412D 00
      RRMSVEC(80,198) = 5.4463D 00
      RRMSVEC(80,199) = 5.4474D 00
      RRMSVEC(80,200) = 5.4551D 00
      RRMSVEC(80,201) = 5.4581D 00
      RRMSVEC(80,202) = 5.4648D 00
      RRMSVEC(80,203) = 5.4679D 00
      RRMSVEC(80,204) = 5.4744D 00
      RRMSVEC(80,205) = 5.4776D 00
      RRMSVEC(80,206) = 5.4837D 00
      RRMSVEC(81,188) = 5.4017D 00
      RRMSVEC(81,190) = 5.4121D 00
      RRMSVEC(81,191) = 5.4169D 00
      RRMSVEC(81,192) = 5.4191D 00
      RRMSVEC(81,193) = 5.4243D 00
      RRMSVEC(81,194) = 5.4259D 00
      RRMSVEC(81,195) = 5.4325D 00
      RRMSVEC(81,196) = 5.4327D 00
      RRMSVEC(81,197) = 5.4388D 00
      RRMSVEC(81,198) = 5.4396D 00
      RRMSVEC(81,199) = 5.4479D 00
      RRMSVEC(81,200) = 5.4491D 00
      RRMSVEC(81,201) = 5.4573D 00
      RRMSVEC(81,202) = 5.4595D 00
      RRMSVEC(81,203) = 5.4666D 00
      RRMSVEC(81,204) = 5.4704D 00
      RRMSVEC(81,205) = 5.4759D 00
      RRMSVEC(81,207) = 5.4853D 00
      RRMSVEC(81,208) = 5.4946D 00
      RRMSVEC(82,182) = 5.3788D 00
      RRMSVEC(82,183) = 5.3869D 00
      RRMSVEC(82,184) = 5.3930D 00
      RRMSVEC(82,185) = 5.3984D 00
      RRMSVEC(82,186) = 5.4027D 00
      RRMSVEC(82,187) = 5.4079D 00
      RRMSVEC(82,188) = 5.4139D 00
      RRMSVEC(82,189) = 5.4177D 00
      RRMSVEC(82,190) = 5.4222D 00
      RRMSVEC(82,191) = 5.4229D 00
      RRMSVEC(82,192) = 5.4300D 00
      RRMSVEC(82,193) = 5.4310D 00
      RRMSVEC(82,194) = 5.4372D 00
      RRMSVEC(82,195) = 5.4389D 00
      RRMSVEC(82,196) = 5.4444D 00
      RRMSVEC(82,197) = 5.4446D 00
      RRMSVEC(82,198) = 5.4524D 00
      RRMSVEC(82,199) = 5.4529D 00
      RRMSVEC(82,200) = 5.4611D 00
      RRMSVEC(82,201) = 5.4629D 00
      RRMSVEC(82,202) = 5.4705D 00
      RRMSVEC(82,203) = 5.4727D 00
      RRMSVEC(82,204) = 5.4803D 00
      RRMSVEC(82,205) = 5.4828D 00
      RRMSVEC(82,206) = 5.4902D 00
      RRMSVEC(82,207) = 5.4943D 00
      RRMSVEC(82,208) = 5.5012D 00
      RRMSVEC(82,209) = 5.5100D 00
      RRMSVEC(82,210) = 5.5208D 00
      RRMSVEC(82,211) = 5.5290D 00
      RRMSVEC(82,212) = 5.5396D 00
      RRMSVEC(82,214) = 5.5577D 00
      RRMSVEC(83,202) = 5.4840D 00
      RRMSVEC(83,203) = 5.4911D 00
      RRMSVEC(83,204) = 5.4934D 00
      RRMSVEC(83,205) = 5.5008D 00
      RRMSVEC(83,206) = 5.5034D 00
      RRMSVEC(83,207) = 5.5103D 00
      RRMSVEC(83,208) = 5.5147D 00
      RRMSVEC(83,209) = 5.5211D 00
      RRMSVEC(83,210) = 5.5300D 00
      RRMSVEC(83,212) = 5.5489D 00
      RRMSVEC(83,213) = 5.5586D 00
      RRMSVEC(84,192) = 5.5220D 00
      RRMSVEC(84,194) = 5.5167D 00
      RRMSVEC(84,196) = 5.5136D 00
      RRMSVEC(84,198) = 5.5146D 00
      RRMSVEC(84,200) = 5.5199D 00
      RRMSVEC(84,202) = 5.5281D 00
      RRMSVEC(84,204) = 5.5378D 00
      RRMSVEC(84,205) = 5.5389D 00
      RRMSVEC(84,206) = 5.5480D 00
      RRMSVEC(84,207) = 5.5501D 00
      RRMSVEC(84,208) = 5.5584D 00
      RRMSVEC(84,209) = 5.5628D 00
      RRMSVEC(84,210) = 5.5704D 00
      RRMSVEC(84,216) = 5.6359D 00
      RRMSVEC(84,218) = 5.6558D 00
      RRMSVEC(86,202) = 5.5521D 00
      RRMSVEC(86,204) = 5.5568D 00
      RRMSVEC(86,205) = 5.5569D 00
      RRMSVEC(86,206) = 5.5640D 00
      RRMSVEC(86,207) = 5.5652D 00
      RRMSVEC(86,208) = 5.5725D 00
      RRMSVEC(86,209) = 5.5743D 00
      RRMSVEC(86,210) = 5.5813D 00
      RRMSVEC(86,211) = 5.5850D 00
      RRMSVEC(86,212) = 5.5915D 00
      RRMSVEC(86,218) = 5.6540D 00
      RRMSVEC(86,219) = 5.6648D 00
      RRMSVEC(86,220) = 5.6731D 00
      RRMSVEC(86,221) = 5.6834D 00
      RRMSVEC(86,222) = 5.6915D 00
      RRMSVEC(87,207) = 5.5720D 00
      RRMSVEC(87,208) = 5.5729D 00
      RRMSVEC(87,209) = 5.5799D 00
      RRMSVEC(87,210) = 5.5818D 00
      RRMSVEC(87,211) = 5.5882D 00
      RRMSVEC(87,212) = 5.5915D 00
      RRMSVEC(87,213) = 5.5977D 00
      RRMSVEC(87,220) = 5.6688D 00
      RRMSVEC(87,221) = 5.6790D 00
      RRMSVEC(87,222) = 5.6890D 00
      RRMSVEC(87,223) = 5.6951D 00
      RRMSVEC(87,224) = 5.7061D 00
      RRMSVEC(87,225) = 5.7112D 00
      RRMSVEC(87,226) = 5.7190D 00
      RRMSVEC(87,227) = 5.7335D 00
      RRMSVEC(87,228) = 5.7399D 00
      RRMSVEC(88,208) = 5.5850D 00
      RRMSVEC(88,209) = 5.5853D 00
      RRMSVEC(88,210) = 5.5917D 00
      RRMSVEC(88,211) = 5.5929D 00
      RRMSVEC(88,212) = 5.5991D 00
      RRMSVEC(88,213) = 5.6020D 00
      RRMSVEC(88,214) = 5.6079D 00
      RRMSVEC(88,220) = 5.6683D 00
      RRMSVEC(88,221) = 5.6795D 00
      RRMSVEC(88,222) = 5.6874D 00
      RRMSVEC(88,223) = 5.6973D 00
      RRMSVEC(88,224) = 5.7046D 00
      RRMSVEC(88,225) = 5.7150D 00
      RRMSVEC(88,226) = 5.7211D 00
      RRMSVEC(88,227) = 5.7283D 00
      RRMSVEC(88,228) = 5.7370D 00
      RRMSVEC(88,229) = 5.7455D 00
      RRMSVEC(88,230) = 5.7551D 00
      RRMSVEC(88,232) = 5.7714D 00
      RRMSVEC(90,227) = 5.7404D 00
      RRMSVEC(90,228) = 5.7488D 00
      RRMSVEC(90,229) = 5.7557D 00
      RRMSVEC(90,230) = 5.7670D 00
      RRMSVEC(90,232) = 5.7848D 00
      RRMSVEC(92,233) = 5.8203D 00
      RRMSVEC(92,234) = 5.8291D 00
      RRMSVEC(92,235) = 5.8337D 00
      RRMSVEC(92,236) = 5.8431D 00
      RRMSVEC(92,238) = 5.8571D 00
      RRMSVEC(94,238) = 5.8535D 00
      RRMSVEC(94,239) = 5.8601D 00
      RRMSVEC(94,240) = 5.8701D 00
      RRMSVEC(94,241) = 5.8748D 00
      RRMSVEC(94,242) = 5.8823D 00
      RRMSVEC(94,244) = 5.8948D 00
      RRMSVEC(95,241) = 5.8928D 00
      RRMSVEC(95,243) = 5.9048D 00
      RRMSVEC(96,242) = 5.8285D 00
      RRMSVEC(96,244) = 5.8429D 00
      RRMSVEC(96,245) = 5.8475D 00
      RRMSVEC(96,246) = 5.8562D 00
      RRMSVEC(96,248) = 5.8687D 00
      
      write(*,*)
      write(*,*) 'RNUCLEUS'
      write(*,*) 'This program defines nuclear data and the radial grid'
      write(*,*) 'Outputfile: isodata'
      write(*,*)
*
*   File  grasp92.iso  is FORMATTED
*
      DEFNAM = 'isodata'
      FORM = 'FORMATTED'
      STATUS = 'NEW'
*
      FILNAM = DEFNAM
*
      CALL OPENFL (22,FILNAM,FORM,STATUS,IERR)
*
      IF (IERR .NE. 0) THEN
         WRITE(istde,*) 'Error when opening isodata'
         STOP
      ENDIF
*
      WRITE(istde,*) 'Enter the atomic number:'

      READ *, Z
      WRITE (22,300) 'Atomic number:'
      WRITE (22,*) Z
*
      WRITE(istde,*) 'Enter the mass number (0 if the'
     &, ' nucleus is to be modelled as a point source:'

      READ *, A

      WRITE (22,300) 'Mass number (integer) :'
      WRITE (22,*) A

      IF (A .EQ. 0.0D 00) THEN

         CPARM = 0.0D 00
         APARM = 0.0D 00
      ELSE
         write(*,*)
         IF (RRMSVEC(int(Z),int(A)) .GT. 0.0D 00) THEN
            RRMS = RRMSVEC(int(Z),int(A))
            WRITE(istde,*) 'Source of rms radius:'
     &           , ' I. Angeli and K.P. Marinova (2013)'       
            WRITE(istde,*) 'The default root mean squared'
     &           , ' radius is ',RRMS,' fm;'
         ELSE
            RRMS =  0.836D 00*A**(1.0D 00/3.0D 00)
     :           +0.570D 00
            WRITE(istde,*) 'Computed rms radius:'
     &           , ' 0.836*A^(1/3) + 0.570'                   
            WRITE(istde,*) 'The default root mean squared'
     &           , ' radius is ',RRMS,' fm;'
         END IF
         TPARM = 2.30D 00
         WRITE(istde,*) 'The default nuclear skin thickness'
     &, ' is   ',TPARM,' fm;'
         WRITE(istde,*) 'Revise these values?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE(istde,*) 'Enter the root mean squared'
     &, ' radius of the nucleus (in fm):'
            READ *, RRMS
            WRITE(istde,*) 'Enter the skin thickness of'
     &, ' the nucleus (in fm):'
            READ *, TPARM
         ENDIF
         APARM = TPARM/(4.0D 00*LOG (3.0D 00))
         CALL GETCPR (RRMS,APARM,CPARM)
      ENDIF
      WRITE (22,300) 'Fermi distribution parameter a:'
      WRITE (22,*) APARM
      WRITE (22,300) 'Fermi distribution parameter c:'
      WRITE (22,*) CPARM
*
      WRITE(istde,*) 'Enter the mass of the neutral'
     &, ' atom (in amu) (0 if the nucleus is to be static):'
      READ *, AMAMU
      IF (AMAMU .NE. 0.0D 00) THEN
          EBIND = 0.D0
         NENEU = NINT (Z)
         IF (EBIND .GT. 0.0D 00) EBIND = -EBIND
         EMNAMU = AMAMU-EMEAMU*DBLE (NENEU)
     :                 -EMEAMU*EBIND/ALFAI**2
      ELSE
         EMNAMU = 0.0D 00
      ENDIF
*
      WRITE (22,300) 'Mass of nucleus (in amu):'
      WRITE (22,*) EMNAMU
*
      WRITE(istde,*) 'Enter the nuclear spin quantum'
     &, ' number (I) (in units of h / 2 pi):'
      READ *, SQN
      WRITE (22,300) 'Nuclear spin (I) (in units of h / 2 pi):'
      WRITE (22,*) SQN
*
      WRITE(istde,*) 'Enter the nuclear dipole moment'
     &, ' (in nuclear magnetons):'
      READ *, DMOMNM
      WRITE (22,300) 'Nuclear dipole moment (in nuclear magnetons):'
      WRITE (22,*) DMOMNM
*
      WRITE(istde,*) 'Enter the nuclear quadrupole'
     &, ' moment (in barns):'
      READ *, QMOMB
      WRITE (22,300) 'Nuclear quadrupole moment (in barns):'
      WRITE (22,*) QMOMB
*     
*     Grid Parameters - Jon Grumer (Lund, 2013)
*      
*     Default values
*
*     If chosen to model nucleus as point source then
      IF (A .EQ. 0) THEN 
         WRITE(istde,*) 'You have chosen to model the nucleus as '
     &,                 'a point source!'
         RNT = EXP (-65.0D 00/16.0D 00) / real(Z)
         H   = 0.5D 00**4
         HP  = 0.d0
         N   = MIN (220,NNNP)
*     Otherwize set normal default values         
      ELSE
         RNT = 2.d-6
         H   = 5.d-2
         HP  = 0.d0
         N   = INT(NNNP) ! default number of grid points from parameters.def
      END IF
*
*     Print default values to screen
*
  101 WRITE(istde,*) '------------------------------------------------'
     &,              '-----------'
      WRITE(istde,*) 'The Grasp2K grid is defined by:'
      WRITE(istde,*) 
      WRITE(istde,*) 'R(i) = RNT*[ exp[ (i-1)*H ] - 1 ]'
      WRITE(istde,*) '  i  = 1, 2, ..., NNNP'
      WRITE(istde,*) 
      WRITE(istde,*) 'The default grid parameters are:'
      WRITE(istde,*) 
      WRITE(istde,*) 'RNT  (first grid point       ) = ', RNT
      WRITE(istde,*) 'H    (grid step-size         ) = ', H
!      WRITE(istde,*) 'HP   (related to linear grid ) = ', HP
      WRITE(istde,*) 'NNNP (max. no. of grid-points) = ', N
      WRITE(istde,*) '------------------------------------------------'
     &,              '-----------'
      WRITE(istde,*) 
      WRITE(istde,*) 'Do you want to revise these values (y/*)?'
*
      READ (*,'(A)') GRIDANSW
*
      IF (GRIDANSW.EQ.'y'.OR.GRIDANSW.EQ.'Y') THEN
*
         WRITE(istde,*)
         WRITE(istde,*) 'NNNP depends on H. A new value of'
         WRITE(istde,*) 'NNNP should satisfy the inequality       '
         WRITE(istde,*) 'NNNP_new >= (0.05/H_new)*590 '
         WRITE(istde,*)


         WRITE(istde,*) 'Do you want to revise RNT (y/*)?'
         READ (*,'(A)'), GRIDANSW2
         IF (GRIDANSW2.EQ.'y'.OR.GRIDANSW2.EQ.'Y') THEN
            WRITE(istde,*) 'Enter RNT (double precision real):'
            READ *, RNT
         ELSE
            WRITE(istde,*) 'Default value of RNT is set'
         END IF
*      
         WRITE(istde,*) 'Do you want to revise H (y/*)?'
         READ (*,'(A)'), GRIDANSW2
         IF (GRIDANSW2.EQ.'y'.OR.GRIDANSW2.EQ.'Y') THEN
            WRITE(istde,*) 'Enter H (double precision real):'
            READ *, H
         ELSE
            WRITE(istde,*) 'Default value of H is set'
         END IF
*      
         WRITE(istde,*) 'Do you want to revise NNNP (y/*)?'
         READ (*,'(A)'), GRIDANSW2
         IF (GRIDANSW2.EQ.'y'.OR.GRIDANSW2.EQ.'Y') THEN
            WRITE(istde,'(a,i6,a)') 'Enter NNNP (integer <=', NNNP,' ):'
            READ *, N
            IF (N.le.(0.05/H)*590) THEN
               WRITE(*,*) 'Relation between NNNP and H not fullfilled'
               STOP
            END IF
         ELSE
            WRITE(istde,*) 'Default value of NNNP is set'
         END IF
*
      ELSE 
         WRITE(istde,*) 'Default grid parameter values are kept!' 
      END IF
*
*     Write grid parameters to isodata
*
      WRITE (22,300) 'RNT:'
      WRITE (22,*) RNT
*
      WRITE (22,300) 'H:'
      WRITE (22,*) H
*
      WRITE (22,300) 'HP:'
      WRITE (22,*) HP
*
      WRITE (22,300) 'NNNP:'
      WRITE (22,*) N
*
      CLOSE (22)
*
      WRITE(istde,*)
      WRITE(istde,*) 'Output file isodata successfully created. '
     &,              'Happy computing!'
      WRITE(istde,*)

      STOP
*
  300 FORMAT (A)
*
      END
