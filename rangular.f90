 
!***********************************************************************
      PROGRAM GENMCP 
!***********************************************************************
!   M o d u l e s 
!-----------------------------------------------
      use const
      use atoms_par
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER NBLK0, NCOUNT1, NCORE, NB
      PARAMETER (NBLK0 = 50) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL ::  DEBUG, RESTRT, YES
      CHARACTER(LEN=128) :: STARTDIR*128, PERMDIR*128, TMPDIR*128
      CHARACTER(LEN=8), DIMENSION(NBLK0):: IDBLK
      INTEGER :: MYID,NPROCS
!-----------------------------------------------------------------------
      write(*,*)
      write(*,*) 'RANGULAR'
      write(*,*) 'This program performs angular integration '
      write(*,*) 'Input file:  rcsf.inp'
      write(*,*) 'Outputfiles: mcp.30, mcp.31, ....'
      write(*,*) 'rangular.log                     '
      write(*,*)
      OPEN(UNIT=739,FILE='rangular.log',STATUS='UNKNOWN')
!=======================================================================
!  Get NDEF
!=======================================================================
      IF (MYID == 0) THEN 
         WRITE (istde,'(A)') ' Full interaction?  (y/n) '
         YES = GETYN() 
         IF (YES) THEN 
            NDEF = 0 
            write(739,'(A)') 'y            ! Full interaction'
         ELSE 
            NDEF = 1 
            write(739,'(A)') 'n            ! Full interaction'
         ENDIF 
      ENDIF 
!=======================================================================
!  Checks and settings... Mostly done in backyard.
!
!    setsum - open the summary file
!    cslh - load header of the csl file
!    setmcp - open and check the  .mcp  files
!    strsum - append a summary of the inputs to the  .sum  file
!    factt - table of logarithms of factorials setup
!=======================================================================
      CALL SETMC 
      CALL CSLH ('rcsf.inp', NCORE, NBLK0, IDBLK) 
      RESTRT = .FALSE. 
      CALL SETMCP2 (MYID, NPROCS, NCORE, IDBLK, 'mcp') 
      IF (NDEF/=0 .AND. MYID==0) CALL STRSUM 
      CALL FACTT 
      CALL genintrk(
!=======================================================================
!     For each block, generate and sort the data
!=======================================================================

      DO NB = 1, NBLOCK 
         NCF = NCFBLK(NB)                        ! This ncf goes to common 
         IF (MYID == 0) THEN 
            WRITE (6, *) 
            WRITE (6, *) 'Block ', NB, ',  ncf = ', NCF 
         ENDIF 
         !*** Load current CSL block. Memories de-allocated in mcp ***
         CALL LODCSH2 (21, NCORE, NB) 
         !*** Open tmp.xx files for block nb ***
         CALL SETTMPGG (nb, 30, 'tmp')
         !*** Generation of MCP coefficients ***
         CALL MCP (NB, RESTRT, MYID, NPROCS, 'mcp') 
      END DO 
      CLOSE(24)                                  ! Summary file 
      CLOSE(739)                                  ! rangular.log
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================
      STOP  
      END PROGRAM GENMCP
