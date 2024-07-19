      PROGRAM NBODY6
*
*             N B O D Y 6++
*             *************
*
*       Regularized AC N-body code with triple & binary collisions.
*       --------------------------------------------------------
*
*       Hermite integration scheme with block-steps (V 4.0.0 April/99).
*       ------------------------------------------------------------------
*
*       Developed by Sverre Aarseth, IOA, Cambridge.
*       ............................................
*       Message Passing Version NBODY6++ for Massively Parallel Systems
*       Developed by Rainer Spurzem, ARI, Heidelberg
*       
*       Hybrid parallelization (GPU, AVX/SSE, OpenMP + MPI) 
*       Developed by Long Wang, KIAA, Peking University
*
      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
      include 'omp_lib.h'
      COMMON/STSTAT/  TINIT,NIR,NIB,NRGL,NKS
#ifdef DEBUG
*     --10/03/14 19:40-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
      COMMON/adt/ adtime,dumptime,dprintt,dtprint,namep
*     --10/03/14 19:41-lwang-end----------------------------------------*
#endif
      EXTERNAL MERGE
      NAMELIST /INNBODY6/ KSTART, TCOMP, TCRTP0,
     &    isernb,iserreg,iserks,adtime,dumptime,dprintt,dtprint,namep
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#define MPIINIT 0
#endif

#if MPIINIT
      character (len=128) env_var
      integer tmp_col, istatus
      !double precision send_dt
      integer cbuf(1), repetitions, smax_0
      double precision next_time

      repetitions = 0

*       Initialize MPI
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

*     Opening on unit 5 the file passed as 1st argument if there is one
      if (iargc() >= 1) then
         call getarg(1, env_var)
         env_var = ADJUSTL(TRIM(env_var))
         INPUT_UNIT=99
         open(unit=99, file=env_var, action='READ')
         write(6,*) 'Reading input from ', env_var
      else
         INPUT_UNIT=5
      end if
      

*     We check if we are in MIMD mode. If so, we must define a
*     environment variable called NB6_PROCS
      !OPEN(999, file='RBD_debug.dat')
      CALL getenv('NB6_PROCS', env_var)
      if (env_var .ne. '') then
*     We switch to MIMD mode
         READ (env_var, '(I6)') isize                                ! Replacing size

         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         call MPI_COMM_SPLIT(MPI_COMM_WORLD, 0, rank, MPI_COMM_NB6
     &        , ierr)

         if (rank .eq. 0) then ! If master process we put the process on the rambody communicator
            tmp_col = 0
         else
            tmp_col = MPI_UNDEFINED
         end if

         call MPI_COMM_SPLIT(MPI_COMM_WORLD, tmp_col, rank, MPI_COMM_RBD
     &        , ierr)
  
         call MPI_COMM_RANK(MPI_COMM_NB6, rank, ierr) ! Replacing rankrr

*     Now communicating greetings with Ramses
         write(6,*) 'NB6 : my rank = ', rank, ' on process ', getpid()
         call flush(6)
         if (rank .eq. 0) then
            write(6,*) 'NB6 : Waiting for greetings from Ramses'
            call flush(6)
            write(6,*) 'RBDCOMM 1'
            call flush(6)
            call MPI_Recv(cbuf, 1, MPI_INTEGER, 1, 1, MPI_COMM_RBD,
     &           MPI_STATUS_IGNORE, ierr)

            if (cbuf(1) .ne. 94111535) then
               write(6,*) 'NB6, ERROR : Greetings not matching Ramses'
               call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
            end if

            write(6,*) 'NB6 : Received = ', cbuf(1)
            call flush(6)

            cbuf(1) = 1180276
            write(6,*) 'RBDCOMM 2'
            call flush(6)
            call MPI_Send(cbuf, 1, MPI_INTEGER, 1, 2, MPI_COMM_RBD,
     &           ierr)


            write(6,*) 'NB6 : Greetings exchanged with Ramses !'
         end if

         write(6,*) ' NB6 : Sync point'
         call flush(6)
         call MPI_Barrier(MPI_COMM_WORLD, ierr)

         MIMD_MODE = .true.
         RBD_ON = .false.
         PREV_T = -1.0

         BAR_SYNC_COUNT=0
      else
*     If we are not in MIMD mode, we just copy the MPI_COMM_NB6 communicator
*     so that it matches MPI_COMM_WORLD
         MPI_COMM_NB6 = MPI_COMM_WORLD
         MIMD_MODE = .false.
      endif
      

      CALL MPI_COMM_GROUP(MPI_COMM_NB6,group,ierr)
      CALL MPI_GROUP_SIZE(group,isize,ierr)
      CALL MPI_GROUP_RANK(group,rank,ierr)
      ibarcount=0
      write(6,11) rank,isize,group
 11   format('MPI-initial: This is rank=',I6,' size=',I6,' group=',I11)
#endif
*
*       Initialize the timer.
      CALL CPUTIM(ttota)

*       Get threads number
#ifdef OMP
!$omp parallel 
      icore=OMP_get_num_threads()
!$omp end parallel
      PRINT*,'RANK: ',rank,' OpenMP Number of Threads: ',icore
#else
      icore = 1
#endif

#ifdef PARALLEL
      call mpi_barrier(MPI_COMM_NB6,ierr)
#endif
      call flush(6)
*
*       Read start/restart indicator & CPU time.
*     19:41-lwang-debug-10/03/14-----------------------------*
      IF(rank.eq.0)THEN
#ifdef DEBUG
      READ (NML=INNBODY6, IOSTAT=IIC, UNIT=5)
      IF (IIC.NE.0) WRITE (6, '("NBODY6D: invalid Namelist format")')
*     READ (5,*)  KSTART, TCOMP, TCRTP0,
*    &    isernb,iserreg,iserks,adtime,dumptime,dprintt,dtprint,namep
#else
      READ (NML=INNBODY6, IOSTAT=IIC, UNIT=5)
      IF (IIC.NE.0) WRITE (6, '("NBODY6: invalid Namelist format")')
*     READ (5,*)  KSTART, TCOMP, TCRTP0,
*    &    isernb,iserreg,iserks
#endif
      END IF
#if MPIINIT
      CALL MPI_BCAST(isernb,1,MPI_INTEGER,0,MPI_COMM_NB6,ierr)
      CALL MPI_BCAST(iserreg,1,MPI_INTEGER,0,MPI_COMM_NB6,ierr)
      CALL MPI_BCAST(iserks,1,MPI_INTEGER,0,MPI_COMM_NB6,ierr)
      CALL MPI_BCAST(KSTART,1,MPI_INTEGER,0,MPI_COMM_NB6,ierr)
      CALL MPI_BCAST(TCOMP,1,MPI_REAL8,0,MPI_COMM_NB6,ierr)
      CALL MPI_BCAST(TCRTP0,1,MPI_REAL8,0,MPI_COMM_NB6,ierr)

#ifdef DEBUG
      CALL MPI_BCAST(adtime,1,MPI_REAL8,0,MPI_COMM_NB6,ierr)
      CALL MPI_BCAST(dumptime,1,MPI_REAL8,0,MPI_COMM_NB6,ierr)
      CALL MPI_BCAST(dprintt,1,MPI_REAL8,0,MPI_COMM_NB6,ierr)
      CALL MPI_BCAST(dtprint,1,MPI_REAL8,0,MPI_COMM_NB6,ierr)
      CALL MPI_BCAST(namep,1,MPI_INTEGER,0,MPI_COMM_NB6,ierr)
#endif
*
      isernb = max(isize,isernb*icore)
      iserreg = max(isize,iserreg*icore)
*      iserks = max(isize,iserks*icore)

      IF(rank.eq.0) THEN
        PRINT*,' MPI: iserreg,isernb,iserks,icore=',
     &        iserreg,isernb,iserks,icore
#ifdef DEBUG
        PRINT*,' DEBUG: adtime,dumptime,dprintt,dtprint,namep =',
     &      adtime,dumptime,dprintt,dtprint,namep
#endif
      END IF
#endif
*
      IF (KSTART.EQ.1) THEN
*
          TCRITP = TCRTP0
          if(rank.eq.0)WRITE (6,*) ' KSTART, TCOMP, TCRITP =',
     &                               KSTART, TCOMP, TCRITP
*       Read input parameters, perform initial setup and obtain output.
          CPU = TCOMP
          CALL START
          call cputim(tt7)
          CALL ADJUST
          call cputim(tt8)
          if(rank.eq.0)ttadj = ttadj + (tt8-tt7)*60.
      ELSE
*       Read previously saved COMMON variables from tape/disc on unit 1.
*       Backup kstart value before call mydump
          KSTART0 = KSTART
          CALL MYDUMP(0,1)
*       Reset kstart to input value
          KSTART = KSTART0
*       
          IF (NDUMP.GE.3) STOP
*       Safety indicator preventing repeated restarts set in routine CHECK.
          CPU = TCOMP
          CPU0 = 0.0 
*       Set IPHASE = -1 for new NLIST in routine INTGRT (Hermite version).
          IPHASE = -1
*
*       Initialize evolution parameters which depend on metallicity.
          IF (KZ(19).GE.3) THEN
              CALL ZCNSTS(ZMET,ZPARS)
          END IF
*
*       Check reading modified restart parameters (KSTART = 3, 4 or 5).
*       Or control output RSp March 2019 (KSTART = 2)
          CALL MODIFY
*
*       Open all other files.
          if(rank.eq.0)CALL FILE_INIT
*
*       Find massive back hole after restart
          IF (KZ(24).EQ.1) call IMBHRESTART

          IF (MIMD_MODE) call rbd_restart

      END IF
*
* (R.Sp.)Set time flag and step number flags for beginning of run
      TINIT = TTOT
      NIR = NSTEPI
      NIB = NSTEPB
      NRGL = NSTEPR
      NKS = NSTEPU
*
      call cputim(tt2)
      if(rank.eq.0)ttinitial = ttinitial + (tt2-ttota)*60.

*     Rambody, setting smax_0 to the initial value of SMAX
      smax_0 = SMAX
      next_time = -1.0
      
*       Advance solutions until next output or change of procedure.
    1 CONTINUE
      
*     Synchronising Rambody
      if (MIMD_MODE) then
         if (TTOT * TSTAR >= next_time) then
            SMAX = smax_0

            !write(6,*) 'Syncing escapers'
            !call rbd_sync_escapers

            if (rank .eq. 0) then
               write(6,*) 'Sending cluster'
               call flush(6)
               call rbd_send_cluster

               write(6,*) 'Sending cluster scale'
               call flush(6)
               call rbd_send_rscale
            
               write(6,*) 'Sending timestep'
               call flush(6)
               call rbd_sync_timestep
            
               write(6,*) 'Waiting for timestep information (k)'
               call flush(6)

               write(6,*) 'RBDCOMM 11'
               call flush(6)
               call MPI_RECV(cbuf, 1, MPI_INTEGER, 1, 11,
     &              MPI_COMM_RBD, MPI_STATUS_IGNORE, ierr)

               repetitions = cbuf(1)
               write(6,*) 'Received from Rambody :'
               write(6,*) ' - Current time : ', TTOT * TSTAR, 'Myrs'
               next_time = (TTOT + repetitions*SMAX) * TSTAR
               write(6,*) ' - Next Sync time : ', next_time, 'Myrs'
               write(6,*) ' - SMAX repetitions : ', repetitions
               write(6,*) ''
               call flush(6)
            end if

            ! Syncing with the other processes
            write(6,*) 'Syncing reps with other processes'
            call MPI_Bcast(next_time, 1, MPI_DOUBLE_PRECISION, 0, 
     &           MPI_COMM_NB6, ierr)

            write(6,*) 'Receiving mesh'
            call flush(6)
            call rbd_sync_mesh
         else
            !repetitions = repetitions - 1
         end if
      end if
      
      call cputim(tt1)
*
*     if(rank.eq.0.and.time.lt.0.01)print*,' Bef INTGRT t=',time
      CALL INTGRT
*     if(rank.eq.0.and.time.lt.0.01)print*,' Aft INTGRT t=',time
*
      call cputim(tt2)
      if(rank.eq.0)ttint = ttint + (tt2-tt1)*60.
*
      IF (IPHASE.EQ.1) THEN
*       Prepare new KS regularization.
      call cputim(tt1)
          CALL KSREG
          CALL FLUSH(6)
      call cputim(tt2)
      if(rank.eq.0)ttks = ttks + (tt2-tt1)*60.
      if(rank.eq.0)ttksinit = ttksinit + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.2) THEN
*       Terminate KS regularization.
      call cputim(tt1)
          CALL KSTERM
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
      ttksterm = ttksterm + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.3) THEN
*       Perform energy check & parameter adjustments and print diagnostics.
          call cputim(tt7)
          CALL ADJUST
          call cputim(tt8)
          ttadj = ttadj + (tt8-tt7)*60.
*
      ELSE IF (IPHASE.EQ.4) THEN
*       Switch to unperturbed three-body regularization.
      call cputim(tt1)
          ISUB = 0 
          CALL TRIPLE(ISUB)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.5) THEN
*       Switch to unperturbed four-body regularization.
      call cputim(tt1)
          ISUB = 0
          CALL QUAD(ISUB)
      call cputim(tt2)
      if(rank.eq.0)ttks = ttks + (tt2-tt1)*60.
*
*       Adopt c.m. approximation for inner binary in hierarchical triple.
      ELSE IF (IPHASE.EQ.6) THEN
      call cputim(tt1)
          CALL MERGE
      call cputim(tt2)
      if(rank.eq.0)ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.7) THEN
*       Restore old binary in hierarchical configuration.
      call cputim(tt1)
          CALL RESET
      call cputim(tt2)
      if(rank.eq.0)ttks = ttks + (tt2-tt1)*60.
*
*       Begin chain regularization.
      ELSE IF (IPHASE.EQ.8) THEN
      call cputim(tt1)
          ISUB = 0
          TIMENB = TTOT
          CALL CHAIN(ISUB,TIMENB)
      call cputim(tt2)
      if(rank.eq.0)ttks = ttks + (tt2-tt1)*60.
      END IF
*
*       Continue integration.
      GO TO 1
*
      END
