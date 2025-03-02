      SUBROUTINE CHINIT(ISUB)
*
*
*       Initialization of chain system.
*       -------------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX3=3*NMX,NMX4=4*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK,ANG(3),FIRR(3),FD(3)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/INCOND/  X4(3,NMX),XDOT4(3,NMX)
      COMMON/ECHAIN/  ECH
      COMMON/SLOW3/   GCRIT,KZ26
*
*
*       Define chain membership.
      call xbpredall
*     --01/03/14 15:15-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      do j=ifirst,ntot
c$$$         write(123+rank,*),j,'n',name(j),'x0',x0(1,j),
c$$$     &        'x',x(1,j),'x0dot',x0dot(1,j),'xdot',xdot(1,j),
c$$$     *        'f',f(1,j),'fdot',fdot(1,j),'t0',t0(j),
c$$$     *        'body',body(j),'time',time
c$$$         call flush(123+rank)
c$$$      end do
*     --01/03/14 15:15-lwang-end----------------------------------------*
      CALL SETSYS
*     --01/03/14 13:06-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      kk= 15967
c$$$      print*,rank,'SETSYS KK',kk,'n',name(kk),'x',x(1,kk),
c$$$     &     'xdot',xdot(1,kk),'fdot',fdot(1,kk),'time',time
c$$$      call flush(6)
c$$$      call mpi_barrier(MPI_COMM_NB6,ierr)
*     --01/03/14 13:06-lwang-end----------------------------------------*
*
*       Initialize c.m. variables.
      DO 2 K = 1,7
          CM(K) = 0.0D0
    2 CONTINUE
*
*       Transform to the local c.m. reference frame.
      DO 4 L = 1,NCH
          J = JLIST(L)
          SIZE(L) = RADIUS(J)
          ISTAR(L) = KSTAR(J)
*       Place the system in first single particle locations.
          CM(7) = CM(7) + M(L)
          DO 3 K = 1,3
              X4(K,L) = X(K,J)
              XDOT4(K,L) = XDOT(K,J)
              CM(K) = CM(K) + M(L)*X4(K,L)
              CM(K+3) = CM(K+3) + M(L)*XDOT4(K,L)
    3     CONTINUE
    4 CONTINUE
*
*       Set c.m. coordinates & velocities of subsystem.
      DO 5 K = 1,6
          CM(K) = CM(K)/CM(7)
    5 CONTINUE
*
*       Specify initial conditions for chain regularization.
      LK = 0
      DO 8 L = 1,NCH
          DO 7 K = 1,3
              LK = LK + 1
              X4(K,L) = X4(K,L) - CM(K)
              XDOT4(K,L) = XDOT4(K,L) - CM(K+3)
              XCH(LK) = X4(K,L)
              VCH(LK) = XDOT4(K,L)
    7     CONTINUE
    8 CONTINUE
*
*       Calculate internal energy and and save in chain energy.
      CALL CONST(XCH,VCH,M,NCH,ENERGY,ANG,GAM)
      ECH = ENERGY
*
*       Set sum of mass products and save separations & RINV for CHLIST.
      SUM = 0.0D0
      RSUM = 0.0D0
      DO 10 L = 1,NCH-1
          DO 9 K = L+1,NCH
              SUM = SUM + M(L)*M(K)
              RLK2 = (X4(1,L) - X4(1,K))**2 + (X4(2,L) - X4(2,K))**2 +
     &                                        (X4(3,L) - X4(3,K))**2
              RSUM = RSUM + SQRT(RLK2)
              RINV(L) = 1.0/SQRT(RLK2)
    9     CONTINUE
   10 CONTINUE
*
*       Reduce RSUM by geometrical factor and check upper limit from IMPACT.
      IF (NCH.EQ.4) RSUM = 0.5*RSUM
      RSUM = MIN(FLOAT(NCH-1)*RSUM/FLOAT(NCH),RMIN)
*
*       Define gravitational radius for initial perturber list.
      RGRAV = SUM/ABS(ENERGY)
*
*       Avoid small value after collision (CHTERM improves perturbers).
      IF (NCH.GT.2) THEN
          RGRAV = MIN(RGRAV,0.5*RSUM)
      END IF
*
*       Set global index of c.m. body and save name (SUBSYS sets NAME = 0).
      IF (TIMEC.GT.0.0D0) ICH0 = ICH
      ICH = JLIST(1)
      NAME0 = NAME(ICH)
*
*       Define subsystem indicator (ISYS = 1, 2, 3 for triple, quad, chain).
      ISYS(NSUB+1) = 3
*
*       Form ghosts and initialize c.m. motion in ICOMP (= JLIST(1)).
      CALL SUBSYS(NCH,CM)
*     --01/03/14 13:06-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      print*,rank,'SUBSYS ICH',ich,'n',name(ich),'x',x(1,ich),
c$$$     &     'xdot',xdot(1,ich),'fdot',fdot(1,ich),'time',time
c$$$      call flush(6)
c$$$      call mpi_barrier(MPI_COMM_NB6,ierr)
*     --01/03/14 13:06-lwang-end----------------------------------------*
*
*       Copy neighbour list for ghost removal.
      NNB = LIST(1,ICH)
      DO 20 L = 2,NNB+1
          JPERT(L-1) = LIST(L,ICH)
   20 CONTINUE
*
*       Check possible switch of reference body on second call from CHAIN.
      IF (TIMEC.GT.0.0D0.AND.ICH.NE.ICH0) THEN
*       Add #ICH to neighbour & perturber lists before removing all ghosts.
          CALL NBREST(ICH0,1,NNB)  
      END IF
*
*       Remove ghosts (saved in JLIST) from neighbour lists of #ICH.
      CALL NBREM(ICH,NCH,NNB)
*
*       Remove ghosts from list of ICOMP (use NTOT as dummy here).
      JPERT(1) = ICOMP
      CALL NBREM(NTOT,NCH,1)
*
*       Initialize perturber list for integration of chain c.m.
      CALL CHLIST(ICH)
*
*       Initialize XC and UC
      CALL XCPRED(0)
      
*       Perform differential F & FDOT corrections due to perturbers.
      DO 25 K = 1,3
          FIRR(K) = 0.0D0
          FD(K) = 0.0
   25 CONTINUE
      CALL CHFIRR(ICH,0,X(1,ICH),XDOT(1,ICH),FIRR,FD)
*     --01/03/14 13:06-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      print*,rank,'CHFIRR ICH',ich,'n',name(ich),'x',x(1,ich),
c$$$     &     'xdot',xdot(1,ich),'fdot',fdot(1,ich),'time',time
c$$$      call flush(6)
c$$$      call mpi_barrier(MPI_COMM_NB6,ierr)
*     --01/03/14 13:06-lwang-end----------------------------------------*
      DO 30 K = 1,3
          F(K,ICH) = F(K,ICH) + 0.5*FIRR(K)
          FDOT(K,ICH) = FDOT(K,ICH) + ONE6*FD(K)
          D1(K,ICH) = D1(K,ICH) + FD(K)
   30 CONTINUE
*     --10/29/13 12:36-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*      print*,rank,'firr',ich,name(ich),firr,f(1:3,ich),fdot(1:3,ich)
*     --10/29/13 12:36-lwang-end----------------------------------------*
*
*       Take maximum integration interval equal to c.m. step.
      TMAX = STEP(ICH)
*
*       Check next treatment time of perturbers.
      CALL TCHAIN(NSUB,TSMIN)
      TMAX = MIN(TMAX,TSMIN)
*
*       Copy binding energy and output & capture option for routine CHAIN.
      CM(8) = E(3)
      KZ26 = KZ(26)
      KZ27 = KZ(27)
      KZ30 = KZ(30)
*       Copy velocity scale factor to VSTAR1.
      VSTAR1 = VSTAR
*
*       Assign new subsystem index and begin chain regularization.
      ISUB = NSUB
      NCHAIN = NCHAIN + 1
*
*       Set phase indicator < 0 to ensure new time-step list in INTGRT.
      IPHASE = -1
*
*     --10/29/13 12:33-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$               do i=1,ifirst-1
c$$$*            if(name(i).eq.16107) then
c$$$            write(102+rank,*) i,name(i),step(i),t0(i),time,body(i),
c$$$     *           gamma((i+1)/2),list(1,i),name((i+1)/2+N)
c$$$c$$$*               write(100+rank,*) i,name(i),x0(1,i),x0dot(1,i),
c$$$c$$$*     *              f(1,i),fdot(1,i),time
c$$$c$$$            end if
c$$$         end do
c$$$         do j=ifirst,ntot
c$$$*            if(name(j).eq.16107) then
c$$$            write(102+rank,*),j,name(j),'x0',(x0(kk,j),kk=1,3),
c$$$     *        'x0dot',(x0dot(kk,j),kk=1,3),'t0',t0(j),
c$$$     &           'step',step(j),'stepr',stepr(j),
c$$$     *        'f',(f(kk,j),kk=1,3),'fdot',(fdot(kk,j),kk=1,3),
c$$$     *       'fi',(fi(kk,j),kk=1,3),'fidot',(fidot(kk,j),kk=1,3),
c$$$     *        'd0',(d0(kk,j),kk=1,3),'d1',(d1(kk,j),kk=1,3),
c$$$     *        'd2',(d2(kk,j),kk=1,3),'d3',(d3(kk,j),kk=1,3),
c$$$     *        'd0r',(d0r(kk,j),kk=1,3),'d1r',(d1r(kk,j),kk=1,3),
c$$$     *        'd2r',(d2r(kk,j),kk=1,3),'d3r',(d3r(kk,j),kk=1,3),
c$$$     *        'body',body(j),'time',time,'list',list(1,j)
c$$$            call flush(100+rank)
c$$$*            end if
c$$$         end do
c$$$         call mpi_barrier(MPI_COMM_NB6,ierr)
c$$$
*     --10/29/13 12:33-lwang-end----------------------------------------*
      RETURN
*
      END
