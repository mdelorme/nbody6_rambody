      SUBROUTINE FPOLY1_MPI(I1,I2,KCASE)
*     
*
*       Force & first derivative.
*       -------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
      REAL*8  A(9),F1(3),F1DOT(3)
      integer inum(maxpe),ista(maxpe)
*
*       Standard case, new c.m. or KS termination (KCASE = 0, 1, 2).
      JLAST = NTOT
*       Reduce loop size for new c.m. polynomial.
      IF (KCASE.EQ.1) JLAST = NTOT - 1
*
*       Loop over all bodies, pair #ICOMP & JCOMP or one single body.
      call cputim(tt998)
      call mpi_barrier(MPI_COMM_NB6,ierr)
      call cputim(tt999)
      if(rank.eq.0)ttbar = ttbar + (tt999-tt998)*60.
      if(rank.eq.0)ibarcount=ibarcount+1
*     
      nl = I2-I1+1
*
      inl = nl/isize
      idiff = nl - isize*inl
      irun = 0
*
      do 1003 ix = 1,isize
      inum(ix)=inl
      if(ix.le.idiff)inum(ix) = inum(ix) + 1
      ista(ix) = irun+1
      irun = irun + inum(ix)
 1003 continue
*
      istart = ista(rank+1) + IFIRST - 1
      iend = ista(rank+1) + inum(rank+1) - 2 + IFIRST
*
      do 40 i = istart,iend
*
*       Initialize forces & first differences for body #I.
      DO 10 K = 1,3
          FI(K,I) = 0.0D0
          FR(K,I) = 0.0D0
          D1(K,I) = 0.0D0
          D1R(K,I) = 0.0D0
   10 CONTINUE
*
*       Obtain force & first derivative by summing over all bodies.
      KDUM = 0
      NNB = LIST(1,I)
*       Set index of first neighbour to be identified in force loop.
      L = 2
      NAMEJ = LIST(L,I)
*
*       Sum over all other bodies.
      DO 30 JDUM = IFIRST,JLAST
          IF (JDUM.EQ.I) GO TO 30
          J = JDUM
          IF (J.GT.N.AND.J.EQ.NAMEJ) THEN
              JPAIR = J - N
*       Use c.m. approximation for unperturbed binary.
              IF (LIST(1,2*JPAIR-1).GT.0) THEN
                  KDUM = 2*JPAIR - 1
                  J = KDUM
              END IF
          END IF
*
   12     DO 15 K = 1,3
              A(K) = X(K,J) - X(K,I)
              A(K+3) = XDOT(K,J) - XDOT(K,I)
   15     CONTINUE
*
          A(7) = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
          A(8) = BODY(J)*A(7)*SQRT(A(7))
          A(9) = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A(7)
*
          DO 20 K = 1,3
              F1(K) = A(K)*A(8)
              F1DOT(K) = (A(K+3) - A(K)*A(9))*A(8)
   20     CONTINUE
*
*       See whether summation index is equal to either component.
          IF (J.EQ.ICOMP.OR.J.EQ.JCOMP) THEN
              IF (KCASE.EQ.1) GO TO 30
*       Note that dominant terms cancel analytically in c.m. polynomial.
          END IF
*
          IF (JDUM.NE.NAMEJ) THEN
              DO 25 K = 1,3
                  FR(K,I) = FR(K,I) + F1(K)
                  D1R(K,I) = D1R(K,I) + F1DOT(K)
   25         CONTINUE
          ELSE
              DO 28 K = 1,3
                  FI(K,I) = FI(K,I) + F1(K)
                  D1(K,I) = D1(K,I) + F1DOT(K)
   28         CONTINUE
*
              IF (J.EQ.KDUM) THEN
                  J = J + 1
                  GO TO 12
              END IF
*
*       Advance the neighbour list until last member has been considered.
              IF (L.LE.NNB) THEN
                  L = L + 1
                  NAMEJ = LIST(L,I)
              END IF
          END IF
   30 CONTINUE
*
*       Check option for interstellar clouds (force & first derivative).
      IF (KZ(13).NE.0.AND.TIME.GT.0.0D0) THEN
          CALL FCLOUD(I,F1,F1DOT,2)
      END IF
*
   40 CONTINUE
*
*        Distribute variables into private vectors again T3D (R.Sp.)
      isend = rank + 1
      if(isend.eq.isize)isend = 0
      irecv = rank - 1
      if(irecv.eq.-1)irecv = isize - 1
*
      do 1001 ir = 0,isize-2
*
      irank = rank - ir
      if(irank.lt.0)irank=irank+isize
*
      istart=ista(irank+1) + IFIRST - 1
      icnt = inum(irank+1)
*
      if(irank.eq.0)irank=isize
      istrec = ista(irank) + IFIRST - 1
      icnt2 = inum(irank)
*
*     print*,' FPOLY1: bef rank,irank=',rank,irank
*     print*,' FPOLY1: bef rank ',rank,' fi(',istart,')=',fi(1,istart)
*     print*,' FPOLY1: bef rank ',rank,' fi(',istrec,')=',fi(1,istrec)
*     print*,' FPOLY1: rank=',rank,' sending ',icnt,' items to ',isend
*     print*,' FPOLY1: rank=',rank,' recving ',icnt2,
*    *  ' items from ',irecv
*
*      call mpi_barrier(MPI_COMM_NB6,ierr)
      CALL MPI_SENDRECV(FI(1,istart),3*icnt,MPI_REAL8,isend,rank,
     *                  FI(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *                  MPI_COMM_NB6,status,ierr)
*      call mpi_barrier(MPI_COMM_NB6,ierr)
      CALL MPI_SENDRECV(FR(1,istart),3*icnt,MPI_REAL8,isend,rank,
     *                  FR(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *                  MPI_COMM_NB6,status,ierr)
*      call mpi_barrier(MPI_COMM_NB6,ierr)
      CALL MPI_SENDRECV(D1(1,istart),3*icnt,MPI_REAL8,isend,rank,
     *                  D1(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *                  MPI_COMM_NB6,status,ierr)
*      call mpi_barrier(MPI_COMM_NB6,ierr)
      CALL MPI_SENDRECV(D1R(1,istart),3*icnt,MPI_REAL8,isend,rank,
     *                  D1R(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *                  MPI_COMM_NB6,status,ierr)
*
      call cputim(tt998)
      call mpi_barrier(MPI_COMM_NB6,ierr)
      call cputim(tt999)
      if(rank.eq.0)ibarcount=ibarcount+1
      if(rank.eq.0)ttbar = ttbar + (tt999-tt998)*60.

 1001 continue
*
*       Check option for external force.
      IF (KZ(14).NE.0.OR.MIMD_MODE) THEN
          CALL XTRNLD(I1,I2,1)
      END IF
*
*       Set total force & first derivative.
      DO 50 I = I1,I2
*
      DO 45 K = 1,3
          F(K,I) = FI(K,I) + FR(K,I)
          FDOT(K,I) = D1(K,I) + D1R(K,I)
   45 CONTINUE
   50 CONTINUE
*
*       Check case of new c.m. (KCASE = 1 with I1 = ICOMP, I2 = JCOMP).
      IF (KCASE.EQ.1) THEN
*       Form total c.m. force & first derivative and also AC components.
          A1 = BODY(ICOMP)
          A2 = BODY(JCOMP)
          DO 60 K = 1,3
              F(K,NTOT) = (A1*F(K,ICOMP) + A2*F(K,JCOMP))/BODY(NTOT)
              FDOT(K,NTOT) = (A1*FDOT(K,ICOMP) + A2*FDOT(K,JCOMP))/
     &                                                        BODY(NTOT)
              FI(K,NTOT) = (A1*FI(K,ICOMP) + A2*FI(K,JCOMP))/BODY(NTOT)
              D1(K,NTOT) = (A1*D1(K,ICOMP) + A2*D1(K,JCOMP))/BODY(NTOT)
              FR(K,NTOT) = (A1*FR(K,ICOMP) + A2*FR(K,JCOMP))/BODY(NTOT)
              D1R(K,NTOT) = (A1*D1R(K,ICOMP) + A2*D1R(K,JCOMP))/
     &                                                        BODY(NTOT)
   60     CONTINUE
          J1 = NTOT
          J2 = NTOT
      ELSE
          J1 = I1
          J2 = I2
      END IF
*
      DO 80 I = J1,J2
          DO 70 K = 1,3
              D0(K,I) = FI(K,I)
              D0R(K,I) = FR(K,I)
              FIDOT(K,I) = D1(K,I)
              FRDOT(K,I) = D1R(K,I)
   70     CONTINUE
   80 CONTINUE
*
      RETURN
*
      END
