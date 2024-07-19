      subroutine rbd_init
      include 'common6.h'
      
      if (KZ(23).LT.1.OR.KZ(23).GT.2) THEN
         write(6,*) 'WARNING : Escapers are not removed (KZ(23) != 1;2)'
         write(6,*) 'RAMBODY Needs this to be set at 1 or 2'
         stop
      end if

      if (KZ(14).NE.0) THEN
         write(6,*) 'WARNING : External perturber set ! (KZ(14) != 0)'
         write(6,*) 'RAMBODY works only with KZ(14) = 0'
         stop
      end if

      RBD_NESCAPERS = 0

      ! TEMP !
      RBD_DBG(1,1) = 1000.0
      RBD_DBG(2,1) = 0.0
      RBD_DBG(3,1) = 0.0

      if (rank .eq. 0) then
         write(6,*) 'RBDCOMM 18'
         call MPI_Recv(RBD_EPSILON, 1, MPI_DOUBLE_PRECISION, 1, 18,
     &        MPI_COMM_RBD, MPI_STATUS_IGNORE, ierr)
      end if
      call MPI_Bcast(RBD_EPSILON, 1, MPI_DOUBLE_PRECISION, 0,
     &     MPI_COMM_NB6, ierr)

      end subroutine


*     Syncing timestep with Ramses, on this side, we only send Ramses the
*     largest block timestep of NBody6
      
      subroutine rbd_sync_timestep
      include 'common6.h'

      real*8 dt_n(1)
      
      ! Sync on this end is very simple : One send to Ramses
      if (rank .eq. 0) then

         dt_n(1) = SMAX * TSTAR

         write(6,*) 'NB6 : Sending timestep to Ramses : ', dt_n
         call flush(6)

         write(6,*) 'RBDCOMM 3'
         call MPI_Send(dt_n, 1, MPI_DOUBLE_PRECISION, 1, 3, 
     &        MPI_COMM_RBD, ierr)


         write(6,*) 'NB6 : Timestep sent, Err = ', ierr
      end if
      end subroutine

      subroutine rbd_send_cluster
      
      include 'common6.h'
      include 'tlist.h'
      include 'galaxy.h'

*     We store everything in a big array
      REAL*8, dimension(7, NMAX) :: buf_send ! Can't be bigger than NMAX (hopefully ?)
      INTEGER :: n_particles
      REAL*8 :: mtot
      INTEGER :: istatus, i, k, i1, i2
      REAL*8, dimension(3) :: xrel, vrel

*     First, sending escapers
      
      write(6,*) 'RBDCOMM 4'
      call flush(6)
      call MPI_Send(RBD_NESCAPERS, 1, MPI_INTEGER, 1, 4, MPI_COMM_RBD,
     &     ierr)
      if (RBD_NESCAPERS > 0) then
         write(6,*) 'RBDCOMM 5'
         call flush(6)
         call MPI_Send(RBD_ESCAPERS, RBD_NESCAPERS*7,
     &        MPI_DOUBLE_PRECISION, 1, 5, MPI_COMM_RBD, ierr)
      end if

      RBD_NESCAPERS = 0

*     For every local process, we get the stars
      mtot = 0.0

*     First, ks pairs (What happens with a pair having a pair ?)
      do i=1, npairs
         call KSRES_OP(i, i1, i2, xrel, vrel, 0)
         
         do k=1, 3
*     Positions
            buf_send(k, i1) = X(k, i1) + BODY(i2)*xrel(k) / BODY(i)
            buf_send(k, i2) = X(k, i2) + BODY(i1)*xrel(k) / BODY(i)
*     Velocities
            buf_send(k+3, i1) = XDOT(k, i1) + BODY(i2)*vrel(k)
     &           / BODY(i)
            buf_send(k+3, i2) = XDOT(k, i2) + BODY(i1)*vrel(k)
     &           / BODY(i)
*     Scaling
            buf_send(k,   i1) = buf_send(k, i1) * RBAR
            buf_send(k,   i2) = buf_send(k, i2) * RBAR
            buf_send(k+3, i1) = buf_send(k+3, i1) * VSTAR
            buf_send(k+3, i2) = buf_send(k+3, i2) * VSTAR
         end do
         
         buf_send(7, i1) = BODY(i1) * ZMBAR
         buf_send(7, i2) = BODY(i2) * ZMBAR
         
         mtot = mtot + buf_send(7, i1)
         mtot = mtot + buf_send(7, i2)
      end do
      
*     Now, single stars
      do i=2*NPAIRS+1, N
         do k=1, 3
            buf_send(k,   i) = X(k, i) * RBAR
            buf_send(k+3, i) = XDOT(k, i) * VSTAR
         end do
         
         buf_send(7, i) = BODY(i) * ZMBAR
         mtot = mtot + buf_send(7, i)
      end do
      
*     Offsetting everything wrt to RDENS
      do i=1,N
         do k=1,3
            buf_send(k, i) = buf_send(k, i) - RDENS(k)
         end do
      end do
      
*     Sending to Ramses
      write(6,*) 'RBDCOMM 6'
      call flush(6)
      call MPI_Send(N, 1, MPI_INTEGER, 1, 6, MPI_COMM_RBD, ierr)
      
      write(6,*) 'RBDCOMM 7'
      call flush(6)
      call MPI_Send(buf_send, N*7, MPI_DOUBLE_PRECISION, 1, 7,
     &     MPI_COMM_RBD, ierr)
      
      write(6,*) 'NB6 : All particles sent ! MTot = ', mtot
      call flush(6)

      
      end subroutine

      subroutine rbd_store_escaper(I)
      include 'common6.h'

      ! Storing escapers to send to Ramses
      RBD_NESCAPERS = RBD_NESCAPERS + 1
      
      do k=1, 3
         RBD_ESCAPERS(k, RBD_NESCAPERS)   = X(K, I) - RDENS(K)
         RBD_ESCAPERS(3+k, RBD_NESCAPERS) = XDOT(K, I)
      end do
      RBD_ESCAPERS(7, RBD_NESCAPERS) = BODY(I)

      end subroutine

      subroutine rbd_send_rscale
      include 'common6.h'
      double precision rsend
      
      write(6,*) 'NB6 : Sending the scale to ramses'
      rsend = RTIDE/RBAR
      call MPI_Send(rsend, 1, MPI_DOUBLE_PRECISION, 1, 8,
     &     MPI_COMM_RBD, ierr)

      end subroutine

      subroutine rbd_sync_mesh
      include 'common6.h'
      include 'galaxy.h'
      include 'tlist.h'

      
      integer np, i, j, i1, i2
      double precision XI(3), XIDOT(3), FMD(3), XG(3), RIJ(3)
      double precision dx, dist
      
      double precision FDIFF(3), FPM(3), FS(3), FM(3), TMP(3)
      double precision FNORM
      double precision r1, r2, enclosed
      
      if (rank .eq. 0) then
         write(6,*) 'NB6 : Waiting for force mesh'
         write(6,*) 'RBDCOMM 16'
         call MPI_Recv(RBD_MESH_SCALE, 1, MPI_DOUBLE_PRECISION, 1, 16,
     &        MPI_COMM_RBD, MPI_STATUS_IGNORE, ierr)

         write(6,*) 'RBDCOMM 9'
         call MPI_Recv(RBD_NMESH, 1, MPI_INTEGER, 1, 9, MPI_COMM_RBD,
     &        MPI_STATUS_IGNORE, ierr)
         
         write(6,*) 'NB6 : Waiting for ', RBD_NMESH, ' points'

         write(6,*) 'RBDCOMM 10'
         call MPI_Recv(RBD_MESH, 6*RBD_NMESH, MPI_DOUBLE_PRECISION, 1, 
     &        10, MPI_COMM_RBD, MPI_STATUS_IGNORE, ierr)
         
         write(6,*) 'RBDCOMM 12'
         call MPI_Recv(RBD_MESH_DESC, 3, MPI_INTEGER, 1, 12,
     &        MPI_COMM_RBD, MPI_STATUS_IGNORE, ierr)

         call MPI_Recv(RBD_DBG, 300, MPI_DOUBLE_PRECISION, 1, 99,
     &        MPI_COMM_RBD, MPI_STATUS_IGNORE, ierr)

         RBD_DBG(:,1) = RBD_DBG(:,1) * RBAR
         RG = RBD_DBG(:,1)
         VG = 0

!     Rescaling
         RBD_MESH_SCALE = RBD_MESH_SCALE / RBAR
         do i=1, RBD_NMESH
            RBD_MESH(1:3,i) = RBD_MESH(1:3,i) / RBAR
            RBD_MESH(4:6,i) = RBD_MESH(4:6,i) / RBAR * TSTAR**2.0
         end do

!     Resetting the force on the cluster
         rbd_cl_force = 0.0
         
!     We define the smoothing length as the distance between two points on the mesh
         !write(6,*) 'Mesh extent:', RBD_MESH(:,2), RBD_MESH(:,RBD_NMESH)
         !write(6,*) 'RBD_EPS = ', RBD_EPSILON
         do i=1, RBD_NMESH
            enclosed = 0.0
            r2=RBD_MESH(1,i)**2+RBD_MESH(2,i)**2+RBD_MESH(3,i)**2
               
            do j=1, N
               r1=X(1,j)**2.0 + X(2,j)**2.0 + X(3,j)**2.0
               RIJ = X(1:3,j) - RBD_MESH(1:3,i)
               if (r1 < r2) then
                  enclosed = enclosed + body(j)
               end if
               
               dist = norm2(RIJ)
               TMP = BODY(j) * RIJ / (dist * (dist**2.0+RBD_EPSILON**2))

               RBD_CL_FORCE(:,i) = RBD_CL_FORCE(:,i) + TMP
            end do
            ! Rescaling to physical units
            RBD_CL_FORCE(:,i) = RBD_CL_FORCE(:,i) * RBAR / TSTAR**2.0
         end do

!     And sending it back
         write(6,*) 'RBDCOMM 17'
         call MPI_Send(RBD_CL_FORCE(:,1:RBD_NMESH), 3*RBD_NMESH, 
     &        MPI_DOUBLE_PRECISION, 1, 17, MPI_COMM_RBD, ierr)
      end if

!     Syncing with other processes
      write(6,*) 'Sharing mesh info with everyone'
      call MPI_BCast(RBD_NMESH, 1, MPI_INTEGER, 0, MPI_COMM_NB6, ierr)
      call MPI_BCast(RBD_MESH_SCALE, 1, MPI_DOUBLE_PRECISION, 0, 
     &     MPI_COMM_NB6, ierr)
      call MPI_BCast(RBD_MESH_DESC, 3, MPI_INTEGER, 0, MPI_COMM_NB6, 
     &     ierr)
      call MPI_BCast(RBD_MESH, 6*RBD_NMESH, MPI_DOUBLE_PRECISION,
     &     0, MPI_COMM_NB6, ierr)
      write(6,*) 'done'

      RBD_ON = .true.

      end subroutine

      subroutine rbd_get_force(XI, FOUT)
      include 'common6.h'

      double precision FOUT(3)
      double precision XI(3)
      integer ii, ij, ik, ngx
      integer iid(8)
      double precision xd, yd, zd
      double precision half_box
      double precision grid_dx
      double precision c00(3), c01(3), c10(3), c11(3)
      double precision c0(3), c1(3)
      
      FOUT = 0.0
      if (RBD_ON .and. norm2(XI) > 0.0) then
         half_box = RBD_MESH_SCALE * 0.5
         ngx = RBD_MESH_DESC(2)
         grid_dx  = RBD_MESH_SCALE / (ngx-1)
         
         ii = (XI(1) + half_box) / grid_dx
         ij = (XI(2) + half_box) / grid_dx
         ik = (XI(3) + half_box) / grid_dx


         ! Extrapolation
         if (ii < 0) ii = 0
         if (ii >= ngx-1) ii=ngx-2
         if (ij < 0) ij = 0
         if (ij >= ngx-1) ij=ngx-2
         if (ik < 0) ik = 0
         if (ik >= ngx-1) ik=ngx-2
         
!     +1 because we start counting at 1
!     and +1 to account for the additional point at the centre of the box so +2
         iid(1) = ii*ngx**2 + ij*ngx + ik + 2 
         iid(2) = iid(1)+1
         iid(3) = iid(1)+ngx
         iid(4) = iid(1)+ngx+1
         iid(5) = iid(1)+ngx**2
         iid(6) = iid(1)+ngx**2+1
         iid(7) = iid(1)+ngx**2+ngx
         iid(8) = iid(1)+ngx**2+ngx+1
         
         xd = (XI(1)-RBD_MESH(1,iid(1))) /
     &        (RBD_MESH(1,iid(8))-RBD_MESH(1,iid(1)))

         yd = (XI(2)-RBD_MESH(2,iid(1))) /
     &        (RBD_MESH(2,iid(8))-RBD_MESH(2,iid(1)))

         zd = (XI(3)-RBD_MESH(3,iid(1))) /
     &        (RBD_MESH(3,iid(8))-RBD_MESH(3,iid(1)))

         ! Trilinear interpolation, for each acceleration direction
         do i=1, 3
            c00(i) = RBD_MESH(3+i,iid(1))*(1-xd)+RBD_MESH(3+i,iid(5))*xd
            c01(i) = RBD_MESH(3+i,iid(2))*(1-xd)+RBD_MESH(3+i,iid(6))*xd
            c10(i) = RBD_MESH(3+i,iid(3))*(1-xd)+RBD_MESH(3+i,iid(7))*xd
            c11(i) = RBD_MESH(3+i,iid(4))*(1-xd)+RBD_MESH(3+i,iid(8))*xd

            c0(i) = c00(i)*(1-yd)+c10(i)*yd
            c1(i) = c01(i)*(1-yd)+c11(i)*yd

            FOUT(i) = c0(i)*(1-zd)+c1(i)*zd
         end do

         if (ISNAN(FOUT(1))) then
!     Problem here ?
            write(6,*) 'ERROR : Force is not a number !'
            
            write(6,*) 'Half box, DX : ', half_box, grid_dx
            write(6,'(a, 3f15.7)') 'LERP diag at position :', XI
            write(6,*) 'Positions of top left, bottom right :',
     &           RBD_MESH(1:3, 2), RBD_MESH(1:3, ngx**3+1)
            write(6,*) 'Positional ids :', ii, ij, ik
            write(6,*) 'Mesh positions/forces :'
            do i=1, 8
               write(6,'(a, 6f15.7)') ' . ', RBD_MESH(:,iid(i))
            end do
            
            write(6,*) 'Forces on both corners of the mesh'
            write(6, '(6f15.7)') RBD_MESH(4:6,2),
     &           RBD_MESH(4:6, ngx**3+1)
            
            write(6,'(a, 6f15.7)') 'Positions/forces of centre point :',
     &           RBD_MESH(:,1)
            

            write(6,*) 'Interp coeffs :', xd, yd, zd

            write(6,*) 'Interpolation on Fx :'
            write(6,*) 'CXX : ', c00(1), c01(1), c10(1), c11(1)
            write(6,*) 'CX : ', c0(1), c1(1)
            write(6,*) 'FX : ', FOUT(1)
         end if
         
      end if
      
      end subroutine

      subroutine rbd_restart
      include 'common6.h'

      if (rank .eq. 0) then
         write(6,*) 'RAMBODY Restart !'

!     Comm tag 13 -> Number of points
         write(6,*) 'RBDCOMM 13'
         call flush(6)
         call MPI_Recv(RBD_NMESH, 1, MPI_INTEGER, 1, 13, MPI_COMM_RBD,
     &        MPI_STATUS_IGNORE, ierr)

         write(6,*) RBD_NMESH, 'points in mesh'
         call flush(6)
         
!     Comm tag 14 -> Mesh buffer
         write(6,*) 'RBDCOMM 14'
         call flush(6)
         call MPI_Recv(RBD_MESH, 6*RBD_NMESH, MPI_DOUBLE_PRECISION, 1,
     &        14, MPI_COMM_RBD, MPI_STATUS_IGNORE, ierr)


!     Comm tag 15 -> Mesh description
         write(6,*) 'RBDCOMM 15'
         call flush(6)
         call MPI_Recv(RBD_MESH_DESC, 3, MPI_INTEGER, 1, 15,
     &        MPI_COMM_RBD, MPI_STATUS_IGNORE, ierr)

         do i=1, RBD_NMESH
            RBD_MESH(1:3,i) = RBD_MESH(1:3,i) / RBAR
            RBD_MESH(4:6,i) = RBD_MESH(4:6,i) / RBAR * TSTAR**2.0
         end do
      end if

      write(6,*) 'Sharing mesh info with everyone (restart)'
      call MPI_BCast(RBD_NMESH, 1, MPI_INTEGER, 0, MPI_COMM_NB6, ierr)
      call MPI_BCast(RBD_MESH_DESC, 3, MPI_INTEGER, 0, MPI_COMM_NB6, 
     &     ierr)
      call MPI_BCast(RBD_MESH, 6*RBD_NMESH, MPI_DOUBLE_PRECISION,
     &     0, MPI_COMM_NB6, ierr)
      write(6,*) 'done'
         
      RBD_ON = .true.
         
      end subroutine
