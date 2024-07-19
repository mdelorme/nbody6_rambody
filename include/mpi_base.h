*     MPI_BASE

      INCLUDE 'mpif.h'
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      LOGICAL MIMD_MODE
      INTEGER MPI_COMM_NB6, MPI_COMM_RBD
      INTEGER icore,isernb,iserreg,iserks
      COMMON/MPIDAT/ group,rank,ierr,isize,status,icore,isernb,iserreg,
     &               iserks, MPI_COMM_NB6, MPI_COMM_RBD, MIMD_MODE
