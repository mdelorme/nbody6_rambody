*       NBODY6 parameters.
*       ------------------
*
*       Choose between small or large run.
      PARAMETER  (NMAX=24576,KMAX=8192,LMAX=400,MMAX=512,
     &            MLD=22,MLR=600,MLV=200,MCL=10,NCMAX=10,NTMAX=10000,
     &            NRBDMAX=100000)
      parameter (maxpe=1024,ithread=8)
*
*
*       ------------------------------------------------------
*       NMAX    Maximum number of single bodies + 3*NBIN + NHI.
*       KMAX    Maximum number of KS solutions.
*       LMAX    Maximum size of neighbour lists.
*       MMAX    Maximum number of merged binaries.
*       MLD     Maximum number of disrupted KS components.
*       MLR     Maximum number of recently regularized KS pairs.
*       MLV     Maximum number of high-velocity particles.
*       MCL     Maximum number of interstellar clouds.
*       NCMAX   Maximum number of chain members (do not change).
*       NTMAX   Maximum number of circularizing binaries.
*       ------------------------------------------------------
*
