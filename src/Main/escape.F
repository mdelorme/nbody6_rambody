      SUBROUTINE ESCAPE
*
*
*       Escaper removal.
*       ----------------
*
      INCLUDE 'common6.h'
#ifdef TT
      include 'tt.h'
#endif
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
      CHARACTER*11  WHICH1
      LOGICAL FIRST,FIRST2
      REAL*8 XLIST(2*NMAX)
#ifdef TT
      LOGICAL RFLAG
#endif
      SAVE FIRST
      DATA FIRST/.true./
*      Write stellar evolution data on esc.11, escbin.31
      REAL*8 TSCLS(20),LUMS(10),GB(10),TM,TN
      REAL*8 M0,M1,RM,LUM,AGE,TEFF,MC,RCC
      REAL*8 M01,M11,RM1,LUM1,AGE1,TEFF1
      REAL*8 M02,M12,RM2,LUM2,AGE2,TEFF2
      REAL*8 MENV,RENV,K2
*
*
*       Adopt twice the tidal radius as escape condition.
      RESC2 = 4.0*RTIDE**2
*       For tidal cutoff check only energy
      IF(KZ(23).EQ.3)THEN
          RESC2 = 0.D0
          ETID = ZMASS/RTIDE
          ZMOLD = ZMASS
      END IF
*
      RTIDE2 = RTIDE**2
      NCORR = 0
      NCRIT1 = 0
      NCRIT2 = 0
      I3HI = 0
      DO 1 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
    1 CONTINUE
*
*       Set the distance (squared) with respect to the density centre.
      I = IFIRST
    5 RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 + 
     &                               (X(3,I) - RDENS(3))**2
      IF (RI2.LT.RTIDE2) NCRIT1 = NCRIT1 + 1
*
      IF(KZ(23).LE.2) THEN
         VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
*       See whether escape is indicated (retain ghost particles).
#ifdef TT
         RFLAG = .false.
         IF (TTMODE.EQ.1) THEN
*     Consider distance along different tensor eigenvector direction for escape criterion
            DO K=1,3
               RTIDEX = X(1,I)*TTEIGENV(1,K) + X(2,I)*TTEIGENV(2,K) 
     &              + X(3,I)*TTEIGENV(3,K)
               IF(ABS(RTIDEX).GT.2*TTRTIDE(K).AND.TTRTIDE(K).GT.0) 
     &              RFLAG = .true.
            END DO
         ELSE
            IF (RI2.GT.RESC2) RFLAG = .true.
         END IF
*       Use energy as escaper criterion also to avoid tidal shock
         IF (RFLAG.AND.VI2.GT.4*PHIDBL(I).AND.RI2.LT.1.0D+10) THEN
#else
         IF (RI2.GT.RESC2.AND.RI2.LT.1.0D+10) THEN
#endif
            RJMIN2 = 1.0D+10
*     Find distance to the nearest neighbour and calculate potential.
            POTI = 0.0D0
            DO 8 J = IFIRST,NTOT
               IF (J.EQ.I) GO TO 8
               RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &              (X(3,I) - X(3,J))**2
               POTI = POTI + BODY(J)/SQRT(RIJ2)
               IF (RIJ2.LT.RJMIN2) THEN
                  RJMIN2 = RIJ2
                  JMIN = J
               END IF
 8          CONTINUE
*     Check escape criterion for tidal case or isolated system.
            EI = 0.5*VI2 - POTI
#ifdef TT
*** FlorentR - add nbody6tt
            IF (KZ(14).EQ.4.OR.KZ(14).EQ.3.OR.KZ(14).EQ.9.OR.MIMD_MODE)
     &           THEN
*** FRenaud
#else
          IF (KZ(14).EQ.4.OR.KZ(14).EQ.3.OR.MIMD_MODE) THEN
#endif
              EI = EI - MP/SQRT(RI2 + AP2)
              IF (BODY(I).EQ.0.0D0) GO TO 30
          END IF
          IF ((KZ(14).GT.0.AND.KZ(14).NE.4).OR.EI.GT.0.0) GO TO 30
          IF (MIMD_MODE) GO TO 30
          IF (NAME(I).LT.0) GO TO 30
         END IF
*
      ELSE
          IF(RI2.LT.1.0D+10) THEN
*       Check tidal cutoff criterion - use last computed phi-value
*       Find distance to the nearest neighbour
              POTI = 0.0D0
              NNB = LIST(1,I)
              RJMIN2 = 1.0D+10
              DO 9 L = 1,NNB
                  J = LIST(L+1,I)
                  RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                          (X(3,I) - X(3,J))**2
                  IF (RIJ2.LT.RJMIN2) THEN
                      RJMIN2 = RIJ2
                      JMIN = J
                  END IF
    9         CONTINUE
              VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
              POTI = -PHIDBL(I)
*       Check escape criterion for tidal case or isolated system.
              EI = 0.5*VI2 - POTI + ETID
*
              IF (EI.GT.0.0.AND.RI2.GT.RTIDE2) THEN
              NCRIT2 = NCRIT2 + 1
              GO TO 30
              END IF
          END IF
      END IF
*
   10 I = I + 1
   12 IF (N.EQ.2) GO TO 13
*
      IF (I.LE.NTOT) GO TO 5
      IF (NCORR.EQ.0) GO TO 25
*
*       Form centre of mass terms.
   13 DO 15 I = 1,N
          DO 14 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)/ZMASS
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)/ZMASS
   14     CONTINUE
   15 CONTINUE
      CMR(4) = SQRT(CMR(1)**2 + CMR(2)**2 + CMR(3)**2)
*
      if(rank.eq.0) then
        ntmp = min(ncorr, 1000)
        
        WRITE (6,18) TTOT, NCORR, N, NSESC, NBESC, ZMASS*ZMBAR, BE(3),
     &       CMR(4),DSQRT(RESC2), ZMASS/FLOAT(N)*ZMBAR, NCRIT1, 
     &       (JLIST(J),J=1,ntmp)
   18 FORMAT (/,' ESCAPE T,NESC =',1P,E13.5,I6,' total N,NSESC,NBESC =',
     & I10,I8,I6,' M* ',1P,E13.5,' total E, RCM, RESC, <m*>, NCRIT1 =',
     &    4E13.5,I10,' JLIST: ',1000I10)
*
      end if
      IF (KZ(23).EQ.4)THEN
          JLAST = MIN(2*NCORR,NMAX)
          if(rank.eq.0) then 
            WRITE (6,20)
     &        (JLIST(J),(XLIST(2*J+KK),KK=-1,0),J=1,ntmp)
   20 FORMAT(' ESCAPE ANGLES ',1P,1000('| ',I10,2E13.5))
          end if
      END IF
*       Check updating of global index for chain c.m.
      IF (NCH.GT.0) THEN
          CALL CHFIND
      END IF
*
*       Set phase indicator < 0 for new NLIST in routine INTGRT (Hermite).
      IPHASE = -1
*
   25 RETURN
*
   30 A2 = (X(1,JMIN) - RDENS(1))**2 + (X(2,JMIN) - RDENS(2))**2 +
     &                                 (X(3,JMIN) - RDENS(3))**2

*       See whether nearest body satisfies escape condition or RIJ > 10*RMIN.
      IF (A2.GT.RESC2.AND.KZ(14).GT.0) GO TO 40
      IF (RJMIN2.GT.100.0*RMIN2) GO TO 40
*
      IF(KZ(23).EQ.3)THEN
          VI2 = XDOT(1,JMIN)**2 + XDOT(2,JMIN)**2 + XDOT(3,JMIN)**2
          EI = 0.5*VI2 - PHIDBL(JMIN) + ETID
          IF(EI.GT.0.D0) GO TO 40
      END IF
*
      A3 = XDOT(1,JMIN)**2 + XDOT(2,JMIN)**2 + XDOT(3,JMIN)**2
      A4 = (XDOT(1,I) - XDOT(1,JMIN))**2 +
     &     (XDOT(2,I) - XDOT(2,JMIN))**2 +
     &     (XDOT(3,I) - XDOT(3,JMIN))**2
      A5 = (BODY(I) + BODY(JMIN))/SQRT(RJMIN2)
*       Check velocity of binary component in case of bound pair.
      A6 = 2.0*A5/A4
      IF (A6.GT.1.0.AND.A3.GT.2.0*ZMASS/SQRT(A2)) GO TO 40
*       Retain escaper if dynamical effect on neighbour is significant.
      IF (A6.GT.0.01) GO TO 10
*
*       Form optional output diagnostics.
   40 X1 = X(1,I) - RDENS(1)
      Y1 = X(2,I) - RDENS(2)
      Z1 = X(3,I) - RDENS(3)
      A3 = Y1/X1
      A4 = Z1/SQRT(X1**2 + Y1**2)
      IF (Y1.GT.0.D0) THEN
          XLIST(2*NCORR+1) = 90.0 + DATAN(A3)*180.0/3.14159
      ELSE
          XLIST(2*NCORR+1) = 270.0 + DATAN(A3)*180.0/3.14159
      END IF
      XLIST(2*NCORR+2) = DATAN(A4)*180.0/3.14159
*       Escape angles with respect to the X-axis and XY-plane.
      RESC = SQRT(RI2)
      STEPI = STEP(I)
      RSI = RS(I)
*
*       Accumulate escaper names and save current name in case I > N.
*     PRINT*,' ESC-Chk I,NM,M,R,n,ncorr,angles=',I,NAME(I),BODY(I),
*    &      SQRT(RI2),N,NCORR
      NCORR = NCORR + 1
      JLIST(NCORR) = NAME(I)
      NAMEI = NAME(I)
      KSTARI = KSTAR(I)
      IF (NAME(I).LT.0) KSTARI = KSTARI + 30
      IF (BODY(I).LE.0.0D0) GO TO 50
*
#ifdef TT
*       Check initialization of tidal tail integration for 3D model.
*** FlorentR - add nbody6tt
      IF (KZ(23).GE.3.AND.(KZ(14).EQ.3.OR.KZ(14).EQ.9.OR.MIMD_MODE))
     &     THEN
*** FRenaud
#else
      IF (KZ(23).GE.3.AND.(KZ(14).EQ.3.OR.MIMD_MODE)) THEN
#endif
         CALL TAIL0(I)
      END IF

      IF (MIMD_MODE) THEN
         CALL rbd_store_escaper(I)
      END IF

*       Obtain binding energy of body #I and update optional tidal radius.
      ZK = 0.5D0*BODY(I)*VI2
#ifdef TT
      IF (KZ(14).GT.0.AND.KZ(14).LT.3) THEN
#else
      IF (KZ(14).GT.0) THEN 
#endif
          CALL XTRNLV(I,I)
          ZK = ZK + HT
          RTIDE = (ZMASS/TIDAL(1))**0.3333
      ELSE IF (KZ(14).EQ.4.OR.KZ(14).EQ.3) THEN
          RTIDE = RTIDE0*ZMASS**0.3333
#ifdef TT
      ELSE IF (KZ(14).EQ.9) THEN
         IF(TTMODE.NE.1) THEN
*     Temporary solution, need to be checked in the future
            RTIDE = 10.*RSCALE
         END IF
*** FRenaud
#endif
      ELSE IF (MIMD_MODE) THEN
         RTIDE = 10.*RSCALE !
      END IF
      EI = ZK - BODY(I)*POTI
*
*       Update tidal radius in case of tidal cutoff (R.Sp.)
      IF(KZ(23).GE.3)THEN
          RTOLD = RTIDE
          RTIDE = RTIDE*(ZMASS/ZMOLD)**0.3333
      END IF
*
*       Correct total energy.
      BE(3) = BE(3) - EI
      ESESC = ESESC + EI
*
*       Update total mass and save energy & number of single escapers.
      IF (I.LE.N) THEN
          ZMASS = ZMASS - BODY(I)
          E(4) = E(4) + EI
          NPOP(4) = NPOP(4) + 1
          NSESC = NSESC + 1
      END IF
*
*       Include optional escape output on unit 11.
      IF (KZ(23).EQ.2.OR.KZ(23).EQ.4.OR.KZ(23).EQ.3) THEN
          TESC = TSCALE*TTOT
          EESC = EI/BODY(I)
*       Distinguish between tidal field and isolated system (ECRIT at RTIDE).
          IF (KZ(14).GT.0.AND.KZ(14).LE.2) THEN
              ECRIT = -1.5*(TIDAL(1)*ZMASS**2)**0.3333
              BESC = ZMBAR*BODY(I)
          ELSE
              BESC = ZMBAR*BODY(I)
          END IF
          VKM = SQRT(VI2)*VSTAR
*
          IF(KZ(19).GT.0.AND.I.LE.N) THEN
              M0 = BODY(I)*ZMBAR
              M1 = M0
              KW = KSTAR(I)
              AGE = TIME*TSTAR - EPOCH(I)
              CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
              CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                    RM,LUM,KW,MC,RCC,MENV,RENV,K2)
              TEFF = 0.25*(LOG10(LUM) - 2.0*LOG10(RM)) + 3.761777537508
          ELSE
              LUM = 0.E0
              RM = 0.E0
              M1 = 0.E0
              TEFF = 0.E0
          END IF
          if(rank.eq.0)then
             if(kstart.eq.1.and.first) then
                first = .false.
                write (11,*) '        ',
     &       'TTOT         ','BODY         ',
     &       'RI           ','VI           ','STEP         ',
     &       'T[Myr]       ','M[M*]      ','EESC      ',
     &       'VI[km/s]     ','K*  ','NAME      ',
     &       'ANGLE PHI    ',' ANGLE THETA ','M1[M*]       ',
     &       'RADIUS[RSun] ',' LUM[LSun]   ',' TEFF        ',
     &       'AGE[Myr]     ',' EPOCH       ',' ZKIN        ',
     &       ' POT         ',' ECRIT       ',' ZMASS       ',
     &       ' TIDAL(1-4)'
             end if
          WRITE (11,500) TTOT,BODY(I),SQRT(RI2),SQRT(VI2),STEP(I), 
     &    TESC,BESC,EESC,VKM,KSTARI,NAMEI,(XLIST(2*NCORR+KK),KK=-1,0),
     &    M1,RM,LUM,TEFF,AGE,EPOCH(I),ZKIN,POT,ECRIT,ZMASS,
     &    (TIDAL(K),K=1,4)
 500  FORMAT(1X,1P,9E13.5,I4,I10,16E13.5)
             call flush(11)
          end if
*     More output of escaping binaries
          IF(I.GT.N) THEN
             IPAIR = I - N
             I1 = 2*IPAIR - 1
             I2 = 2*IPAIR
             IF(KZ(19).GT.0) THEN
                 M01 = BODY(I1)*ZMBAR
                 M11 = M01
                 KW1 = KSTAR(I1)
                 AGE1 = TIME*TSTAR - EPOCH(I1)
                 CALL star(KW1,M01,M11,TM,TN,TSCLS,LUMS,GB,ZPARS)
                 CALL hrdiag(M01,AGE1,M11,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                       RM1,LUM1,KW1,MC,RCC,MENV,RENV,K2)
                 TEFF1=0.25*(LOG10(LUM1)-2.0*LOG10(RM1))+3.761777537508
                 M02 = BODY(I2)*ZMBAR
                 M12 = M02
                 KW2 = KSTAR(I2)
                 AGE2 = TIME*TSTAR - EPOCH(I2)
                 CALL star(KW2,M02,M12,TM,TN,TSCLS,LUMS,GB,ZPARS) 
                 CALL hrdiag(M02,AGE2,M12,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                       RM2,LUM2,KW2,MC,RCC,MENV,RENV,K2)
                 TEFF2=0.25*(LOG10(LUM2)-2.0*LOG10(RM2))+3.761777537508
             ELSE
                 LUM1 = 0.E0
                 RM1 = 0.E0
                 M11 = 0.E0
                 TEFF1 = 0.E0
                 LUM2 = 0.E0
                 RM2 = 0.E0
                 M12 = 0.E0
                 TEFF2 = 0.E0
             END IF
             if(rank.eq.0)then
             if(kstart.eq.1.and.first2) then
                first2 = .false.
                write (31,502)
 502  FORMAT(1X,'TTOT         ','BODY1,2      ',13X,
     &          'NAME1,2  ',10X,'KW1,2  ','M1,2 [Sol]   ',13X,
     &          'RAD1,2[Sol]  ',13X,'LUM1,2[Sol]  ',13X,
     &          'TEFF1,2 [K]  ',13X,'AGE1,2[Myr]  ',13X,
     &          'EPOCH1,2 [Myr]')
             end if
          WRITE (31,501) TTOT,BODY(I1),BODY(I2),NAME(I1),NAME(I2),
     &    KW1,KW2,M11,M12,RM1,RM2,LUM1,LUM2,TEFF1,TEFF2,AGE1,AGE2,
     &    EPOCH(I1),EPOCH(I2)
 501  FORMAT(1X,1P,3E13.5,2I10,2I4,12E13.5)
             call flush(31)
             end if
          END IF
      END IF
*
*       Reduce particle number & total membership and check NNBMAX.
   50 N = N - 1
      NTOT = NTOT - 1
      NNBMAX = MIN(N/2,NNBMAX)
      ZNBMAX = 0.9*FLOAT(NNBMAX)
      ZNBMIN = MAX(0.01*FLOAT(NNBMAX),1.0)
*       Ensure NNBOPT is smaller than NNBMAX
      NNBOPT = MIN(NNBMAX-1,NNBOPT)
      NNBOPT = MAX(NNBOPT,1)
*       Set indicator for removal of c.m. or KS components.
      KCASE = 1
*
*       Update COMMON arrays to remove the escaper and correct F & FDOT.
      CALL REMOVE(I,1)
*
*       Delete escaper from neighbour lists and reduce higher locations.
   60 DO 150 J = 1,NTOT
          NNB = LIST(1,J)
          IF (NNB.EQ.0) GO TO 150
          L = 2
   70     IF (LIST(L,J).NE.I) GO TO 130
*
*       Move up the remaining list members and reduce membership by one.
          DO 80 K = L,NNB
              LIST(K,J) = LIST(K+1,J)
   80     CONTINUE
          LIST(1,J) = LIST(1,J) - 1
*       Reduce the steps to minimize error effect (do not allow DT < 0).
*         STEP(J) = MAX(0.5D0*STEP(J),TIME - T0(J))
*         STEPR(J) = MAX(0.5D0*STEPR(J),TIME - T0R(J))
          IF (LIST(1,J).GT.0) GO TO 130
*
*       Add a distant body as neighbour if list only contains escaper.
          K = IFIRST - 1
  100     K = K + 1
          RJK2 = (X(1,J) - X(1,K))**2 + (X(2,J) - X(2,K))**2 +
     &                                  (X(3,J) - X(3,K))**2
          IF (RJK2.LT.RESC2.AND.K.LT.N) GO TO 100
          LIST(1,J) = 1
          LIST(2,J) = K
          GO TO 150
*
*       Reduce higher particle locations by one.
  130     IF (LIST(L,J).GT.I) LIST(L,J) = LIST(L,J) - 1
          L = L + 1
          IF (L.LE.LIST(1,J) + 1) GO TO 70
  150 CONTINUE
*
*       Update list of old KS components (remove #I and rename > I).
      NNB = LISTR(1)
      DO 170 L = 2,NNB+1
          IF (LISTR(L).EQ.I) THEN
*       Remove both components of pair and reduce membership by two.
              J = 2*KVEC(L-1)
              DO 165 K = J,NNB-1
                  LISTR(K) = LISTR(K+2)
  165         CONTINUE
              LISTR(1) = LISTR(1) - 2
          END IF
  170 CONTINUE
*
*       Reduce higher particle locations by one (separate loop for pairs).
      DO 180 L = 2,NNB+1
          IF (LISTR(L).GT.I) LISTR(L) = LISTR(L) - 1
  180 CONTINUE
*
*       Update list of high velocity particles (remove #I and rename > I).
      NNB = LISTV(1)
      DO 190 L = 2,NNB+1
          IF (LISTV(L).EQ.I) THEN
*       Remove escaper and reduce the membership.
              DO 185 K = L,NNB
                  LISTV(K) = LISTV(K+1)
  185         CONTINUE
              LISTV(1) = LISTV(1) - 1
          END IF
*       Reduce higher particle locations by one (or three for c.m.).
          IF (LISTV(L).GT.I) THEN
              LISTV(L) = LISTV(L) - 1
              IF (I.GT.N + 1) LISTV(L) = LISTV(L) - 2
          END IF
  190 CONTINUE
*
*       Check special case of second KS component removal.
      IF (KCASE.GT.1) GO TO 205
*
*       See whether the escaper is a single particle or c.m.
      IF (I.LE.N + 1) GO TO 12
*
*       Prepare removal of regularized pair.
      IPAIR = I - N - 1
*
*       Skip correction if ghost is also merged binary (NAMEI = 0 below).
      IF (NAMEI.NE.0.AND.BODY(2*IPAIR-1).GT.0.0D0) THEN
          ZMB = BODY(2*IPAIR-1) + BODY(2*IPAIR)
          ZM2 = BODY(2*IPAIR)
          EB = H(IPAIR)*BODY(2*IPAIR-1)*BODY(2*IPAIR)/ZMB
          SEMI = -0.5*ZMB/H(IPAIR)
          ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(SEMI*ZMB)
          ECC = SQRT(ECC2)
          PMIN  = SEMI*(1.0 - ECC)
          RATIO = MAX(RADIUS(2*IPAIR-1),RADIUS(2*IPAIR))/PMIN
          I1 = 2*IPAIR - 1
          I2 = 2*IPAIR
          NAME1 = NAME(I1)
          NAME2 = NAME(I2)
          KSTAR1 = KSTAR(2*IPAIR-1)
          PB = DAYS*SEMI*SQRT(SEMI/ZMB)
      ELSE
*       Obtain two-body elements of ghost binary and update energies.
          EB = 0.0D0
          RATIO = 0.0
          BODYCM = CM(3,IM) + CM(4,IM)
          SEMI = -0.5*BODYCM/H(JPAIR)
          ZMU = CM(3,IM)*CM(4,IM)/BODYCM
          PB = DAYS*SEMI*SQRT(SEMI/BODYCM)
          ECC2 = (1.0-R(JPAIR)/SEMI)**2 + TDOT2(JPAIR)**2/(BODYCM*SEMI)
          ECC = SQRT(ECC2)
          EB1 = ZMU*H(JPAIR)
          BE(3) = BE(3) - EB1
          EMESC = EMESC + EB1
          EMERGE = EMERGE - EB1
          DEB = DEB + EB1
      END IF
*
*       Update total energy (ECOLL with EB < -1 & #27 > 0 affects BINOUT).
      BE(3) = BE(3) - EB
      EBESC = EBESC + EB
*
*       Check optional diagnostics for hierarchical systems.
      IF (NAMEI.LE.0.AND.(KZ(18).EQ.1.OR.KZ(18).EQ.3)) THEN
          IPHASE = 7
          CALL HIARCH(IPAIR)
      END IF
*
*       Specify binary type (0: standard; -1: merged binary).
      M = 0
      IF (NAMEI.LE.0) M = -1
*
*       Distinguish between actual and ghost binary (mergers come later).
      IF (BODY(2*IPAIR-1).GT.0.0D0) THEN
*
*       Include rare case of higher-order system (4, 5 or 6 members).
          IF (NAMEI.LT.-2*NZERO) THEN
              JM = 1
              DO 195 K = 1,NMERGE
                  IF (NAMEM(K).EQ.NAMEI) JM = K
  195         CONTINUE
              I3HJ = N
              DO 196 J = IFIRST,NTOT
                  IF (NAME(J).EQ.NAMEG(JM)) I3HJ = J
  196         CONTINUE
*       Identify quartet, quintuplet or sextuplet.
              IF (NAME2.LE.NZERO.AND.NAME(I3HJ).LE.NZERO) THEN
                  WHICH1 = '  QUARTET  '
*       Include both types of quintuplet: [[B,S],B] and [[B,B],S].
              ELSE IF (NAME2.LE.NZERO.OR.NAME(I3HJ).LE.NZERO) THEN
                  WHICH1 = ' QUINTUPLET'
              ELSE
                  WHICH1 = ' SEXTUPLET '
              END IF
              ZM1 = CM(1,JM)*ZMBAR
              ZM2 = CM(2,JM)*ZMBAR
              EB3 = CM(1,JM)*CM(2,JM)*HM(JM)/(CM(1,JM) + CM(2,JM))
              SEMI3 = -0.5*(CM(1,JM) + CM(2,JM))/HM(JM)
              PB3 = DAYS*SEMI3*SQRT(SEMI3/(CM(1,JM) + CM(2,JM)))
              if(rank.eq.0)
     &        WRITE (6,199)  WHICH1, NAME(2*IPAIR-1), NAME(2*IPAIR),
     &                       NAME(I3HJ), ZM1, ZM2, EB3, SEMI3, PB3
  199         FORMAT (/,A11,' ESCAPE    NM =',3I6,'  M =',2F5.1,
     &                    '  EB3 =',F8.4,'  A3 =',1P,E8.1,'  P3 =',E8.1)
          END IF
*
          if(rank.eq.0)then
          VI = DSQRT(VI2)
          RI = DSQRT(RI2)
          NICM = JLIST(NCORR)
          WRITE (6,200)  IPAIR,N,NCORR,NAME1,NAME2,NICM,
     &          KSTAR(I1),KSTAR(I2),KSTARI,LIST(2,I2),
     &          BODY(2*IPAIR-1)*ZMBAR,BODY(2*IPAIR)*ZMBAR,
     &          SEMI,ECC,EB,PB,RATIO,VI,RI,
     &          EI,RADIUS(I1)*SU,RADIUS(I2)*SU
  200     FORMAT (/,' BINARY ESCAPE  IP,N,NCORR =',I5,I10,I5,
     &                '  NM1,2,S=',3I10,'  KW1,2,S,LS=',4I4,
     &                '  M1,2* =',1P,2E13.5,
     &                '  a,e,eb,pb =',4E13.5,' R*/PM =',E13.5,
     &                '  RI,VI =',2E13.5,'  EI =',E13.5,
     &                '  RAD1,2[*]',2E13.5)
          end if
      ELSE
          if(rank.eq.0)
     &    WRITE (6,202)  IPAIR,NAME1,NAME2,LIST(2,2*IPAIR),KSTAR(JM),
     &                   CM(3,IM)*ZMBAR,CM(4,IM)*ZMBAR,EB1,ECC,SEMI,PB
  202     FORMAT (/,' QUAD ESCAPE    KS =',I5,'  NM =',2I6,
     &                '  K* =',2I3,'  M =',2F5.1,'  EB =',F8.4,
     &                '  E =',F7.3,'  A =',1P,E8.1,'  P =',E8.1)
          EB = EB + EB1
          EI = 0.0
      END IF
*
*       Accumulate escaping binary energies and increase the counter.
      IF (LIST(2,2*IPAIR).EQ.-1) THEN
          E(5) = E(5) + EB
          E(6) = E(6) + EI
          NPOP(5) = NPOP(5) + 1
      ELSE
          E(7) = E(7) + EB
          E(8) = E(8) + EI
          NPOP(6) = NPOP(6) + 1
      END IF
*
      IF (M.EQ.0) THEN
          NBESC = NBESC + 1
      ELSE
          NMESC = NMESC + 1
      END IF
*
*       Reduce particle number, pair index & single particle index. 
      N = N - 1
      NPAIRS = NPAIRS - 1
      IFIRST = 2*NPAIRS + 1
*       Initialize NNTB correctly
      IF (NPAIRS.EQ.0) NNTB = 0 
*
*       Move up all tables of regularized pairs below IPAIR.
      IF (IPAIR.LE.NPAIRS) THEN
          CALL REMOVE(IPAIR,2)
      END IF
*
*       Increase index for removing KS components.
  205 KCASE = KCASE + 1
      IF (KCASE.LE.3) THEN
*       Remove COMMON arrays of the second component before the first.
          I = 2*IPAIR + 2 - KCASE
          NTOT = NTOT - 1
*       Reduce NTOT by 3 and N by 2 when KS pair escapes.
          CALL REMOVE(I,3)
          GO TO 60
      END IF
*
*       Check selection of possible ghost in higher-order system.
      IF (I3HI.GT.0) THEN
          I = 0
          DO 208 J = IFIRST,NTOT
              IF (NAME(J).EQ.NAMEG(JM)) THEN
                  I = J
              END IF
  208     CONTINUE
          IF (I.GT.0) THEN
              I3HI = 0
              GO TO 50
          END IF
      END IF
*
*       Include the case of escaping merger.
      IF (NAMEI.GE.0) GO TO 250
*
*       Locate current position in the merger table (standard case).
      IM = 0
      DO 210 K = 1,NMERGE
          IF (NAMEM(K).EQ.NAMEI) IM = K
  210 CONTINUE
*       Skip on failed detection just in case.
      IF (IM.EQ.0) GO TO 250
*
*       Include case of higher-order system (outer single or binary).
      DEB = 0.0
      IF (NAMEI.LT.0) THEN
*       Determine the ghost index.
          I3HI = 0
          DO 215 J = IFIRST,NTOT
              IF (NAME(J).EQ.NAMEG(IM)) I3HI = J
  215     CONTINUE
*       Specify nominal escape distance for ghost removal.
          X(1,I3HI) = 1.0D+04
*       Define possible KS pair index for quadruple system correction.
          JPAIR = I3HI - N
      END IF
*
*       Consider current ghost unless deeper hierarchy is present.
      JM = IM
      IF (I3HI.GT.0) THEN
*       Search for c.m. name one level below current (NAMEI < 0).
          IF (NAMEI.LT.-2*NZERO) THEN
              DO 220 K = 1,NMERGE
                  IF (NAMEM(K).EQ.NAMEI + 2*NZERO) JM = K
  220         CONTINUE
*       Use previous merger index to look for binary ghost at earlier level.
              I30 = I3HI
              DO 225 J = IFIRST,NTOT
                  IF (NAME(J).EQ.NAMEG(JM)) I3HI = J
  225         CONTINUE
              IF (I3HI.GT.N) THEN
                  JPAIR = I3HI - N
              ELSE
*       Adopt original solution on failure to identify binary.
                  I3HI = I30
              END IF
*       Set nominal escape distance of any new ghost (I3HI <= N possible).
              X(1,I3HI) = 1.0D+04
          END IF
      END IF
*
*       Copy merger energy to respective energy bins (primordial or new).
      ZMU = CM(1,JM)*CM(2,JM)/(CM(1,JM) + CM(2,JM))
      IF (MIN(NAME1,NAMEG(JM)).LE.2*NBIN0) THEN
          E(5) = E(5) + ZMU*HM(JM) + DEB
      ELSE
          E(7) = E(7) + ZMU*HM(JM) + DEB
      END IF
*
*       Reduce membership if IM is last (otherwise done in RESET).
      IF (IM.EQ.NMERGE) THEN
          NMERGE = NMERGE - 1
      END IF
*
*       Identify merged ghost particle (single body or binary c.m.).
      JCOMP = -1
      DO 230 J = IFIRST,NTOT
          IF (BODY(J).EQ.0.0D0.AND.NAME(J).EQ.NAMEG(IM)) JCOMP = J
  230 CONTINUE
*       Include search over lower level on failed identification.
      IF (JCOMP.EQ.-1.AND.I3HI.GT.0) THEN
          DO 232 J = IFIRST,NTOT
              IF (BODY(J).EQ.0.0D0.AND.NAME(J).EQ.NAMEG(JM)) JCOMP = J
  232     CONTINUE
      END IF
*
*       Skip if correct ghost not identified (note I3HI # JCOMP if JM # IM).
      IF (JCOMP.GT.0) THEN
          I = I3HI
*       Form two-body elements and period of inner binary.
          NAMEI = 0
          NPOP(7) = NPOP(7) + 1
          BODYCM = CM(1,JM) + CM(2,JM)
          SEMI0 = -0.5*BODYCM/HM(JM)
          ZMU = CM(1,JM)*CM(2,JM)/BODYCM
          EB = ZMU*HM(JM)
          PB = DAYS*SEMI0*SQRT(SEMI0/BODYCM)
          BE(3) = BE(3) - EB
          EMESC = EMESC + EB
          EMERGE = EMERGE - EB
*       Include extra correction for mergers between binary pairs.
          IF (JM.LT.IM) THEN
              ZMU2 = CM(1,IM)*CM(2,IM)/(CM(1,IM) + CM(2,IM))
              EB2 = ZMU2*HM(IM)
              BE(3) = BE(3) - EB2
              EMESC = EMESC + EB2
              EMERGE = EMERGE - EB2
*       Copy smaller index for QUAD binary output.
              IM = JM
          END IF
          RJ = 0.0
          TD2 = 0.0
          DO 235 K = 1,4
              RJ = RJ + UM(K,JM)**2
              TD2 = TD2 + 2.0*UM(K,JM)*UMDOT(K,JM)
  235     CONTINUE
          ECC2 = (1.0 - RJ/SEMI0)**2 + TD2**2/(BODYCM*SEMI0)
          ECC0 = SQRT(ECC2)
          Q = ZM2/BODYCM
          XFAC = (1.0 + Q)*(1.0 + ECC)/SQRT(1.0 - ECC)
          PCRIT = 2.8*XFAC**0.4*SEMI0
*
          if(rank.eq.0)
     &    WRITE (6,240)  NAME1, NAMEG(JM), KSTAR1, KSTAR(JCOMP),
     &                   KSTARM(JM), CM(1,JM)*ZMBAR, CM(2,JM)*ZMBAR,
     &                   EB, ECC0, PMIN/PCRIT, SEMI0, PB
  240     FORMAT (/,' HIARCH ESCAPE    NM =',2I6,'  K* =',3I3,
     &              '  M =',2F5.1,'  EB =',F8.4,'  E =',F7.3,
     &              '  PM/PC =',1P,E8.1,'  A =',E8.1,'  P =',E8.1)
*
*       Remove the ghost particle (NAME = 0 & EB = 0 for second binary).
          GO TO 50
      END IF
*
  250 I = IPAIR + N
      GO TO 12
*
      END

