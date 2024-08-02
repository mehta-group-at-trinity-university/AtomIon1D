c234567890
      program CoupledAdiabaticBoxAvg
      implicit none
      integer iband, iBox, iE, ic,ii,iR,is,ithresh,jj,kk,kl,ku,MatrixDim
      integer LegPoints,xNumPoints
      integer Order,ishift,NumRmatch,NumNewSectors
      integer NumInpChannels,NumChannels,NumOpenChannels,NumSectorsTot
      integer NumOpenL, NumOpenR, NumTotStates,NumSectors,NumOpen
      integer NumDataPoints,NumBoxes,NumE
      integer TotalAngMomentum
      integer, allocatable :: MEMap(:,:),NumberL(:),NumberR(:),NumSectorsBox(:)
      double precision, allocatable :: Thresholds(:),leff(:),xPoints(:)

      double precision rStart,rEnd,xStart,xEnd,xStep,rmsError,rmsErrorP
      double precision rmsErrorMax,rmsErrorPMax,WhereStart,BoxAvg,CrossTot
      double precision, allocatable :: xSector(:),logderiv(:),xSectorTot(:)
      double precision, allocatable :: xData(:),xSave(:)
      double precision, allocatable :: VQmat(:,:,:),Pmat(:,:,:)
      double precision, allocatable :: ub(:,:),up(:,:)
      double precision, allocatable :: PFit(:,:,:),VQFit(:,:,:)
      double precision, allocatable :: PExponents(:,:,:),VQExponents(:,:,:)

      integer i,j,k,ix,i1,i2,HalfBandWidth,lwork,info
      integer NumClosedChannels
      integer iEnergy,NumEnergySteps
      integer xDim,xOpenDim,xClosedDim,LeadDim
      integer ClosedMatrixDim,OpenMatrixDim
      integer, allocatable :: xBounds(:)
      integer, allocatable :: ipiv(:)
      double precision Tol,Einc
      double precision TotalMemory,ddot,Mass,kelvinPERau
      double precision m1,m2,m3,ma,mi,muai
      double precision Energy,Einput,RMatch,RMatch1,f,fp,g,gp,Pi,Wronskian
      double precision, allocatable :: work(:), Egrid(:)
      double precision, allocatable :: xLeg(:),wLeg(:)
      double precision, allocatable :: u(:,:,:),ux(:,:,:),uxx(:,:,:)
      double precision, allocatable :: uTemp(:),uxTemp(:)
      double precision, allocatable :: RData(:),SData(:,:,:,:),HData(:,:,:,:,:)
      double precision, allocatable :: Scc(:,:),Sco(:,:),Soc(:,:),Soo(:,:)
      double precision, allocatable :: Hcc(:,:),Hco(:,:),Hoc(:,:),Hoo(:,:)
      double precision, allocatable :: Gcc(:,:),Gco(:,:),Goc(:,:),Goo(:,:)
      double precision, allocatable :: b(:)

      double precision, allocatable :: psi_in(:,:),solution(:,:),deriv(:)

      double precision, allocatable :: Tmatrix(:,:),Kmatrix(:,:)
      double precision, allocatable :: TmatrixAvg(:,:),KmatrixAvg(:,:)

      double precision, allocatable :: CrossSections(:,:)
      double precision, allocatable :: BoxAvgPartial(:,:)
      double complex, allocatable :: Smatrix(:,:), SmatrixAvg(:,:),SmatrixOld(:,:), TimeDelay(:,:)
      double complex, allocatable :: Test(:,:)
      character*64 LegendreFile,SFile,HFile,PFile,QFile,FitFile

      double precision muK
      double precision sec,time,secp,timep
      double precision step

      double precision EinputAux
      double precision Ri,Rf,Xi,Xf,StepX,R(1:100000000)
      double precision Emin,Emax,wlenght,rstep,XP
      integer NPoints,npwave,NPointsW,NPS
      double precision URmin,wlenghti,rstepi,dumb,ascat,phaseshift,Parity
      integer ncount,ndumb,P
      character(len=20), external :: str
      character*64 GridName

C     TIME INICIALIZATION
      call cpu_time(time)
      SEC = time
      secp = time


      Pi = 3.1415926535897932385d0
      kelvinPERau=315774.65d0
      muK = 3.166815355d-12     !jpdincao ! muK=1.d-6/kelvinPERau

c     read in Gauss-Legendre info
      read(5,*)
      read(5,1002) LegendreFile
      read(5,*)

      read(5,*)
      read(5,*) LegPoints, WhereStart
      read(5,*)

      read(5,*)
      read(5,*) TotalAngMomentum,Emin,Emax,NumE,ithresh,Parity
      read(5,*)

      read(5,*)
      read(5,*) ma,mi
      read(5,*)

      ma = ma/mi
      mi = 1d0
      muai=mi*ma/(mi+ma)        ! ion-atom reduced mass
      Mass = dsqrt(ma/mi)        ! Hyperradial reduced mass

      read(5,*)
      read(5,*) NumChannels,NumOpenChannels,NumBoxes
      read(5,*)

      read(5,*)
      allocate(NumberL(NumBoxes),NumberR(NumBoxes))
      do i = 1,NumBoxes
         read(5,*) NumberL(i),NumberR(i)
      enddo
c      write(6,*) "allocating SmatrixOld, and TimeDelay to be of size", NumberR(NumBoxes)
      allocate(SmatrixOld(NumberR(NumBoxes),NumberR(NumBoxes)))
      allocate(TimeDelay(NumberR(NumBoxes),NumberR(NumBoxes)))
      read(5,*)
      read(5,*)
      read(5,*)NumRmatch,NumNewSectors
c     allocate(leff(NumberL(NumBoxes)),Thresholds(NumberR(NumBoxes)))
      allocate(leff(NumberL(NumBoxes)),Thresholds(NumberL(NumBoxes)))
      read(5,*)
      read(5,*)  
      do i = 1,NumberR(NumBoxes)
         read(5,*) leff(i),Thresholds(i)
      enddo
      Emin = Emin + Thresholds(ithresh)  ! Emin and Emax are read in as collision energies, but should be absolute energies.
      Emax = Emax + Thresholds(ithresh)
      read(5,*)
      read(5,*)
      read(5,1002) PFile
      write(6,*) "PFile = ", PFile
      read(5,*)
      read(5,*)
      read(5,1002) QFile  
      read(5,*)
      read(5,*)
      read(5,1002) FitFile

c     Thresholds(i) = 2B Energies from the Extrapolation Code
      open(18,file=FitFile,status='old')
      read(18,*) NumTotStates,NumDataPoints
      close(18)
      open(19,file='FitLeff.data',status='old') ! Contains only thresholds?
      read(19,*)
      read(19,*)
      read(19,*)
      
      do i=1,NumChannels!NumTotStates
         read(19,*) ndumb, Thresholds(i),dumb,dumb,dumb
      enddo
      close(19)
      allocate(Egrid(NumE))
      call GridMaker(Egrid, NumE, Emin, Emax, 'linear')
      Energy = Egrid(1)!Einput + Thresholds(ithresh)

c     Checking WhereStart
      open(18,file=FitFile,status='old')
      read(18,*)NumTotStates,NumDataPoints
      close(18)

      allocate(VQmat(NumDataPoints,Numchannels,Numchannels))
      allocate(xdata(NumDataPoints))
      open(300,file=QFile,status='old')
      do i=1,NumDataPoints
         read(300,*)xdata(i)
         do j=1,Numchannels
            read(300,11)(VQmat(i,j,k),k=1,Numchannels)
         enddo
         do j=Numchannels+1,NumTotStates
            read(300,11)
         enddo
      enddo
      close(300)
 11   format(1P,100e22.12)
 301  format(100(e18.12,1x))

      if (xdata(NumDataPoints).ne.WhereStart) then
         WhereStart = xdata(NumDataPoints)
      endif

      deallocate(VQmat,xdata)

      read(5,*)
      read(5,*)
      read(5,*) Ri, Rf
      read(5,*)
      read(5,*)
      read(5,*) NPoints
      read(5,*)
      read(5,*)
      read(5,*) npwave

c     RADIAL GRID
      
      Xi = Ri**(1.0d0/3.0d0)
      Xf = WhereStart**(1.0d0/3.0d0)
      StepX = (Xf-Xi)/dfloat(NPoints-1)
      do iR = 1,NPoints
         R(iR) = (Xi+(iR-1)*StepX)**3
      enddo

      Pi = dacos(-1.d0)
      wlenght = 2.d0*Pi/dsqrt(2.d0*Mass*(Emax-Thresholds(1)))
      rstep = wlenght/dfloat(npwave) 

      NPointsW = nint((Rf-WhereStart)/rstep)+1
      do iR = 1,NPointsW
         R(iR+NPoints) = WhereStart+(iR)*rstep
      enddo

      if (R(NPoints)-R(Npoints-1).gt.R(NPoints+1)-R(Npoints)) then
         write(*,*)'Step size bigger than the extrapolation step'
         write(*,*)'... Increase NPoints ...'
         write(*,*)'Program STOPED'
         stop
      endif   

      NumSectorsTot = NPoints+NPointsW

      allocate(xSectorTot(NumSectorsTot))
      do iR = 1,NumSectorsTot
         xSectorTot(iR) = R(iR)
      enddo

      allocate(NumSectorsBox(NumBoxes))
      XP = 1.d0/(NumBoxes-1)
      NPS = int(XP*NumSectorsTot)
      do k = 1,NumBoxes
         NumSectorsBox(k) = NPS
         if (k.eq.(NumBoxes-1)) NumSectorsBox(k) = 0.5*NPS
         if (k.eq.(NumBoxes))   NumSectorsBox(k) = 2
      enddo


      GridName = 'grid.used.data'
      open(unit=17,file=GridName)
      write(17,*)
      write(17,*)'----------------------'      
      write(17,1705)URmin
      write(17,1700)NumSectorsTot
      write(17,1704)NPoints,NPointsW
      write(17,1701)wlenghti,wlenght
      write(17,1703)rstepi,rstep
      write(17,1702)npwave
      write(17,*)'----------------------'      
      do k = 1,NumBoxes
      write(17,171)k,NumSectorsBox(k)   
      enddo
      write(17,*)'----------------------'      
      do iR = 1,NumSectorsTot
      write(17,*)iR,xSectorTot(iR)
      enddo
      write(17,*)'----------------------'      
      close(17)

 1700 format(1x,'NumSectorsTot =',I8)
 1701 format(1x,'Wave Lenght =',f14.8',',f14.8)
 1702 format(1x,'Num. Points per Wave Lenght =',I3)
 1703 format(1x,'Step =',f14.8',',f14.8)
 1704 format(1x,'NPoints =',I6,',',I6)
 1705 format(1x,'Potential Minimum =',f14.8)
 171  format(1x,'Box(',I3,')=',I8)

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      allocate(xLeg(LegPoints),wLeg(LegPoints))
      call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)
      
      allocate(VQExponents(Numchannels,Numchannels,3))
      allocate(VQFit(Numchannels,Numchannels,3))
      allocate(PExponents(Numchannels,Numchannels,3))
      allocate(PFit(Numchannels,Numchannels,3))
      
      open(18,file=FitFile,status='old') ! Contains fit parameters for VQ and for P?
      read(18,*) NumTotStates,NumDataPoints
      rmsErrorMax=0.d0
      rmsErrorPMax=0.d0
      do j=1,NumChannels
         do k=j,NumChannels
            read(18,*) jj,kk,(VQExponents(jj,kk,ii),ii=1,3)
            read(18,*) (VQFit(jj,kk,ii),ii=1,3),rmsError
            IF(j.NE.k) THEN
               read(18,*) jj,kk,(PExponents(jj,kk,ii),ii=1,3)
               read(18,*) (PFit(jj,kk,ii),ii=1,3),rmsErrorP
            ENDIF
            IF(jj.NE.j .OR. kk.NE.k) THEN
               write(6,*) 'index mismatch!',j,jj,k,kk
            ENDIF
            DO i=1,3
               ! enforce symmetry in VQ and antisymmetry in P
               VQFit(k,j,i)=VQFit(j,k,i)
               VQExponents(k,j,i)=VQExponents(j,k,i)
               PFit(k,j,i)=-PFit(j,k,i)
               PExponents(k,j,i)=PExponents(j,k,i)
               IF(j.EQ.k) THEN
                  PFit(k,k,i)=0.d0
                  PExponents(k,k,i)=1
               ENDIF
            ENDDO
            IF(rmsError .GT. rmsErrorMax) THEN
               rmsErrorMax=dabs(rmsError)
            ENDIF
            IF(rmsErrorP.GT. rmsErrorPMax) THEN
               rmsErrorPMax=dabs(rmsErrorP)
            ENDIF
         ENDDO
         do k=Numchannels+1,NumTotStates
            read(18,*)
            read(18,*)
            IF(j.NE.k) THEN
               read(18,*)
               read(18,*)
            ENDIF
         enddo
      ENDDO
      close(18)

      allocate(xData(NumDataPoints))
      allocate(VQmat(NumDataPoints,Numchannels,Numchannels))
      allocate(Pmat(NumDataPoints,Numchannels,Numchannels))
      
      call potential(NumChannels,NumDataPoints,NumTotStates,
     >     PFile,QFile,xdata,VQmat,Pmat) ! read in the P, VQ matrix data to be interpolated later 


      open(unit = 25, file = "Observables"//trim(str(nint(Parity)))//".dat")
      open(unit = 26, file = "Kmatrix"//trim(str(nint(Parity)))//".dat")
      open(unit = 27, file = "Smatrix"//trim(str(nint(Parity)))//".dat")
      open(unit = 28, file = "TimeDelay"//trim(str(nint(Parity)))//".dat")
      
      do ie=1,NumE
         energy = Egrid(iE)

         ishift = 0
         do ibox = 1,NumBoxes-1

            NumSectors = NumSectorsBox(ibox)
            NumOpenL = NumberL(ibox)
            NumOpenR = NumberR(ibox)
            NumOpen = NumOpenL + NumOpenR
            MatrixDim = Numchannels*NumSectors*4
            kl = 6*Numchannels - 1
            ku = kl
            iband = 3*kl + 1
            allocate(xSector(NumSectors+1))
            
            if (ibox .gt. 1) ishift =  ishift + NumSectorsBox(ibox-1)
            do is=1,NumSectors+1
               xSector(is) = xSectorTot(ishift+is)
            enddo
            Rstart = xSector(1)
            Rmatch = xSector(NumSectors+1)
            
            allocate(MEMap(Numchannels*NumSectors,6))
            call GenMEMap(Numchannels,Numsectors,MEMap)
            
            allocate(Gcc(iband,MatrixDim))
            call CalcClosed(Numchannels,NumSectors,MatrixDim,kl,iband,MEMap,Mass,energy,
     >           xsector,LegPoints,xLeg,wLeg,NumDataPoints,xdata,VQmat,Pmat,PFit,VQFit,
     >           WhereStart,PExponents,VQExponents,Gcc)
            allocate(Gco(MatrixDim,NumOpen),Goo(NumOpen,NumOpen))
            call CalcOpen(NumOpenL,NumOpenR,Numchannels,NumSectors,MatrixDim,MEMap,Mass,
     >           energy,xsector,LegPoints,xLeg,wLeg,NumDataPoints,xdata,VQmat,Pmat,
     >           PFit,VQFit,WhereStart,PExponents,VQExponents,Gco,Goo)
            allocate(Goc(NumOpen,MatrixDim))
            
            do i=1,NumOpen
               do j=1,MatrixDim
                  Goc(i,j) = Gco(j,i)
               enddo
            enddo
            
c     Lapack Banded Linear Equation Solver
            
            allocate(ipiv(MatrixDim))
            
            call dgbsv(MatrixDim,kl,ku,NumOpen,Gcc,iband,ipiv,Gco,
     >           MatrixDim,info)
            
            call dgemm('N','N',NumOpen,NumOpen,MatrixDim,-1.0d0,Goc,NumOpen,
     >           Gco,MatrixDim,1.0d0,Goo,NumOpen)
            
            lwork = 3*NumOpen
            allocate(work(lwork)) 
            allocate(logderiv(NumOpen))
            
            call dsyev('V','U',NumOpen,Goo,NumOpen,logderiv,work,lwork,info)

            deallocate(xSector,work,ipiv)
            deallocate(MEMap,Gcc,Gco,Goc)
            
c     starting matching section
            if (ibox .eq. 1)then
               allocate(psi_in(NumOpenR,NumOpenR),b(NumOpenR))
               psi_in = Goo
               b = logderiv
               deallocate(Goo,logderiv)
            else
!     write(6,*)'starting matching section'

               allocate(solution(NumOpenR,NumOpenR))
               allocate(deriv(NumOpenR))
               call BoxMatching(NumOpenL,NumOpenR,psi_in,b,Goo,logderiv,
     >              solution,deriv)

               deallocate(psi_in,b,Goo,logderiv)
               allocate(psi_in(NumOpenR,NumOpenR),b(NumOpenR))
               psi_in = solution
               b = deriv
               deallocate(solution,deriv)
            endif

         enddo
         
         NumOpenR = NumberR(NumBoxes) ! Added this 5/28/2024 NPM        
         allocate(BoxAvgPartial(NumberR(NumBoxes),NumberR(NumBoxes)))
         allocate(TmatrixAvg(NumOpenR,NumOpenR))
         allocate(KmatrixAvg(NumOpenR,NumOpenR))
         allocate(SmatrixAvg(NumOpenR,NumOpenR))
         allocate(test(NumOpenR,NumOpenR))
         
         BoxAvg = 0.0d0
         BoxAvgPartial = 0.0d0
         TmatrixAvg = 0.0d0
         KmatrixAvg = 0.0d0
         SmatrixAvg = 0d0
         TimeDelay = 0d0
         
         ! Now do the final box, averaging over the NumRmatch sectors. (NPM)
         do ir=1,NumRmatch

            NumSectors = NumSectorsBox(NumBoxes) + (ir-1)*NumNewSectors
            allocate(xSector(NumSectors+1))
            if(ir .eq. 1 .and. NumBoxes .eq. 1) ishift = 0 
            if(ir .eq. 1 .and. NumBoxes .gt. 1) ishift = ishift + NumSectorsBox(NumBoxes-1)
            do is=1,NumSectors+1
               xSector(is) = xSectorTot(ishift+is)
            enddo
            
            Rstart = xSector(1)
            Rmatch = xSector(NumSectors+1)

            NumOpenL = NumberL(NumBoxes)
            NumOpenR = NumberR(NumBoxes)
            NumOpen = NumOpenL + NumOpenR
            MatrixDim = Numchannels*NumSectors*4
            kl = 6*Numchannels - 1
            ku = kl
            iband = 3*kl + 1
            
            allocate(MEMap(Numchannels*NumSectors,6))
            call GenMEMap(Numchannels,Numsectors,MEMap)
            
            allocate(Gcc(iband,MatrixDim))
            
!     write(6,*)'starting CalcClosed'
            call CalcClosed(Numchannels,NumSectors,MatrixDim,kl,iband,MEMap,Mass,energy,
     >           xsector,LegPoints,xLeg,wLeg,NumDataPoints,xdata,VQmat,Pmat,PFit,VQFit,
     >           WhereStart,PExponents,VQExponents,Gcc)
            
            allocate(Gco(MatrixDim,NumOpen),Goo(NumOpen,NumOpen))
            
            call CalcOpen(NumOpenL,NumOpenR,Numchannels,NumSectors,MatrixDim,MEMap,Mass,
     >           energy,xsector,LegPoints,xLeg,wLeg,NumDataPoints,xdata,VQmat,Pmat,
     >           PFit,VQFit,WhereStart,PExponents,VQExponents,Gco,Goo)
           
            allocate(Goc(NumOpen,MatrixDim))
            
            do i=1,NumOpen
               do j=1,MatrixDim
                  Goc(i,j) = Gco(j,i)
               enddo
            enddo
            
c     Lapack Banded Linear Equation Solver
            
            allocate(ipiv(MatrixDim))
            
            call dgbsv(MatrixDim,kl,ku,NumOpen,Gcc,iband,ipiv,Gco,
     >           MatrixDim,info)
            
            call dgemm('N','N',NumOpen,NumOpen,MatrixDim,-1.0d0,Goc,NumOpen,
     >           Gco,MatrixDim,1.0d0,Goo,NumOpen)
            
            lwork = 3*NumOpen
            allocate(work(lwork))
            allocate(logderiv(NumOpen))
            
            call dsyev('V','U',NumOpen,Goo,NumOpen,logderiv,work,lwork,info)
            
            deallocate(xSector,work,ipiv)
            deallocate(MEMap,Gcc,Gco,Goc)
            
c     starting matching section
!     write(6,*)'starting matching section'
            
            allocate(solution(NumOpenR,NumOpenR))
            allocate(deriv(NumOpenR))
            call BoxMatching(NumOpenL,NumOpenR,psi_in,b,Goo,logderiv,
     >           solution,deriv)

            deallocate(Goo,logderiv)

c     calculate S-matrix and scattering cross sections

            allocate(Smatrix(NumOpenR,NumOpenR),Tmatrix(NumOpenR,NumOpenR),Kmatrix(NumOpenR,NumOpenR))
            allocate(CrossSections(NumOpenR,NumOpenR))

            call CalcSmatrix1D(NumOpenR,leff,Thresholds,TotalAngMomentum,Mass,Energy,RMatch,
     >           solution,deriv,Smatrix,Tmatrix,CrossSections,Kmatrix,Parity)

            do i = 1,NumOpenR
               do j = 1,NumOpenR
                  BoxAvgPartial(i,j) = BoxAvgPartial(i,j) + CrossSections(i,j)
                  SmatrixAvg(i,j) = SmatrixAvg(i,j) + Smatrix(i,j)
                  TmatrixAvg(i,j) = TmatrixAvg(i,j) + Tmatrix(i,j)
                  KmatrixAvg(i,j) = KmatrixAvg(i,j) + Kmatrix(i,j)
               enddo
            enddo

            do ic=1,ithresh-1
               Crosstot = Crosstot + CrossSections(ic,ithresh)
            enddo
            Boxavg = Boxavg + Crosstot
                        
            deallocate (Smatrix,Tmatrix,CrossSections,Kmatrix)
            deallocate(solution,deriv)
            
         enddo
         
         do i = 1,NumberR(NumBoxes)
            do j = 1,NumberR(NumBoxes)
               BoxAvgPartial(i,j) = BoxAvgPartial(i,j)/dfloat(NumRmatch)
               SmatrixAvg(i,j) = SmatrixAvg(i,j)/dfloat(NumRmatch)
               TmatrixAvg(i,j) = TmatrixAvg(i,j)/dfloat(NumRmatch)
               KmatrixAvg(i,j) = KmatrixAvg(i,j)/dfloat(NumRmatch)
            enddo
         enddo
         
         Timedelay = (0d0,0d0)
         if(iE.gt.1) then
            
            Timedelay = (SmatrixAvg - SmatrixOld)/(Egrid(iE)-Egrid(iE-1))            
            call zgemm('N', 'C', NumOpenR, NumOpenR, NumOpenR, -(0d0,1d0), SmatrixAvg,
     .           NumOpenR, TimeDelay, NumOpenR, (0d0,0d0), TimeDelay, NumOpenR)
            do i = 1, NumOpenR
               do j = 1, NumOpenR
                  test(i,j) = 0.5d0*(TimeDelay(i,j) + conjg(TimeDelay(j,i)))
               enddo
            enddo
            TimeDelay = test
            
c            call zgemm('C', 'N', NumOpenR, NumOpenR, NumOpenR, (1d0,0d0), SmatrixOld,
c     .           NumOpenR, SmatrixOld, NumOpenR, (0d0,0d0), test, NumOpenR) ! Check the unitarity of S
c     write(6,*) test
c     call zprintmatrix(test,NumOpenR,NumOpenR,6)
         endif
         SmatrixOld = SmatrixAvg

         BoxAvg = BoxAvg/float(NumRmatch)

         call cpu_time(timep)
         P = nint((1d0 + (-1)**ithresh*Parity)/2)
         if(P.eq.1) then
            ascat = 1d0/(dsqrt(2.0d0*Mass*(Energy-Thresholds(ithresh)))*KmatrixAvg(ithresh,ithresh)) !a=1/(k tandelta) for P=1
         else if(P.eq.0) then
            ascat = -KmatrixAvg(ithresh,ithresh)/(dsqrt(2.0d0*Mass*(Energy-Thresholds(ithresh)))) ! a = -k/tandel for P=0
         endif
         phaseshift = atan(KmatrixAvg(ithresh,ithresh))

         write(6,21) 'min. to go = ',((timep-secp)*float(NumE-iE))/(60d0), (energy-Thresholds(ithresh)), ascat, phaseshift,
     .        BoxAvgPartial(1,1)     
         write(25,20) energy, ascat, ((BoxAvgPartial(i,j), i=1,NumOpenR), j=1,NumOpenR)
         write(26,20) energy, ((KmatrixAvg(i,j), i=1,NumOpenR), j=1,NumOpenR)
         write(27,20) energy, ((dreal(conjg(SmatrixAvg(i,j))*SmatrixAvg(i,j)), i=1,NumOpenR), j=1,NumOpenR)
         write(28,20) energy, ((TimeDelay(i,j), i=1,NumOpenR), j=1,NumOpenR)

         call cpu_time(timep)
         secp = timep 

         deallocate(BoxAvgPartial,TmatrixAvg,KmatrixAvg,SmatrixAvg,test)
         deallocate(psi_in,b)

      enddo

 1234 continue

      deallocate(xData,VQmat,Pmat,Pfit,VQfit)
      deallocate(VQexponents,Pexponents)
      deallocate(xLeg,wLeg)

 20   format(1P,100e16.8)
 21   format(a13,f7.2,100e16.8)
 28   format(f10.2,1P,100e14.5)
 407  format(10(e16.8,1x))
 1002 format(a64)

      call cpu_time(time)
      WRITE(*,*)'>>> EXECUTION TIME =',time-SEC

      stop
      end

c******************************************************************
      subroutine GetGaussFactors(File,Points,x,w)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine retrieves the points and weights for
c      Gaussian Quadrature from a given file
c
c     Variables:
c      File		name of file containing points and 
c			 weights
c      Points		number of points to be used for 
c			quadrature
c      x		array of nodes
c      w		array of weights
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      integer Points
      double precision x(Points),w(Points)
      character*64 File
      
      integer i,tempPoints
      
      open(unit=7,file=File(1:index(File,' ')-1))

       do i = 1,18
        read(7,*)
       enddo
       read(7,*) tempPoints
       do while (tempPoints .ne. Points)
        do i = 1,tempPoints
         read(7,*)
        enddo
        read(7,*) tempPoints
       enddo
 
       do i = 1,Points
        read(7,*) x(i),w(i)
       enddo
      
      close(unit=7)

      return
      end
      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine BesselBasePair(k,r,l,f,fp,g,gp)

      double precision l,k,r,f,fp,g,gp

      integer lInt,info
      double precision x,a,b,alpha,alphaComp,Pi
      double precision pj,py,mj,my
      double precision, allocatable :: j(:),y(:)

      Pi = 3.1415926535897932385d0

      x = k*r
      a = dsqrt(r)
      b = 1.0d0/a
      lInt = nint(l+0.5d0)
      alpha = l+0.5d0-dfloat(lInt)
      alphaComp = dabs(alpha-1.0d0)

      if (lInt .eq. 0) then

       allocate(j(2),y(2))

       call RJBESL(x,alpha,2,j,info)
       call RYBESL(x,alpha,2,y,info)
       pj = j(1)
       py = y(1)

       call RJBESL(x,alphaComp,2,j,info)
       call RYBESL(x,alphaComp,2,y,info)
       mj = dcos(Pi*alphaComp)*j(1)-dsin(Pi*alphaComp)*y(1)
       my = dsin(Pi*alphaComp)*j(1)+dcos(Pi*alphaComp)*y(1)

       f = a*pj
       g = a*py

       fp = b*(x*mj-l*pj)
       gp = b*(x*my-l*py)

       deallocate(j,y)

      else

       allocate(j(lInt+1),y(lInt+1))

       call RJBESL(x,alpha,lInt+1,j,info)
       call RYBESL(x,alpha,lInt+1,y,info)

       f = a*j(lInt+1)
       g = a*y(lInt+1)

       fp = b*(x*j(lInt)-l*j(lInt+1))
       gp = b*(x*y(lInt)-l*y(lInt+1))

       deallocate(j,y)

      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE RJBESL(X, ALPHA, NB, B, NCALC)
C---------------------------------------------------------------------
C This routine calculates Bessel functions J sub(N+ALPHA) (X)
C   for non-negative argument X, and non-negative order N+ALPHA.
C
C
C  Explanation of variables in the calling sequence.
C
C   X     - working precision non-negative real argument for which
C           J's are to be calculated.
C   ALPHA - working precision fractional part of order for which
C           J's or exponentially scaled J'r (J*exp(X)) are
C           to be calculated.  0 <= ALPHA < 1.0.
C   NB  - integer number of functions to be calculated, NB > 0.
C           The first function calculated is of order ALPHA, and the
C           last is of order (NB - 1 + ALPHA).
C   B  - working precision output vector of length NB.  If RJBESL
C           terminates normally (NCALC=NB), the vector B contains the
C           functions J/ALPHA/(X) through J/NB-1+ALPHA/(X), or the
C           corresponding exponentially scaled functions.
C   NCALC - integer output variable indicating possible errors.
C           Before using the vector B, the user should check that
C           NCALC=NB, i.e., all orders have been calculated to
C           the desired accuracy.  See Error Returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C  Explanation of machine-dependent constants
C
C   it     = Number of bits in the mantissa of a working precision
C            variable
C   NSIG   = Decimal significance desired.  Should be set to
C            INT(LOG10(2)*it+1).  Setting NSIG lower will result
C            in decreased accuracy while setting NSIG higher will
C            increase CPU time without increasing accuracy.  The
C            truncation error is limited to a relative error of
C            T=.5*10**(-NSIG).
C   ENTEN  = 10.0 ** K, where K is the largest integer such that
C            ENTEN is machine-representable in working precision
C   ENSIG  = 10.0 ** NSIG
C   RTNSIG = 10.0 ** (-K) for the smallest integer K such that
C            K .GE. NSIG/4
C   ENMTEN = Smallest ABS(X) such that X/4 does not underflow
C   XLARGE = Upper limit on the magnitude of X.  If ABS(X)=N,
C            then at least N iterations of the backward recursion
C            will be executed.  The value of 10.0 ** 4 is used on
C            every machine.
C
C
C     Approximate values for some important machines are:
C
C
C                            it    NSIG    ENTEN       ENSIG
C
C   CRAY-1        (S.P.)     48     15    1.0E+2465   1.0E+15
C   Cyber 180/855
C     under NOS   (S.P.)     48     15    1.0E+322    1.0E+15
C   IEEE (IBM/XT,
C     SUN, etc.)  (S.P.)     24      8    1.0E+38     1.0E+8
C   IEEE (IBM/XT,
C     SUN, etc.)  (D.P.)     53     16    1.0D+308    1.0D+16
C   IBM 3033      (D.P.)     14      5    1.0D+75     1.0D+5
C   VAX           (S.P.)     24      8    1.0E+38     1.0E+8
C   VAX D-Format  (D.P.)     56     17    1.0D+38     1.0D+17
C   VAX G-Format  (D.P.)     53     16    1.0D+307    1.0D+16
C
C
C                           RTNSIG      ENMTEN      XLARGE
C
C   CRAY-1        (S.P.)    1.0E-4    1.84E-2466   1.0E+4
C   Cyber 180/855
C     under NOS   (S.P.)    1.0E-4    1.25E-293    1.0E+4
C   IEEE (IBM/XT,
C     SUN, etc.)  (S.P.)    1.0E-2    4.70E-38     1.0E+4
C   IEEE (IBM/XT,
C     SUN, etc.)  (D.P.)    1.0E-4    8.90D-308    1.0D+4
C   IBM 3033      (D.P.)    1.0E-2    2.16D-78     1.0D+4
C   VAX           (S.P.)    1.0E-2    1.17E-38     1.0E+4
C   VAX D-Format  (D.P.)    1.0E-5    1.17D-38     1.0D+4
C   VAX G-Format  (D.P.)    1.0E-4    2.22D-308    1.0D+4
C
C*******************************************************************
C*******************************************************************
C
C  Error returns
C
C    In case of an error,  NCALC .NE. NB, and not all J's are
C    calculated to the desired accuracy.
C
C    NCALC .LT. 0:  An argument is out of range. For example,
C       NBES .LE. 0, ALPHA .LT. 0 or .GT. 1, or X is too large.
C       In this case, B(1) is set to zero, the remainder of the
C       B-vector is not calculated, and NCALC is set to
C       MIN(NB,0)-1 so that NCALC .NE. NB.
C
C    NB .GT. NCALC .GT. 0: Not all requested function values could
C       be calculated accurately.  This usually occurs because NB is
C       much larger than ABS(X).  In this case, B(N) is calculated
C       to the desired accuracy for N .LE. NCALC, but precision
C       is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish
C       for N .GT. NCALC (because it is too small to be represented),
C       and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
C       significant figures of B(N) can be trusted.
C
C
C  Intrinsic and other functions required are:
C
C     ABS, AINT, COS, DBLE, GAMMA (or DGAMMA), INT, MAX, MIN,
C
C     REAL, SIN, SQRT
C
C
C  Acknowledgement
C
C   This program is based on a program written by David J. Sookne
C   (2) that computes values of the Bessel functions J or I of real
C   argument and integer order.  Modifications include the restriction
C   of the computation to the J Bessel function of non-negative real
C   argument, the extension of the computation to arbitrary positive
C   order, and the elimination of most underflow.
C
C  References: "A Note on Backward Recurrence Algorithms," Olver,
C               F. W. J., and Sookne, D. J., Math. Comp. 26, 1972,
C               pp 941-947.
C
C              "Bessel Functions of Real Argument and Integer Order,"
C               Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp
C               125-132.
C
C  Latest modification: March 19, 1990
C
C  Author: W. J. Cody
C          Applied Mathematics Division
C          Argonne National Laboratory
C          Argonne, IL  60439
C
C---------------------------------------------------------------------
      INTEGER I,J,K,L,M,MAGX,N,NB,NBMX,NCALC,NEND,NSTART
CS    REAL               GAMMA,
      DOUBLE PRECISION  DGAMMA,
     1 ALPHA,ALPEM,ALP2EM,B,CAPP,CAPQ,CONV,EIGHTH,EM,EN,ENMTEN,ENSIG,
     2 ENTEN,FACT,FOUR,FUNC,GNU,HALF,HALFX,ONE,ONE30,P,PI2,PLAST,
     3 POLD,PSAVE,PSAVEL,RTNSIG,S,SUM,T,T1,TEMPA,TEMPB,TEMPC,TEST,
     4 THREE,THREE5,TOVER,TWO,TWOFIV,TWOPI1,TWOPI2,X,XC,XIN,XK,XLARGE,
     5 XM,VCOS,VSIN,Z,ZERO
      DIMENSION B(NB), FACT(25)
C---------------------------------------------------------------------
C  Mathematical constants
C
C   PI2    - 2 / PI
C   TWOPI1 - first few significant digits of 2 * PI
C   TWOPI2 - (2*PI - TWOPI) to working precision, i.e.,
C            TWOPI1 + TWOPI2 = 2 * PI to extra precision.
C---------------------------------------------------------------------
CS    DATA PI2, TWOPI1, TWOPI2 /0.636619772367581343075535E0,6.28125E0,
CS   1 1.935307179586476925286767E-3/
CS    DATA ZERO, EIGHTH, HALF, ONE /0.0E0,0.125E0,0.5E0,1.0E0/
CS    DATA TWO, THREE, FOUR, TWOFIV /2.0E0,3.0E0,4.0E0,25.0E0/
CS    DATA ONE30, THREE5 /130.0E0,35.0E0/
      DATA PI2, TWOPI1, TWOPI2 /0.636619772367581343075535D0,6.28125D0,
     1 1.935307179586476925286767D-3/
      DATA ZERO, EIGHTH, HALF, ONE /0.0D0,0.125D0,0.5D0,1.0D0/
      DATA TWO, THREE, FOUR, TWOFIV /2.0D0,3.0D0,4.0D0,25.0D0/
      DATA ONE30, THREE5 /130.0D0,35.0D0/
C---------------------------------------------------------------------
C  Machine-dependent parameters
C---------------------------------------------------------------------
CS    DATA ENTEN, ENSIG, RTNSIG /1.0E38,1.0E8,1.0E-2/
CS    DATA ENMTEN, XLARGE /1.2E-37,1.0E4/
      DATA ENTEN, ENSIG, RTNSIG /1.0D308,1.0D16,1.0D-4/
      DATA ENMTEN, XLARGE /8.9D-308,1.0D4/
C---------------------------------------------------------------------
C     Factorial(N)
C---------------------------------------------------------------------
CS    DATA FACT /1.0E0,1.0E0,2.0E0,6.0E0,24.0E0,1.2E2,7.2E2,5.04E3,
CS   1 4.032E4,3.6288E5,3.6288E6,3.99168E7,4.790016E8,6.2270208E9,
CS   2 8.71782912E10,1.307674368E12,2.0922789888E13,3.55687428096E14,
CS   3 6.402373705728E15,1.21645100408832E17,2.43290200817664E18,
CS   4 5.109094217170944E19,1.12400072777760768E21,
CS   5 2.585201673888497664E22,6.2044840173323943936E23/
      DATA FACT /1.0D0,1.0D0,2.0D0,6.0D0,24.0D0,1.2D2,7.2D2,5.04D3,
     1 4.032D4,3.6288D5,3.6288D6,3.99168D7,4.790016D8,6.2270208D9,
     2 8.71782912D10,1.307674368D12,2.0922789888D13,3.55687428096D14,
     3 6.402373705728D15,1.21645100408832D17,2.43290200817664D18,
     4 5.109094217170944D19,1.12400072777760768D21,
     5 2.585201673888497664D22,6.2044840173323943936D23/
C---------------------------------------------------------------------
C Statement functions for conversion and the gamma function.
C---------------------------------------------------------------------
CS    CONV(I) = REAL(I)
CS    FUNC(X) = GAMMA(X)
      CONV(I) = DBLE(I)
      FUNC(X) = DGAMMA(X)
C---------------------------------------------------------------------
C Check for out of range arguments.
C---------------------------------------------------------------------
      MAGX = INT(X)
      IF ((NB.GT.0) .AND. (X.GE.ZERO) .AND. (X.LE.XLARGE) 
     1       .AND. (ALPHA.GE.ZERO) .AND. (ALPHA.LT.ONE))  
     2   THEN
C---------------------------------------------------------------------
C Initialize result array to zero.
C---------------------------------------------------------------------
            NCALC = NB
            DO 20 I=1,NB
              B(I) = ZERO
   20       CONTINUE
C---------------------------------------------------------------------
C Branch to use 2-term ascending series for small X and asymptotic
C form for large X when NB is not too large.
C---------------------------------------------------------------------
            IF (X.LT.RTNSIG) THEN
C---------------------------------------------------------------------
C Two-term ascending series for small X.
C---------------------------------------------------------------------
               TEMPA = ONE
               ALPEM = ONE + ALPHA
               HALFX = ZERO
               IF (X.GT.ENMTEN) HALFX = HALF*X
               IF (ALPHA.NE.ZERO)
     1            TEMPA = HALFX**ALPHA/(ALPHA*FUNC(ALPHA))
               TEMPB = ZERO
               IF ((X+ONE).GT.ONE) TEMPB = -HALFX*HALFX
               B(1) = TEMPA + TEMPA*TEMPB/ALPEM
               IF ((X.NE.ZERO) .AND. (B(1).EQ.ZERO)) NCALC = 0
               IF (NB .NE. 1) THEN
                  IF (X .LE. ZERO) THEN
                        DO 30 N=2,NB
                          B(N) = ZERO
   30                   CONTINUE
                     ELSE
C---------------------------------------------------------------------
C Calculate higher order functions.
C---------------------------------------------------------------------
                        TEMPC = HALFX
                        TOVER = (ENMTEN+ENMTEN)/X
                        IF (TEMPB.NE.ZERO) TOVER = ENMTEN/TEMPB
                        DO 50 N=2,NB
                          TEMPA = TEMPA/ALPEM
                          ALPEM = ALPEM + ONE
                          TEMPA = TEMPA*TEMPC
                          IF (TEMPA.LE.TOVER*ALPEM) TEMPA = ZERO
                          B(N) = TEMPA + TEMPA*TEMPB/ALPEM
                          IF ((B(N).EQ.ZERO) .AND. (NCALC.GT.N))
     1                       NCALC = N-1
   50                   CONTINUE
                  END IF
               END IF
            ELSE IF ((X.GT.TWOFIV) .AND. (NB.LE.MAGX+1)) THEN
C---------------------------------------------------------------------
C Asymptotic series for X .GT. 21.0.
C---------------------------------------------------------------------
               XC = SQRT(PI2/X)
               XIN = (EIGHTH/X)**2
               M = 11
               IF (X.GE.THREE5) M = 8
               IF (X.GE.ONE30) M = 4
               XM = FOUR*CONV(M)
C---------------------------------------------------------------------
C Argument reduction for SIN and COS routines.
C---------------------------------------------------------------------
               T = AINT(X/(TWOPI1+TWOPI2)+HALF)
               Z = ((X-T*TWOPI1)-T*TWOPI2) - (ALPHA+HALF)/PI2
               VSIN = SIN(Z)
               VCOS = COS(Z)
               GNU = ALPHA + ALPHA
               DO 80 I=1,2
                 S = ((XM-ONE)-GNU)*((XM-ONE)+GNU)*XIN*HALF
                 T = (GNU-(XM-THREE))*(GNU+(XM-THREE))
                 CAPP = S*T/FACT(2*M+1)
                 T1 = (GNU-(XM+ONE))*(GNU+(XM+ONE))
                 CAPQ = S*T1/FACT(2*M+2)
                 XK = XM
                 K = M + M
                 T1 = T
                 DO 70 J=2,M
                   XK = XK - FOUR
                   S = ((XK-ONE)-GNU)*((XK-ONE)+GNU)
                   T = (GNU-(XK-THREE))*(GNU+(XK-THREE))
                   CAPP = (CAPP+ONE/FACT(K-1))*S*T*XIN
                   CAPQ = (CAPQ+ONE/FACT(K))*S*T1*XIN
                   K = K - 2
                   T1 = T
   70            CONTINUE
                 CAPP = CAPP + ONE
                 CAPQ = (CAPQ+ONE)*(GNU*GNU-ONE)*(EIGHTH/X)
                 B(I) = XC*(CAPP*VCOS-CAPQ*VSIN)
                 IF (NB.EQ.1) GO TO 300
                 T = VSIN
                 VSIN = -VCOS
                 VCOS = T
                 GNU = GNU + TWO
   80         CONTINUE
C---------------------------------------------------------------------
C If  NB .GT. 2, compute J(X,ORDER+I)  I = 2, NB-1
C---------------------------------------------------------------------
               IF (NB .GT. 2) THEN
                  GNU = ALPHA + ALPHA + TWO
                  DO 90 J=3,NB
                    B(J) = GNU*B(J-1)/X - B(J-2)
                    GNU = GNU + TWO
   90             CONTINUE
               END IF
C---------------------------------------------------------------------
C Use recurrence to generate results.  First initialize the
C calculation of P*S.
C---------------------------------------------------------------------
            ELSE
               NBMX = NB - MAGX
               N = MAGX + 1
               EN = CONV(N+N) + (ALPHA+ALPHA)
               PLAST = ONE
               P = EN/X
C---------------------------------------------------------------------
C Calculate general significance test.
C---------------------------------------------------------------------
               TEST = ENSIG + ENSIG
               IF (NBMX .GE. 3) THEN
C---------------------------------------------------------------------
C Calculate P*S until N = NB-1.  Check for possible overflow.
C---------------------------------------------------------------------
                  TOVER = ENTEN/ENSIG
                  NSTART = MAGX + 2
                  NEND = NB - 1
                  EN = CONV(NSTART+NSTART) - TWO + (ALPHA+ALPHA)
                  DO 130 K=NSTART,NEND
                     N = K
                     EN = EN + TWO
                     POLD = PLAST
                     PLAST = P
                     P = EN*PLAST/X - POLD
                     IF (P.GT.TOVER) THEN
C---------------------------------------------------------------------
C To avoid overflow, divide P*S by TOVER.  Calculate P*S until
C ABS(P) .GT. 1.
C---------------------------------------------------------------------
                        TOVER = ENTEN
                        P = P/TOVER
                        PLAST = PLAST/TOVER
                        PSAVE = P
                        PSAVEL = PLAST
                        NSTART = N + 1
  100                   N = N + 1
                           EN = EN + TWO
                           POLD = PLAST
                           PLAST = P
                           P = EN*PLAST/X - POLD
                        IF (P.LE.ONE) GO TO 100
                        TEMPB = EN/X
C---------------------------------------------------------------------
C Calculate backward test and find NCALC, the highest N such that
C the test is passed.
C---------------------------------------------------------------------
                        TEST = POLD*PLAST*(HALF-HALF/(TEMPB*TEMPB))
                        TEST = TEST/ENSIG
                        P = PLAST*TOVER
                        N = N - 1
                        EN = EN - TWO
                        NEND = MIN(NB,N)
                        DO 110 L=NSTART,NEND
                           POLD = PSAVEL
                           PSAVEL = PSAVE
                           PSAVE = EN*PSAVEL/X - POLD
                           IF (PSAVE*PSAVEL.GT.TEST) THEN
                              NCALC = L - 1
                              GO TO 190
                           END IF
  110                   CONTINUE
                        NCALC = NEND
                        GO TO 190
                     END IF
  130             CONTINUE
                  N = NEND
                  EN = CONV(N+N) + (ALPHA+ALPHA)
C---------------------------------------------------------------------
C Calculate special significance test for NBMX .GT. 2.
C---------------------------------------------------------------------
                  TEST = MAX(TEST,SQRT(PLAST*ENSIG)*SQRT(P+P))
               END IF
C---------------------------------------------------------------------
C Calculate P*S until significance test passes.
C---------------------------------------------------------------------
  140          N = N + 1
                  EN = EN + TWO
                  POLD = PLAST
                  PLAST = P
                  P = EN*PLAST/X - POLD
               IF (P.LT.TEST) GO TO 140
C---------------------------------------------------------------------
C Initialize the backward recursion and the normalization sum.
C---------------------------------------------------------------------
  190          N = N + 1
               EN = EN + TWO
               TEMPB = ZERO
               TEMPA = ONE/P
               M = 2*N - 4*(N/2)
               SUM = ZERO
               EM = CONV(N/2)
               ALPEM = (EM-ONE) + ALPHA
               ALP2EM = (EM+EM) + ALPHA
               IF (M .NE. 0) SUM = TEMPA*ALPEM*ALP2EM/EM
               NEND = N - NB
               IF (NEND .GT. 0) THEN
C---------------------------------------------------------------------
C Recur backward via difference equation, calculating (but not
C storing) B(N), until N = NB.
C---------------------------------------------------------------------
                  DO 200 L=1,NEND
                     N = N - 1
                     EN = EN - TWO
                     TEMPC = TEMPB
                     TEMPB = TEMPA
                     TEMPA = (EN*TEMPB)/X - TEMPC
                     M = 2 - M
                     IF (M .NE. 0) THEN
                        EM = EM - ONE
                        ALP2EM = (EM+EM) + ALPHA
                        IF (N.EQ.1) GO TO 210
                        ALPEM = (EM-ONE) + ALPHA
                        IF (ALPEM.EQ.ZERO) ALPEM = ONE
                        SUM = (SUM+TEMPA*ALP2EM)*ALPEM/EM
                     END IF
  200             CONTINUE
               END IF
C---------------------------------------------------------------------
C Store B(NB).
C---------------------------------------------------------------------
  210          B(N) = TEMPA
               IF (NEND .GE. 0) THEN
                  IF (NB .LE. 1) THEN
                        ALP2EM = ALPHA
                        IF ((ALPHA+ONE).EQ.ONE) ALP2EM = ONE
                        SUM = SUM + B(1)*ALP2EM
                        GO TO 250
                     ELSE
C---------------------------------------------------------------------
C Calculate and store B(NB-1).
C---------------------------------------------------------------------
                        N = N - 1
                        EN = EN - TWO
                        B(N) = (EN*TEMPA)/X - TEMPB
                        IF (N.EQ.1) GO TO 240
                        M = 2 - M
                        IF (M .NE. 0) THEN
                           EM = EM - ONE
                           ALP2EM = (EM+EM) + ALPHA
                           ALPEM = (EM-ONE) + ALPHA
                           IF (ALPEM.EQ.ZERO) ALPEM = ONE
                           SUM = (SUM+B(N)*ALP2EM)*ALPEM/EM
                        END IF
                  END IF
               END IF
               NEND = N - 2
               IF (NEND .NE. 0) THEN
C---------------------------------------------------------------------
C Calculate via difference equation and store B(N), until N = 2.
C---------------------------------------------------------------------
                  DO 230 L=1,NEND
                     N = N - 1
                     EN = EN - TWO
                     B(N) = (EN*B(N+1))/X - B(N+2)
                     M = 2 - M
                     IF (M .NE. 0) THEN
                        EM = EM - ONE
                        ALP2EM = (EM+EM) + ALPHA
                        ALPEM = (EM-ONE) + ALPHA
                        IF (ALPEM.EQ.ZERO) ALPEM = ONE
                        SUM = (SUM+B(N)*ALP2EM)*ALPEM/EM
                     END IF
  230             CONTINUE
               END IF
C---------------------------------------------------------------------
C Calculate B(1).
C---------------------------------------------------------------------
               B(1) = TWO*(ALPHA+ONE)*B(2)/X - B(3)
  240          EM = EM - ONE
               ALP2EM = (EM+EM) + ALPHA
               IF (ALP2EM.EQ.ZERO) ALP2EM = ONE
               SUM = SUM + B(1)*ALP2EM
C---------------------------------------------------------------------
C Normalize.  Divide all B(N) by sum.
C---------------------------------------------------------------------
  250          IF ((ALPHA+ONE).NE.ONE)
     1              SUM = SUM*FUNC(ALPHA)*(X*HALF)**(-ALPHA)
               TEMPA = ENMTEN
               IF (SUM.GT.ONE) TEMPA = TEMPA*SUM
               DO 260 N=1,NB
                 IF (ABS(B(N)).LT.TEMPA) B(N) = ZERO
                 B(N) = B(N)/SUM
  260          CONTINUE
            END IF
C---------------------------------------------------------------------
C Error return -- X, NB, or ALPHA is out of range.
C---------------------------------------------------------------------
         ELSE
            B(1) = ZERO
            NCALC = MIN(NB,0) - 1
      END IF
C---------------------------------------------------------------------
C Exit
C---------------------------------------------------------------------
  300 RETURN
C ---------- Last line of RJBESL ----------
      END

CS    REAL FUNCTION GAMMA(X)
      DOUBLE PRECISION FUNCTION DGAMMA(X)
C----------------------------------------------------------------------
C
C This routine calculates the GAMMA function for a real argument X.
C   Computation is based on an algorithm outlined in reference 1.
C   The program uses rational functions that approximate the GAMMA
C   function to at least 20 significant decimal digits.  Coefficients
C   for the approximation over the interval (1,2) are unpublished.
C   Those for the approximation for X .GE. 12 are from reference 2.
C   The accuracy achieved depends on the arithmetic system, the
C   compiler, the intrinsic functions, and proper selection of the
C   machine-dependent constants.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C beta   - radix for the floating-point representation
C maxexp - the smallest positive power of beta that overflows
C XBIG   - the largest argument for which GAMMA(X) is representable
C          in the machine, i.e., the solution to the equation
C                  GAMMA(XBIG) = beta**maxexp
C XINF   - the largest machine representable floating-point number;
C          approximately beta**maxexp
C EPS    - the smallest positive floating-point number such that
C          1.0+EPS .GT. 1.0
C XMININ - the smallest positive floating-point number such that
C          1/XMININ is machine representable
C
C     Approximate values for some important machines are:
C
C                            beta       maxexp        XBIG
C
C CRAY-1         (S.P.)        2         8191        966.961
C Cyber 180/855
C   under NOS    (S.P.)        2         1070        177.803
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)        2          128        35.040
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)        2         1024        171.624
C IBM 3033       (D.P.)       16           63        57.574
C VAX D-Format   (D.P.)        2          127        34.844
C VAX G-Format   (D.P.)        2         1023        171.489
C
C                            XINF         EPS        XMININ
C
C CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
C Cyber 180/855
C   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
C IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
C VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
C VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns the value XINF for singularities or
C     when overflow would occur.  The computation is believed
C     to be free of underflow and overflow.
C
C
C  Intrinsic functions required are:
C
C     INT, DBLE, EXP, LOG, REAL, SIN
C
C
C References: "An Overview of Software Development for Special
C              Functions", W. J. Cody, Lecture Notes in Mathematics,
C              506, Numerical Analysis Dundee, 1975, G. A. Watson
C              (ed.), Springer Verlag, Berlin, 1976.
C
C              Computer Approximations, Hart, Et. Al., Wiley and
C              sons, New York, 1968.
C
C  Latest modification: October 12, 1989
C
C  Authors: W. J. Cody and L. Stoltz
C           Applied Mathematics Division
C           Argonne National Laboratory
C           Argonne, IL 60439
C
C----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY
CS    REAL 
      DOUBLE PRECISION 
     1    C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE,
     2    TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
C----------------------------------------------------------------------
C  Mathematical constants
C----------------------------------------------------------------------
CS    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/,
CS   1     SQRTPI/0.9189385332046727417803297E0/,
CS   2     PI/3.1415926535897932384626434E0/
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
     1     SQRTPI/0.9189385332046727417803297D0/,
     2     PI/3.1415926535897932384626434D0/
C----------------------------------------------------------------------
C  Machine dependent parameters
C----------------------------------------------------------------------
CS    DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,
CS   1     XINF/3.4E38/
      DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
     1     XINF/1.79D308/
C----------------------------------------------------------------------
C  Numerator and denominator coefficients for rational minimax
C     approximation over (1,2).
C----------------------------------------------------------------------
CS    DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,
CS   1       -3.79804256470945635097577E+2,6.29331155312818442661052E+2,
CS   2       8.66966202790413211295064E+2,-3.14512729688483675254357E+4,
CS   3       -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
CS    DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,
CS   1      -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,
CS   2        2.25381184209801510330112E+4,4.75584627752788110767815E+3,
CS   3      -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
      DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
     1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
     2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
     3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
      DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
     1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
     2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
     3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
C----------------------------------------------------------------------
C  Coefficients for minimax approximation over (12, INF).
C----------------------------------------------------------------------
CS    DATA C/-1.910444077728E-03,8.4171387781295E-04,
CS   1     -5.952379913043012E-04,7.93650793500350248E-04,
CS   2     -2.777777777777681622553E-03,8.333333333333333331554247E-02,
CS   3      5.7083835261E-03/
      DATA C/-1.910444077728D-03,8.4171387781295D-04,
     1     -5.952379913043012D-04,7.93650793500350248D-04,
     2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
     3      5.7083835261D-03/
C----------------------------------------------------------------------
C  Statement functions for conversion between integer and float
C----------------------------------------------------------------------
CS    CONV(I) = REAL(I)
      CONV(I) = DBLE(I)
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .LE. ZERO) THEN
C----------------------------------------------------------------------
C  Argument is negative
C----------------------------------------------------------------------
            Y = -X
            Y1 = AINT(Y)
            RES = Y - Y1
            IF (RES .NE. ZERO) THEN
                  IF (Y1 .NE. AINT(Y1*HALF)*TWO) PARITY = .TRUE.
                  FACT = -PI / SIN(PI*RES)
                  Y = Y + ONE
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C----------------------------------------------------------------------
C  Argument is positive
C----------------------------------------------------------------------
      IF (Y .LT. EPS) THEN
C----------------------------------------------------------------------
C  Argument .LT. EPS
C----------------------------------------------------------------------
            IF (Y .GE. XMININ) THEN
                  RES = ONE / Y
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
         ELSE IF (Y .LT. TWELVE) THEN
            Y1 = Y
            IF (Y .LT. ONE) THEN
C----------------------------------------------------------------------
C  0.0 .LT. argument .LT. 1.0
C----------------------------------------------------------------------
                  Z = Y
                  Y = Y + ONE
               ELSE
C----------------------------------------------------------------------
C  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
C----------------------------------------------------------------------
                  N = INT(Y) - 1
                  Y = Y - CONV(N)
                  Z = Y - ONE
            END IF
C----------------------------------------------------------------------
C  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
C----------------------------------------------------------------------
            XNUM = ZERO
            XDEN = ONE
            DO 260 I = 1, 8
               XNUM = (XNUM + P(I)) * Z
               XDEN = XDEN * Z + Q(I)
  260       CONTINUE
            RES = XNUM / XDEN + ONE
            IF (Y1 .LT. Y) THEN
C----------------------------------------------------------------------
C  Adjust result for case  0.0 .LT. argument .LT. 1.0
C----------------------------------------------------------------------
                  RES = RES / Y1
               ELSE IF (Y1 .GT. Y) THEN
C----------------------------------------------------------------------
C  Adjust result for case  2.0 .LT. argument .LT. 12.0
C----------------------------------------------------------------------
                  DO 290 I = 1, N
                     RES = RES * Y
                     Y = Y + ONE
  290             CONTINUE
            END IF
         ELSE
C----------------------------------------------------------------------
C  Evaluate for argument .GE. 12.0,
C----------------------------------------------------------------------
            IF (Y .LE. XBIG) THEN
                  YSQ = Y * Y
                  SUM = C(7)
                  DO 350 I = 1, 6
                     SUM = SUM / YSQ + C(I)
  350             CONTINUE
                  SUM = SUM/Y - Y + SQRTPI
                  SUM = SUM + (Y-HALF)*LOG(Y)
                  RES = EXP(SUM)
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C----------------------------------------------------------------------
C  Final adjustments and return
C----------------------------------------------------------------------
      IF (PARITY) RES = -RES
      IF (FACT .NE. ONE) RES = FACT / RES
CS900 GAMMA = RES
  900 DGAMMA = RES
      RETURN
C ---------- Last line of GAMMA ----------
      END


      SUBROUTINE RYBESL(X,ALPHA,NB,BY,NCALC)
C----------------------------------------------------------------------
C
C  This routine calculates Bessel functions Y SUB(N+ALPHA) (X)
C  for non-negative argument X, and non-negative order N+ALPHA.
C
C
C Explanation of variables in the calling sequence
C
C X     - Working precision non-negative real argument for which
C         Y's are to be calculated.
C ALPHA - Working precision fractional part of order for which
C         Y's are to be calculated.  0 .LE. ALPHA .LT. 1.0.
C NB    - Integer number of functions to be calculated, NB .GT. 0.
C         The first function calculated is of order ALPHA, and the 
C         last is of order (NB - 1 + ALPHA).
C BY    - Working precision output vector of length NB.  If the
C         routine terminates normally (NCALC=NB), the vector BY
C         contains the functions Y(ALPHA,X), ... , Y(NB-1+ALPHA,X),
C         If (0 .LT. NCALC .LT. NB), BY(I) contains correct function
C         values for I .LE. NCALC, and contains the ratios
C         Y(ALPHA+I-1,X)/Y(ALPHA+I-2,X) for the rest of the array.
C NCALC - Integer output variable indicating possible errors.
C         Before using the vector BY, the user should check that 
C         NCALC=NB, i.e., all orders have been calculated to
C         the desired accuracy.  See error returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   beta   = Radix for the floating-point system
C   p      = Number of significant base-beta digits in the
C            significand of a floating-point number
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C   EPS    = beta ** (-p)
C   DEL    = Machine number below which sin(x)/x = 1; approximately
C            SQRT(EPS).
C   XMIN   = Smallest acceptable argument for RBESY; approximately
C            max(2*beta**minexp,2/XINF), rounded up
C   XINF   = Largest positive machine number; approximately
C            beta**maxexp
C   THRESH = Lower bound for use of the asymptotic form; approximately
C            AINT(-LOG10(EPS/2.0))+1.0
C   XLARGE = Upper bound on X; approximately 1/DEL, because the sine
C            and cosine functions have lost about half of their 
C            precision at that point.
C
C
C     Approximate values for some important machines are:
C
C                        beta    p     minexp      maxexp      EPS
C
C  CRAY-1        (S.P.)    2    48     -8193        8191    3.55E-15
C  Cyber 180/185 
C    under NOS   (S.P.)    2    48      -975        1070    3.55E-15
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)    2    24      -126         128    5.96E-8
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)    2    53     -1022        1024    1.11D-16
C  IBM 3033      (D.P.)   16    14       -65          63    1.39D-17
C  VAX           (S.P.)    2    24      -128         127    5.96E-8
C  VAX D-Format  (D.P.)    2    56      -128         127    1.39D-17
C  VAX G-Format  (D.P.)    2    53     -1024        1023    1.11D-16
C
C
C                         DEL      XMIN      XINF     THRESH  XLARGE
C
C CRAY-1        (S.P.)  5.0E-8  3.67E-2466 5.45E+2465  15.0E0  2.0E7
C Cyber 180/855
C   under NOS   (S.P.)  5.0E-8  6.28E-294  1.26E+322   15.0E0  2.0E7
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)  1.0E-4  2.36E-38   3.40E+38     8.0E0  1.0E4
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)  1.0D-8  4.46D-308  1.79D+308   16.0D0  1.0D8
C IBM 3033      (D.P.)  1.0D-8  2.77D-76   7.23D+75    17.0D0  1.0D8
C VAX           (S.P.)  1.0E-4  1.18E-38   1.70E+38     8.0E0  1.0E4
C VAX D-Format  (D.P.)  1.0D-9  1.18D-38   1.70D+38    17.0D0  1.0D9
C VAX G-Format  (D.P.)  1.0D-8  2.23D-308  8.98D+307   16.0D0  1.0D8
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  In case of an error, NCALC .NE. NB, and not all Y's are
C  calculated to the desired accuracy.
C
C  NCALC .LT. -1:  An argument is out of range. For example,
C       NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.
C       XMAX.  In this case, BY(1) = 0.0, the remainder of the
C       BY-vector is not calculated, and NCALC is set to
C       MIN0(NB,0)-2  so that NCALC .NE. NB.
C  NCALC = -1:  Y(ALPHA,X) .GE. XINF.  The requested function
C       values are set to 0.0.
C  1 .LT. NCALC .LT. NB: Not all requested function values could
C       be calculated accurately.  BY(I) contains correct function
C       values for I .LE. NCALC, and and the remaining NB-NCALC
C       array elements contain 0.0.
C
C
C Intrinsic functions required are:
C
C     DBLE, EXP, INT, MAX, MIN, REAL, SQRT
C
C
C Acknowledgement
C
C  This program draws heavily on Temme's Algol program for Y(a,x)
C  and Y(a+1,x) and on Campbell's programs for Y_nu(x).  Temme's
C  scheme is used for  x < THRESH, and Campbell's scheme is used
C  in the asymptotic region.  Segments of code from both sources
C  have been translated into Fortran 77, merged, and heavily modified.
C  Modifications include parameterization of machine dependencies,
C  use of a new approximation for ln(gamma(x)), and built-in
C  protection against over/underflow.
C
C References: "Bessel functions J_nu(x) and Y_nu(x) of real
C              order and real argument," Campbell, J. B.,
C              Comp. Phy. Comm. 18, 1979, pp. 133-142.
C
C             "On the numerical evaluation of the ordinary
C              Bessel function of the second kind," Temme,
C              N. M., J. Comput. Phys. 21, 1976, pp. 343-350.
C
C  Latest modification: March 19, 1990
C
C  Modified by: W. J. Cody
C               Applied Mathematics Division
C               Argonne National Laboratory
C               Argonne, IL  60439
C
C----------------------------------------------------------------------
      INTEGER I,K,NA,NB,NCALC
CS    REAL
      DOUBLE PRECISION
     1  ALFA,ALPHA,AYE,B,BY,C,CH,COSMU,D,DEL,DEN,DDIV,DIV,DMU,D1,D2,
     2  E,EIGHT,EN,ENU,EN1,EPS,EVEN,EX,F,FIVPI,G,GAMMA,H,HALF,ODD,
     3  ONBPI,ONE,ONE5,P,PA,PA1,PI,PIBY2,PIM5,Q,QA,QA1,Q0,R,S,SINMU,
     4  SQ2BPI,TEN9,TERM,THREE,THRESH,TWO,TWOBYX,X,XINF,XLARGE,XMIN,
     5  XNA,X2,YA,YA1,ZERO
      DIMENSION BY(NB),CH(21)
C----------------------------------------------------------------------
C  Mathematical constants
C    FIVPI = 5*PI
C    PIM5 = 5*PI - 15
C    ONBPI = 1/PI
C    PIBY2 = PI/2
C    SQ2BPI = SQUARE ROOT OF 2/PI
C----------------------------------------------------------------------
CS    DATA ZERO,HALF,ONE,TWO,THREE/0.0E0,0.5E0,1.0E0,2.0E0,3.0E0/
CS    DATA EIGHT,ONE5,TEN9/8.0E0,15.0E0,1.9E1/
CS    DATA FIVPI,PIBY2/1.5707963267948966192E1,1.5707963267948966192E0/
CS    DATA PI,SQ2BPI/3.1415926535897932385E0,7.9788456080286535588E-1/
CS    DATA PIM5,ONBPI/7.0796326794896619231E-1,3.1830988618379067154E-1/
      DATA ZERO,HALF,ONE,TWO,THREE/0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/
      DATA EIGHT,ONE5,TEN9/8.0D0,15.0D0,1.9D1/
      DATA FIVPI,PIBY2/1.5707963267948966192D1,1.5707963267948966192D0/
      DATA PI,SQ2BPI/3.1415926535897932385D0,7.9788456080286535588D-1/
      DATA PIM5,ONBPI/7.0796326794896619231D-1,3.1830988618379067154D-1/
C----------------------------------------------------------------------
C  Machine-dependent constants
C----------------------------------------------------------------------
CS    DATA DEL,XMIN,XINF,EPS/1.0E-4,2.36E-38,3.40E38,5.96E-8/
CS    DATA THRESH,XLARGE/8.0E0,1.0E4/
      DATA DEL,XMIN,XINF,EPS/1.0D-8,4.46D-308,1.79D308,1.11D-16/
      DATA THRESH,XLARGE/16.0D0,1.0D8/
C----------------------------------------------------------------------
C  Coefficients for Chebyshev polynomial expansion of 
C         1/gamma(1-x), abs(x) .le. .5
C----------------------------------------------------------------------
CS    DATA CH/-0.67735241822398840964E-23,-0.61455180116049879894E-22,
CS   1         0.29017595056104745456E-20, 0.13639417919073099464E-18,
CS   2         0.23826220476859635824E-17,-0.90642907957550702534E-17,
CS   3        -0.14943667065169001769E-14,-0.33919078305362211264E-13,
CS   4        -0.17023776642512729175E-12, 0.91609750938768647911E-11,
CS   5         0.24230957900482704055E-09, 0.17451364971382984243E-08,
CS   6        -0.33126119768180852711E-07,-0.86592079961391259661E-06,
CS   7        -0.49717367041957398581E-05, 0.76309597585908126618E-04,
CS   8         0.12719271366545622927E-02, 0.17063050710955562222E-02,
CS   9        -0.76852840844786673690E-01,-0.28387654227602353814E+00,
CS   A         0.92187029365045265648E+00/
      DATA CH/-0.67735241822398840964D-23,-0.61455180116049879894D-22,
     1         0.29017595056104745456D-20, 0.13639417919073099464D-18,
     2         0.23826220476859635824D-17,-0.90642907957550702534D-17,
     3        -0.14943667065169001769D-14,-0.33919078305362211264D-13,
     4        -0.17023776642512729175D-12, 0.91609750938768647911D-11,
     5         0.24230957900482704055D-09, 0.17451364971382984243D-08,
     6        -0.33126119768180852711D-07,-0.86592079961391259661D-06,
     7        -0.49717367041957398581D-05, 0.76309597585908126618D-04,
     8         0.12719271366545622927D-02, 0.17063050710955562222D-02,
     9        -0.76852840844786673690D-01,-0.28387654227602353814D+00,
     A         0.92187029365045265648D+00/
C----------------------------------------------------------------------
      EX = X
      ENU = ALPHA
      IF ((NB .GT. 0) .AND. (X .GE. XMIN) .AND. (EX .LT. XLARGE)
     1       .AND. (ENU .GE. ZERO) .AND. (ENU .LT. ONE))  THEN
            XNA = AINT(ENU+HALF)
            NA = INT(XNA)
            IF (NA .EQ. 1) ENU = ENU - XNA
            IF (ENU .EQ. -HALF) THEN
                  P = SQ2BPI/SQRT(EX)
                  YA = P * SIN(EX)
                  YA1 = -P * COS(EX)
               ELSE IF (EX .LT. THREE) THEN
C----------------------------------------------------------------------
C  Use Temme's scheme for small X
C----------------------------------------------------------------------
                  B = EX * HALF
                  D = -LOG(B)
                  F = ENU * D
                  E = B**(-ENU)
                  IF (ABS(ENU) .LT. DEL) THEN
                        C = ONBPI
                     ELSE
                        C = ENU / SIN(ENU*PI)
                  END IF
C----------------------------------------------------------------------
C  Computation of sinh(f)/f
C----------------------------------------------------------------------
                  IF (ABS(F) .LT. ONE) THEN
                        X2 = F*F
                        EN = TEN9
                        S = ONE
                        DO 80 I = 1, 9
                           S = S*X2/EN/(EN-ONE)+ONE
                           EN = EN - TWO
   80                   CONTINUE
                     ELSE 
                        S = (E - ONE/E) * HALF / F
                  END IF
C----------------------------------------------------------------------
C  Computation of 1/gamma(1-a) using Chebyshev polynomials
C----------------------------------------------------------------------
                  X2 = ENU*ENU*EIGHT
                  AYE = CH(1)
                  EVEN = ZERO
                  ALFA = CH(2)
                  ODD = ZERO
                  DO 40 I = 3, 19, 2
                     EVEN = -(AYE+AYE+EVEN)
                     AYE = -EVEN*X2 - AYE + CH(I)
                     ODD = -(ALFA+ALFA+ODD)
                     ALFA = -ODD*X2 - ALFA + CH(I+1)
   40             CONTINUE
                  EVEN = (EVEN*HALF+AYE)*X2 - AYE + CH(21)
                  ODD = (ODD+ALFA)*TWO
                  GAMMA = ODD*ENU + EVEN
C----------------------------------------------------------------------
C  End of computation of 1/gamma(1-a)
C----------------------------------------------------------------------
                  G = E * GAMMA
                  E = (E + ONE/E) * HALF
                  F = TWO*C*(ODD*E+EVEN*S*D)
                  E = ENU*ENU
                  P = G*C
                  Q = ONBPI / G
                  C = ENU*PIBY2
                  IF (ABS(C) .LT. DEL) THEN
                        R = ONE
                     ELSE 
                        R = SIN(C)/C
                  END IF
                  R = PI*C*R*R
                  C = ONE
                  D = - B*B
                  H = ZERO
                  YA = F + R*Q
                  YA1 = P
                  EN = ZERO
  100             EN = EN + ONE
                  IF (ABS(G/(ONE+ABS(YA)))
     1                      + ABS(H/(ONE+ABS(YA1))) .GT. EPS) THEN
                        F = (F*EN+P+Q)/(EN*EN-E)
                        C = C * D/EN
                        P = P/(EN-ENU)
                        Q = Q/(EN+ENU)
                        G = C*(F+R*Q)
                        H = C*P - EN*G
                        YA = YA + G
                        YA1 = YA1+H
                        GO TO 100
                  END IF
                  YA = -YA
                  YA1 = -YA1/B
               ELSE IF (EX .LT. THRESH) THEN
C----------------------------------------------------------------------
C  Use Temme's scheme for moderate X
C----------------------------------------------------------------------
                  C = (HALF-ENU)*(HALF+ENU)
                  B = EX + EX
                  E = (EX*ONBPI*COS(ENU*PI)/EPS)
                  E = E*E
                  P = ONE
                  Q = -EX
                  R = ONE + EX*EX
                  S = R
                  EN = TWO
  200             IF (R*EN*EN .LT. E) THEN
                        EN1 = EN+ONE
                        D = (EN-ONE+C/EN)/S
                        P = (EN+EN-P*D)/EN1
                        Q = (-B+Q*D)/EN1
                        S = P*P + Q*Q
                        R = R*S
                        EN = EN1
                        GO TO 200
                  END IF
                  F = P/S
                  P = F
                  G = -Q/S
                  Q = G
  220             EN = EN - ONE  
                  IF (EN .GT. ZERO) THEN
                        R = EN1*(TWO-P)-TWO
                        S = B + EN1*Q
                        D = (EN-ONE+C/EN)/(R*R+S*S)
                        P = D*R
                        Q = D*S
                        E = F + ONE
                        F = P*E - G*Q
                        G = Q*E + P*G
                        EN1 = EN
                        GO TO 220
                  END IF
                  F = ONE + F
                  D = F*F + G*G
                  PA = F/D
                  QA = -G/D
                  D = ENU + HALF -P
                  Q = Q + EX
                  PA1 = (PA*Q-QA*D)/EX
                  QA1 = (QA*Q+PA*D)/EX
                  B = EX - PIBY2*(ENU+HALF)
                  C = COS(B)
                  S = SIN(B)
                  D = SQ2BPI/SQRT(EX)
                  YA = D*(PA*S+QA*C)
                  YA1 = D*(QA1*S-PA1*C)
               ELSE
C----------------------------------------------------------------------
C  Use Campbell's asymptotic scheme.
C----------------------------------------------------------------------
                  NA = 0
                  D1 = AINT(EX/FIVPI)
                  I = INT(D1)
                  DMU = ((EX-ONE5*D1)-D1*PIM5)-(ALPHA+HALF)*PIBY2
                  IF (I-2*(I/2) .EQ. 0) THEN
                        COSMU = COS(DMU)
                        SINMU = SIN(DMU)
                     ELSE
                        COSMU = -COS(DMU)
                        SINMU = -SIN(DMU)
                  END IF
                  DDIV = EIGHT * EX
                  DMU = ALPHA
                  DEN = SQRT(EX)
                  DO 350 K = 1, 2
                     P = COSMU
                     COSMU = SINMU
                     SINMU = -P
                     D1 = (TWO*DMU-ONE)*(TWO*DMU+ONE)
                     D2 = ZERO
                     DIV = DDIV
                     P = ZERO
                     Q = ZERO
                     Q0 = D1/DIV
                     TERM = Q0
                     DO 310 I = 2, 20
                        D2 = D2 + EIGHT
                        D1 = D1 - D2
                        DIV = DIV + DDIV
                        TERM = -TERM*D1/DIV
                        P = P + TERM
                        D2 = D2 + EIGHT
                        D1 = D1 - D2
                        DIV = DIV + DDIV
                        TERM = TERM*D1/DIV
                        Q = Q + TERM
                        IF (ABS(TERM) .LE. EPS) GO TO 320
  310                CONTINUE
  320                P = P + ONE
                     Q = Q + Q0
                     IF (K .EQ. 1) THEN
                           YA = SQ2BPI * (P*COSMU-Q*SINMU) / DEN
                        ELSE
                           YA1 = SQ2BPI * (P*COSMU-Q*SINMU) / DEN
                     END IF
                     DMU = DMU + ONE
  350             CONTINUE
            END IF
            IF (NA .EQ. 1) THEN
               H = TWO*(ENU+ONE)/EX
               IF (H .GT. ONE) THEN
                  IF (ABS(YA1) .GT. XINF/H) THEN
                     H = ZERO
                     YA = ZERO
                  END IF
               END IF
               H = H*YA1 - YA
               YA = YA1
               YA1 = H
            END IF
C----------------------------------------------------------------------
C  Now have first one or two Y's
C----------------------------------------------------------------------
            BY(1) = YA
            BY(2) = YA1
            IF (YA1 .EQ. ZERO) THEN
                  NCALC = 1
               ELSE
                  AYE = ONE + ALPHA
                  TWOBYX = TWO/EX
                  NCALC = 2
                  DO 400 I = 3, NB
                     IF (TWOBYX .LT. ONE) THEN
                           IF (ABS(BY(I-1))*TWOBYX .GE. XINF/AYE)
     1                                                     GO TO 450
                        ELSE
                           IF (ABS(BY(I-1)) .GE. XINF/AYE/TWOBYX )
     1                                                     GO TO 450
                     END IF
                     BY(I) = TWOBYX*AYE*BY(I-1) - BY(I-2) 
                     AYE = AYE + ONE
                     NCALC = NCALC + 1
  400             CONTINUE
            END IF
  450       DO 460 I = NCALC+1, NB
               BY(I) = ZERO
  460       CONTINUE
         ELSE
            BY(1) = ZERO
            NCALC = MIN(NB,0) - 1
      END IF
  900 RETURN
C---------- Last line of RYBESL ----------
      END



      subroutine potential(Numchannels,NumDataPoints,NumTotStates,
     >    PFile,QFile,xdata,VQmat,Pmat)

      integer Numchannels,NumDataPoints,NumTotStates

      double precision VQmat(NumDataPoints,Numchannels,Numchannels)
      double precision xdata(NumDataPoints)
      double precision Pmat(NumDataPoints,NumChannels,NumChannels)
      character*64 PFile,QFile

      open(200,file=PFile,status='old')
      open(300,file=QFile,status='old')


      Pmat = 0.d0
      VQmat = 0.d0

      do i=1,NumDataPoints
       read(200,*) xdata(i)
       read(300,*)
       do j=1,Numchannels
        read(200,11)(Pmat(i,j,k),k=1,Numchannels)
        read(300,11)(VQmat(i,j,k),k=1,Numchannels)
       enddo
       do j=Numchannels+1,NumTotStates
        read(200,11)
        read(300,11)
       enddo
      enddo

c      zero the couplings to see what happens
c$$$      do k = 1, NumDataPoints
c$$$         do i = 1, NumChannels
c$$$            do j = 1, NumChannels
c$$$               if(i.ne.j) then
c$$$                  Pmat(k,i,j) = 0d0
c$$$                  VQmat(k,i,j) = 0d0
c$$$               endif
c$$$            enddo
c$$$         enddo
c$$$      enddo
                  
cc      XXX
c       do i = 1,NumDataPoints
c       do j = 1,NumChannels
c       do k = j+1,NumChannels
c       j=1; k=3;
c       VQmat(i,j,k) = 0.d0
c       VQmat(i,k,j) = 0.d0
c       Pmat(i,j,k) = 0.d0
c       Pmat(i,k,j) = 0.d0
cc       VQmat(i,j,k) = VQmat(i,k,j) 
cc        Pmat(i,j,k) = -Pmat(i,k,j) 
c       enddo
c       enddo
c       enddo

cc     Filtering noise (jpdincao)
c      do i=1,NumDataPoints
c       do j=1,Numchannels
c       do k=1,Numchannels
c           if (dabs(Pmat(i,j,k)).le.1.d-14) Pmat(i,j,k) = 0.d0 !sign(1.d-14,Pmat(i,j,k))
c           if (dabs(VQmat(i,j,k)).le.1.d-14) VQmat(i,j,k) = 0.d0 !sign(1.d-14,VQmat(i,j,k))
c       enddo
c       enddo
c      enddo


      close(200)
      close(300)
c     3001 format(100(e20.12,1x))
 11   format(1P,100e22.12)
 3001 format(100(e18.12,1x))
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function u(n,x)

      integer n
      double precision x

      if (n .eq. 1) u = x**2-1.25d0*x**3-0.5d0*x**4+0.75d0*x**5
      if (n .eq. 2) u = (x**2-x**3-x**4+x**5)/4.0d0
      if (n .eq. 3) u = 1-2.0d0*x**2+x**4
      if (n .eq. 4) u = x-2.0d0*x**3+x**5
      if (n .eq. 5) u = x**2+1.25d0*x**3-0.5d0*x**4-0.75d0*x**5
      if (n .eq. 6) u = (-x**2-x**3+x**4+x**5)/4.0d0

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine BasicOverlapMatrix(BasicS)

      double precision BasicS(6,6)

      BasicS(1,1) = 1046.0d0/3465.0d0
      BasicS(1,2) = 38.0d0/1155.0d0
      BasicS(1,3) = 8.0d0/63.0d0
      BasicS(1,4) = -32.0d0/693.0d0
      BasicS(1,5) = 131.0d0/3465.0d0
      BasicS(1,6) = -29.0d0/3465.0d0

      BasicS(2,2) = 16.0d0/3465.0d0
      BasicS(2,3) = 8.0d0/315.0d0
      BasicS(2,4) = -8.0d0/1155.0d0
      BasicS(2,5) = 29.0d0/3465.0d0
      BasicS(2,6) = -2.0d0/1155.0d0

      BasicS(3,3) = 256.0d0/315.0d0
      BasicS(3,4) = 0.0d0
      BasicS(3,5) = 8.0d0/63.0d0
      BasicS(3,6) = -8.0d0/315.0d0

      BasicS(4,4) = 256.0d0/3465.0d0
      BasicS(4,5) = 32.0d0/693.0d0
      BasicS(4,6) = -8.0d0/1155.0d0

      BasicS(5,5) = 1046.0d0/3465.0d0
      BasicS(5,6) = -38.0d0/1155.0d0

      BasicS(6,6) = 16.0d0/3465.0d0

      do i = 2,6
       do j = 1,i-1
        BasicS(i,j) = BasicS(j,i)
       enddo
      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine BasicKineticMatrix(BasicT)

c -1/2 factor is included:  this is (-1/2)* \int_{-1}^{1}{ui(x) uj''(x)}

      double precision BasicT(6,6)

      BasicT(1,1) = 139.0d0/210.0d0
      BasicT(1,2) = 223.0d0/420.0d0
      BasicT(1,3) = -64.0d0/105.0d0
      BasicT(1,4) = 4.0d0/21.0d0
      BasicT(1,5) = -11.0d0/210.0d0
      BasicT(1,6) = -1.0d0/140.0d0

      BasicT(2,2) = 2.0d0/45.0d0
      BasicT(2,3) = -4.0d0/105.0d0
      BasicT(2,4) = -4.0d0/315.0d0
      BasicT(2,5) = 1.0d0/140.0d0
      BasicT(2,6) = -1.0d0/126.0d0

      BasicT(3,3) = 128.0d0/105.0d0
      BasicT(3,4) = 0.0d0
      BasicT(3,5) = -64.0d0/105.0d0
      BasicT(3,6) = 4.0d0/105.0d0

      BasicT(4,4) = 128.0d0/315.0d0
      BasicT(4,5) = -4.0d0/21.0d0
      BasicT(4,6) = -4.0d0/315.0d0

      BasicT(5,5) = 139.0d0/210.0d0
      BasicT(5,6) = -223.0d0/420.0d0

      BasicT(6,6) = 2.0d0/45.0d0

      do i = 2,6
       do j = 1,i-1
        BasicT(i,j) = BasicT(j,i)
       enddo
      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine GenMEMap(Numchannels,NumSectors,MEMap)

      integer NumSectors,Numchannels
      integer MEMap(Numchannels*Numsectors,6)

c form matrix element map

      icnt = 1
      do k=1,Numsectors
       do j=1,6
        do i=1,NumChannels

         index = Numchannels*(k-1) + i

         if (k .eq. 1 .and. j .eq. 1)then
          MEMap(index,j) = 0
         elseif (k .eq. Numsectors .and. j .eq. 5)then
          MEMap(index,j) = 0
         else
          MEMap(index,j) = icnt
          icnt = icnt + 1
         endif

         enddo
        enddo
        icnt = icnt - 2*Numchannels
       enddo

1901  format(6(i4,1x))

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine CalcClosed(Numchannels,NumSectors,MatrixDim,kl,iband,MEMap,Mass,energy,
     >   xsector,LegPoints,xLeg,wLeg,NumDataPoints,xdata,VQmat,Pmat,PFit,VQFit,
     >   WhereStart,PExponents,VQExponents,Gcc)

      integer NumSectors,Numchannels,MatrixDim,kl,iband,LegPoints
      integer MEMap(NumSectors*NumChannels,6)
      double precision Mass,xSector(*),xLeg(*),wLeg(*),energy
      double precision Gcc(iband,MatrixDim)

      double precision Scale(6),IntScale,BasicT(6,6)
      double precision BasicVQ(6,6),BasicS(6,6),BasicP(6,6)
      double precision ScaledZero,x(LegPoints),xL(LegPoints),xR(LegPoints)
      double precision u(6,LegPoints),VQInterp(LegPoints)
      double precision up(6,LegPoints),PInterp(LegPoints)
      double precision dPInterp(LegPoints),PLInterp(LegPoints)
      double precision PRInterp(LegPoints)
      double precision VQmat(NumDataPoints,Numchannels,Numchannels)
      double precision Pmat(NumDataPoints,Numchannels,Numchannels)
      double precision xdata(NumDataPoints),ydata(NumDataPoints)

      double precision PFit(Numchannels,Numchannels,3),VQFit(Numchannels,Numchannels,3)
      double precision WhereStart,vvj,ppj,psign
      double precision PExponents(Numchannels,Numchannels,3)
      double precision VQExponents(Numchannels,Numchannels,3)
      double precision, external :: AsymptoticP, AsymptoticVQ

      call BasicKineticMatrix(BasicT)

      call CalcLocalBasis(LegPoints,xLeg,u,up)

      call BasicOverlapMatrix(BasicS)

      ku = kl
      kb = kl + ku + 1

      Gcc = 0.0d0
      Scale = 1.0d0

      knc = 1
      kjc = 1

c index one sector at a time
      do n = 1,NumSectors

       IntScale = (xSector(n+1)-xSector(n))/2.0d0
       ScaledZero = (xSector(n+1)+xSector(n))/2.0d0

       do nc = 1,Numchannels
          
c     integrate potentials in each block
          
          do jc = 1,Numchannels
             psign = sign(1d0,AsymptoticP(nc,jc,xdata(NumDataPoints))*Pmat(NumDataPoints,nc,jc))
c             write(6,*) "psign = ", psign
             do ik=1,NumDataPoints
                ydata(ik)=VQmat(ik,nc,jc)
             enddo
             
             do j = 1,LegPoints
                x(j) = IntScale*xLeg(j)+ScaledZero
                
                IF(x(j) .ge. WhereStart) THEN
                   
c-----------------------------------------------------------------------
c     comment out this original section in favor of exact asymptotic forms
                   vvj=VQFit(nc,jc,1)/x(j)**VQExponents(nc,jc,1)
                   vvj=vvj+VQFit(nc,jc,2)/x(j)**VQExponents(nc,jc,2)
                   vvj=vvj+VQFit(nc,jc,3)/x(j)**VQExponents(nc,jc,3)
                   VQInterp(j)=vvj
                   
                   ppj=PFit(nc,jc,1)/x(j)**PExponents(nc,jc,1)
                   ppj=ppj+PFit(nc,jc,2)/x(j)**PExponents(nc,jc,2)
                   ppj=ppj+PFit(nc,jc,3)/x(j)**PExponents(nc,jc,3)
                   PInterp(j)=ppj
                   
c     Exact asymptotic forms here:
c$$$                   vvj = AsymptoticVQ(Mass,nc,jc,x(j))
c$$$                   VQInterp(j) = vvj
c$$$                   ppj = AsymptoticP(nc,jc,x(j))
c$$$                   PInterp(j) = ppj*psign
c-----------------------------------------------------------------------
                endif
             enddo
             
             IF(x(1) .LT. WhereStart) THEN
                call Akima(NumDataPoints,xdata,ydata,LegPoints,x,VQInterp)
             ENDIF
             
c     c       Filtering noise (jpdincao)
c     do j = 1,LegPoints
c     if (dabs(VQInterp(j)).le.1.d-14) VQInterp(j) = 0.d0 !sign(1.d-14,VQInterp(j))
c     enddo
             
c     c       Cutoff Long-range Q-coupling (jpdincao)
c     if (nc.eq.1.or.nc.eq.3) then
c     if (jc.eq.1.or.jc.eq.3) then
c     do j = 1,LegPoints
c     if (x(j).ge.WhereStart) VQInterp(j) = 0.d0 
c     enddo
c     endif
c     endif
             
             do i = 1,6
                do ip = 1,6
                   BasicVQ(i,ip) = 0.0d0
                   do j = 1,LegPoints
                      BasicVQ(i,ip) = BasicVQ(i,ip) + wLeg(j)*u(i,j)*VQInterp(j)*u(ip,j)
                   enddo
                enddo
             enddo
             
             do ik=1,NumDataPoints
                ydata(ik)=Pmat(ik,nc,jc)
             enddo
             
             IF(x(1) .LT. WhereStart) THEN
                call akima(NumDataPoints,xdata,ydata,LegPoints,x,Pinterp)
             ENDIF
             
c     c       Filtering noise (jpdincao)
c     do j = 1,LegPoints
c     if (dabs(PInterp(j)).le.1.d-14) PInterp(j) = 0.d0 !sign(1.d-14,PInterp(j))
c     enddo   

c     c       Cutoff Long-range P-coupling (jpdincao)
c     if (nc.eq.1.or.nc.eq.3) then
c     if (jc.eq.1.or.jc.eq.3) then
c     do j = 1,LegPoints
c           if (x(j).ge.WhereStart) PInterp(j) = 0.d0 
c        enddo
c        endif
c        endif

             do i = 1,6
                do ip = 1,6
                   BasicP(i,ip) = 0.0d0
                   do j = 1,LegPoints
                      BasicP(i,ip) = BasicP(i,ip) + wLeg(j)*PInterp(j)*
     >                     (u(i,j)*up(ip,j) - up(i,j)*u(ip,j))
                   enddo
                enddo
             enddo
             
             if (n .ne. 1) then
                Scale(2) = (xSector(n+1)-xSector(n))/(xSector(n)-xSector(n-1))
             else
                Scale(2) = 1.0d0
             endif
             
             do i = 1,6
                do ip = 1,6
                   
                   if ((MEMap(knc,i) .ne. 0) .AND. (MEMap(kjc,ip) .ne. 0)) then
                      
                      imap = MEMap(knc,i)
                      jmap = MEMap(kjc,ip)
                      
                      if (knc .eq. kjc)then
                         
c     diagonal blocks
                         
                         Gcc(kb+imap-jmap,jmap) = Gcc(kb+imap-jmap,jmap) + 
     >                        2.0d0*mass*Scale(i)*Scale(ip)*
     >                        (energy*Intscale*BasicS(i,ip) -
     >                        (BasicT(i,ip)/(Mass*IntScale)+BasicVQ(i,ip)*IntScale
     >                        + BasicP(i,ip)))
                         
                      else 
                         
c     off-diagonal blocks
                         
                         Gcc(kb+imap-jmap,jmap) = Gcc(kb+imap-jmap,jmap) - 
     >                        2.0d0*mass*Scale(i)*Scale(ip)*(BasicVQ(i,ip)*IntScale
     >                        + BasicP(i,ip))
                         
                      endif
                   endif
                enddo
             enddo
             
             kjc = kjc + 1
          enddo
          
          kjc = Numchannels*(n-1) + 1
          knc = knc + 1
       enddo
       
       knc = Numchannels*n + 1
       kjc = knc
      enddo

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine CalcLocalBasis(LegPoints,xLeg,u,up)

      integer LegPoints
      double precision xLeg(LegPoints),u(6,LegPoints)
      double precision up(6,LegPoints)

      integer n
      double precision x,x2,x3,x4,x5

      do n = 1,LegPoints

       x = xLeg(n)
       x2 = x*x
       x3 = x2*x
       x4 = x3*x
       x5 = x4*x

       u(1,n) = x2 - 1.25d0*x3 - 0.5d0*x4 + 0.75d0*x5
       up(1,n) = 2.0d0*x - 3.750*x2 - 2*x3 + 3.750d0*x4
c       upp(1,n) = 2.0d0 - 7.5d0*x - 6.0d0*x2 + 15.0d0*x3

       u(2,n) = 0.25d0*(x2 - x3 - x4 + x5)
       up(2,n) = 0.5d0*x - 0.75d0*x2 - x3 + 1.25d0*x4
c       upp(2,n) = 0.5d0 - 1.5d0*x - 3.0d0*x2 + 5.0d0*x3

       u(3,n) = 1.0d0 - 2.0d0*x2 + x4
       up(3,n) =  -4.0d0*x + 4.0d0*x3
c       upp(3,n) =  -4.0d0 + 12.0d0*x2

       u(4,n) = x - 2.0d0*x3 + x5
       up(4,n) = 1.0d0 - 6.0d0*x2 + 5.0d0*x4
c       upp(4,n) =  -12.0d0*x + 20.0d0*x3

       u(5,n) = x2 + 1.25d0*x3 - 0.5d0*x4 - 0.75d0*x5
       up(5,n) = 2.0d0*x + 3.75d0*x2 - 2.0d0*x3 - 3.75d0*x4
c       upp(5,n) = 2.0d0 + 7.5d0*x - 6.0d0*x2 - 15.0d0*x3

       u(6,n) = 0.25d0*( -x2 - x3 + x4 + x5)
       up(6,n) =  -0.5d0*x - 0.75d0*x2 + x3 + 1.25d0*x4
c       upp(6,n) =  -0.5d0 - 1.5d0*x + 3.0d0*x2 + 5.0d0*x3

      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine CalcOpen(NumOpenL,NumOpenR,Numchannels,NumSectors,MatrixDim,MEMap,Mass,
     >  energy,xsector,LegPoints,xLeg,wLeg,NumDataPoints,xdata,VQmat,Pmat,
     >  PFit,VQFit,WhereStart,PExponents,VQExponents,Gco,Goo)

      integer NumSectors,NumOpenL,NumOpenR,Numchannels,LegPoints,NumDataPoints
      integer MEMap(NumSectors*Numchannels,6),MatrixDim
      double precision Mass,xSector(*),energy
      double precision Gco(MatrixDim,NumOpenL+NumOpenR)
      double precision Goo(NumOpenL+NumOpenR,NumOpenL+NumOpenR)
      double precision xLeg(LegPoints),wLeg(LegPoints)

      double precision Scale(6),IntScale,BasicT(6,6)
      double precision BasicP(6,6),BasicVQ(6,6),BasicS(6,6)
      double precision ScaledZero,x(LegPoints),up(6,LegPoints)
      double precision u(6,LegPoints),VQInterp(LegPoints)
      double precision PInterp(LegPoints)
      double precision VQmat(NumDataPoints,NumChannels,NumChannels)
      double precision Pmat(NumDataPoints,NumChannels,NumChannels)
      double precision xdata(NumDataPoints)
      double precision ydata(NumDataPoints)

      double precision PFit(NumChannels,NumChannels,3),VQFit(NumChannels,NumChannels,3)
      double precision WhereStart,vvj,ppj,Dpj,psign
      double precision PExponents(NumChannels,NumChannels,3)
      double precision VQExponents(NumChannels,NumChannels,3)
      double precision, external :: AsymptoticP, AsymptoticVQ
      
      call BasicKineticMatrix(BasicT)

      call CalcLocalBasis(LegPoints,xLeg,u,up)

      call BasicOverlapMatrix(BasicS)

      NumOpen = NumOpenL + NumOpenR

      Gco = 0.0d0
      Goo = 0.0d0

      Scale = 1.0d0

c open channels on left boundary

      ip = 1
      n  = 1
      
      do nc = 1,NumOpenL
         IntScale = (xSector(n+1)-xSector(n))/2.0d0 ! a_n in Burke's thesis
         ScaledZero = (xSector(n+1)+xSector(n))/2.0d0 ! d_n in Burke's thesis
         
c     integrate potentials in each block
         
         do jc = 1,NumChannels
            psign = sign(1d0,AsymptoticP(nc,jc,xdata(NumDataPoints))*Pmat(NumDataPoints,nc,jc))
            do ik=1,NumDataPoints
               ydata(ik) = VQmat(ik,nc,jc)
            enddo
            
            do j = 1,LegPoints
               x(j) = IntScale*xLeg(j)+ScaledZero
               
c-----------------------------------------------------------------------
c     comment out this original section in favor of exact asymptotic forms
               vvj=VQFit(nc,jc,1)/x(j)**VQExponents(nc,jc,1)
               vvj=vvj+VQFit(nc,jc,2)/x(j)**VQExponents(nc,jc,2)
               vvj=vvj+VQFit(nc,jc,3)/x(j)**VQExponents(nc,jc,3)
               VQInterp(j)=vvj
               
               ppj=PFit(nc,jc,1)/x(j)**PExponents(nc,jc,1)
               ppj=ppj+PFit(nc,jc,2)/x(j)**PExponents(nc,jc,2)
               ppj=ppj+PFit(nc,jc,3)/x(j)**PExponents(nc,jc,3)
               PInterp(j)=ppj
c$$$c     Exact asymptotic forms:
c$$$               vvj = AsymptoticVQ(Mass,nc,jc,x(j))
c$$$               VQInterp(j) = vvj
c$$$               ppj = AsymptoticP(nc,jc,x(j))
c$$$               PInterp(j) = ppj*psign
               
            enddo
            ! If x < WhereStart, then just interpolate the known points instead of using the extrapolation fits.
            IF(x(1) .LT. WhereStart) THEN
               call Akima(NumDataPoints,xdata,ydata,LegPoints,x,VQInterp) ! do the interpolation for this particular matrix element only (nc,jc)? (npmehta)
            ENDIF
            
c     c       Filtering noise (jpdincao)
c     do j = 1,LegPoints
c           if (dabs(VQInterp(j)).le.1.d-14) VQInterp(j) = 0.d0 !sign(1.d-14,VQInterp(j))
c        enddo   

cc       Cutoff Long-range Q-coupling (jpdincao)
c        if (nc.eq.1.or.nc.eq.3) then
c        if (jc.eq.1.or.jc.eq.3) then
c        do j = 1,LegPoints
c           if (x(j).ge.WhereStart) VQInterp(j) = 0.d0 
c        enddo
c        endif
c        endif

        do i = 1,6
         BasicVQ(i,ip) = 0.0d0
         do j = 1,LegPoints
          BasicVQ(i,ip) = BasicVQ(i,ip) + wLeg(j)*u(i,j)*VQInterp(j)*u(ip,j)
         enddo
       enddo

        do ik=1,NumDataPoints
         ydata(ik)=Pmat(ik,nc,jc)
        enddo

	 IF(x(1) .LT. WhereStart) THEN
         call akima(NumDataPoints,xdata,ydata,LegPoints,x,Pinterp)
	 ENDIF


cc       Filtering noise (jpdincao)
c        do j = 1,LegPoints
c           if (dabs(PInterp(j)).le.1.d-14) PInterp(j) = 0.d0 !sign(1.d-14,PInterp(j))
c        enddo   

cc       Cutoff Long-range P-coupling (jpdincao)
c        if (nc.eq.1.or.nc.eq.3) then
c        if (jc.eq.1.or.jc.eq.3) then
c        do j = 1,LegPoints
c           if (x(j).ge.WhereStart) PInterp(j) = 0.d0 
c        enddo
c        endif
c        endif

       do i = 1,6
         BasicP(i,ip) = 0.0d0
         do j = 1,LegPoints
          BasicP(i,ip) = BasicP(i,ip) + wLeg(j)*PInterp(j)*
     >        (u(ip,j)*up(i,j) - up(ip,j)*u(i,j))
         enddo
       enddo
       
       do i = 1,6
          irow = MEMap(jc,i)
          
c     diagonal blocks
          if (nc .eq. jc)then
             if (i .ne. ip)then
                Gco(irow,nc) = 2.0d0*mass*Scale(i)*Scale(ip)*
     >               (energy*Intscale*BasicS(i,ip) -
     >               (BasicT(i,ip)/(Mass*IntScale) +
     >               BasicVQ(i,ip)*IntScale
     >               + BasicP(i,ip)))
             else
                if (jc .le. NumOpenL)then
                   Goo(jc,nc) =  2.0d0*mass*Scale(i)*Scale(ip)*
     >                  (energy*Intscale*BasicS(i,ip) -
     >                  (BasicT(i,ip)/(Mass*IntScale) +
     >                  BasicVQ(i,ip)*IntScale
     >                  + BasicP(i,ip))) 
                endif
             endif
             
c     bloch operator
             if (i .eq. 2)Gco(irow,nc)=Gco(irow,nc) + 1.0d0/Intscale
             
c     off-diagonal blocks
          else
             if (i .ne. ip)then
                Gco(irow,nc) = -2.0d0*mass*Scale(i)*Scale(ip)*
     >               (BasicVQ(i,ip)*IntScale + BasicP(i,ip))
             else
                if (jc .le. NumOpenL)then
                   Goo(jc,nc) =  -2.0d0*mass*Scale(i)*Scale(ip)*
     >                  (BasicVQ(i,ip)*IntScale + BasicP(i,ip))
                endif
             endif
          endif
          
       enddo
       
      enddo
      enddo
      
c     now do open channels on right boundary (Gco)
      
      ip = 5
      n = Numsectors
      
      do nc = NumOpenL+1,NumOpen
         IntScale = (xSector(n+1)-xSector(n))/2.0d0
         ScaledZero = (xSector(n+1)+xSector(n))/2.0d0
         knv = nc - NumOpenL
         
c     integrate potentials in each block
         
         do jc = 1,Numchannels
            kjcs = jc  + Numchannels*(Numsectors-1)
            kjv = jc 
            
            do ik=1,NumDataPoints
               ydata(ik) = VQmat(ik,kjv,knv)
            enddo
            
            do j = 1,LegPoints
               x(j) = IntScale*xLeg(j)+ScaledZero
               
               vvj=VQFit(kjv,knv,1)/x(j)**VQExponents(kjv,knv,1)
               vvj=vvj+VQFit(kjv,knv,2)/x(j)**VQExponents(kjv,knv,2)
               vvj=vvj+VQFit(kjv,knv,3)/x(j)**VQExponents(kjv,knv,3)
               VQInterp(j)=vvj
               
               ppj=PFit(kjv,knv,1)/x(j)**PExponents(kjv,knv,1)
               ppj=ppj+PFit(kjv,knv,2)/x(j)**PExponents(kjv,knv,2)
               ppj=ppj+PFit(kjv,knv,3)/x(j)**PExponents(kjv,knv,3)
               PInterp(j)=ppj
            enddo
            
            IF(x(1) .LT. WhereStart) THEN
               call Akima(NumDataPoints,xdata,ydata,LegPoints,x,VQInterp)
            ENDIF
            
            do i = 1,6
               BasicVQ(i,ip) = 0.0d0
               do j = 1,LegPoints
                  BasicVQ(i,ip) = BasicVQ(i,ip) + wLeg(j)*u(i,j)*VQInterp(j)*u(ip,j)
               enddo
            enddo
            
            do ik=1,NumDataPoints
               ydata(ik)=Pmat(ik,kjv,knv)
            enddo
            
            
            IF(x(1) .LT. WhereStart) THEN
               call akima(NumDataPoints,xdata,ydata,LegPoints,x,Pinterp)
            ENDIF
            
            do i = 1,6
               BasicP(i,ip) = 0.0d0
               do j = 1,LegPoints
                  BasicP(i,ip) = BasicP(i,ip) + wLeg(j)*PInterp(j)*
     >                 (u(i,j)*up(ip,j) - up(i,j)*u(ip,j))
               enddo
            enddo
            
            If(n .ne. 1)Scale(2) = (xSector(n+1)-xSector(n))/(xSector(n)-xSector(n-1))
            
            do i = 1,6
               irow = MEMap(kjcs,i)
               
c     diagonal blocks
               if (nc .eq. jc+NumOpenL)then
                  if (i .ne. ip)then
                     Gco(irow,nc) = 2.0d0*mass*Scale(i)*Scale(ip)*
     >                    (energy*Intscale*BasicS(i,ip) -
     >                    (BasicT(i,ip)/(Mass*IntScale)+BasicVQ(i,ip)*IntScale
     >                    + BasicP(i,ip)))
                  else 
                     if (jc .le. NumOpenR)then
                        Goo(jc+NumOpenL,nc) =  2.0d0*mass*Scale(i)*Scale(ip)*
     >                       (energy*Intscale*BasicS(i,ip) -
     >                       (BasicT(i,ip)/(Mass*IntScale)+BasicVQ(i,ip)*IntScale
     >                       + BasicP(i,ip))) 
                     endif
                  endif
                  
c     bloch operator
                  if (i .eq. 6)Gco(irow,nc)=Gco(irow,nc) - 1.0d0/Intscale
                  
                  
c     off-diagonal blocks
               else
                  if (i .ne. ip)then
                     Gco(irow,nc) = -2.0d0*mass*Scale(i)*Scale(ip)*
     >                    (BasicVQ(i,ip)*IntScale+ BasicP(i,ip))
                  else
                     if (jc .le. NumOpenR)then
                        Goo(jc+NumOpenL,nc) =  -2.0d0*mass*Scale(i)*Scale(ip)*
     >                       (BasicVQ(i,ip)*IntScale + BasicP(i,ip))
                     endif
                  endif
               endif
               
            enddo
            
         enddo
      enddo
      
      return
      end
      
c     ***********************************************************************

      subroutine CalcSmatrix(NumOpenR,leff,Thresholds,TotalAngMomentum,Mass,
     >                       Energy,RMatch,solution,deriv,Smatrix,Tmatrix,CrossSections)
      integer TotalAngMomentum
      integer NumOpenR,ipiv(NumOpenR)
      double precision leff(NumOpenR),Thresholds(NumOpenR),Mass,RMatch
      double precision solution(NumOpenR,NumOpenR)
      double precision deriv(NumOpenR),cmPERau,secPERau,convertCGS
      double precision convertCGS_K3,convertCGS_RX
      double precision CrossSections(NumOpenR,NumOpenR)
      double precision kVector(NumOpenR),f,fp,g,gp,Wronskian,pi,Energy
      double precision Imat(NumOpenR,NumOpenR)
      double precision Jmat(NumOpenR,NumOpenR)
      double precision enorm(NumOpenR),x
      double precision Tmatrix(NumOpenR,NumOpenR),TmatrixAux(NumOpenR,NumOpenR),ImatAux(NumOpenR,NumOpenR)
      double complex Smatrix(NumOpenR,NumOpenR)
      double complex IJmat(NumOpenR,NumOpenR)

      pi = dacos(-1.0d0)
      cmPERau = 5.29178d-09
      secPERau = 2.4189d-17
      convertCGS_K3 = cmPERau**6/secPERau
      convertCGS_RX = cmPERau**3/secPERau

      do i=1,NumOpenR

       kVector(i) = dsqrt(2.0d0*Mass*(Energy-Thresholds(i)))
       enorm(i) = dsqrt(2.0d0*kVector(i)/pi)

c     call BesselBasePair(kVector(i),RMatch,leff(i),f,fp,g,gp)
 
      x = kVector(i)*Rmatch
      call sphbes(leff(i),x,f,g,fp,gp)
      fp = enorm(i)*(f + x*fp)
      f = enorm(i)*Rmatch*f
      gp = enorm(i)*(g + x*gp)
      g = enorm(i)*Rmatch*g

       Wronskian = 1.0d0/(g*fp-gp*f)
       !write(6,*) 'Wronskian check:    ',2.0d0*Wronskian/Pi
36    format(5(e12.6,1x))
 
       !write(6,*) 'Bessel funcs:'
       !write(6,*) i,f,fp
       !write(6,*) i,g,gp
       !write(6,*)

       do j = 1,NumOpenR
        IMat(i,j) = -(g*deriv(j)+gp)*Solution(i,j)*Wronskian
        JMat(i,j) = -(f*deriv(j)+fp)*Solution(i,j)*Wronskian
       enddo
 
      enddo 

c     calculate S-matrix

      IJMat = (0.0d0,0.0d0)
 
       do i = 1,NumOpenR
        do j = 1,NumOpenR
         IJMat(i,j)   = IMat(j,i)-(0.0d0,1.0d0)*JMat(j,i)
         Smatrix(i,j) = IMat(j,i)+(0.0d0,1.0d0)*JMat(j,i)
        enddo
       enddo
 
       call zgesv(NumOpenR,NumOpenR,IJMat,NumOpenR,ipiv,Smatrix,NumOpenR,info)

       do i = 1,NumOpenR
        Smatrix(i,i) = Smatrix(i,i) - (1.0d0,0.0d0)
       enddo


c     calculate T-matrix

       do i = 1,NumOpenR
        do j = 1,NumOpenR
         ImatAux(i,j) = Imat(j,i)
         TmatrixAux(i,j) = Jmat(j,i)
        enddo
       enddo

       call dgesv(NumOpenR,NumOpenR,ImatAux,NumOpenR,ipiv,TmatrixAux,NumOpenR,info)

       do i = 1,NumOpenR
        do j = 1,NumOpenR
         Tmatrix(i,j) = TmatrixAux(j,i)
        enddo
       enddo

 
c i -> final state
c j -> initial state
 
       do i = 1,NumOpenR
          do j = 1,NumOpenR

c     For 1D two-body channels
             
             
c       For Recombination 
c        CrossSections(i,j) = convertCGS_K3*dfloat(2*TotalAngMomentum+1)*
c     .                       192.d0*pi*pi/(Mass*kVector(j)**4)*
c     .                       dreal(dconjg(Smatrix(i,j))*Smatrix(i,j))
c       For Relaxation
c       CrossSections(i,j) = convertCGS_RX*dfloat(2*TotalAngMomentum+1)*
c     .                        pi/(Mass*kVector(j)**1)*
c     .                        dreal(dconjg(Smatrix(i,j))*Smatrix(i,j))
        enddo
       enddo 

       return
       end
c-------------------------------------------------------------------------------------------------------      
      subroutine CalcSmatrix1D(NumOpenR,leff,Thresholds,TotalAngMomentum,Mass,
     >     Energy,RMatch,solution,deriv,Smatrix,Tmatrix,CrossSections,Kmat,Parity)
      implicit none
      integer TotalAngMomentum,i,info,j
      integer NumOpenR,ipiv(NumOpenR)
      double precision leff(NumOpenR),Thresholds(NumOpenR),Mass,RMatch,P,Parity
      double precision solution(NumOpenR,NumOpenR)
      double precision deriv(NumOpenR),cmPERau,secPERau,convertCGS
      double precision convertCGS_K3,convertCGS_RX
      double precision CrossSections(NumOpenR,NumOpenR) !For this 1D Calcualtion, let CrossSections be simply the K-matrix
      double precision kVector(NumOpenR),f,fp,g,gp,Wronskian,pi,Energy
      double precision Imat(NumOpenR,NumOpenR)
      double precision ImatTemp(NumOpenR,NumOpenR)
      double precision Jmat(NumOpenR,NumOpenR)
      double precision Kmat(NumOpenR,NumOpenR)
      double precision enorm(NumOpenR),x
      double precision Tmatrix(NumOpenR,NumOpenR),TmatrixAux(NumOpenR,NumOpenR),ImatAux(NumOpenR,NumOpenR)
      double complex Smatrix(NumOpenR,NumOpenR)
      double complex IJmat(NumOpenR,NumOpenR)
      double precision, external :: kdelta
      
      pi = dacos(-1.0d0)
      cmPERau = 5.29178d-09
      secPERau = 2.4189d-17
      convertCGS_K3 = cmPERau**6/secPERau
      convertCGS_RX = cmPERau**3/secPERau
      
      do i=1,NumOpenR
         
         kVector(i) = dsqrt(2.0d0*Mass*(Energy-Thresholds(i)))
c     call BesselBasePair(kVector(i),RMatch,leff(i),f,fp,g,gp)
         

C     In 1D we switch the roles of f and g so look carefully at the argument order of the following call
c     Note that leff(i) should be zero for all two-body channels in 1D.
c     Make sure this is the case in the input file
c         enorm(i) = dsqrt(2.0d0*kVector(i)/pi)


c$$$         x = kVector(i)*Rmatch
c$$$         call sphbes(leff(i),x,g,f,gp,fp)
c     Now g ~ sin(kx) and f ~ -cos(kx) so we'll need to multiply f by a negative sign.         
c$$$         fp = -enorm(i)*(f + x*fp)
c$$$         f = -enorm(i)*Rmatch*f
c$$$         gp = enorm(i)*(g + x*gp)
c$$$         g = enorm(i)*Rmatch*g
c         write(6,*) "leff(",i,") = ", leff(i)
          x = kVector(i)*Rmatch
c$$$          call sphbes(leff(i),x,f,g,fp,gp)
c$$$          fp = enorm(i)*(f + x*fp)
c$$$          f = enorm(i)*Rmatch*f
c$$$          gp = enorm(i)*(g + x*gp)
c$$$          g = enorm(i)*Rmatch*g
          P = (1d0 + Parity*(-1d0)**i)/2d0
          enorm(i) = dsqrt(Mass/(2d0*pi*kVector(i)))
          f = enorm(i)*sin(x + P*0.5d0*pi)
          fp = enorm(i)*kVector(i)*cos(x + P*0.5d0*pi)
          g = -enorm(i)*cos(x + P*0.5d0*pi)
          gp = enorm(i)*kVector(i)*sin(x + P*0.5d0*pi)
c          f = enorm(i)*cos(x - 0.25d0*pi)
c          fp = -enorm(i)*kVector(i)*sin(x - 0.25d0*pi)
c          g = enorm(i)*sin(x - 0.25d0*pi)
c          gp = enorm(i)*kVector(i)*cos(x - 0.25d0*pi)
         
          Wronskian = 1.0d0/(g*fp-gp*f) ! W(g,f)
     
!write(6,*) 'Wronskian check:    ',2.0d0*Wronskian/Pi
 36      format(5(e12.6,1x))
         
!     write(6,*) 'Bessel funcs:'
!     write(6,*) i,f,fp
!     write(6,*) i,g,gp
!     write(6,*)
c----------------------------------------------------------------------------------------------------
c     Here, the solution matrix F = f I - g J and we want to match to the form M = f - g K
c     Therefore M = F I^(-1) and K = J I^(-1)
c     Each element of J is computed W(F,f) = -W(g,f) J ==> J = -W(F,f)/W(g,f)
c     Similarly,     W(F,g) = W(f,g) I ==> I = W(F,g)/W(f,g) = -W(F,g)/W(g,f)
         do j = 1,NumOpenR
            IMat(i,j) = -(g*deriv(j)+gp)*Solution(i,j)*Wronskian
            JMat(i,j) = -(f*deriv(j)+fp)*Solution(i,j)*Wronskian
         enddo
         
      enddo
c     Compute the K-matrix
      ImatTemp = Imat      
      call SqrMatInv(ImatTemp,NumOpenR)
      Kmat = MATMUL(Jmat,ImatTemp)
c--------------------------------------------------------      
c     calculate S-matrix

      IJMat = (0.0d0,0.0d0)
 
       do i = 1,NumOpenR
        do j = 1,NumOpenR
         IJMat(i,j)   = IMat(j,i)-(0.0d0,1.0d0)*JMat(j,i)
         Smatrix(i,j) = IMat(j,i)+(0.0d0,1.0d0)*JMat(j,i)
        enddo
       enddo
       
       call zgesv(NumOpenR,NumOpenR,IJMat,NumOpenR,ipiv,Smatrix,NumOpenR,info)

c     Subtract 1 from the diagonal elements
c       do i = 1,NumOpenR
c        Smatrix(i,i) = Smatrix(i,i) - (1.0d0,0.0d0)
c       enddo

c     calculate T-matrix

       do i = 1,NumOpenR
        do j = 1,NumOpenR
         ImatAux(i,j) = Imat(j,i)
         TmatrixAux(i,j) = Jmat(j,i)
        enddo
       enddo

       call dgesv(NumOpenR,NumOpenR,ImatAux,NumOpenR,ipiv,TmatrixAux,NumOpenR,info)

       do i = 1,NumOpenR
        do j = 1,NumOpenR
         Tmatrix(i,j) = TmatrixAux(j,i)
        enddo
       enddo
       
c i -> final state
c j -> initial state
c       CrossSections = Kmat
       do i = 1,NumOpenR
          do j = 1,NumOpenR
!     In 1D cross section is:
c             CrossSections(i,j) = kvector(i)/kvector(j) *  dreal(dconjg(Smatrix(i,j))*Smatrix(i,j))
c     CrossSections(i,j) =  dreal(dconjg(Smatrix(i,j)-kdelta(i,j))*(Smatrix(i,j)-kdelta(i,j)))/4d0
             CrossSections(i,j) =  kvector(i)/kvector(j) * dreal(dconjg(Smatrix(i,j)-kdelta(i,j))*(Smatrix(i,j)-kdelta(i,j)))/4d0
c             CrossSections(i,j) =  dreal(dconjg(Smatrix(i,j))*Smatrix(i,j))
             
c       For Recombination 
c        CrossSections(i,j) = convertCGS_K3*dfloat(2*TotalAngMomentum+1)*
c     .                       192.d0*pi*pi/(Mass*kVector(j)**4)*
c     .                       dreal(dconjg(Smatrix(i,j))*Smatrix(i,j))
c       For Relaxation
c       CrossSections(i,j) = convertCGS_RX*dfloat(2*TotalAngMomentum+1)*
c     .                        pi/(Mass*kVector(j)**1)*
c     .                        dreal(dconjg(Smatrix(i,j))*Smatrix(i,j))
        enddo
       enddo 

       return
       end      

c************************************************************************************
      subroutine BoxMatching(NumOpenL,NumOpenR,psi_in,b,Goo,logderiv,solution,deriv)
      integer NumOpenChannels,NumOpen,NumOpenR,ierr
      integer Matz,NumCoef,NumOpenTot
      double precision b(NumOpenL),logderiv(NumOpenL+NumOpenR)
      double precision psi_in(NumOpenL,NumOpenL)
      double precision Goo(NumOpenL+NumOpenR,NumOpenL+NumOpenR)
      double precision left(2*NumOpenL+NumOpenR,2*NumOpenL+NumOpenR)
      double precision right(2*NumOpenL+NumOpenR,2*NumOpenL+NumOpenR)
      double precision areal(2*NumOpenL+NumOpenR)
      double precision aimag(2*NumOpenL+NumOpenR)
      double precision cpsi(2*NumOpenL+NumOpenR,2*NumOpenL+NumOpenR)
      double precision denom(2*NumOpenL+NumOpenR)
      double precision sum,solution(NumOpenR,NumOpenR)
      double precision coef(2*NumOpenL+NumOpenR,NumOpenR)
      double precision deriv(NumOpenR)

      NumOpenTot = NumOpenL+NumOpenR

      do i=1,NumOpenL
        do j=1,NumOpenL
         left(i,j) = psi_in(i,j)
         right(i,j) = 0.0d0
        enddo
       enddo
 
       do i=1,NumOpenL
        do j=1,NumOpenTot
         left(i,NumOpenL+j) = -Goo(i,j)
         right(i,NumOpenL+j) = 0.0d0
        enddo
       enddo
 
       do i=1,NumOpenL
        do j=1,NumOpenL
         left(i+NumOpenL,j) = -psi_in(i,j)*b(j)
         right(i+NumOpenL,j) = 0.0d0 
        enddo
       enddo
 
       do i=1,NumOpenL
        do j=1,NumOpenTot
         left(i+NumOpenL,NumOpenL+j) = -Goo(i,j)*logderiv(j)
         right(i+NumOpenL,NumOpenL+j) = 0.0d0
        enddo
       enddo
 
       do i=1,NumOpenR
        do j=1,NumOpenTot
         left(i+2*NumOpenL,j) = 0.0d0
         right(i+2*NumOpenL,j) = 0.0d0
        enddo
       enddo
 
       do i=1,NumOpenR
        do j=1,NumOpenTot
         left(i+2*NumOpenL,NumOpenL+j) =
     >       -Goo(i+NumOpenL,j)*logderiv(j)
         right(i+2*NumOpenL,NumOpenL+j) = Goo(i+NumOpenL,j)
        enddo
       enddo 

       NumCoef = NumOpenL + NumOpenTot
       matz = 1
       call rgg(NumCoef,NumCoef,left,right,areal,aimag,denom,matz,cpsi,ierr)
       !write(6,*)
       !write(6,*)'return from rgg, ierr=',ierr
       !write(6,*)

        iflag = 1
        do i=1,NumCoef
         if (abs(denom(i)) .gt. 1.0d-8) then  ! Select the non-infinite eigenvalues (NPM)
          deriv(iflag) = -areal(i)/denom(i) ! The normal log-derivative = -b
          do j=1,NumCoef
           coef(j,iflag) = cpsi(j,i) ! The D coefficients in Burke's thesis
          enddo
          iflag = iflag+1
         endif
        enddo   

        do i=1,NumOpenR
         do j=1,NumOpenR
          solution(i,j) = 0.0d0
          do k=1,NumOpentot
           solution(i,j) = solution(i,j) +
     >           Goo(i+NumOpenL,k)*coef(k+NumOpenL,j)
          enddo
         enddo
        enddo
 
       do i=1,NumOpenR
        sum = 0.0d0
        do j=1,NumOpenR
         sum = sum + solution(j,i)*solution(j,i)
        enddo
        sum = dsqrt(sum)
        do j=1,NumOpenR
         solution(j,i) = solution(j,i)/sum
        enddo
       enddo 

       return 
       end


c*************************************************************
      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      implicit real*8(a-h,o-z)
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=7,NUSE2=8)
CU    USES chebev
      double precision xx,c1(7),c2(8),chebev
      SAVE c1,c2
      DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4,
     *-3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,
     *-4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/

      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
      gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END
c*************************************************************
      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)
      implicit real*8 (a-h,o-z)
      INTEGER MAXIT
      double precision rj,rjp,ry,ryp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.d-16,FPMIN=1.d-30,MAXIT=200000,XMIN=2.d0,
     *PI=3.141592653589793d0)
CU    USES beschb
      INTEGER i,isign,l,nl
      DOUBLE PRECISION a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,
     *f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,
     *r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,
     *temp,w,x2,xi,xi2,xmu,xmu2

      if(x.le.0..or.xnu.lt.0.) then
         write(6,*) 'bad arguments in bessjy. Press enter to continue'
         read(5,*)
      endif
      if(x.lt.XMIN)then
        nl=int(xnu+.5d0)
      else
        nl=max(0,int(xnu-x+1.5d0))
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI
      isign=1
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=b-d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b-1.d0/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=c*d
        h=del*h
        if(d.lt.0.d0)isign=-isign
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
      write(*,*)
      write(*,*)'x =',x
      write(6,*) 'x too large in bessjy; try asymptotic expansion. Press enter to continue'
      read(5,*)
1     continue
      rjl=isign*FPMIN
      rjpl=h*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do 12 l=nl,1,-1
        rjtemp=fact*rjl+rjpl
        fact=fact-xi
        rjpl=fact*rjtemp-rjl
        rjl=rjtemp
12    continue
      if(rjl.eq.0.d0)rjl=EPS
      f=rjpl/rjl
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
        e=exp(e)
        p=e/(gampl*PI)
        q=1.d0/(e*PI*gammi)
        pimu2=0.5d0*pimu
        if(abs(pimu2).lt.EPS)then
          fact3=1.d0
        else
          fact3=sin(pimu2)/pimu2
        endif
        r=PI*pimu2*fact3*fact3
        c=1.d0
        d=-x2*x2
        sum=ff+r*q
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*(ff+r*q)
          sum=sum+del
          del1=c*p-i*del
          sum1=sum1+del1
          if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2
13      continue
        write(6,*) 'bessy series failed to converge. Press enter to continue'
        read(5,*)
2       continue
        rymu=-sum
        ry1=-sum1*xi2
        rymup=xmu*xi*rymu-ry1
        rjmu=w/(rymup-f*rymu)
      else
        a=.25d0-xmu2
        p=-.5d0*xi
        q=1.d0
        br=2.d0*x
        bi=2.d0
        fact=a*xi/(p*p+q*q)
        cr=br+q*fact
        ci=bi+p*fact
        den=br*br+bi*bi
        dr=br/den
        di=-bi/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        do 14 i=2,MAXIT
          a=a+2*(i-1)
          bi=bi+2.d0
          dr=a*dr+br
          di=a*di+bi
          if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
          fact=a/(cr*cr+ci*ci)
          cr=br+cr*fact
          ci=bi-ci*fact
          if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
          den=dr*dr+di*di
          dr=dr/den
          di=-di/den
          dlr=cr*dr-ci*di
          dli=cr*di+ci*dr
          temp=p*dlr-q*dli
          q=p*dli+q*dlr
          p=temp
          if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
14      continue
        write(6,*) 'cf2 failed in bessjy'
        write(6,*) "Press enter to continue."
        read(5,*)
3       continue
        gam=(p-f)/q
        rjmu=sqrt(w/((p-f)*gam+q))
        rjmu=sign(rjmu,rjl)
        rymu=rjmu*gam
        rymup=rymu*(p+q/gam)
        ry1=xmu*xi*rymu-rymup
      endif
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
      do 15 i=1,nl
        rytemp=(xmu+i)*xi2*ry1-rymu
        rymu=ry1
        ry1=rytemp
15    continue
      ry=rymu
      ryp=xnu*xi*rymu-ry1

 123  return
      END
c**************************************************************
      Double Precision FUNCTION chebev(a,b,c,m,x)
      implicit real*8(a-h,o-z)
      INTEGER m
      double precision a,b,x,c(m)
      INTEGER j
      double precision d,dd,sv,y,y2

      if ((x-a)*(x-b).gt.0.) then
         write(6,*) "error in chebev."
         write(6,*) "Press enter to continue."
         read(5,*)
      endif
      d=0.
      dd=0.
      y=(2.*x-a-b)/(b-a)
      y2=2.*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5*c(1)
      return
      END
c**************************************************************
      SUBROUTINE sphbes(n,x,sj,sy,sjp,syp)
      implicit real*8(a-h,o-z)
      double precision sj,sjp,sy,syp,x,n
CU    USES bessjy
      double precision factor,order,rj,rjp,ry,ryp,RTPIO2

      pi = acos(-1.0d0)
      rtpio2 = dsqrt(0.5d0*pi)
      if(n.lt.0.d0.or.x.le.0.d0) then
         write(6,*) "error in sphbes."
         write(6,*) "Press enter to continue."
         read(5,*)
      endif
      order=n+0.5d0
      if (x.le.1.86d5) then
      call bessjy(x,order,rj,ry,rjp,ryp)
      else
      call BesselAsymNew(order,x,rj,ry,rjp,ryp)
      endif
      factor=RTPIO2/sqrt(x)
      sj=factor*rj
      sy=factor*ry
      sjp=factor*rjp-sj/(2.d0*x)
      syp=factor*ryp-sy/(2.d0*x)
      return
      END




c**************************************************************

      subroutine BesselAsymNew(v,x,sj,sy,sjp,syp)
      integer k,b,kstop
      double precision v,x,sj,sy,sjp,syp,u
      double precision P,Q,R,S,pi,factor,chi
      double precision denomP,denomQ,denomP2,denomQ2
      double precision Pfunction,Qfunction,sign
      double precision coefR,coefS,factrl

c     asymptotic expansions from Abramowitz and Stegun, pg. 364
c     following assumes v = n+1/2 where n=0,1,2,......

      pi = dacos(-1.0d0)
      factor = dsqrt(2.0d0/(pi*x))
      chi = x - pi*(0.5d0*v + 0.25d0)
      u = 4.0d0*v*v

      P = 0.0d0
      R = 0.0d0
      do k=0,int(v/2.d0-0.25d0)
         sign = (-1.0d0)**k
         denomP = (2.0d0*x)**(2*k)
         denomP2 = factrl(2*k)
         coefR = (u-1.0d0+16.0d0*dfloat(k*k))/
     .       (u-(4.0d0*dfloat(k)-1.0d0)*(4.0d0*dfloat(k)-1.0d0))
         Pfunction = 1.0d0
         do b=0,4*k-1
          Pfunction = Pfunction*(v-0.5d0+float(2*k-b))
         enddo
      P = P + sign*Pfunction/(denomP*denomP2)
      R = R + sign*coefR*Pfunction/(denomP*denomP2)
      enddo

      Q = 0.0d0
      S = 0.0d0
      do k=0,int(v/2.d0-0.75d0)
         sign = (-1.0d0)**k
         denomQ = denomP*(2.0d0*x)
         denomQ2 = factrl(2*k+1)
         coefS = (u+4.0d0*(2.0d0*dfloat(k)+1.0d0)*(2.0d0*dfloat(k)+1.0d0)-
     .           1.0d0)/(u-(4.0d0*dfloat(k)+1.0d0)*(4.0d0*dfloat(k)+1.0d0))
         Qfunction = 1.0d0
         do b=0,4*k+1
         Qfunction = Qfunction*(v+0.5d0+float(2*k-b))
         enddo       
      Q = Q + sign*Qfunction/(denomQ*denomQ2)
      S = S + sign*coefS*Qfunction/(denomQ*denomQ2)
      enddo


      if (v.eq.0.5d0) then
         R = 1.d0
         S = 1.d0/(2.d0*x)
      endif

      sj = factor*(P*dcos(chi) - Q*dsin(chi)) 
      sy = factor*(P*dsin(chi) + Q*dcos(chi))
      sjp=-factor*(R*dsin(chi) - S*dcos(chi))
      syp= factor*(R*dcos(chi) - S*dsin(chi))

      return
      end

c**************************************************************
      subroutine BesselAsym(v,x,sj,sy,sjp,syp)
      integer k,b
      double precision v,x,sj,sy,sjp,syp,u
      double precision P,Q,R,S,pi,factor,chi
      double precision denomP,denomQ,denomP2,denomQ2
      double precision Pfunction,Qfunction,sign
      double precision coefR,coefS,factrl

c asymptotic expansions from Abramowitz and Stegun, pg. 364
c following assumes v = n+1/2 where n=0,1,2,......

      pi = dacos(-1.0d0)
      factor = dsqrt(2.0d0/(pi*x))
      chi = x - pi*(0.5d0*v + 0.25d0)
      u = 4.0d0*v*v

      P = 0.0d0
      Q = 0.0d0
      R = 0.0d0
      S = 0.0d0
      do k=0,nint(v-0.5d0)
       !write(6,*)'Doing Loop'
       sign = (-1.0d0)**k
       denomP = (2.0d0*x)**(2*k)
       denomQ = denomP*(2.0d0*x)
       denomP2 = factrl(2*k)
       denomQ2 = factrl(2*k+1)
       coefR = (u-1.0d0+16.0d0*float(k*k))/
     > (u-(4.0d0*float(k)-1.0d0)*(4.0d0*float(k)-1.0d0))
       coefS = (u+4.0d0*(2.0d0*float(k)+1.0d0)*(2.0d0*float(k)+1.0d0)-
     >          1.0d0)/(u-(4.0d0*float(k)+1.0d0)*(4.0d0*float(k)+1.0d0))
       !write(6,*)'CoefS',coefR,coefS,k

       Pfunction = 1.0d0
       do b=0,4*k-1
        Pfunction = Pfunction*(v-0.5d0+float(2*k-b))
       enddo
       Qfunction = 1.0d0
       do b=0,4*k+1
        Qfunction = Qfunction*(v+0.5d0+float(2*k-b))
       enddo       
     
       P = P + sign*Pfunction/(denomP*denomP2)
       Q = Q + sign*Qfunction/(denomQ*denomQ2)
       R = R + sign*coefR*Pfunction/(denomP*denomP2)
       S = S + sign*coefS*Qfunction/(denomQ*denomQ2)
      enddo

      !write(6,*)'R,S',R,S*2.0d0*x
      sj = factor*(P*dcos(chi) - Q*dsin(chi)) 
      sy = factor*(P*dsin(chi) + Q*dcos(chi))
      sjp=-factor*(R*dsin(chi) - S*dcos(chi))
      syp= factor*(R*dcos(chi) - S*dsin(chi))
      
      return
      end
c**************************************************************
      double precision FUNCTION factrl(n)
      INTEGER n
C     USES gammln
      INTEGER j,ntop
      double precision a(33),gammln
      SAVE ntop,a
      DATA ntop,a(1)/0,1./
      if (n.lt.0) then
         write(6,*) "negative factorial in factrl. "
         write(6,*) "Press enter to continue."
         read(5,*)
      else if (n.le.ntop) then
         factrl=a(n+1)
      else if (n.le.32) then
         do 11 j=ntop+1,n
            a(j+1)=j*a(j)
11      continue
        ntop=n
        factrl=a(n+1)
      else
         factrl=exp(gammln(float(n)+1.0d0))
      endif
      return
      END
c**************************************************************
      double precision FUNCTION gammln(xx)
      double precision xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
c**************************************************************      
      subroutine SqrMatInv(A, N)
      implicit none
      integer N,info,lwk
      integer, allocatable :: ipiv(:)
      double precision, allocatable :: work(:)
      double precision A(N,N)
      allocate(ipiv(N))
      call dgetrf(N, N, A, N, ipiv, info)
      allocate(work(1))
      lwk = -1
      call dgetri(N, A, N, ipiv, work, lwk, info)
      lwk = work(1)
      deallocate(work)
      allocate(work(lwk))
      call dgetri(N, A, N, ipiv, work, lwk, info)
      deallocate(ipiv,work)
      end
c**************************************************************      
      SUBROUTINE GridMaker(grid,numpts,E1,E2,scale)
      implicit none
      DOUBLE PRECISION grid(numpts)
      DOUBLE PRECISION E1,E2,LE1,LE2,DE,LDE
      INTEGER numpts, iE
      CHARACTER(LEN=*), INTENT(IN) :: scale
!--------------------------------------------
! Linear grid:
!--------------------------------------------
      grid(1)=E1
      IF((scale.EQ."linear").and.(numpts.gt.1)) THEN
      DE=(E2-E1)/DBLE(numpts-1)
      DO iE=1,numpts
         grid(iE) = E1 + (iE-1)*DE
      ENDDO
      ENDIF
!--------------------------------------------
! Log grid:
!--------------------------------------------
      IF((scale.EQ."log").and.(numpts.gt.1)) THEN
      LE1=dlog(E1)
      LE2=dlog(E2)

      LDE=(LE2-LE1)/DBLE(numpts-1d0)
      DO iE=1,numpts
        grid(iE) = dexp(LE1 + (iE-1)*LDE)
        !        write(6,*) LE1, LE2, LDE, grid(iE)
      ENDDO
      ENDIF
!--------------------------------------------
! NegLog grid:  (use this for a log spacing of negative numbers)
!--------------------------------------------
       IF((scale.EQ."neglog").and.(numpts.gt.1)) THEN
          LE1=dlog(-E2)
          LE2=dlog(-E1)
          
          LDE=(LE2-LE1)/DBLE(numpts-1d0)
          DO iE=1,numpts
             grid(iE) = -dexp(LE2 - (iE-1)*LDE)
!        write(6,*) iE, grid(iE)
          ENDDO
       ENDIF
!--------------------------------------------
! Quadratic grid:
!--------------------------------------------
       IF((scale.EQ."quadratic").and.(numpts.gt.1)) THEN
          DE=(E2-E1)
          DO iE=1,numpts
             grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**2*DE
          ENDDO
       ENDIF
!--------------------------------------------
! Cubic grid:
!--------------------------------------------
       IF((scale.EQ."cubic").and.(numpts.gt.1)) THEN
          DE=(E2-E1)
          DO iE=1,numpts
             grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**3*DE
          ENDDO
       ENDIF
!--------------------------------------------
! quartic grid:
!--------------------------------------------
       IF((scale.EQ."quartic").and.(numpts.gt.1)) THEN
          DE=(E2-E1)
          DO iE=1,numpts
             grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**4*DE
          ENDDO
       ENDIF
!--------------------------------------------
! sqrroot grid:
!--------------------------------------------
       IF((scale.EQ."sqrroot").and.(numpts.gt.1)) THEN
          DE=(E2-E1)
          DO iE=1,numpts
             grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**0.5d0*DE
          ENDDO
       ENDIF
!--------------------------------------------
! cuberoot grid:
!--------------------------------------------
       IF((scale.EQ."cuberoot").and.(numpts.gt.1)) THEN
          DE=(E2-E1)
          DO iE=1,numpts
             grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**0.33333333333333333333333333333d0*DE
          ENDDO
       ENDIF
       
       END SUBROUTINE GridMaker

      double precision function kdelta(mch,nch)
      implicit none
      integer mch,nch
      if (mch.eq.nch) then
         kdelta = 1.0d0
      else 
         kdelta = 0.0d0
      endif
      return
      end function kdelta
      
      character(len=20) function str(k)
!     "Convert an integer to string."
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
      end function str
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! This is the asymptotic form of the potential and the Q-matrix for states that go collision channels (exculded are molecular ion channels)
      double precision function AsymptoticVQ(mu,m,n,R)
      implicit none
      integer m,n
      double precision R, VQ, U, mu
      double precision, external :: kdelta

c     The asymptotic form of Qtilde
      VQ = kdelta(m,n)*dble(2 + n**2 + n*(3 + iabs(n-1)))
     .     -kdelta(m,n-4)*sqrt(dble(n*(n-1)*(n-2)*(n-3)))
     .     -kdelta(m,n+4)*sqrt(dble((n+1)*(n+2)*(n+3)*(n+4)))

      VQ = VQ*0.25d0/(2d0*mu*R**2)

c     Add the asymptotic form of Un(R) and the -1/4/(2 mu R^2) from removal of 1st deriv terms
      if(m.eq.n) then
         VQ = VQ + dble(n) + 0.5d0 - dble(1 + 2*n*(n+1))*0.25d0/(2d0*mu*R*R) - 0.25d0/(2d0*mu*R*R) 
      endif
      AsymptoticVQ = VQ
      
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function AsymptoticP(m,n,R)
      implicit none

      integer m,n
      double precision R
      double precision, external :: kdelta
      if(m.ne.n) then
         AsymptoticP = kdelta(m,n-2)*0.5d0*sqrt(dble(n*n-1)) - kdelta(m,n+2)*0.5d0*sqrt(dble(n+1)*(n+2))
         AsymptoticP = AsymptoticP/R
      else
         AsymptoticP = 0d0
      endif
      
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      SUBROUTINE printmatrix(M,nr,nc,file)
      IMPLICIT NONE
      INTEGER nr,nc,file,j,k
      DOUBLE PRECISION M(nr,nc)
      
      DO j = 1,nr
         WRITE(file,40) (M(j,k), k = 1,nc)
      ENDDO
      
 20   FORMAT(1P,100D16.8)
 30   FORMAT(100D14.4)
 40   FORMAT(100F20.10)
      
      END SUBROUTINE printmatrix

            
      SUBROUTINE zprintmatrix(M,nr,nc,file)
      IMPLICIT NONE
      INTEGER nr,nc,file,j,k
      DOUBLE complex M(nr,nc)
      
      DO j = 1,nr
         WRITE(file,40) (M(j,k), k = 1,nc)
      ENDDO
      
 20   FORMAT(1P,100D16.8)
 30   FORMAT(100D14.4)
 40   FORMAT(100F20.10)
      
      END SUBROUTINE zprintmatrix
      
