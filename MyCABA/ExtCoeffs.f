      Program ExtCoeffs

      implicit none

      integer i,j,k,iR,NumBound,ii,jj,nn
      integer NumStates,RSteps,NumRPoints
      integer SVDFlag,NumSVDSectors,NumSVDPoints
      integer ncount,NumFitTerms,NumFitPoints
      double precision m1,m2,m3,mu,Cfactor,ang
      double precision RFirst,RLast,RiSVD,RfSVD
      double precision FitRangeMin,FitRangeMax,FitRangeMinVar,FitRangeMaxVar
      double precision, allocatable:: U(:,:),R(:),Pmat(:,:,:),VQmat(:,:,:)
      double precision, allocatable:: thresh(:),leff(:),Ropt(:),lambda(:)
      double precision, allocatable:: leffmin(:),Rmin(:),xmin(:),xleff(:)
      double precision, allocatable:: xbind(:),bindmin(:)

      double precision, allocatable:: xPoints(:),yPoints(:)
      double precision, allocatable:: VQFit(:,:,:),PmatFit(:,:,:)
      double precision, allocatable:: VQExponents(:,:,:),PmatExp(:,:,:)
      double precision, allocatable:: kExponents(:),FitCoeffs(:),Errors(:)
      double precision rmsError

      NumFitTerms = 3
      open(18,file='Fit.data')
      open(19,file='FitLeff.data')

c     readin R grid info
      open(777,file='infoExtCoefs.data',status='old')
      read(777,*)
      read(777,*) m1,m2,m3
      read(777,*)
      read(777,*) NumStates,SVDFlag
      read(777,*)
      read(777,*) RiSVD,RfSVD,NumSVDSectors,NumSVDPoints
      read(777,*)
      read(777,*) RFirst,RLast,RSteps

      if (SVDFlag.ne.0) then
         NumRPoints = NumSVDSectors*NumSVDPoints+RSteps
      else
         NumRPoints = RSteps
      endif

c     readin potentials and couplings
      allocate(R(NumRPoints),thresh(NumStates))
      allocate(U(NumRPoints,NumStates),Pmat(NumRPoints,NumStates,NumStates))
      allocate(VQmat(NumRPoints,NumStates,NumStates))

      open(100,File='fort.100')
      open(99,File='fort.99')
      open(98,File='fort.98')

      do iR = 1,NumRPoints
         read(100,11) R(iR),(U(iR,j),j=1,NumStates)
         if (SVDFlag.ne.0) then
            if (iR.gt.NumRPoints-RSteps) then
            read(98,*) R(iR)
            read(99,*) R(iR)
            do j=1,NumStates
               read(98,20)(Pmat(iR,j,k),k=1,NumStates)
               read(99,20)(VQmat(iR,j,k),k=1,NumStates)
            enddo
            endif
         else   
         read(98,*) R(iR)
         read(99,*) R(iR)
         do j=1,NumStates
            read(98,20)(Pmat(iR,j,k),k=1,NumStates)
            read(99,20)(VQmat(iR,j,k),k=1,NumStates)
         enddo
         endif   
      enddo

      close(100)
      close(98)
      close(99)

      mu = dsqrt(m1*m2*m3/(m1+m2+m3))
      Cfactor = 8.0d0*mu

      Pmat = -0.5d0/mu*Pmat
      VQmat = -0.5d0/mu*VQmat
      
      do iR = 1,NumRPoints   
         write(120,20) (Pmat(iR,j,6),j=1,NumStates)
         write(130,20) (VQmat(iR,j,6),j=1,NumStates)
      enddo

      do i = 1,NumStates
         do iR = 1,NumRPoints   
            VQmat(iR,i,i) = U(iR,i)+VQmat(iR,i,i)
         enddo   
      enddo  

      do iR = 1,NumRPoints   
         write(110,20) (VQmat(iR,j,j),j=1,NumStates)
      enddo

   11 format(f12.3,1P,600e20.12)
   20 format(1P,600e20.12)
   22 format(600e20.12)  

      stop
      end





