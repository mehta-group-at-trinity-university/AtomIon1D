      Program FitProg
!implicit real*8 (a-h,o-z)
      implicit none
      double precision ma,mi
      integer ncount, MaxNumChannels,MaxNumDataPts
      double precision, parameter :: Rdelay=10.d0
      integer, parameter :: NumFitTerms=3
!     This is for the atom-ion problem so use the Yb, Li masses
      integer mrdim,NumDataPoints
      integer*2 ifa,ifb,mfa,mfb,onoff, i, j, k, ii, iR

      double precision, allocatable :: VQmat(:,:,:), Pmat(:,:,:)
      double precision, allocatable :: VQExponents(:,:,:),PmatExp(:,:,:)
      double precision, allocatable :: xdata(:), ydata(:), xPoints(:), yPoints(:)
      double precision, allocatable :: VQFit(:,:,:), PmatFit(:,:,:)
      double precision, allocatable :: Errors(:), Leff(:)
      double precision, allocatable :: FitRangeMaxOpt(:),FitRangeMinOpt(:)
      double precision Mass, massfactor, ang, dummy, Cfactor
      double precision FitCoeffs(NumFitTerms)
      integer NumFitPoints
      double precision rmsError,FitRangeMin,FitRangeMax
      double precision kExponents(NumFitTerms)
      double precision FitRangeMinVar,aux
      character*64 DummyChar

c****************************************************************
c     Note: All open channels should precede any closed channels in *
c     the matrix setup.                                             *
c****************************************************************
c      open(100,file='DiabaticEnergies-LR=0-0.dat',status='old')
c      open(200,file='Pmat-LR=0-0.dat',status='old')
c      open(300,file='VQmat-LR=0-0.dat',status='old')
      open(100,file='DiabaticEnergies.dat',status='old')
      open(200,file='Pmat.dat',status='old')
      open(300,file='VQmat.dat',status='old')

      read(100,*) DummyChar, MaxNumDataPts, MaxNumChannels

      allocate(VQmat(MaxNumDataPts,MaxNumChannels,MaxNumChannels))
      allocate(Pmat(MaxNumDataPts,MaxNumChannels,MaxNumChannels))
      allocate(xdata(MaxNumDataPts),ydata(MaxNumDataPts))
      allocate(VQFit(MaxNumChannels,MaxNumChannels,NumFitTerms))
      allocate(PmatFit(MaxNumChannels,MaxNumChannels,NumFitTerms))
      allocate(xPoints(MaxNumDataPts),yPoints(MaxNumDataPts))
      allocate(Errors(MaxNumDataPts))
      allocate(VQExponents(MaxNumChannels,MaxNumChannels,NumFitTerms))
      allocate(FitRangeMaxOpt(MaxNumChannels),FitRangeMinOpt(MaxNumChannels))
      allocate(PmatExp(MaxNumChannels,MaxNumChannels,NumFitTerms))
      allocate(Leff(MaxNumChannels))
      
      open(18, file = 'Fit.data')
      open(19, file = 'FitLeff.data')
      write(18,*) MaxNumChannels,MaxNumDataPts
c     mass = (m1*m2*m3*m4/(m1+m2+m3+m4))**(1.d0/3.d0)
      mi=170.9363230d0
      ma=6.015122d0
      ma = ma/mi
      mi = 1d0
      Mass = dsqrt(ma)
      Cfactor = 8.0d0*Mass
      
      massfactor = -0.5d0/Mass
      mrdim = MaxNumChannels
      write(6,*) "Mass = ", Mass
      write(6,*) "Cfactor = ", Cfactor
      write(6,*) "massfactor = ", massfactor

      write(6,*) "MaxNumDataPts = ", MaxNumDataPts, "mrdim = ", mrdim

      read(100,*)
      do i=1,MaxNumDataPts
         read(200,*) xdata(i)
         read(300,*)
         do j=1,mrdim
            read(200,11)(Pmat(i,j,k),k=1,mrdim)
            read(300,11)(VQmat(i,j,k),k=1,mrdim)
         enddo
         !write(6,*) xdata(i), Pmat(i,1,2)
      enddo

 11   format(1P,100e22.12)
 20   format(1P,100e20.12)

      FitRangeMax = xdata(MaxNumDataPts) 
      NumDataPoints = MaxNumDataPts
      
      Leff = 0d0

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FitRangeMin = FitRangeMax - Rdelay
c     write(19,*)'FitRangeMin =   OPTIMIZED'
      write(19,*)'FitRangeMin =', FitRangeMin
      write(19,*)'FitRangeMax =',FitRangeMax
      write(19,*)


c$$$      do i=1,NumDataPoints
c$$$         write(3000,*)xdata(i)
c$$$         write(3005,*)xdata(i)
c$$$         do j=1,mrdim
c$$$            write(3000,3001)(VQmat(i,j,k),k=1,mrdim)
c$$$            write(3005,3001)(Pmat(i,j,k),k=1,mrdim)
c$$$         enddo
c$$$      enddo
 3001 format(100(e20.12,1x))

      open(unit = 777, file = "QuickPmat.dat")
      open(unit = 888, file = "QuickVQmat.dat")
      do i=1,MaxNumDataPts
         write(777,11)xdata(i),((Pmat(i,j,k),j=k,3),k=1,3)
         write(888,11)xdata(i),((VQmat(i,j,k),j=k,3),k=1,3)
      enddo   

      
      DO j=1,mrdim

         ncount = 0 
         FitRangeMinVar = FitRangeMax-Rdelay

 123     NumFitPoints=0
         DO iR=1,NumDataPoints
            IF (xdata(iR).LT.FitRangeMax.AND.xdata(iR).GE.FitRangeMinVar) THEN
               NumFitPoints=NumFitPoints+1
               xPoints(NumFitPoints)=xdata(iR)
               yPoints(NumFitPoints)=VQmat(iR,j,j)
            ENDIF
         ENDDO

         kExponents(1)=0.0d0
         kExponents(2)=2.0d0
         kExponents(3)=3.0d0
c     Uu-(1/2mu)*Quu : (BC-BC)
c     Uu-(1/2mu)*Quu : (CC-CC)
         
         CALL FitTail(NumFitTerms,kExponents,xPoints,yPoints,Errors,
     .        NumFitPoints,MaxNumDataPts,rmsError,FitCoeffs)

         DO ii=1,NumFitTerms
            VQFit(j,j,ii)=FitCoeffs(ii)
            VQExponents(j,j,ii)=kExponents(ii)
         ENDDO

      ENDDO

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      FitRangeMinOpt = FitRangeMin

      DO j=1,mrdim
         DO k=j,mrdim
            NumFitPoints=0
            DO iR=1,NumDataPoints
               IF(xdata(iR).LT.FitRangeMax .AND. xdata(iR).GE.FitRangeMinOpt(j))THEN
                  NumFitPoints=NumFitPoints+1
                  xPoints(NumFitPoints)=xdata(iR)
                  yPoints(NumFitPoints)=VQmat(iR,j,k)
               ENDIF
            ENDDO

c     Uu-(1/2mu)*Quu 
            IF (j.EQ.k) THEN    ! Diagonal initial guesses
               kExponents(1)=0.0d0
               kExponents(2)=2.0d0
               kExponents(3)=3.0d0
c     DO WE WANT TO HAVE A DIFFERENT FIT FOR THE "Bound States?"  If so we need to add a conditional statement here               
c     (CC-CC)
c     kExponents(1)=-2.0d0 ! atom-ion bound channels go as R^2 at long range so fit to a harmonic potential
c     kExponents(2)=0.0d0
c     kExponents(3)=-3.0d0
            ENDIF
            
c     Quv : 
            IF (j.NE.k) THEN
               kExponents(1)=2.0d0
               kExponents(2)=3.0d0
               kExponents(3)=4.0d0
            ENDIF


            CALL FitTail(NumFitTerms,kExponents,xPoints,yPoints,Errors,
     >           NumFitPoints,MaxNumDataPts,rmsError,FitCoeffs)
            
            DO ii=1,NumFitTerms
               VQFit(j,k,ii)=FitCoeffs(ii)
               VQExponents(j,k,ii)=kExponents(ii)
               VQFit(k,j,ii)=FitCoeffs(ii)
               VQExponents(k,j,ii)=kExponents(ii)
            ENDDO

            if (j .eq. k)then
               ang = 0.5d0*(-1.0d0+dsqrt(1.0d0+Cfactor*abs(VQFit(j,k,2))))
               write(19,3313)j,VQFit(j,k,1),(dabs(VQFit(j,k,2))-Leff(j)*(Leff(j)+1.d0)/2.d0/Mass),ang,FitRangeMinOpt(j)
            endif
            write(18,5151) j,k,(VQExponents(j,k,ii),ii=1,NumFitTerms)
            write(18,5152) (VQFit(j,k,ii),ii=1,NumFitTerms),rmsError
 5151       format(2(1x,i2),3(3x,f8.2), ' **VQfit parameters')
c     5152   format(3(1x,1pd14.7),5x,1pd10.3)
 5152       format(3(1x,1pe20.12),5x,1pe10.3)

c     ** Next fit the P-matrix elements when j.ne.k
            IF(j.NE.k) THEN

               NumFitPoints=0
               DO iR=1,NumDataPoints
                  IF(xdata(iR).LT.FitRangeMax .AND. xdata(iR).GE.FitRangeMinOpt(j)) THEN
                     NumFitPoints=NumFitPoints+1
                     xPoints(NumFitPoints)=xdata(iR)
                     yPoints(NumFitPoints)=Pmat(iR,j,k)
                  ENDIF
               ENDDO

               kExponents(1)=1.0d0
               kExponents(2)=2.0d0
               kExponents(3)=3.0d0

               CALL FitTail(NumFitTerms,kExponents,xPoints,yPoints,Errors,
     >              NumFitPoints,MaxNumDataPts,rmsError,FitCoeffs)


               DO ii=1,NumFitTerms
                  PmatFit(j,k,ii)=FitCoeffs(ii)
                  PmatExp(j,k,ii)=kExponents(ii)
                  PmatFit(k,j,ii)=-FitCoeffs(ii)
                  PmatExp(k,j,ii)=kExponents(ii)
               ENDDO
               write(18,5153) j,k,(kExponents(ii),ii=1,NumFitTerms)
               write(18,5152) (FitCoeffs(ii),ii=1,NumFitTerms),rmsError
 5153          format(2(1x,i2),3(3x,f8.2), ' **PmatFit parameters')
            ENDIF

         ENDDO
      ENDDO
      
      close(100)
      close(200)
      close(300)
      close(400)

 331  format(4(e16.8,1x))
 3313 format(I4,1X,(e21.13),1X,(e21.13),1X,f14.8,1X,f12.4)
c33   13 format(I4,1X,e16.8,1X,e16.8,1X,f16.8,1X,f16.8)

      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     --  This fitting routine is written in Fortran 90
      Subroutine FitTail(nTerms,kExp,xPoints,yPoints,Errors,nData,MaxPts,rmsError,cc)
      implicit real*8(a-h,o-z)
      integer nData,nTerms,MaxPts
      double precision cc(nTerms),Errors(MaxPts)
      double precision xPoints(MaxPts),yPoints(MaxPts),rmsError
      double precision r0,fitfuncs,dfloat
      double precision DELik,qmax,factor,sum1,sum2,error1,error2
      double precision kExp(nTerms)
      double precision, allocatable :: xMat(:,:),aVec(:)
      double precision, allocatable :: Aux(:,:),SaveMatrix(:,:),SecondVec(:)
      double precision, allocatable :: SolutionVec(:),SaveVec(:)
      integer n,i,j,k,Row,Col,nCols,MasterCol

c     c    external fitfuncs

      n = nTerms
      r0=xPoints(nData/2)

      nCols=2*n+1
      allocate(xMat(n,n),aVec(n),Aux(n,nCols),SaveMatrix(n,n),SecondVec(n))
      allocate(SolutionVec(n),SaveVec(n))

C     --  now set up the matrix xMat and store it in the left-most part of Aux:

      do kp=1,n
         do k=1,n
            xMat(kp,k)=0.d0
            do i=1,nData
               xMat(kp,k) = xMat(kp,k)+fitfuncs(kExp(kp),xPoints(i),r0)*
     >              fitfuncs(kExp(k),xPoints(i),r0)
            enddo
            Aux(kp,k) = xMat(kp,k)
         enddo
         Aux(kp,n+1) = 0.d0
         do i=1,nData
            Aux(kp,n+1)=Aux(kp,n+1)+yPoints(i)*fitfuncs(kExp(kp),xPoints(i),r0)
         enddo
         SaveVec(kp)=Aux(kp,n+1)
      enddo

      do i=1,n
         do k=1,n
            DELik=0.d0
            if(i.eq.k) DELik=1.d0
            Aux(i,k+n+1)=DELik
         enddo
      enddo

 155  format(4(1x,1pd14.7))

      SaveMatrix=xMat

C     -- Now solve the necessary linear system by performing 
C     ---- row operations on this matrix Aux:
C     -- The algorithm I will adopt here is a simple one without "pivoting",
C     ---  "Gaussian elimination with backsubstitution" (see Numerical Recipes)
      
C     -- Step 1: stabilize by rescaling matrices so that the maximum element of xMat 
C     ----  in any given row is unity.

      DO Row=1,n
         qmax=0.d0
         DO Col=1,n
            IF(dabs(Aux(Row,Col)).GT.qmax) THEN
               qmax = dabs(Aux(Row,Col))
            ENDIF
         ENDDO

         DO Col=1,nCols
            Aux(Row,Col) = Aux(Row,Col)/qmax
         ENDDO
      ENDDO

C     -- Now perform the triangularization of the matrix:

      DO MasterCol=1,n-1
         IF(dabs(Aux(MasterCol,MasterCol)).LT.1.d-20) then
            write(6,*) 'Division by zero problem -- need to implement pivoting'
            write(6,*) MasterCol,Aux(MasterCol,MasterCol)
            STOP
         ENDIF	  
         
         DO Row=MasterCol+1,n
            factor = Aux(Row,MasterCol)/Aux(MasterCol,MasterCol)
            DO Col=1,nCols
               Aux(Row,Col) = Aux(Row,Col)-factor*Aux(MasterCol,Col)
            ENDDO
         ENDDO
      ENDDO

C     -- Now it should be triangular, so we can "back-substitute".
C     ---- First get the desired solution vector:

      aVec(n) = Aux(n,n+1)/Aux(n,n)
      DO Row = n-1,1,-1
         IF(dabs(Aux(Row,Row)).LT.1.d-20) then
            write(6,*) '**Division by zero problem** '
            write(6,*) Row,Aux(Row,Row)
            STOP
         ENDIF
         aVec(Row) = Aux(Row,n+1) 
         DO i=Row+1,n
            aVec(Row) = aVec(Row) - Aux(Row,i)*aVec(i)
         ENDDO
         aVec(Row) = aVec(Row)/Aux(Row,Row)
      ENDDO

      DO Col=1,n
         cc(Col)=r0**kExp(Col)*aVec(Col)
      ENDDO

C     -- Next cycle through the modified column vectors of the unit matrix to
C     -- construct the inverse of the original matrix in case it is desired
C     -- To save storage, overwrite the original matrix.

C     -- Compare the final accuracy of this calculation by testing the ability of
C     --- the polynomial to reproduce the data points:

      rmsError=0.d0
      DO Row=1,nData
         sum1=0.d0
         DO i=1,n
            sum1=sum1+cc(i)*fitfuncs(kExp(i),xPoints(Row),r0)/r0**kExp(i)
         ENDDO
         error1=sum1-yPoints(Row)
         Errors(Row)=error1
         rmsError = rmsError + error1**2
      ENDDO	
      rmsError = dsqrt(rmsError/dfloat(nData))
      deallocate(aVec,xMat,Aux,SaveMatrix,SecondVec,SolutionVec,SaveVec)
      RETURN
      END
      
      double precision function fitfuncs(kExp,x,r0)
      double precision kExp
      double precision x,r0
      fitfuncs = (r0/x)**kExp
      return
      end
      
