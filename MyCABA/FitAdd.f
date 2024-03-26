      Program FitProg
      implicit real*8 (a-h,o-z)
      double precision m1,m2,m3

      Parameter(MaxNumDataPts=140)
      integer NumDataPoints

      double precision K3(MaxNumDataPts),K3x(MaxNumDataPts)
      double precision xdata(-50:MaxNumDataPts)
      integer NumFitPoints,ir,jr,NumFitTerms
      double precision xPoints(MaxNumDataPts),yPoints(MaxNumDataPts)
      double precision Errors(MaxNumDataPts)
      double precision rmsError,FitRangeMin,FitRangeMax
      double precision, allocatable :: kExponents(:),FitCoeffs(:)
      double precision muK,K3Coeff,K3Ext

      muK = 3.166815355d-12

      do i=1,MaxNumDataPts
      read(25,*)xdata(i),K3(i),K3x(i)
      enddo
      
      K3 = 1.d0/K3

      Nadd = 40
      do i=1,Nadd
      xdata(i-Nadd) = xdata(1)*1.d-2*muK*(10.d0**(dble(i-1)*0.05d0))/muK
      enddo

      write(*,*)'NumFitTerms'
      read(*,*)NumFitTerms

      write(*,*)'FitRangeMin,FitRangeMax'
      read(*,*)FitRangeMin,FitRangeMax

      allocate(kExponents(NumFitTerms),FitCoeffs(NumFitTerms))

      NumFitPoints=0
      DO ir=1,MaxNumDataPts
        IF(xdata(ir).LT.FitRangeMax.AND.xdata(ir).GE.FitRangeMin)THEN
        NumFitPoints=NumFitPoints+1
        xPoints(NumFitPoints)=xdata(ir)
        yPoints(NumFitPoints)=K3(ir)
        ENDIF
      ENDDO

      do ir=1,NumFitTerms
      kExponents(ir)=-2.0d0*(ir-1)
      enddo

      CALL FitTail(NumFitTerms,kExponents,xPoints,yPoints,Errors,
     .     NumFitPoints,MaxNumDataPts,rmsError,FitCoeffs)
      
      write(*,*)(FitCoeffs(ir),ir=1,NumFitTerms)

      open(100,file='ExtraPol.data')

      write(100,*)'FitRangeMin =',FitRangeMin
      write(100,*)'FitRangeMax =',FitRangeMax
      write(100,*)
      do ir=1,NumFitTerms
      write(100,*)kExponents(ir),FitCoeffs(ir)
      enddo

      DO ir=1-Nadd,MaxNumDataPts
        IF(xdata(ir).LT.FitRangeMax.AND.xdata(ir).GE.0.d0)THEN
        K3Ext = 0.d0
        do jr=1,NumFitTerms 
           K3Ext = K3Ext+FitCoeffs(jr)*xdata(ir)**(-kExponents(jr))
        enddo   
        write(250,*)xdata(ir),1.d0/K3Ext
        ELSE
        write(250,*)xdata(ir),1.d0/K3(ir)
        ENDIF
      ENDDO

      DO ir=1-Nadd,MaxNumDataPts
         if (ir.le.0) then
         write(251,*)xdata(ir),0.d0   
         else
         write(251,*)xdata(ir),K3x(ir)
         endif
      ENDDO
   
      stop
      end


c -- This fitting routine is written in Fortran 90
	Subroutine FitTail(nTerms,kExp,xPoints,yPoints,Errors,
     >                     nData,MaxPts,rmsError,cc)
        implicit real*8(a-h,o-z)
        integer nData,nTerms,MaxPts
        double precision cc(nTerms),Errors(MaxPts)
        double precision xPoints(MaxPts),yPoints(MaxPts),rmsError
        double precision r0,fitfuncs,dfloat
        double precision DELik,qmax,factor,sum1,sum2,error1,error2
        double precision kExp(3)

        double precision, allocatable :: xMat(:,:),aVec(:)
        double precision, allocatable :: Aux(:,:),SaveMatrix(:,:),SecondVec(:)
        double precision, allocatable :: SolutionVec(:),SaveVec(:)
        integer n,i,j,k,Row,Col,nCols,MasterCol

cc	external fitfuncs

        n = nTerms

        r0=xPoints(nData/2)

        nCols=2*n+1
        allocate(xMat(n,n),aVec(n),Aux(n,nCols),SaveMatrix(n,n),SecondVec(n))
        allocate(SolutionVec(n),SaveVec(n))

C -- now set up the matrix xMat and store it in the left-most part of Aux:

        do kp=1,n
         do k=1,n
         xMat(kp,k)=0.d0
          do i=1,nData
           xMat(kp,k) = xMat(kp,k)+fitfuncs(kExp(kp),xPoints(i),r0)*
     >                  fitfuncs(kExp(k),xPoints(i),r0)
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

 155	format(4(1x,1pd14.7))

        SaveMatrix=xMat

C -- Now solve the necessary linear system by performing 
C ---- row operations on this matrix Aux:
C -- The algorithm I will adopt here is a simple one without "pivoting",
C ---  "Gaussian elimination with backsubstitution" (see Numerical Recipes)
 
C -- Step 1: stabilize by rescaling matrices so that the maximum element of xMat 
C ----  in any given row is unity.

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

C -- Now perform the triangularization of the matrix:

        DO MasterCol=1,n-1
         IF(dabs(Aux(MasterCol,MasterCol)).LT.1.d-20) then
            write(6,*) 'Division by zero problem'
            write(6,*) '-- need to implement pivoting'
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

C -- Now it should be triangular, so we can "back-substitute".
C ---- First get the desired solution vector:

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

C -- Next cycle through the modified column vectors of the unit matrix to
C -- construct the inverse of the original matrix in case it is desired
C -- To save storage, overwrite the original matrix.

C -- Compare the final accuracy of this calculation by testing the ability of
C --- the polynomial to reproduce the data points:

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

