      Program Thermal_Avg
      implicit real*8(a-h,o-z)
      parameter (ngauss=85000)
      dimension x(ngauss),w(ngauss)
      integer NumIntPoints
      double precision K3A(1:200,1:200),energy(1:200)
      double precision K3AInt(1:200,ngauss)
      double precision xa(1:200),ya(1:200),k3,Pi,Ecoef,Lcoef,muK,m40K
      integer ScaL(1:200)
      character*64 A(1:200),B(1:300)

      Pi = dacos(-1.d0)

      muK = 3.1668d-12
      m40K = 7.28497736805d4
      Ecoef = 1.d0/m40K/(65.d0)**2
      Lcoef = 65.d0


c     Total Number of points
      Nadd = 0
      Np = 180+Nadd

      energy = 0.d0
      K3A = 0.d0
      iA = 1
      do iE = 1,Np
        read(250,*)e,k3 
        energy(iE) = e
        K3A(iA,iE) = k3
      enddo
      do iE = 1,Nadd
        energy(iE) = energy(Nadd+1)*1.d-2*(10.d0**(float(iE-1)*.05d0)) 
        K3A(iA,iE) = K3A(iA,Nadd+1)
      enddo  

      energy = Ecoef*energy
      k3a = k3a*Lcoef

      do iE = 1,Np
         write(100,*)energy(iE),K3A(iA,iE)
      enddo   

      xinf = energy(1)
      xsup = energy(Np)
      NumIntPoints = ngauss

!     ** abcissas and weights **
      x = 0.d0; w = 0.d0;
      alpha = 0.d0; beta = 0.d0;
!     call gaussjacobi(x,w,NumIntPoints,alpha,beta)
      do ni = 1,ngauss
      read(111,*)nix,x(ni),w(ni)
      enddo
      x = (xsup-xinf)/2.d0*x+(xsup+xinf)/2.d0

c     ** Interpolating **
      K3AInt = 0.d0
      write(*,*)'Interpolating'   

      iA = 1
      do ni = 1,NumIntPoints
            xx = x(ni)   
            ntry = 5
            nstart = 0
            do ii= 1,Np 
               if (energy(ii)/xx.le.1.d0) then
               nstart = nstart+1 
               endif
            enddo    
            if (nstart.lt.ntry) nstart = 1
               if ((nstart-1+ntry).le.Np) then
	       ncount = 0
	       do nns = nstart,nstart-1+ntry
      	       ncount = ncount+1
	       xa(ncount) = energy(nns)
	       ya(ncount) = K3a(iA,nns)
	       enddo
            else
	       ncount = 0
	       do nns = Np-ntry+1,Np
      	       ncount = ncount+1
	       xa(ncount) = energy(nns)
	       ya(ncount) = K3a(iA,nns)
	       enddo
            endif
            call polint(xa,ya,ntry,xx,y,dy)
            K3AInt(iA,ni) = y
         enddo

!c     Interpolated Result (just to check)
!      do iA = 1,1 !1,Nk
!         do ni = 1,ngauss
!           write(300,*)x(ni),K3AInt(iA,ni)
!         enddo  
!      enddo  


      write(*,*)'Calculating Thermal Average'   
      ncount = 0
      iA = 1
      do iT = Nadd-1,Np-1 
         T = energy(1)*(10.d0**(float(iT)*.05d0)) 
         kb = 1.d0
         xint = 0.d0
         do ni = 1,NumIntPoints
         xx = x(ni)   
         aux = (2.d0/dsqrt(Pi)/((kb*T)**1.5d0))*K3AInt(iA,ni)*(xx**0.5d0)*
     .          dexp(-xx/(kb*T)) 
         xint = xint+w(ni)*aux
         enddo
         xint = (xsup-xinf)/2.d0*xint
!        write(*,*)iT,dabs(xint-K3A(iA,1))/K3A(iA,1)
         if (dabs(xint-K3A(iA,1))/K3A(iA,1).lt.1.e-3.or.ncount.ne.0) then
            write(255,*)T,xint
            ncount = ncount+1
         else   
            write(255,*)T,K3A(iA,1)
         endif   


      enddo            

      stop
      end

c     ********** rotina de integracao gauss **********

      subroutine gauss(ff,xin,ngauss,dfake)
      implicit real*8 (a-h,o-z)
      dimension ff(-2:2,200),t(0:2),xk(0:2)
      c = 1.d0/dsqrt(0.175d0)
      d = dfake
      sas = 65.d0/c
      xk(0) = 128.d0
      xk(1) = 80.5d0+sas
      xk(2) = 80.5d0-sas
      do 1604 li = 0,2
      t(li) = 0.d0
 1604 continue
      do 1605 ki = -2,2
      do 1605 kj = 1,ngauss
      l = iabs(ki)
      t(l) = t(l)+ff(ki,kj)
 191  format(1x,i3,e45.32)
 1605 continue
      e = 0.d0
      do 333 kk = 0,2
      e = e+xk(kk)*t(kk)
 1331 format(1x,i3,3e25.17)
 333  continue
      e = e*d/450.d0
      xin = e 
      return
      end


      SUBROUTINE POLCOF(XA,YA,N,COF)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NPMAX=100)
      DIMENSION COF(N),XA(N),YA(N)
      DIMENSION X(NPMAX),Y(NPMAX)
C     USES POLINT
      DO 11 J=1,N
        X(J)=XA(J)
        Y(J)=YA(J)
11    CONTINUE
      DO 14 J=1,N
        CALL POLINT(X,Y,N+1-J,0.D0,COF(J),DY)
        XMIN=1.E38
        K=0
        DO 12 I=1,N+1-J
          IF (DABS(X(I)).LT.XMIN)THEN
            XMIN=DABS(X(I))
            K=I
          ENDIF
          IF(X(I).NE.0.D0) Y(I)=(Y(I)-COF(J))/X(I)
12      CONTINUE
        DO 13 I=K+1,N+1-J
          Y(I-1)=Y(I)
          X(I-1)=X(I)
13      CONTINUE
14    CONTINUE
      RETURN
      END

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NPMAX=100)
      DIMENSION XA(N),YA(N)
      DIMENSION C(NPMAX),D(NPMAX)
      NS=1
      DIF=DABS(X-XA(1))
      DO 11 I=1,N
        DIFT=DABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF (DEN.EQ.0.D0) PAUSE 'FAILURE IN POLINT'
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

