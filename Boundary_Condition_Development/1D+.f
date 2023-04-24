c     This program uses the newer "CalcBasisFuncsBP" routine to calculate the 
c     bspline basis functions capable of arbitrary Bethe-Peierls boundary conditions
      program ONED

      implicit none
      integer LegPoints,xNumPoints,LeadDim,MatrixDim,HalfBandWidth,NLeft,NRight,Left,Right,Order,xDim
      integer iX,k,NumSteps
      integer NumStates,NumStatesAsk
      double precision Shift,testval
      double precision, allocatable :: Psi(:,:),Energies(:),LUFac(:,:)

      integer info,j,RSteps,Uin,Pin,Qin,Out
      integer, allocatable :: xBounds(:)
      double precision xMin,xMax,ans,a,b,m1,m2,mu,Dcoef,Rscale,xLast,firstX,lastX,xInitial,Stepx,phirecon,rpart,ipart,aLeft,aRight,x
      double precision, allocatable :: xLeg(:),wLeg(:),xPoints(:),Alphar(:),Alphai(:),Beta(:),Vl(:,:),Vr(:,:),Work(:),eEnergies(:)
      double precision, allocatable :: u(:,:,:),ux(:,:,:),uxx(:,:,:),S(:,:),H(:,:),RecPsi(:),ExactEnergies(:)
      double precision, allocatable :: R(:),UE(:,:),P(:,:,:),Q(:,:,:),PotKnots(:)
      double precision, allocatable :: Ubcoef(:,:),Pbcoef(:,:,:),Qbcoef(:,:,:)
      character*64 LegendreFile
      integer*1 n
      double precision Pi


      Pi = acos(-1d0)

c     read in info
c     States
      read(5,*) ! 5 says to read in from the command line, * says take any type of input, # says look for line called that and then you can format it
      read(5,*) NumStates

c     BC
      read(5,*)
      read(5,*)
      read(5,*) NLeft,NRight

c     Gauss-Legendre info
      read(5,*)
      read(5,*)
      read(5,1002) LegendreFile
      read(5,*)
      read(5,*)
      read(5,*) LegPoints

c     Order & Shift
      read(5,*)
      read(5,*)
      read(5,*) Order,Shift

c     Masses
      read(5,*)
      read(5,*)
      read(5,*) m1,m2
      mu = (m1*m2)/(m1+m2)

c     Potential Coefficients
      read(5,*)
      read(5,*)
      read(5,*) Dcoef,Rscale,aLeft,aRight !aLeft = -1/a_l; similar for aRight

c     Grid Information
      read(5,*)
      read(5,*)
      read(5,*) xNumPoints,xLast,RSteps

c     Input/Output files
      read(5,*)
      read(5,*)
      read(5,*) Uin,Pin,Qin,Out

      xMin = 0.0d0
      xMax = xLast
      Left = NLeft              ! boundary condition at 0
      Right = NRight            ! boundary condition at xLast

      write(6,*) 'mass 1'
      write(6,*) m1
      write(6,*) 'mass 2'
      write(6,*) m2
      write(6,*) 'mu'
      write(6,*) mu

      allocate(xLeg(LegPoints),wLeg(LegPoints))
      call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

      xDim = xNumPoints+Order-3
      if(Left.eq.2)xDim = xDim+1
      if(Right.eq.2)xDim = xDim+1

      MatrixDim = xDim*xDim
      HalfBandWidth = xDim*Order+Order
      LeadDim = xDim+2*Order

      allocate(xPoints(xNumPoints),u(LegPoints,xNumPoints,xDim),ux(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim))
      allocate(xBounds(xNumPoints+2*Order))
      allocate(Psi(MatrixDim,NumStates),Energies(NumStates),LUFac(LeadDim,xDim))
      allocate(Alphar(xDim),Alphai(xDim),Beta(xDim),Vl(xDim,xDim),Vr(xDim,xDim),Work(25000))
c      allocate(S(xDim+2*Order,xDim),H(xDim+2*Order,xDim),eEnergies(xDim))
      allocate(S(xDim,xDim),H(xDim,xDim),eEnergies(xDim),ExactEnergies(xDim)) !dggev
      allocate(R(RSteps),UE(NumStates,Rsteps),P(NumStates,NumStates,Rsteps),Q(NumStates,NumStates,RSteps))
      allocate(PotKnots(RSteps+Order),Ubcoef(NumStates,RSteps),Pbcoef(NumStates,NumStates,RSteps))
      allocate(Qbcoef(NumStates,NumStates,RSteps))

c      call readCouplings(Uin,Pin,Qin,RSteps,R,UE,P,Q,NumStates)
c      call setupInterp(R,Order,PotKnots,RSteps,P,Q,UE,Ubcoef,Pbcoef,Qbcoef,NumStates)
c      call setupInterp(xPoints,Order,PotKnots,RSteps,P,Q,UE,Ubcoef,Pbcoef,Qbcoef,NumStates)

      call GridMaker(xMin,xMax,xNumPoints,xPoints,RSteps,R)
c      call CheckBasisPhi(xMin,xMax,Left,Right,aLeft,aRight,xDim,xNumPoints,xPoints,0,Order,300)
c      call CheckBasisPhi(xMin,xMax,Left,Right,aLeft,aRight,xDim,xNumPoints,xPoints,1,Order,301)
c      STOP
      call CalcBasisFuncsBP(Left,Right,aLeft,aRight,Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,0,u)
      call CalcBasisFuncsBP(Left,Right,aLeft,aRight,Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,1,ux)
      call CalcBasisFuncsBP(Left,Right,aLeft,aRight,Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,2,uxx)
      call CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,xDim,xNumPoints,u,xBounds,HalfBandWidth,MatrixDim,S)
      call CalcHamiltonian(Order,Dcoef,Rscale,mu,xPoints,LegPoints,xLeg,wLeg,xNumPoints,xDim,xBounds,
     >     u,ux,uxx,MatrixDim,HalfBandWidth,H,PotKnots,RSteps,NumStates,Ubcoef,Qbcoef,R,Left,Right,aLeft,aRight)

      call dggev('N','V',xDim,H,xDim,S,xDim,Alphar,Alphai,Beta,Vl,xDim,Vr,xDim,Work,25000,info)
c      call MyLargeDsband(NumStates,Shift,Energies,Psi,MatrixDim,H,S,LeadDim,HalfBandWidth,LUFac,Order,xDim)

      eEnergies = -1*(Alphar/Beta) !*27.2113959819d0 !for H
      call EIGSRT(eEnergies,Vr,xDim,xDim)
      eEnergies = -1*eEnergies

      write(Out,*) 'Calculated eEnergies'
      do k = 1,xDim
         write(Out,*) eEnergies(k)
      enddo

c      write(910,*) 'Predicted eEnergies'
c      do k = 1,xDim
c         n = n+1
c         write(910,*) (n**2)*(Pi**2)/(2*mu)
c      enddo
c      n = 0
c      write(915,*) 'eEnergy Differences'
c      do k = 1,xDim
c         n = n+1
c         write(915,*) (eEnergies(k) - (n**2)*(Pi**2)/(2*mu))
c      enddo

      NumSteps=1000
      Stepx=(xMax-xMin)/(NumSteps-1)
      do iX = 1, NumSteps
         x=xMin + (iX - 1)*Stepx
         write(200,11) x,phirecon(x,1,Vr,Left,Right,aLeft,aRight,xDim,xDim,xNumPoints,xPoints,Order,0), !wavefunction
     >        phirecon(x,1,Vr,Left,Right,aLeft,aRight,xDim,xDim,xNumPoints,xPoints,Order,1), !wavefunction derivative
     >        phirecon(x,1,Vr,Left,Right,aLeft,aRight,xDim,xDim,xNumPoints,xPoints,Order,1) / 
     >        phirecon(x,1,Vr,Left,Right,aLeft,aRight,xDim,xDim,xNumPoints,xPoints,Order,0) !wavefunction log derivative
         
       enddo

      do k = 1,NumStates
         write(Out,*) Energies(k) !*27.2113959819d0 !for H
      enddo

      deallocate(xLeg,wLeg,xPoints,xBounds,S,H,u,ux,uxx,ExactEnergies)
      deallocate(Psi,Energies,LUfac)
      deallocate(Alphar,Alphai,Beta,Vl,Vr,Work)

 11   format(f16.12,1p,100e20.12)
 1002 format(a64)

      end program

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine GridMaker(xMin,xMax,xNumPoints,xPoints,RSteps,R)

      implicit none
      integer xNumPoints,i,k
      double precision xMin,xMax,xdelta1,xdelta2,xdelta3,Pi
      double precision xPoints(xNumPoints)

      integer RSteps
      double precision R(RSteps)

      !Pi = 3.14159265358979323846d0

c     even grid
      xdelta1 = (xMax-xMin)/(xNumPoints-1)
      do i = 0,xNumPoints-1
         xPoints(i+1) = (xMin + i*xdelta1)
      enddo

c      xPoints = R

cprints the grid
c      do i = 1,xNumPoints
c         write(700,*) i,xPoints(i)
c      enddo

      return
      end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine TestGridMaker(xMin,xMax,xNumPoints,xPoints)

      implicit none
      integer xNumPoints,i
      double precision xMin,xMax,xdelta,xPoints(xNumPoints)

      xdelta = (xMax-xMin)/(xNumPoints-1)
      do i = 0,xNumPoints-1
         xPoints(i+1) = (xMin + i*xdelta)
      enddo

      write(300,*) xPoints

      end subroutine TestGridMaker


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,xDim,xNumPoints,u,xBounds,HalfBandWidth,MatrixDim,S)

      implicit none
      integer Order,LegPoints,XDim,xNumPoints,HalfBandWidth,MatrixDim
      integer xBounds(xNumPoints+2*Order)
      double precision xPoints(*),xLeg(LegPoints),wLeg(LegPoints)
      double precision u(LegPoints,xNumPoints,xDim)
c      double precision S(xDim+2*Order,xDim)
      double precision S(xDim,xDim) !dggev

      integer ix,ixp,kx,lx,Row,NewRow,Col,k,i1,i1p,ix2,i2,ix2p,i2p,k0,kp0,j
      integer, allocatable :: kxMin(:,:),kxMax(:,:)
      double precision ax,bx,tempS,xScaledZero
      double precision, allocatable :: tempSM(:,:),xScale(:)

      allocate(kxMin(xDim,xDim),kxMax(xDim,xDim),tempSM(xDim,xDim),xScale(xNumPoints))

      S = 0.0d0

      do kx = 1,xNumPoints-1
         ax = xPoints(kx)
         bx = xPoints(kx+1)
         xScale(kx) = 0.5d0*(bx-ax)
      enddo

      do ix = 1,xDim
         do ixp = 1,xDim
            kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
            kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
         enddo
      enddo

      tempSM = 0.0d0
      do ix = 1,xDim
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            do kx = kxMin(ixp,ix),kxMax(ixp,ix)
               tempS = 0.0d0
               do lx = 1,LegPoints
                  tempS = tempS + wLeg(lx)*u(lx,kx,ix)*u(lx,kx,ixp)
               enddo
c               tempSM(ixp,ix) = tempSM(ixp,ix) + xScale(kx)*tempS
               S(ixp,ix) = S(ixp,ix) + xScale(kx)*tempS !dggev
            enddo
         enddo
      enddo

c      do ix = 1,xDim
c         do ixp = max(1,ix-Order),min(xDim,ix+Order)
c            Row = 2*Order+1+ixp-ix
c            S(Row,ix) = tempSM(ixp,ix)
c         enddo
c      enddo

c      do ix = 1,xDim
c         do ixp = 1,xDim
c            write(20,11,advance='NO') S(ix,ixp)
c         enddo
c         write(20,*)
c      enddo

      deallocate(kxMin,kxMax,tempSM,xScale)
 11   format(f16.12,1p,100e20.12)
      return
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine CalcHamiltonian(Order,Dcoef,Rscale,mu,xPoints,LegPoints,xLeg,wLeg,xNumPoints,xDim,
     >     xBounds,u,ux,uxx,MatrixDim,HalfBandWidth,H,PotKnots,RSteps,NumStates,Ubcoef,Qbcoef,R,Left,Right,kLeft,kRight)
      
      implicit none
      integer Order,LegPoints,xDim,xNumPoints,k,kx,ix,lx,ixp,MatrixDim,HalfBandWidth
      integer xBounds(xNumPoints+2*Order),i1,i1p,ix2,ix2p,i2,i2p,Row,Col,NewRow,kx0,kxp0,kp
      integer, allocatable :: kxMin(:,:),kxMax(:,:)
      double precision Pi,hbar,xScaledZero,tempPot,Vpot,Dcoef,Rscale,mu,tempH
      double precision xPoints(xNumPoints),xLeg(LegPoints),wLeg(LegPoints),ax,bx
      double precision, allocatable :: xScale(:)
      double precision u(LegPoints,xNumPoints,xDim),ux(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim)
      double precision x(LegPoints,xNumPoints-1),tempHM(xDim,xDim)
c      double precision H(xDim+2*Order,xDim)
      double precision H(xDim,xDim) !dggev

      integer NumStates,RSteps,Left,Right
      double precision PotKnots(NumStates+Order),Ubcoef(NumStates,RSteps),Qbcoef(NumStates,NumStates,RSteps),R(RSteps)
      double precision kLeft,kRight

      double precision, external :: BasisPhi

      !Pi = 3.14159265358979323846d0
      H = 0.0d0
      allocate(kxMin(xDim,xDim),kxMax(xDim,xDim),xScale(xNumPoints))

      do kx = 1,xNumPoints-1 !# points to calculate bspline OVERLAP integral on; each overlapping segment has the 10 legendre point/guassian quadriture calculated on it 
         ax = xPoints(kx)
         bx = xPoints(kx+1)
         xScale(kx) = 0.5d0*(bx-ax) !b-a/2
         xScaledZero = 0.5d0*(bx+ax) !middle point in new,rescaled gaussian quadriture point
         do lx = 1,LegPoints !calculates all 10 points within given range
            x(lx,kx) = xScale(kx)*xLeg(lx) + xScaledZero !stretching out the x and offsetting it
         enddo
      enddo

      do ix = 1,xDim !determines section of overlap between b-splines -DO NOT MODIFY
         do ixp = 1,xDim
            kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp)) !min kx and max kx needed for a matrix elements containing a pair
            kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
         enddo
      enddo

!u is basis
!ux is derivative
!uxx = ddx^2

      tempHM = 0.0d0
      do ix = 1,xDim
         kx0 = xBounds(ix)-1 !bra index: leftmost section for bra basis element
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            kxp0 = xBounds(ixp)-1 !leftmost section of the ket
            do kx = kxMin(ixp,ix),kxMax(ixp,ix)
               tempH = 0.0d0
               do lx = 1,LegPoints !actually doing the Gaussian quadriture in each section
                  tempH = tempH + wLeg(lx)*(1.d0/(2.d0*mu)*ux(lx,kx,ix)*ux(lx,kx,ixp))
!     >                 +0.5d0*mu*x(lx,kx)**2  !+ (-10d0)/dcosh(x(lx,kx)/1d0)**2  !!you can see the **2 part, that is x
!     >                 *u(lx,kx,ixp)*u(lx,kx,ix))
c                  tempH = tempH + wLeg(lx)*(1.d0/(2.d0*mu)*ux(lx,kx,ix)*ux(lx,kx,ixp)
c     >                 +Vpot(x(lx,kx),Rscale,Order,PotKnots,RSteps,NumStates,Ubcoef(1,:))
c     >                 *u(lx,kx,ixp)*u(lx,kx,ix))
c                  tempH = tempH + wLeg(lx)*u(lx,kx,ix)*((-1.d0/(2.d0*mu))*uxx(lx,kx,ixp)
c     >                 +(Vpot(x(lx,kx),Rscale,Order,PotKnots,RSteps,NumStates,Ubcoef(1,:))
c     >                 )*u(lx,kx,ixp))
               enddo
c               tempHM(ixp,ix) = tempHM(ixp,ix) + xScale(kx)*tempH
               H(ixp,ix) = H(ixp,ix) + xScale(kx)*tempH !dggev !!bras and kets got flipped, but doesn't matter because it's symmetriccal
            enddo
            if ((ix.eq.1).and.(ixp.eq.1).and.(Left.eq.3)) then !if left index is 1, right index is 1, and the log derivative is 0
               H(ixp,ix) = H(ixp,ix) + (0.5d0/mu)*kLeft  ! this adds a pieces from the surface term kd(i,1)*kd(j,1)*ui(a)*uj(a)*(uj'(a)/uj(a))
            endif
            if ((ix.eq.xDim).and.(ixp.eq.xDim).and.(Right.eq.3)) then
               H(ixp,ix) = H(ixp,ix) - (0.5d0/mu)*kRight
            endif
c            H(ixp,ix) = H(ixp,ix) + 1.d0/(2.d0*mu)*BasisPhi(0.d0,Left,Right,kLeft,kRight,Order,xDim,xPoints,xNumPoints,0,ix)*
c     >                  BasisPhi(0.d0,Left,Right,kLeft,kRight,Order,xDim,xPoints,xNumPoints,0,ixp)*kLeft
c            H(ixp,ix) = H(ixp,ix) - 1.d0/(2.d0*mu)*BasisPhi(1.d2,Left,Right,kLeft,kRight,Order,xDim,xPoints,xNumPoints,0,ix)*
c     >                  BasisPhi(1.d2,Left,Right,kLeft,kRight,Order,xDim,xPoints,xNumPoints,0,ixp)*kRight
         enddo
      enddo

c      do ix = 1,xDim
c         do ixp = max(1,ix-Order),min(xDim,ix+Order)
c            Row = 2*Order+1+ixp-ix
c            H(Row,ix) = tempHM(ixp,ix)
c         enddo
c      enddo

c      do ix = 1,xDim
c         do ixp = 1,xDim
c            write(21,11,advance='NO') H(ix,ixp)
c         enddo
c         write(21,*)
c      enddo

      deallocate(kxMin,kxMax,xScale)
 11   format(f16.12,1p,100e20.12)
      return
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function Vpot(r,Rscale,Order,KnotSeq,InterpPoints,NumChan,Coef)

      implicit none
      integer Order,InterpPoints,NumChan
      double precision r,Rscale,x,Pi,KnotSeq(InterpPoints+Order),Coef(InterpPoints),dbsval

      !Pi = 3.1415926535897932385d0

c      x = r/Rscale
c      Vpot = -1.d0/x

c      Vpot = dbsval(r,Order,KnotSeq,InterpPoints,Coef)
      Vpot = 0.0d0


      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine GLint(a,b,LegPoints,xLeg,wLeg,ans)

      implicit none
      integer LegPoints,i
      double precision a,b,ans,dx,xsquared,xScaledZero
      double precision xLeg(LegPoints),wLeg(LegPoints),xScale

      xScaledZero = 0.5d0*(b+a)
      xScale = 0.5d0*(b-a)
      do i = 1,LegPoints
         dx = xScale*xLeg(i)
         ans = ans + wLeg(i)*xsquared(xScaledZero+dx)
      enddo

      print*,ans

      return
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function xsquared(x)
      double precision x

      xsquared = x**2

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine CalcBasisFuncsBP(Left,Right,kLeft,kRight,Order,xPoints,LegPoints,
     >     xLeg,MatrixDim,xBounds,xNumPoints,Deriv,u)

      integer Left,Right,Order,LegPoints,MatrixDim,xBounds(*),
     >     xNumPoints,Deriv
      double precision xPoints(*),xLeg(*)
      double precision u(LegPoints,xNumPoints,MatrixDim)

      integer i,k,l,Count
      integer, allocatable :: t(:)
      double precision MyBSpline
      double precision x,ax,bx,xIntScale,xScaledZero
      double precision kLeft,kRight,constLeft,constRight,lc1n,lc2n,rc1n,rc2n

      allocate(t(xNumPoints+2*Order))

      do i = 1,Order
       t(i) = 1
      enddo
      do i = 1,xNumPoints
       t(i+Order) = i
      enddo
      do i = 1,Order
       t(i+Order+xNumPoints) = xNumPoints
      enddo

      select case (Left)
      case (0:1)
         select case (Right)
         case (0:1)
            do i = 2,xNumPoints+2*Order-1
               xBounds(i-1) = t(i)
            enddo
         case (2)
            do i = 2,xNumPoints+2*Order
               xBounds(i-1) = t(i)
            enddo
         case(3)
            do i = 2,xNumPoints+2*Order-1
               xBounds(i-1) = t(i)
            enddo
         end select
      case (2)
         select case (Right)
         case (0:1)
            do i = 1,xNumPoints+2*Order-1
               xBounds(i) = t(i)
            enddo
         case (2)
            do i = 1,xNumPoints+2*Order
               xBounds(i) = t(i)
            enddo
         case (3)
            do i = 1,xNumPoints+2*Order-1
               xBounds(i) = t(i)
            enddo
         end select
      case (3)
         select case (Right)
         case (0:1)
            do i = 2,xNumPoints+2*Order-1
               xBounds(i-1) = t(i)
            enddo
         case (2)
            do i = 2,xNumPoints+2*Order
               xBounds(i-1) = t(i)
            enddo
         case(3)
            do i = 2,xNumPoints+2*Order-1
               xBounds(i-1) = t(i)
            enddo
         end select
      end select
      
      deallocate(t)

      Count = 1
      select case (Left)
      case (0)
         do k = 1,xNumPoints-1
            ax = xPoints(k)
            bx = xPoints(k+1)
            xScale = 0.5d0*(bx-ax)
            xScaledZero = 0.5d0*(bx+ax)
            do l = 1,LegPoints
               x = xScale*xLeg(l)+xScaledZero
               u(l,k,Count) = MYBSpline(Order,Deriv,xNumPoints,xPoints,2,x)
            enddo
         enddo
         Count = Count + 1
      case (1)
         do k = 1,xNumPoints-1
            ax = xPoints(k)
            bx = xPoints(k+1)
            xScale = 0.5d0*(bx-ax)
            xScaledZero = 0.5d0*(bx+ax)
            do l = 1,LegPoints
               x = xScale*xLeg(l)+xScaledZero
               u(l,k,Count) = MYBSpline(Order,Deriv,xNumPoints,xPoints,1,x)+
     >              MYBSpline(Order,Deriv,xNumPoints,xPoints,2,x)
            enddo
         enddo
         Count = Count + 1
      case(2)
         do i = 1,2
            do k = 1,xNumPoints-1
               ax = xPoints(k)
               bx = xPoints(k+1)
               xScale = 0.5d0*(bx-ax)
               xScaledZero = 0.5d0*(bx+ax)
               do l = 1,LegPoints
                  x = xScale*xLeg(l)+xScaledZero
                  u(l,k,Count) = MYBSpline(Order,Deriv,xNumPoints,xPoints,i,x)
               enddo
            enddo
            Count = Count + 1
         enddo
      case(3)
c         constLeft = MYBSpline(Order,1,xNumPoints,xPoints,2,xPoints(1)) - MYBSpline(Order,0,xNumPoints,xPoints,2,xPoints(1))*kLeft /!!!!change this block and the constRight block to match the conventions used in notes
c     >        (MYBSpline(Order,0,xNumPoints,xPoints,1,xPoints(1))*kLeft - 
c     >        MYBSpline(Order,1,xNumPoints,xPoints,1,xPoints(1)))
         constLeft = (MYBSpline(Order,1,xNumPoints,xPoints,2,xPoints(1))
     >       + MYBSpline(Order,0,xNumPoints,xPoints,2,xPoints(1))*kLeft)/
     >        (MYBSpline(Order,0,xNumPoints,xPoints,1,xPoints(1))*kLeft + 
     >        MYBSpline(Order,1,xNumPoints,xPoints,1,xPoints(1)))
         do k = 1,xNumPoints-1
            ax = xPoints(k)
            bx = xPoints(k+1)
            xScale = 0.5d0*(bx-ax)
            xScaledZero = 0.5d0*(bx+ax)
            do l = 1,LegPoints
               x = xScale*xLeg(l)+xScaledZero
               u(l,k,Count) = MYBSpline(Order,Deriv,xNumPoints,xPoints,1,x)+
     >                        MYBSpline(Order,Deriv,xNumPoints,xPoints,2,x)/constLeft
            enddo
         enddo
         Count = Count + 1
      end select

      do i = 3,xNumPoints+Order-3
         do k = 1,xNumPoints-1
            ax = xPoints(k)
            bx = xPoints(k+1)
            xIntScale = 0.5d0*(bx-ax)
            xScaledZero = 0.5d0*(bx+ax)
            do l = 1,LegPoints
               x = xIntScale*xLeg(l)+xScaledZero
               u(l,k,Count) = MYBSpline(Order,Deriv,xNumPoints,xPoints,i,x)
            enddo
         enddo
         Count = Count + 1
      enddo

      select case (Right)
      case (0)
         do k = 1,xNumPoints-1
            ax = xPoints(k)
            bx = xPoints(k+1)
            xScale = 0.5d0*(bx-ax)
            xScaledZero = 0.5d0*(bx+ax)
            do l = 1,LegPoints
               x = xScale*xLeg(l)+xScaledZero
               u(l,k,Count) = MYBSpline(Order,Deriv,xNumPoints,xPoints,
     >              xNumPoints+Order-2,x)
            enddo
         enddo
      case (1)
         do k = 1,xNumPoints-1
            ax = xPoints(k)
            bx = xPoints(k+1)
            xScale = 0.5d0*(bx-ax)
            xScaledZero = 0.5d0*(bx+ax)
            do l = 1,LegPoints
               x = xScale*xLeg(l)+xScaledZero
               u(l,k,Count) = MYBSpline(Order,Deriv,xNumPoints,xPoints,
     >              xNumPoints+Order-2,x)+
     >              MYBSpline(Order,Deriv,xNumPoints,xPoints,
     >              xNumPoints+Order-1,x)
            enddo
         enddo
      case(2)
         do i = xNumPoints+Order-2,xNumPoints+Order-1
            do k = 1,xNumPoints-1
               ax = xPoints(k)
               bx = xPoints(k+1)
               xScale = 0.5d0*(bx-ax)
               xScaledZero = 0.5d0*(bx+ax)
               do l = 1,LegPoints
                  x = xScale*xLeg(l)+xScaledZero
                  u(l,k,Count) = MYBSpline(Order,Deriv,xNumPoints,xPoints,i,x)
               enddo
            enddo
            Count = Count + 1
         enddo
      case(3)
         constRight = (MYBSpline(Order,1,xNumPoints,xPoints,xNumPoints+Order-2,xPoints(xNumPoints))
     >        -MYBSpline(Order,1,xNumPoints,xPoints,xNumPoints+Order-2,xPoints(xNumPoints))*kRight) / ( 
     >        MYBSpline(Order,1,xNumPoints,xPoints,xNumPoints+Order-1,xPoints(xNumPoints))*kRight - 
     >        MYBSpline(Order,0,xNumPoints,xPoints,xNumPoints+Order-1,xPoints(xNumPoints)))
         do k = 1,xNumPoints-1
            ax = xPoints(k)
            bx = xPoints(k+1)
            xScale = 0.5d0*(bx-ax)
            xScaledZero = 0.5d0*(bx+ax)
            do l = 1,LegPoints
               x = xScale*xLeg(l)+xScaledZero
               u(l,k,Count) = MYBSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-2,x)/constRight+
     >              MYBSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-1,x)
            enddo
         enddo
      end select

      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function MyBSpline(Order,Deriv,xNumPoints,
     >     xPoints,n,x)
      
      integer Order,Deriv,xNumPoints,n
      double precision xPoints(*),x
      
      integer i
      double precision bvalue
      double precision, allocatable :: t(:),b(:)
      
      allocate(t(xNumPoints+2*Order))
      allocate(b(xNumPoints+Order))
      
      do i = 1,Order
         t(i) = xPoints(1)
      enddo
      do i = 1,xNumPoints
         t(i+Order) = xPoints(i)
      enddo
      do i = 1,Order
         t(i+Order+xNumPoints) = xPoints(xNumPoints)
      enddo

      do i = 1,xNumPoints+Order
         b(i) = 0.0d0
      enddo
      b(n) = 1.0d0
      
      MyBSpline = bvalue(t,b,n,Order+1,x,Deriv)
c      print*, MyBSpline
      deallocate(t)
      deallocate(b)
      
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function bvalue ( t, bcoef, n, k, x, jderiv )
c  from  * a practical guide to splines *  by c. de boor    
calls  interv
c
calculates value at  x  of  jderiv-th derivative of spline from b-repr.
c  the spline is taken to be continuous from the right, EXCEPT at the
c  rightmost knot, where it is taken to be continuous from the left.
c
c******  i n p u t ******
c  t, bcoef, n, k......forms the b-representation of the spline  f  to
c        be evaluated. specifically,
c  t.....knot sequence, of length  n+k, assumed nondecreasing.
c  bcoef.....b-coefficient sequence, of length  n .
c  n.....length of  bcoef  and dimension of spline(k,t),
c        a s s u m e d  positive .
c  k.....order of the spline .
c
c  w a r n i n g . . .   the restriction  k .le. kmax (=20)  is imposed
c        arbitrarily by the dimension statement for  aj, dl, dr  below,
c        but is  n o w h e r e  c h e c k e d  for.
c
c  x.....the point at which to evaluate .
c  jderiv.....integer giving the order of the derivative to be evaluated
c        a s s u m e d  to be zero or positive.
c
c******  o u t p u t  ******
c  bvalue.....the value of the (jderiv)-th derivative of  f  at  x .
c
c******  m e t h o d  ******
c     The nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo-
c  cated with the aid of  interv . The  k  b-coeffs of  f  relevant for
c  this interval are then obtained from  bcoef (or taken to be zero if
c  not explicitly available) and are then differenced  jderiv  times to
c  obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
c  Precisely, with  j = jderiv, we have from x.(12) of the text that
c
c     (d**j)f  =  sum ( bcoef(.,j)*b(.,k-j,t) )
c
c  where
c                   / bcoef(.),                     ,  j .eq. 0
c                   /
c    bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1)
c                   / ----------------------------- ,  j .gt. 0
c                   /    (t(.+k-j) - t(.))/(k-j)
c
c     Then, we use repeatedly the fact that
c
c    sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
c  with
c                 (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1)
c    a(.,x)  =    ---------------------------------------
c                 (x - t(.))      + (t(.+m-1) - x)
c
c  to write  (d**j)f(x)  eventually as a linear combination of b-splines
c  of order  1 , and the coefficient for  b(i,1,t)(x)  must then be the
c  desired number  (d**j)f(x). (see x.(17)-(19) of text).
c
      integer jderiv,k,n,   i,ilo,imk,j,jc,jcmin,jcmax,jj,kmax,kmj,km1
     *                     ,mflag,nmi,jdrvp1
      parameter (kmax = 20)
C     double precision bcoef(n),t(1),x,   aj(20),dl(20),dr(20),fkmj
      double precision bcoef(*),t(*),x,   aj(kmax),dl(kmax),dr(kmax),fkmj
c      dimension t(n+k)
c  former fortran standard made it impossible to specify the length of  t
c  precisely without the introduction of otherwise superfluous addition-
c  al arguments.
      bvalue = 0.0d0
      if (jderiv .ge. k)                go to 99
c
c  *** Find  i   s.t.   1 .le. i .lt. n+k   and   t(i) .lt. t(i+1)   and
c      t(i) .le. x .lt. t(i+1) . If no such i can be found,  x  lies
c      outside the support of  the spline  f , hence  bvalue = 0.
c      (The asymmetry in this choice of  i  makes  f  rightcontinuous, except
c      at  t(n+k) where it is leftcontinuous.)
      call interv ( t, n+k, x, i, mflag )
      if (mflag .ne. 0)                 go to 99
c  *** if k = 1 (and jderiv = 0), bvalue = bcoef(i).
      km1 = k - 1
      if (km1 .gt. 0)                   go to 1
      bvalue = bcoef(i)
                                        go to 99
c
c  *** store the k b-spline coefficients relevant for the knot interval
c     (t(i),t(i+1)) in aj(1),...,aj(k) and compute dl(j) = x - t(i+1-j),
c     dr(j) = t(i+j) - x, j=1,...,k-1 . set any of the aj not obtainable
c     from input to zero. set any t.s not obtainable equal to t(1) or
c     to t(n+k) appropriately.
    1 jcmin = 1
      imk = i - k
      if (imk .ge. 0)                   go to 8
      jcmin = 1 - imk
      do j=1,i
         dl(j) = x - t(i+1-j)
      enddo
      do j=i,km1
         aj(k-j) = 0.0d0
         dl(j) = dl(i)
      enddo
                                        go to 10
    8 do j=1,km1
         dl(j) = x - t(i+1-j)
      enddo
c
   10 jcmax = k
      nmi = n - i
      if (nmi .ge. 0)                   go to 18
      jcmax = k + nmi
      do j=1,jcmax
         dr(j) = t(i+j) - x
      enddo
      do j=jcmax,km1
         aj(j+1) = 0.0d0
         dr(j) = dr(jcmax)
      enddo
                                        go to 20
   18 do j=1,km1
         dr(j) = t(i+j) - x
      enddo
c
   20 do jc=jcmin,jcmax
         aj(jc) = bcoef(imk + jc)
      enddo
c
c               *** difference the coefficients  jderiv  times.
      if (jderiv .eq. 0)                go to 30
      do j=1,jderiv
         kmj = k-j
         fkmj = float(kmj)
         ilo = kmj
         do jj=1,kmj
            aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
            ilo = ilo - 1
         enddo
      enddo
c
c  *** compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
c     given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv).
   30 if (jderiv .eq. km1)              go to 39
      jdrvp1 = jderiv + 1     
      do j=jdrvp1,km1
         kmj = k-j
         ilo = kmj
         do jj=1,kmj
            aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
            ilo = ilo - 1
         enddo
      enddo
   39 bvalue = aj(1)
c
   99                                   return
      end
      subroutine interv ( xt, lxt, x, left, mflag )
c  from  * a practical guide to splines *  by C. de Boor    
computes  left = max( i :  xt(i) .lt. xt(lxt) .and.  xt(i) .le. x )  .
c
c******  i n p u t  ******
c  xt.....a double precision sequence, of length  lxt , assumed to be nondecreasing
c  lxt.....number of terms in the sequence  xt .
c  x.....the point whose location with respect to the sequence  xt  is
c        to be determined.
c
c******  o u t p u t  ******
c  left, mflag.....both integers, whose value is
c
c   1     -1      if               x .lt.  xt(1)
c   i      0      if   xt(i)  .le. x .lt. xt(i+1)
c   i      0      if   xt(i)  .lt. x .eq. xt(i+1) .eq. xt(lxt)
c   i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
c
c        In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
c        indicates that  x  lies outside the CLOSED interval
c        xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
c        intervals is due to the decision to make all pp functions cont-
c        inuous from the right, but, by returning  mflag = 0  even if
C        x = xt(lxt), there is the option of having the computed pp function
c        continuous from the left at  xt(lxt) .
c
c******  m e t h o d  ******
c  The program is designed to be efficient in the common situation that
c  it is called repeatedly, with  x  taken from an increasing or decrea-
c  sing sequence. This will happen, e.g., when a pp function is to be
c  graphed. The first guess for  left  is therefore taken to be the val-
c  ue returned at the previous call and stored in the  l o c a l  varia-
c  ble  ilo . A first check ascertains that  ilo .lt. lxt (this is nec-
c  essary since the present call may have nothing to do with the previ-
c  ous call). Then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
c  ilo  and are done after just three comparisons.
c     Otherwise, we repeatedly double the difference  istep = ihi - ilo
c  while also moving  ilo  and  ihi  in the direction of  x , until
c                      xt(ilo) .le. x .lt. xt(ihi) ,
c  after which we use bisection to get, in addition, ilo+1 = ihi .
c  left = ilo  is then returned.
c
      integer left,lxt,mflag,   ihi,ilo,istep,middle
      double precision x,xt(lxt)
      data ilo /1/
      save ilo  
      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
         if (x .ge. xt(lxt))            go to 110
         if (lxt .le. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt
c
   20 if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100
c
c              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
      istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)                go to 35
         if (x .ge. xt(ilo))            go to 50
         istep = istep*2
                                        go to 31
   35 ilo = 1
      if (x .lt. xt(1))                 go to 90
                                        go to 50
c              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         if (ihi .ge. lxt)              go to 45
         if (x .lt. xt(ihi))            go to 50
         istep = istep*2
                                        go to 41
   45 if (x .ge. xt(lxt))               go to 110
      ihi = lxt
c
c           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)              go to 100
c     note. it is assumed that middle = ilo in case ihi = ilo+1 .
      if (x .lt. xt(middle))            go to 53
         ilo = middle
                                        go to 50
   53    ihi = middle
                                        go to 50
c**** set output and return.
   90 mflag = -1
      left = 1
                                        return
  100 mflag = 0
      left = ilo
                                        return
  110 mflag = 1
	  if (x .eq. xt(lxt)) mflag = 0
      left = lxt
  111 if (left .eq. 1)                  return
	  left = left - 1
	  if (xt(left) .lt. xt(lxt))        return
										go to 111
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) The corresponding approximate eigenvectors;
c
c      (2) An orthonormal (Lanczos) basis for the associated approximate
c          invariant subspace;
c
c      (3) Both.
c
c  Matrices A and B are stored in LAPACK-style band form.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  (Lanczos) basis is always computed.  There is an additional storage cost 
c  of n*nev if both are requested (in this case a separate array Z must be 
c  supplied).
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are called Ritz values and Ritz vectors respectively.  They are referred 
c  to as such in the comments that follow.  The computed orthonormal basis 
c  for the invariant subspace corresponding to these Ritz values is referred 
c  to as a Lanczos basis.
c
c \Usage
c   call dsband
c      ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, N, AB, MB, LDA, 
c        RFAC, KL, KU, WHICH, BMAT, NEV, TOL, RESID, NCV, V, 
c        LDV, IPARAM, WORKD, WORKL, LWORKL, IWORK, INFO )
c
c \Arguments
c
c  RVEC    Logical (INPUT)
c          Specifies whether Ritz vectors corresponding to the Ritz value 
c          approximations to the eigenproblem A*z = lambda*B*z are computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute the associated Ritz vectors. 
c
c  HOWMNY  Character*1  (INPUT) 
c          Specifies how many Ritz vectors are wanted and the form of Z
c          the matrix of Ritz vectors. See remark 1 below.
c          = 'A': compute all Ritz vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' , SELECT is not referenced.
c
c  D       Double precision array of dimension NEV.  (OUTPUT)
c          On exit, D contains the Ritz value approximations to the
c          eigenvalues of A*z = lambda*B*z. The values are returned
c          in ascending order. If IPARAM(7) = 3,4,5 then D represents
c          the Ritz values of OP computed by dsaupd transformed to
c          those of the original eigensystem A*z = lambda*B*z. If 
c          IPARAM(7) = 1,2 then the Ritz values of OP are the same 
c          as the those of A*z = lambda*B*z.
c
c  Z       Double precision N by NEV array if HOWMNY = 'A'.  (OUTPUT)
c          On exit, Z contains the B-orthonormal Ritz vectors of the
c          eigensystem A*z = lambda*B*z corresponding to the Ritz
c          value approximations.
c
c          If  RVEC = .FALSE. then Z is not referenced.
c          NOTE: The array Z may be set equal to first NEV columns of the 
c          Lanczos basis array V computed by DSAUPD.
c
c  LDZ     Integer.  (INPUT) 
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
c
c  SIGMA   Double precision  (INPUT)
c          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
c          IPARAM(7) = 1 or 2.
c 
c  N       Integer.  (INPUT) 
c          Dimension of the eigenproblem.  
c 
c  AB      Double precision array of dimension LDA by N. (INPUT)
c          The matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c
c  MB      Double precision array of dimension LDA by N. (INPUT)
c          The matrix M in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set. 
c          The j-th column of M is stored in the j-th column of the
c          array AB as follows:
c          MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c          Not referenced if IPARAM(7) = 1
c
c  LDA     Integer. (INPUT)
c          Leading dimension of AB, MB, RFAC.
c
c  RFAC    Double precision array of LDA by N. (WORKSPACE/OUTPUT)
c          RFAC is used to store the LU factors of MB when IPARAM(7) = 2 
c          is invoked.  It is used to store the LU factors of
c          (A-sigma*M) when IPARAM(7) = 3,4,5 is invoked.
c          It is not referenced when IPARAM(7) = 1.
c
c  KL      Integer. (INPUT)
c          Max(number of subdiagonals of A, number of subdiagonals of M)
c
c  KU      Integer. (OUTPUT)
c          Max(number of superdiagonals of A, number of superdiagonals of M)
c
c  WHICH   Character*2.  (INPUT)
c          When IPARAM(7)= 1 or 2,  WHICH can be set to any one of
c          the following.
c  
c            'LM' -> want the NEV eigenvalues of largest magnitude.
c            'SM' -> want the NEV eigenvalues of smallest magnitude.
c            'LA' -> want the NEV eigenvalues of largest REAL part.
c            'SA' -> want the NEV eigenvalues of smallest REAL part.
c            'BE' -> Compute NEV eigenvalues, half from each end of the 
c                    spectrum.  When NEV is odd, compute one more from 
c                    the high end than from the low end. 
c
c          When IPARAM(7) = 3, 4, or 5,  WHICH should be set to 'LM' only. 
c          
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x

c  NEV     Integer. (INPUT)
c          Number of eigenvalues of OP to be computed.
c   
c  TOL     Double precision scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value 
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
c          If TOL .LE. 0. is passed a default is set:
c          DEFAULT = DLAMCH('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine DLAMCH).
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT:
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c          Represents the dimension of the Lanczos basis constructed
c          by dsaupd for OP.
c
c  V       Double precision array N by NCV.  (OUTPUT)
c          Upon INPUT: the NCV columns of V contain the Lanczos basis 
c                      vectors as constructed by dsaupd for OP.
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns 
c                       represent the Ritz vectors that span the desired 
c                       invariant subspace.
c          NOTE: The array Z may be set equal to first NEV columns of the 
c          Lanczos basis vector array V computed by dsaupd. In this case
c          if RVEC=.TRUE., the first NCONV=IPARAM(5) columns of V contain
c          the desired Ritz vectors. 
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: 
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          It is set to 1 in this subroutine.  The user do not need
c          to set this parameter.
c          ------------------------------------------------------------
c          ISHIFT = 1: exact shifts with respect to the reduced 
c                      tridiagonal matrix T.  This is equivalent to 
c                      restarting the iteration with a starting vector 
c                      that is a linear combination of Ritz vectors 
c                      associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = No longer referenced. 
c
c          IPARAM(3) = MXITER
c          On INPUT:  max number of Arnoldi update iterations allowed.
c          On OUTPUT: actual number of Arnoldi update iterations taken.
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" eigenvalues.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used. 
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4,5; See under \Description of dsband for the 
c          five modes available.
c
c          IPARAM(8) = NP
c          Not referenced.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.
c
c WORKD    Double precision work array of length at least 3*n. (WORKSPACE)
c
c WORKL    Double precision work array of length LWORKL.  (WORKSPACE)
c
c LWORKL   Integer.  (INPUT)
c          LWORKL must be at least NCV**2 + 8*NCV.
c
c IWORK    Integer array of dimension at least N. (WORKSPACE)
c          Used when IPARAM(7)=2,3,4,5 to store the pivot information in the 
c          factorization of M or (A-SIGMA*M).
c            
c INFO     Integer.  (INPUT/OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)  
c                returns the number of wanted converged Ritz values.
c          =  3: No shifts could be applied during a cycle of the 
c                Implicitly restarted Arnoldi iteration. One possibility 
c                is to increase the size of NCV relative to NEV. 
c                See remark 4 in DSAUPD.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Informational error from LAPACK routine dsteqr.
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: NEV and WHICH = 'BE' are incompatible.
c          = -13: HOWMNY must be one of 'A' or 'P'
c          = -14: DSAUPD did not find any eigenvalues to sufficient
c                 accuracy.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current
c                   Arnoldi factorization.
c
c\Routines called:
c     dsaupd  ARPACK reverse communication interface routine.
c     dseupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     dgbtrf  LAPACK band matrix factorization routine.
c     dgbtrs  LAPACK band linear system solve routine. 
c     dlacpy  LAPACK matrix copy routine.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the dot product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dgbmv   Level 2 BLAS that computes the band matrix vector product.
c
c\Remarks
c  1. The converged Ritz values are always returned in increasing 
c     (algebraic) order.
c
c  2. Currently only HOWMNY = 'A' is implemented. It is included at this
c     stage for the user who wants to incorporate it.
c
c\Author    
c     Danny Sorensen
c     Richard Lehoucq
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c$$$      subroutine MyDsband(select,d,z,ldz,sigma,n,ab,mb,lda,rfac,ldrfac,k,nev,
c$$$     >     tol,resid,ncv,v,ldv,iparam,workd,workl,lworkl,iwork,info)
c$$$
c$$$      character        which*2, bmat, howmny
c$$$      integer          n, lda, ldrfac, k, nev, ncv, ldv, ldz, lworkl, info  
c$$$      Double precision tol, sigma
c$$$      logical          rvec
c$$$
c$$$      integer          iparam(*), iwork(*)
c$$$      logical          select(*)
c$$$      Double precision d(*), resid(*), v(ldv,*), z(ldz,*), ab(lda,n), mb(lda,n), rfac(ldrfac,n), workd(*), workl(*)
c$$$
c$$$      integer          ipntr(14)
c$$$
c$$$      integer          ido, i, j, Row, Col, type, ierr,ix,ixp
c$$$
c$$$      Double precision one, zero
c$$$      parameter        (one = 1.0d0, zero = 0.0d0)
c$$$
c$$$      Double precision ddot, dnrm2, dlapy2
c$$$      external         ddot, dcopy, dgbmv, dgbtrf, dgbtrs, dnrm2, dlapy2, dlacpy
c$$$
c$$$c     iparam(3) : Max number of Arnoldi iterations
c$$$      iparam(3) = 300
c$$$      iparam(7) = 3
c$$$      rvec = .TRUE.
c$$$      howmny = 'A'
c$$$      which = 'LA'
c$$$      bmat = 'G'
c$$$      type = 4 
c$$$      ido = 0
c$$$      iparam(1) = 1
c$$$
c$$$      rfac = 0.0d0
c$$$      do i = 1,n
c$$$         do j = 1,min(i+k,n)
c$$$            Row = k+1+i-j
c$$$            Col = j
c$$$            rfac(k+Row,Col) = ab(Row,Col) - sigma*mb(Row,Col)
c$$$         enddo
c$$$         do j = max(1,i-k),i-1
c$$$            Row = 2*k+1
c$$$            Col = j
c$$$            rfac(Row+i-j,j) = rfac(Row+j-i,i)
c$$$c            rfac(Row+j-i,i) = rfac(Row+i-j,j)
c$$$         enddo
c$$$      enddo
c$$$
c$$$      do ix = 1,ldrfac
c$$$         do ixp = 1,n
c$$$            write(22,11,advance='NO') rfac(ix,ixp)
c$$$         enddo
c$$$         write(22,*)
c$$$      enddo
c$$$ 11   format(f16.12,1p,100e20.12)
c$$$
c$$$      call dgbtrf(n,n,k,k,rfac,ldrfac,iwork,ierr)
c$$$      if ( ierr .ne. 0 )  then
c$$$         print*, ' '
c$$$         print*, '_SBAND: Error with _gbtrf:',ierr
c$$$         print*, ' '
c$$$         go to 9000
c$$$      end if
c$$$
c$$$  90  continue 
c$$$
c$$$      call dsaupd(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
c$$$
c$$$      if (ido .eq. -1) then
c$$$         call dsbmv('U',n,k,1.0d0,mb,lda,workd(ipntr(1)),1,0.0d0,workd(ipntr(2)),1)
c$$$         call dgbtrs('Notranspose',n,k,k,1,rfac,ldrfac,iwork,workd(ipntr(2)),n,ierr)
c$$$         if (ierr .ne. 0) then
c$$$            print*, ' ' 
c$$$            print*, '_SBAND: Error with _gbtrs.'
c$$$            print*, ' ' 
c$$$            go to 9000
c$$$         end if
c$$$      else if (ido .eq. 1) then
c$$$         call dcopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
c$$$         call dgbtrs('Notranspose',n,k,k,1,rfac,ldrfac,iwork,workd(ipntr(2)),n,ierr)
c$$$         if (ierr .ne. 0) then 
c$$$            print*, ' '
c$$$            print*, '_SBAND: Error with _gbtrs.' 
c$$$            print*, ' '
c$$$            go to 9000
c$$$         end if
c$$$      else if (ido .eq. 2) then
c$$$         call dsbmv('U',n,k,1.0d0,mb,lda,workd(ipntr(1)),1,0.0d0,workd(ipntr(2)),1)
c$$$      else 
c$$$         if ( info .lt. 0) then
c$$$            print *, ' '
c$$$            print *, ' Error with _saupd info = ',info
c$$$            print *, ' Check the documentation of _saupd '
c$$$            print *, ' '
c$$$            go to 9000
c$$$         else 
c$$$            if ( info .eq. 1) then
c$$$c     print *, ' '
c$$$c     print *, ' Maximum number of iterations reached.'
c$$$c     print *, ' '
c$$$            else if ( info .eq. 3) then
c$$$               print *, ' '
c$$$               print *, ' No shifts could be applied during implicit'
c$$$               print *, 'Arnoldi update, try increasing NCV.'
c$$$               print *, ' '
c$$$            end if
c$$$            if (iparam(5) .gt. 0) then
c$$$               call dseupd(rvec,'A',select,d,z,ldz,sigma,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
c$$$               if ( info .ne. 0) then
c$$$                  print *, ' ' 
c$$$                  print *, ' Error with _neupd = ', info
c$$$                  print *, ' Check the documentation of _neupd '
c$$$                  print *, ' ' 
c$$$                  go to 9000
c$$$               endif
c$$$            endif
c$$$         endif
c$$$         go to 9000
c$$$      endif
c$$$
c$$$      go to 90 
c$$$ 9000 continue
c$$$
c$$$      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dsband( rvec, howmny, select, d, z, ldz, sigma, 
     &           n, ab, mb, lda, rfac, kl, ku, which, bmat, nev, 
     &           tol, resid, ncv, v, ldv, iparam, workd, workl, 
     &           lworkl, iwork, info)
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c 
      character        which*2, bmat, howmny
      integer          n, lda, kl, ku, nev, ncv, ldv,
     &                 ldz, lworkl, info  
      Double precision
     &                 tol, sigma
      logical          rvec
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer          iparam(*), iwork(*)
      logical          select(*)
      Double precision
     &                 d(*), resid(*), v(ldv,*), z(ldz,*),
     &                 ab(lda,*), mb(lda,*), rfac(lda,*), 
     &                 workd(*), workl(*)
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer          ipntr(14)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer          ido, i, j, type, imid, itop, ibot, ierr
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &                  one, zero
      parameter        (one = 1.0, zero = 0.0)
c
c
c     %-----------------------------%
c     | LAPACK & BLAS routines used |
c     %-----------------------------%
c
      Double precision
     &                 ddot, dnrm2, dlapy2
      external         ddot, dcopy, dgbmv, dgbtrf, 
     &                 dgbtrs, dnrm2, dlapy2, dlacpy
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c     
c     %----------------------------------------------------------------%
c     | Set type of the problem to be solved. Check consistency        |
c     | between BMAT and IPARAM(7).                                    |
c     | type = 1 --> Solving standard problem in regular mode.         |
c     | type = 2 --> Solving standard problem in shift-invert mode.    | 
c     | type = 3 --> Solving generalized problem in regular mode.      |
c     | type = 4 --> Solving generalized problem in shift-invert mode. |
c     | type = 5 --> Solving generalized problem in Buckling mode.     |
c     | type = 6 --> Solving generalized problem in Cayley mode.       | 
c     %----------------------------------------------------------------%
c
      if ( iparam(7) .eq. 1 ) then
         type = 1
      else if ( iparam(7) .eq. 3 .and. bmat .eq. 'I') then
         type = 2
      else if ( iparam(7) .eq. 2 ) then
         type = 3
      else if ( iparam(7) .eq. 3 .and. bmat .eq. 'G') then
         type = 4 
      else if ( iparam(7) .eq. 4 ) then
         type = 5
      else if ( iparam(7) .eq. 5 ) then 
         type = 6
      else
         print*, ' '
         print*, 'BMAT is inconsistent with IPARAM(7).'
         print*, ' ' 
         go to 9000
      end if
c
c     %------------------------%
c     | Initialize the reverse |
c     | communication flag.    |
c     %------------------------%
c
      ido   = 0
c
c     %----------------%
c     | Exact shift is |
c     | used.          | 
c     %----------------%
c 
      iparam(1) = 1
c
c     %-----------------------------------%
c     | Both matrices A and M are stored  |
c     | between rows itop and ibot.  Imid |
c     | is the index of the row that      |
c     | stores the diagonal elements.     |
c     %-----------------------------------%
c
      itop = kl + 1
      imid = kl + ku + 1
      ibot = 2*kl + ku + 1
c
      if ( type .eq. 2 .or. type .eq. 6 .and. bmat .eq. 'I' ) then
c
c         %----------------------------------%
c         | Solving a standard eigenvalue    |
c         | problem in shift-invert or       |
c         | Cayley mode. Factor (A-sigma*I). |
c         %----------------------------------%
c
          call dlacpy ('A', ibot, n, ab, lda, rfac, lda )
          do 10 j = 1, n
             rfac(imid,j) =  ab(imid,j) - sigma
  10      continue
          call dgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr )
          if (ierr .ne. 0) then
             print*, ' ' 
             print*, ' _SBAND: Error with _gbtrf. '
             print*, ' '
             go to  9000
          end if
c
      else if ( type .eq. 3 ) then
c
c        %----------------------------------------------%
c        | Solving generalized eigenvalue problem in    |
c        | regular mode. Copy M to rfac and Call LAPACK |
c        | routine dgbtrf to factor M.                  |
c        %----------------------------------------------%
c
         call dlacpy ('A', ibot, n, mb, lda, rfac, lda )
         call dgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr) 
         if (ierr .ne. 0) then 
             print*, ' '
             print*,'_SBAND:  Error with _gbtrf.'
             print*, ' ' 
             go to 9000 
         end if
c
      else if ( type .eq. 4 .or. type .eq. 5 .or. type .eq. 6 
     &         .and. bmat .eq. 'G' ) then
c
c        %-------------------------------------------%
c        | Solving generalized eigenvalue problem in |
c        | shift-invert, Buckling, or Cayley mode.   |
c        %-------------------------------------------%
c 
c        %-------------------------------------%
c        | Construct and factor (A - sigma*M). |
c        %-------------------------------------%
c
         do 60 j = 1,n
            do 50 i = itop, ibot 
               rfac(i,j) = ab(i,j) - sigma*mb(i,j)
  50        continue
  60     continue
c
         call dgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr)
         if ( ierr .ne. 0 )  then
             print*, ' '
             print*, '_SBAND: Error with _gbtrf.'
             print*, ' '
             go to 9000
         end if
c 
      end if 
c
c     %--------------------------------------------%
c     |  M A I N   L O O P (reverse communication) |
c     %--------------------------------------------%
c
  90  continue 
c
      call dsaupd ( ido, bmat, n, which, nev, tol, resid, ncv,
     &              v, ldv, iparam, ipntr, workd, workl, lworkl,
     &              info )
c
      if (ido .eq. -1) then
c
         if ( type .eq. 1) then
c
c           %----------------------------%
c           | Perform  y <--- OP*x = A*x |
c           %----------------------------%
c
            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
c
         else if ( type .eq. 2 ) then
c
c           %----------------------------------%
c           |             Perform              |
c           | y <--- OP*x = inv[A-sigma*I]*x   |
c           | to force the starting vector     |
c           | into the range of OP.            |
c           %----------------------------------%
c
            call dcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                    iwork, workd(ipntr(2)), n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' ' 
               print*, ' _SBAND: Error with _bgtrs. '
               print*, ' '
               go to 9000
            end if
c
         else if ( type .eq. 3 ) then
c
c           %-----------------------------------%
c           | Perform  y <--- OP*x = inv[M]*A*x |
c           | to force the starting vector into | 
c           | the range of OP.                  |
c           %-----------------------------------%
c
            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                  lda, workd(ipntr(1)), 1, zero, 
     &                  workd(ipntr(2)), 1)
            call dcopy(n, workd(ipntr(2)), 1, workd(ipntr(1)), 1)
            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                    iwork, workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print*, ' '
               print*, '_SBAND: Error with sbgtrs.'
               print*, ' '
               go to 9000
            end if
c
         else if ( type .eq. 4 ) then
c
c           %-----------------------------------------%
c           | Perform y <-- OP*x                      |
c           |           = inv[A-SIGMA*M]*M            | 
c           | to force the starting vector into the   |
c           | range of OP.                            |
c           %-----------------------------------------%
c
            call dgbmv('Notranspose', n, n, kl, ku, one, mb(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                   iwork, workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print*, ' ' 
               print*, '_SBAND: Error with _gbtrs.'
               print*, ' ' 
               go to 9000
            end if
c
         else if ( type .eq. 5) then
c
c           %---------------------------------------% 
c           | Perform y <-- OP*x                    |
c           |    = inv[A-SIGMA*M]*A                 |
c           | to force the starting vector into the |
c           | range of OP.                          |
c           %---------------------------------------%
c
            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
            call dgbtrs('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                   iwork, workd(ipntr(2)), n, ierr)
c
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' _SBAND: Error with _gbtrs. '
               print*, ' '
               go to 9000
            end if
c
         else if ( type .eq. 6 ) then
c
c           %---------------------------------------%
c           | Perform y <-- OP*x                    |
c           | = (inv[A-SIGMA*M])*(A+SIGMA*M)*x      | 
c           | to force the starting vector into the |
c           | range of OP.                          | 
c           %---------------------------------------%
c
            if ( bmat .eq. 'G' ) then
               call dgbmv('Notranspose', n, n, kl, ku, one, 
     &                    ab(itop,1), lda, workd(ipntr(1)), 1, 
     &                    zero, workd(ipntr(2)), 1)
               call dgbmv('Notranspose', n, n, kl, ku, sigma, 
     &                    mb(itop,1), lda, workd(ipntr(1)), 1, 
     &                    one, workd(ipntr(2)), 1)
            else 
               call dcopy(n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
               call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                    lda, workd(ipntr(1)), 1, sigma, 
     &                    workd(ipntr(2)), 1)
            end if   
c
            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                   iwork, workd(ipntr(2)), n, ierr)
c
            if (ierr .ne. 0) then 
               print*, ' '
               print*, '_SBAND: Error with _gbtrs.' 
               print*, ' '
               go to 9000
            end if
c
         end if
c
      else if (ido .eq. 1) then
c
         if ( type .eq. 1) then
c
c           %----------------------------%
c           | Perform  y <--- OP*x = A*x |
c           %----------------------------%
c
            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
c
         else if ( type .eq. 2) then
c
c              %----------------------------------%
c              |             Perform              |
c              | y <--- OP*x = inv[A-sigma*I]*x.  |
c              %----------------------------------%
c
               call dcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
               call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                       iwork, workd(ipntr(2)), n, ierr)
               if ( ierr .ne. 0 ) then
                  print*, ' '
                  print*, '_SBAND: Error with _gbtrs.' 
                  print*, ' '
                  go to 9000
               end if
c
         else if ( type .eq. 3 ) then
c
c           %-----------------------------------%
c           | Perform  y <--- OP*x = inv[M]*A*x |
c           %-----------------------------------%
c
            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                  lda, workd(ipntr(1)), 1, zero, 
     &                  workd(ipntr(2)), 1)
            call dcopy(n, workd(ipntr(2)), 1, workd(ipntr(1)), 1) 
            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                    iwork, workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print*, ' '
               print*, '_SBAND: error with _bgtrs.'
               print*, ' ' 
               go to 9000
            end if
c
         else if ( type .eq. 4 ) then
c
c           %-------------------------------------%
c           | Perform y <-- inv(A-sigma*M)*(M*x). |
c           | (M*x) has been computed and stored  |
c           | in workd(ipntr(3)).                 |           
c           %-------------------------------------%
c
            call dcopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                    iwork, workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then 
               print*, ' '
               print*, '_SBAND: Error with _gbtrs.' 
               print*, ' '
               go to 9000
            end if
c 
         else if ( type .eq. 5 ) then
c
c           %-------------------------------% 
c           | Perform y <-- OP*x            |
c           |    = inv[A-SIGMA*M]*A*x       |
c           | B*x = A*x has been computed   |
c           | and saved in workd(ipntr(3)). |
c           %-------------------------------%
c
            call dcopy (n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
            call dgbtrs('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                   iwork, workd(ipntr(2)), n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' _SBAND: Error with _gbtrs. '
               print*, ' '
               go to 9000
            end if
c
         else if ( type .eq. 6) then
c
c           %---------------------------------%
c           | Perform y <-- OP*x              |
c           | = inv[A-SIGMA*M]*(A+SIGMA*M)*x. | 
c           | (M*x) has been saved in         |
c           | workd(ipntr(3)).                |
c           %---------------------------------%
c
            if ( bmat .eq. 'G' ) then
               call dgbmv('Notranspose', n, n, kl, ku, one, 
     &                    ab(itop,1), lda, workd(ipntr(1)), 1, 
     &                    zero, workd(ipntr(2)), 1)
               call daxpy( n, sigma, workd(ipntr(3)), 1, 
     &                    workd(ipntr(2)), 1 )
            else 
               call dcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
               call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                    lda, workd(ipntr(1)), 1, sigma, 
     &                    workd(ipntr(2)), 1)
            end if
            call dgbtrs('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                   iwork, workd(ipntr(2)), n, ierr)
c
         end if
c
      else if (ido .eq. 2) then
c
c        %----------------------------------%
c        |        Perform y <-- B*x         | 
c        | Note when Buckling mode is used, |
c        | B = A, otherwise B=M.            | 
c        %----------------------------------%
c
         if (type .eq. 5) then
c
c           %---------------------%
c           | Buckling Mode, B=A. |
c           %---------------------%
c
            call dgbmv('Notranspose', n, n, kl, ku, one, 
     &                ab(itop,1), lda, workd(ipntr(1)), 1, 
     &                zero, workd(ipntr(2)), 1)
         else

            call dgbmv('Notranspose', n, n, kl, ku, one, 
     &                mb(itop,1), lda, workd(ipntr(1)), 1, 
     &                zero, workd(ipntr(2)), 1)
         end if
c
      else 
c
c        %-----------------------------------------%
c        | Either we have convergence, or there is | 
c        | error.                                  |
c        %-----------------------------------------%
c
         if ( info .lt. 0) then
c
c           %--------------------------%
c           | Error message, check the |
c           | documentation in DSAUPD  |
c           %--------------------------%
c
            print *, ' '
            print *, ' Error with _saupd info = ',info
            print *, ' Check the documentation of _saupd '
            print *, ' '
            go to 9000
c
         else 
c
            if ( info .eq. 1) then
               print *, ' '
               print *, ' Maximum number of iterations reached.'
               print *, ' '
            else if ( info .eq. 3) then
               print *, ' '
               print *, ' No shifts could be applied during implicit',
     &                  ' Arnoldi update, try increasing NCV.'
               print *, ' '
            end if
c
            if (iparam(5) .gt. 0) then
c
               call dseupd ( rvec, 'A', select, d, z, ldz, sigma, 
     &                  bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &                  iparam, ipntr, workd, workl, lworkl, info )            
c
               if ( info .ne. 0) then
c 
c                 %------------------------------------%
c                 | Check the documentation of dneupd. |
c                 %------------------------------------%
c
                  print *, ' ' 
                  print *, ' Error with _neupd = ', info
                  print *, ' Check the documentation of _neupd '
                  print *, ' ' 
                  go to 9000
c 
               end if
c
            end if
c
         end if
c
         go to 9000
c
      end if
c
c     %----------------------------------------%
c     | L O O P  B A C K to call DSAUPD again. |
c     %----------------------------------------%
c
      go to 90 
c
 9000 continue
c
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine MyLargeDsband(NumStates,Shift,Energies,Psi,MatrixDim,H,S,LeadDim,HalfBandWidth,LUFac,Order,xDim)

      implicit none
      integer NumStates,ncv,lworkl,info,MatrixDim,LeadDim,HalfBandWidth,i,j,iparam(11),kl,ku,Order,xDim
      integer, allocatable :: iwork(:)
      logical, allocatable :: Select(:)
      double precision Shift,Tol,Energies(NumStates),Psi(MatrixDim,NumStates),LUFac(LeadDim,MatrixDim)
      double precision H(HalfBandWidth+1,MatrixDim),S(HalfBandWidth+1,MatrixDim)
      double precision, allocatable :: workd(:),workl(:),Residuals(:)
      double precision, allocatable :: TempPsi(:,:),TempEnergies(:),V(:,:)

      Tol = -1.0d0
      ncv = 2*NumStates
      lworkl = ncv*ncv+8*ncv
      allocate(Select(ncv),iwork(xDim),workd(3*xDim),workl(lworkl),Residuals(xDim))
      allocate(TempPsi(xDim,NumStates),TempEnergies(NumStates),V(xDim,ncv))
      info = 0
      iparam(3) = 700 !o - 300
      iparam(7) = 3

c      call MyDsband(Select,TempEnergies,TempPsi,MatrixDim,Shift,xDim,H,S,xDim+2*Order,LUFac,LeadDim,
c     >     Order,NumStates,Tol,Residuals,ncv,TempPsi,MatrixDim,iparam,workd,workl,lworkl,iwork,info)

      call dsband(.true.,'A',Select,TempEnergies,TempPsi,xDim,Shift, 
     &           xDim,H,S,xDim+2*Order,LUFac,Order,Order,'LA','G',NumStates, 
     &           Tol,Residuals,ncv,V,xDim,iparam, workd,workl, 
     &           lworkl,iwork,info)

      do i = 1,NumStates
         Energies(i) = TempEnergies(i)
         do j = 1,xDim
            Psi(j,i) = TempPsi(j,i)
         enddo
      enddo

c      call CalcEigenErrors(info,MatrixDim,H,S,HalfBandWidth,NumStates,Psi,Energies,NumStates)

      deallocate(Select,iwork,workd,workl,Residuals,TempPsi,TempEnergies)
      return
      end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine CalcEigenErrors(info,MatrixDim,H,S,HalfBandWidth,NumStates,Psi,Energies,MaxNumStates)

      integer info,MatrixDim,LeadDim,HalfBandWidth,NumStates,MaxNumStates
      double precision H(HalfBandWidth+1,MatrixDim),S(HalfBandWidth+1,MatrixDim)
      double precision Psi(MatrixDim,MaxNumStates),Energies(MaxNumStates,2)

      integer j
      double precision dnrm2
      double precision, allocatable :: HPsi(:),SPsi(:)

      if ( info .eq. 0) then

c Compute the residual norm: ||  A*x - lambda*x ||

       allocate(HPsi(MatrixDim))
       allocate(SPsi(MatrixDim))
       do j = 1,NumStates
        call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,H,HalfBandWidth+1,Psi(1,j),1,0.0d0,HPsi,1)
        call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,Psi(1,j),1,0.0d0,SPsi,1)
        call daxpy(MatrixDim,-Energies(j,1),SPsi,1,HPsi,1)
        Energies(j,2) = dnrm2(MatrixDim,HPsi,1)
        Energies(j,2) = Energies(j,2)/dabs(Energies(j,1))
       enddo
       deallocate(HPsi)
       deallocate(SPsi)
      else
       !write(6,*) ' '
       !write(6,*) ' Error with _sband, info= ', info
       !write(6,*) ' Check the documentation of _sband '
       !write(6,*) ' '
      end if

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function  BasisPhi(x,Left,Right,kLeft,kRight,Order,RDim,xPoints,
     >     xNumPoints,Deriv,Count)
      implicit none
      integer Left,Right,Order,LegPoints,RDim,xNumPoints,Deriv
      double precision xPoints(*)
      integer i,k,l,Count
      double precision MYBSpline
      double precision x,ax,bx,xIntScale,xScaledZero
      double precision kLeft,kRight,constLeft,constRight


c      print*,'Count = ',Count,'Left = ',Left,'Right = ',Right
      select case (Left)
      case (0)
         if(Count.eq.1) then
            BasisPhi = MYBSpline(Order,Deriv,xNumPoints,xPoints,2,x)
         endif
         if(Count.gt.1) then
            BasisPhi = MYBSpline(Order,Deriv,xNumPoints,xPoints,Count+1,x)
         endif
      case (1)
         if(Count.eq.1) then
            BasisPhi = MYBSpline(Order,Deriv,xNumPoints,xPoints,1,x) + 
     >           MYBSpline(Order,Deriv,xNumPoints,xPoints,2,x)
         endif
         if(Count.gt.1) then
            BasisPhi = MYBSpline(Order,Deriv,xNumPoints,xPoints,Count+1,x)
         endif
      case (2)
         BasisPhi = MYBSpline(Order,Deriv,xNumPoints,xPoints,Count,x)
      case (3)
         constLeft = -(kLeft*MYBSpline(Order,0,xNumPoints,xPoints,2,xPoints(1)) -
     >        MYBSpline(Order,1,xNumPoints,xPoints,2,xPoints(1))) / (
     >        MYBSpline(Order,0,xNumPoints,xPoints,1,xPoints(1))*kLeft - 
     >        MYBSpline(Order,1,xNumPoints,xPoints,1,xPoints(1)))
         if(Count.eq.1) then
            BasisPhi = MYBSpline(Order,Deriv,xNumPoints,xPoints,1,x)+ !c1=1, c2=constLeft
     >           MYBSpline(Order,Deriv,xNumPoints,xPoints,2,x)/constLeft
c$$$            BasisPhi = constLeft*MYBSpline(Order,Deriv,xNumPoints,xPoints,1,x)+ !c1=1, c2=constLeft
c$$$     >           MYBSpline(Order,Deriv,xNumPoints,xPoints,2,x)
         endif
         if(count.gt.1) then
            BasisPhi = MYBSpline(Order,Deriv,xNumPoints,xPoints,Count+1,x)
         endif
      end select
      
      
      select case (Right)
      case (0)
         if(Count.eq.RDim) then
            BasisPhi = MYBSpline(Order,Deriv,xNumPoints,xPoints,
     >           xNumPoints+Order-2,x)
         endif
      case (1)
         if(Count.eq.RDim) then
            BasisPhi = MYBSpline(Order,Deriv,xNumPoints,xPoints,
     >           xNumPoints+Order-2,x)+
     >           MYBSpline(Order,Deriv,xNumPoints,xPoints,
     >           xNumPoints+Order-1,x)
         endif
      case (2)
         if(Count.eq.RDim) then
            BasisPhi = MYBSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-1,x)
         endif
      case (3)
         constRight = -(kRight*MYBSpline(Order,0,xNumPoints,xPoints,xNumPoints+Order-2,xPoints(xNumPoints)) -
     >        MYBSpline(Order,1,xNumPoints,xPoints,xNumPoints+Order-2,xPoints(xNumPoints))) / ( !c2/c1
     >        MYBSpline(Order,0,xNumPoints,xPoints,xNumPoints+Order-1,xPoints(xNumPoints))*kRight - 
     >        MYBSpline(Order,1,xNumPoints,xPoints,xNumPoints+Order-1,xPoints(xNumPoints)))
         if(count.eq.RDim) then
            BasisPhi = MYBSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-2,x)/constRight +
     >           MYBSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-1,x)
c$$$            BasisPhi = MYBSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-2,x)+
c$$$     >           constRight*MYBSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-1,x)
         endif
      end select


      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine CheckBasisPhi(RMin,RMax,Left,Right,aLeft,aRight,RDim,RNumPoints,RPoints,Deriv,order,file)
      implicit none
      double precision, external :: BasisPhi
      integer MatrixDim,RDim,nch,beta,i,RNumPoints,Left,Right,Deriv,order,file,ix
      double precision R,RMin,RMax,RPoints(RNumPoints)
      double precision aLeft,aRight

      do ix=1,RDim
         R=RMin
         do while (R.le.RMax)
            write(file,*) R, BasisPhi(R,Left,Right,aLeft,aRight,order,RDim,RPoints,RNumPoints,Deriv,ix)      
            R = R+0.0001d0
         enddo
         write(file,*) ' '
      enddo

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c$$$      double precision function CalcPsi(eVec,x,NumIt,jchoose,Order,Deriv,xNumPoints,xPoints)
c$$$
c$$$      implicit none
c$$$      integer jchoose,Order,Deriv,xNumPoints,NumIt,i
c$$$      double precision MyBSpline,x,xPoints(*),eVec(NumIt,jchoose)
c$$$
c$$$      do i = 1,NumIt
c$$$         CalcPsi = CalcPsi + eVec(i,jchoose)*MyBSpline(Order,Deriv,xNumPoints,xPoints,1,x)
c$$$      enddo
c$$$
c$$$      return
c$$$      end function

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function phirecon(R,beta,evec,left,right,kLeft,kRight,RDim,MatrixDim,RNumPoints,RPoints,order,Deriv)
      double precision, external :: BasisPhi
      integer MatrixDim,RDim,nch,beta,i,RNumPoints,left,right,order,Deriv
      double precision R,evec(MatrixDim,MatrixDim),RPoints(RNumPoints),kLeft,kRight
      phirecon= 0.0d0
      do i = 1,RDim
         phirecon = phirecon + evec(i,beta)*BasisPhi(R,left,right,kLeft,kRight,order,RDim,RPoints,
     >        RNumPoints,Deriv,i)
      enddo
      return
      end function phirecon

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE EIGSRT(D,V,N,NP)
      implicit none
      integer N,NP,I,J,K
      double precision D(NP),V(NP,NP),P
      DO 13 I=1,N-1
        K=I 
        P=D(I) 
        DO 11 J=I+1,N 
          IF(D(J).GE.P)THEN 
            K=J 
            P=D(J) 
          ENDIF 
11      CONTINUE 
        IF(K.NE.I)THEN 
          D(K)=D(I) 
          D(I)=P 
          DO 12 J=1,N 
            P=V(J,I) 
            V(J,I)=V(J,K) 
            V(J,K)=P 
12        CONTINUE 
        ENDIF 
13    CONTINUE 
      RETURN 
      END 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine setupInterp(R,PotInterpOrder,Potknots,PotInterpPoints,P,Q,UCurves,Ubcoef,Pbcoef,Qbcoef,NumChan)

      implicit none
      integer PotInterpPoints, NumChan,PotInterpOrder
      double precision UCurves(NumChan,PotInterpPoints),P(NumChan,NumChan,PotInterpPoints) ! 
      double precision Q(NumChan,NumChan,PotInterpPoints),R(PotInterpPoints) ! 
      double precision Potknots(PotInterpPoints+PotInterpOrder), Ubcoef(NumChan,PotInterpPoints) ! 
      double precision Pbcoef(NumChan,NumChan,PotInterpPoints),Qbcoef(NumChan,NumChan,PotInterpPoints) ! 
      integer nch, mch

      call dbsnak(PotInterpPoints, R, PotInterpOrder, Potknots)

      do nch=1, NumChan
         call dbsint(PotInterpPoints,R,UCurves(nch,:),PotInterpOrder,Potknots,Ubcoef(nch,:)) ! 
         do mch=1,NumChan
            call dbsint(PotInterpPoints,R,P(nch,mch,:),PotInterpOrder,Potknots,Pbcoef(nch,mch,:)) ! 
            call dbsint(PotInterpPoints,R,Q(nch,mch,:),PotInterpOrder,Potknots,Qbcoef(nch,mch,:)) ! 
         enddo
      enddo

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine readCouplings(fileinU,fileinP,fileinQ,numR,R,U,P,Q,NumChan)

      implicit none
      integer fileinU,fileinP,fileinQ,numR,i,j,k,NumChan
      double precision R(numR),U(NumChan,numR),P(NumChan,NumChan,numR),Q(NumChan,NumChan,numR)

      open(fileinU)
      open(fileinP)
      open(fileinQ)
      do i = 1,numR
         read(fileinU,20) R(i),(U(k,i),k=1,NumChan)
         do j = 1,NumChan
            read(fileinP,20) (P(j,k,i),k=1,NumChan)
            read(fileinQ,20) (Q(j,k,i),k=1,NumChan)
         enddo
      enddo
      close(fileinU)
      close(fileinP)
      close(fileinQ)

 11   format(f16.12,1P,100e20.12)
 20   format(1P,100e20.12)

      end subroutine readCouplings
