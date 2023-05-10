module BasisSets
  implicit none
  double precision, allocatable :: lu(:,:,:),luxx(:,:,:)  ! left point
  double precision, allocatable :: cu(:,:,:),cuxx(:,:,:)  ! center point
  double precision, allocatable :: ru(:,:,:),ruxx(:,:,:)  ! right point
  double precision, allocatable :: pu(:,:,:),puxx(:,:,:)  ! primitive set

contains
  subroutine AllocateBasisSet(Left,Right,xNumPoints,LegPoints,Order,u,uxx)
    implicit none
    integer Left,Right,LegPoints,xNumPoints,Order,xDim
    double precision, allocatable :: u(:,:,:),uxx(:,:,:)
    xDim = xNumPoints+Order-3
    if (Left .eq. 2) xDim = xDim + 1
    if (Right .eq. 2) xDim = xDim + 1
    allocate(u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim))    
  end subroutine AllocateBasisSet

  !derived data type: set/basis <- have allocated arrays for XBounds
  !<- XDim
  !<- u, uxx (maybe just ux -- not a priority)

end module BasisSets


program HHL1DHyperspherical
  use BasisSets
  implicit none
  integer LegPoints,xNumPoints
  integer NumStates,PsiFlag,Order,Left,Right
  integer RSteps,CouplingFlag,CalcNewBasisFunc
  double precision alpha,mass,Shift,Shift2,NumStateInc,mi,ma,phiai,mgamma
  double precision RLeft,RRight,RDerivDelt,DD,L
  DOUBLE PRECISION RFirst,RLast,XFirst,XLast,StepX,StepR
  double precision xMin,xMax
  double precision, allocatable :: R(:)
  double precision, allocatable :: xPoints(:)

  logical, allocatable :: Select(:)

  integer iparam(11),ncv,info
  integer i,j,k,iR,NumFirst,NumBound
  integer LeadDim,MatrixDim,HalfBandWidth
  integer xDim, Nbs
  integer, allocatable :: iwork(:)
  integer, allocatable :: xBounds(:)
  double precision Tol,RChange
  double precision TotalMemory
  double precision mu, mu12, mu123, r0diatom, dDiatom, etaOVERpi, Pi
  double precision kLeft, kRight

  double precision, allocatable :: LUFac(:,:),workl(:)
  double precision, allocatable :: workd(:),Residuals(:)
  double precision, allocatable :: xLeg(:),wLeg(:)
!  double precision, allocatable :: u(:,:,:),uxx(:,:,:)
  double precision, allocatable :: S(:,:),H(:,:)
  double precision, allocatable :: lPsi(:,:),mPsi(:,:),rPsi(:,:),Energies(:,:)
  double precision, allocatable :: P(:,:),Q(:,:),dP(:,:)
  double precision ur(1:50000),acoef,bcoef,diff
  double precision sec,time,Rinitial,secp,timep,Rvalue, rbc, C4,lho
  double precision hbar, phase2B, amu,omega,Rstar, dum
  character*64 LegendreFile
  common /Rvalue/ Rvalue      
  hbar = 1.054571817d-34
  amu = 1.660539d-27
  !     read in number of energies and states to print
  read(5,*)
  read(5,*) NumStates,PsiFlag,CouplingFlag
  write(6,*) NumStates,PsiFlag,CouplingFlag

  !     read in Gauss-Legendre info
  read(5,*)
  read(5,*)
  read(5,1002) LegendreFile
  write(6,1002) LegendreFile
  read(5,*)
  read(5,*)
  read(5,*) LegPoints
  write(6,*) LegPoints,' LegPoints'

  !     read in boundary conditions
  read(5,*)
  read(5,*)
  read(5,*) Shift,Shift2,Order,Left,Right
  print*, 'Shift,Shift2, Order, Left, Right'
  print*, Shift,Shift2,Order,Left,Right

  !     read in potential parameters
  read(5,*)
  read(5,*)
  read(5,*) alpha,mi,ma,DD,L
  write(6,*) alpha,mi,ma,DD,L
  mi = mi*amu
  ma = ma*amu

  mu12=mi*ma/(mi+ma)  ! ion-atom reduced mass
  mu=dsqrt(ma/mi)
  Pi=dacos(-1.d0)
  write(6,*) 'mi = ', mi, 'ma = ', ma, 'mu12 = ', mu12, 'mu = ', mu
  mgamma = mu
  phiai = datan(mgamma)     !same as Seths beta !!! arctan mu, d?atan
  write(6,*) "phiai = ", phiai

  !     read in grid information
  read(5,*)
  read(5,*)
  read(5,*) xNumPoints, omega, Nbs, C4, phase2B
  write(6,*) xNumPoints, omega, Nbs, C4, phase2B
  omega = 2d0*Pi*omega
  lho = dsqrt(hbar/mi/omega)
  Rstar = dsqrt(mu12*C4/hbar**2)/lho
  C4 = C4/(hbar*omega*lho**4)
  !rbc = Rstar*sqrt(2d0)/((dble(Nbs) + 1d0)*Pi - phase2B)
  rbc = 0.22 * Rstar
  
  write(6,*) "C4 = ", C4
  write(6,*) "Rstar = ", Rstar
  write(6,*) "rbc = ", rbc

  ! Re-define in oscillator units

  ma = ma/mi
  mi = 1d0
  mu12=mi*ma/(mi+ma)  ! ion-atom reduced mass
  !      mu=dsqrt(mi*ma)

  read(5,*)
  read(5,*)
  read(5,*) RSteps,RDerivDelt,RFirst,RLast
  write(6,*) RSteps,RDerivDelt,RFirst,RLast

  !     c	XFirst=dsqrt(RFirst)
  !     c	XLast=dsqrt(RLast)
  !     c	XFirst=RFirst**(1.d0/3.d0)
  !     c	XLast=RLast**(1.d0/3.d0)

  !     SET UP A LOG-GRID IN THE HYPERRADIUS
!  XFirst = dlog10(RFirst)
!  XLast = dlog10(RLast)
!  StepX=(XLast-XFirst)/(RSteps-1.d0)

  allocate(R(RSteps))
  StepR = (RLast-RFirst)/(dble(RSteps-1))
  do i = 1,RSteps
     !     read(5,*) R(i)
     !     R(i)= (XFirst+(i-1)*StepX)**3 ! Cubic Grid
     R(i) = (RFirst+(i-1)*StepR) !Linear Grid
 !    R(i) = 10.d0**(XFirst+(i-1)*StepX) ! Log Grid
  enddo

  !      if (mod(xNumPoints,2) .ne. 0) then
  !         write(6,*) 'xNumPoints not divisible by 2'
  !         xNumPoints = (xNumPoints/2)*2
  !         write(6,*) '   truncated to ',xNumPoints
  !      endif

  allocate(xLeg(LegPoints),wLeg(LegPoints))
  call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  xDim = xNumPoints+Order-3
  if (Left .eq. 2) xDim = xDim + 1
  if (Right .eq. 2) xDim = xDim + 1

  MatrixDim = xDim
  HalfBandWidth = Order
  LeadDim = 3*HalfBandWidth+1

  TotalMemory = 2.0d0*(HalfBandWidth+1)*MatrixDim ! S, H
  TotalMemory = TotalMemory + 2.0d0*LegPoints*(Order+2)*xDim ! x splines
  TotalMemory = TotalMemory + 2.0d0*NumStates*NumStates ! P and Q matrices
  TotalMemory = TotalMemory + LeadDim*MatrixDim ! LUFac
  TotalMemory = TotalMemory + 4.0d0*NumStates*MatrixDim ! channel functions
  TotalMemory = TotalMemory + 4*xDim*xDim ! (CalcHamiltonian)
  TotalMemory = TotalMemory + LegPoints**2*xNumPoints ! (CalcHamiltonian)
  TotalMemory = 8.0d0*TotalMemory/(1024.0d0*1024.0d0)

  write(6,*)
  write(6,*) 'MatrixDim ',MatrixDim
  write(6,*) 'HalfBandWidth ',HalfBandWidth
  write(6,*) 'Approximate peak memory usage (in Mb) ',TotalMemory
  write(6,*)

  allocate(xPoints(xNumPoints))
  allocate(xBounds(xNumPoints+2*Order))



  !  allocate(u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim))
  !u(:,:,1) = alpha*pu(:,:,1) + pu(:,:,2)

!   write(6,*) "Left and Right:"
!   write(6,*) Left, Right
!   write(6,*) "Shape of u:"
!   write(6,*) shape(u)
!   write(6,*) "Shape of pu:"
!   write(6,*) shape(pu)
!   stop

  allocate(S(HalfBandWidth+1,MatrixDim),H(HalfBandWidth+1,MatrixDim))
  allocate(P(NumStates,NumStates),Q(NumStates,NumStates),dP(NumStates,NumStates))

  ncv = 2*NumStates
  LeadDim = 3*HalfBandWidth+1
  allocate(iwork(MatrixDim))
  allocate(Select(ncv))
  allocate(LUFac(LeadDim,MatrixDim))
  allocate(workl(ncv*ncv+8*ncv))
  allocate(workd(3*MatrixDim))
  allocate(lPsi(MatrixDim,ncv),mPsi(MatrixDim,ncv),rPsi(MatrixDim,ncv))
  allocate(Residuals(MatrixDim))
  allocate(Energies(ncv,2))
  info = 0
  iR=1
  CalcNewBasisFunc=1
  Tol=1e-20

  NumBound=0

  call AllocateBasisSet(2,2,xNumPoints,LegPoints,Order,pu,puxx)

  if (CouplingFlag .ne. 0) then

   call AllocateBasisSet(Left,Right,xNumPoints,LegPoints,Order,lu,luxx)
   call AllocateBasisSet(Left,Right,xNumPoints,LegPoints,Order,cu,cuxx)
   call AllocateBasisSet(Left,Right,xNumPoints,LegPoints,Order,ru,ruxx)

   write(6,*) 'Left, right, and center basis splines allocated.'

  end if

  RChange=100.d0
  do iR = RSteps,1,-1!RSteps
     NumFirst=NumStates
     if (R(iR).gt.RChange) then
        NumFirst=NumBound
     endif
     NumStateInc=NumStates-NumFirst
     !c----------------------------------------------------------------------------------------
     !     must move this block inside the loop over iR if the grid is adaptive

     print*, 'calling GridMaker'
     xMin = phiai + asin(dsqrt(mu/(1d0+mu**2))* rbc/R(iR)) ! theta_bc -- only need to call gridmaker once for each value of 
     xMax = Pi + phiai - asin(dsqrt(mu/(1d0+mu**2))* rbc/R(iR))

     print*, 'iR = ', iR
     ! Double check the scaling of R to be consistent with oscillator units in this boundary condition.
     kLeft = R(iR)/rbc * sqrt(((1+mu**2)/mu)-(rbc**2)/R(iR)**2)*COTAN(-(1/rbc)+phiai) !But this rbc is different from rbc
     kRight = -kLeft

     if(xMax.le.xMin) then
        write(6,*) "minimum hyperradius too small."
        stop
     endif
     call GridMakerIA(mu,mu12,phiai,Rstar,R(iR),rbc,xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
     if(CalcNewBasisFunc.eq.1) then !!Call primitive basis calcs for the R(iR)
        ! each triade of R (+,-,0) will have a different physical set -- 
        print*, 'done... Calculating Primitive Basis Functions'
        call CalcBasisFuncsBP(2,2, kLeft, kRight, Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,0,pu)
        call CalcBasisFuncsBP(2,2, kLeft, kRight, Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,2,puxx)
        ! two more calcbasisfuncs for kLeft +- delta R
     endif
     print*, 'done... Calculating overlap matrix'
     !     must move this block inside the loop if the grid is adaptive
     !----------------------------------------------------------------------------------------
     call CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,xDim,xNumPoints,pu,xBounds,HalfBandWidth,S)


!******** CONSTRUCTION OF CENTER BASIS SETS ********

      write(6,*) 'Constructing Center Basis Set'

      kLeft = R(iR)/rbc * sqrt(((1+mu**2)/mu)-(rbc**2)/R(iR)**2)*COTAN(-(1/rbc)+phiai) !But this rbc is different from rbc
      kRight = -kLeft

      call CalcBasisFuncsBP(Left,Right, kLeft, kRight, Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,0,cu)
      call CalcBasisFuncsBP(Left,Right, kLeft, kRight, Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,2,cuxx)
      
        call CalcHamiltonian(alpha,R(iR),mu,mi,phiai,C4,L,Order,xPoints,&
             LegPoints,xLeg,wLeg,xDim,xNumPoints,cu,cuxx,xBounds,HalfBandWidth,H)
        !            call MyLargeDsband(NumFirst,Shift2,NumStateInc,Energies,rPsi,
        !     >           Shift,MatrixDim,H,S,LUFac,LeadDim,HalfBandWidth,NumStates)
        !            call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,NumStates,lPsi,rPsi)
        call MyDsband(Select,Energies,rPsi,MatrixDim,Shift,MatrixDim,H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,&
             NumStates,Tol,Residuals,ncv,rPsi,MatrixDim,iparam,workd,workl,ncv*ncv+8*ncv,iwork,info)
        call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,lPsi,rPsi)
        call CalcEigenErrors(info,iparam,MatrixDim,H,HalfBandWidth+1,S,HalfBandWidth,NumStates,rPsi,Energies,ncv)


     if (CouplingFlag .ne. 0) then
        ! Inside this IF block we must have:
        ! 1) Construction of the physical basis set at the left and right point
        ! 2) ???

!******** CONSTRUCTION OF LEFT BASIS SETS ********

      write(6,*) 'Constructing Left Basis Set'

      RLeft = 0.99d0*R(iR) !R(iR)-RDerivDelt < might be better

      xMin = phiai + asin(dsqrt(mu/(1d0+mu**2))* rbc/RLeft) ! theta_bc -- only need to call gridmaker once for each value of 
      xMax = Pi + phiai - asin(dsqrt(mu/(1d0+mu**2))* rbc/RLeft)

      kLeft = RLeft/rbc * sqrt(((1+mu**2)/mu)-(rbc**2)/RLeft**2)*COTAN(-(1/rbc)+phiai) !But this rbc is different from rbc
      kRight = -kLeft

      !!Recompute kleft and kright at each triade
      
      call CalcBasisFuncsBP(Left,Right, kLeft, kRight, Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,0,lu)
      call CalcBasisFuncsBP(Left,Right, kLeft, kRight, Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,2,luxx)
      
        call CalcHamiltonian(alpha,RLeft,mu,mi,phiai,C4,L,Order,xPoints,&
             LegPoints,xLeg,wLeg,xDim,xNumPoints,lu,luxx,xBounds,HalfBandWidth,H)
        call MyDsband(Select,Energies,rPsi,MatrixDim,Shift,MatrixDim,H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,&
             NumStates,Tol,Residuals,ncv,rPsi,MatrixDim,iparam,workd,workl,ncv*ncv+8*ncv,iwork,info)
        call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,lPsi,rPsi)
        call CalcEigenErrors(info,iparam,MatrixDim,H,HalfBandWidth+1,S,HalfBandWidth,NumStates,rPsi,Energies,ncv)

!******** CONSTRUCTION OF RIGHT BASIS SETS ********

      write(6,*) 'Constructing Right Basis Set'

      RRight = 1.01d0*R(iR) !R(iR)+RDerivDelt

      xMin = phiai + asin(dsqrt(mu/(1d0+mu**2))* rbc/RRight) ! theta_bc -- only need to call gridmaker once for each value of 
      xMax = Pi + phiai - asin(dsqrt(mu/(1d0+mu**2))* rbc/RRight)

      kLeft = RRight/rbc * sqrt(((1+mu**2)/mu)-(rbc**2)/RRight**2)*COTAN(-(1/rbc)+phiai) !But this rbc is different from rbc
      kRight = -kLeft

      call CalcBasisFuncsBP(Left,Right, kLeft, kRight, Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,0,ru)
      call CalcBasisFuncsBP(Left,Right, kLeft, kRight, Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,2,ruxx)

         call CalcHamiltonian(alpha,RRight,mu,mi,phiai,C4,L,Order,xPoints,&
             LegPoints,xLeg,wLeg,xDim,xNumPoints,ru,ruxx,xBounds,HalfBandWidth,H)
        call MyDsband(Select,Energies,rPsi,MatrixDim,Shift,MatrixDim,H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,&
             NumStates,Tol,Residuals,ncv,rPsi,MatrixDim,iparam,workd,workl,ncv*ncv+8*ncv,iwork,info)
        call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,lPsi,rPsi)
        call CalcEigenErrors(info,iparam,MatrixDim,H,HalfBandWidth+1,S,HalfBandWidth,NumStates,rPsi,Energies,ncv)

     endif

!     write(6,*) 'Calculating Hamiltonian on primitive set with R = ', R(iR)
!     call CalcHamiltonian(alpha,R(iR),mu,mi,phiai,C4,L,Order,xPoints,LegPoints,&
!     xLeg,wLeg,xDim,xNumPoints,pu,puxx,xBounds,HalfBandWidth,H)

!!     !write(6,*) 'Calling MyLargeDsband'
!!$c     call MyLargeDsband(NumFirst,Shift2,NumStateInc,Energies,mPsi,
!!$c     >        Shift,MatrixDim,H,S,LUFac,LeadDim,HalfBandWidth,NumStates)
!!$c     write(6,*) 'done...'
!!$c     if (CouplingFlag .ne. 0) then
!!$c     if(iR.gt.1) then
!!$c     write(6,*) 'Calling FixPhase'
!!$c     call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,NumStates,rPsi,mPsi)
!!$c     endif
!!$c     endif

     call MyDsband(Select,Energies,mPsi,MatrixDim,Shift,MatrixDim,&
          H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,NumStates,&
          Tol,Residuals,ncv,mPsi,MatrixDim,iparam,workd,workl,&
          ncv*ncv+8*ncv,iwork,info)
     if (CouplingFlag .ne. 0) call FixPhase(NumStates,HalfBandWidth,&
          MatrixDim,S,ncv,rPsi,mPsi)

     call CalcEigenErrors(info,iparam,MatrixDim,H,HalfBandWidth+1,S,&
          HalfBandWidth,NumStates,mPsi,Energies,ncv)
     write(6,*) 'writing the energies'
     write(200,20) R(iR),(Energies(i,1), i = 1,min(NumStates,iparam(5)))

     write(6,*)
     write(6,*) 'RMid = ', R(iR)
     do i = 1,min(NumStates,iparam(5))
        write(6,*) 'Energy(',i,') = ',Energies(i,1),'  Error = ', Energies(i,2)
     enddo


!    Adjusting Shift
     ur(iR) = Energies(1,1)
     Shift = -200d0
     if(ur(iR).lt.0d0) then
        Shift = ur(iR)*10d0
     endif         

     if (CouplingFlag .ne. 0) then
        call CalcCoupling(NumStates,HalfBandWidth,MatrixDim,RDerivDelt,lPsi,mPsi,rPsi,S,P,Q,dP)

        write(101,*) R(iR)
        write(102,*) R(iR)
        write(103,*) R(iR)
        do i = 1,min(NumStates,iparam(5))
           write(101,20) (P(i,j), j = 1,min(NumStates,iparam(5)))
           write(102,20) (Q(i,j), j = 1,min(NumStates,iparam(5)))
           write(103,20) (dP(i,j), j = 1,min(NumStates,iparam(5)))
        enddo
     endif
     !         write(400,20) R(iR),R(iR)**3.0d0*(Energies(2,1)-Q(2,2))
     if (PsiFlag .ne. 0) then
        do i = 1,xNumPoints
           write(97,*) xPoints(i)
        enddo
        do i = 1,MatrixDim
           write(999+iR,20) (mPsi(i,j), j = 1,NumStates)
        enddo
        close(unit=999+iR)
     endif

  enddo

  deallocate(S,H)
  deallocate(Energies)
  deallocate(iwork)
  deallocate(Select)
  deallocate(LUFac)
  deallocate(workl)
  deallocate(workd)
  deallocate(lPsi,mPsi,rPsi)
  deallocate(Residuals)
  deallocate(P,Q,dP)
  deallocate(xPoints)
  deallocate(xLeg,wLeg)
  deallocate(xBounds)
!  deallocate(u,uxx)

  deallocate(R)

10 format(1P,100e25.15)
20 format(1P,100e16.8)
1002 format(a64)

  stop
end program HHL1DHyperspherical

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,xDim,&
     xNumPoints,u,xBounds,HalfBandWidth,S)
  implicit none
  integer Order,LegPoints,xDim,xNumPoints,xBounds(xNumPoints+2*Order),HalfBandWidth
  double precision xPoints(*),xLeg(*),wLeg(*)
  double precision S(HalfBandWidth+1,xDim)
  double precision u(LegPoints,xNumPoints,xDim)

  integer ix,ixp,kx,lx
  integer i1,i1p
  integer Row,NewRow,Col
  integer, allocatable :: kxMin(:,:),kxMax(:,:)
  double precision a,b,m
  double precision xTempS
  double precision ax,bx
  double precision, allocatable :: xIntScale(:),xS(:,:)

  allocate(xIntScale(xNumPoints),xS(xDim,xDim))
  allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))

  S = 0.0d0

  do kx = 1,xNumPoints-1
     ax = xPoints(kx)
     bx = xPoints(kx+1)
     xIntScale(kx) = 0.5d0*(bx-ax)
  enddo

  !      do ix=1,xNumPoints+2*Order
  !         print*, ix, xBounds(ix)
  !      enddo

  do ix = 1,xDim
     do ixp = 1,xDim
        kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
        kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
     enddo
  enddo

  do ix = 1,xDim
     do ixp = max(1,ix-Order),min(xDim,ix+Order)
        xS(ixp,ix) = 0.0d0
        do kx = kxMin(ixp,ix),kxMax(ixp,ix)
           xTempS = 0.0d0
           do lx = 1,LegPoints
              a = wLeg(lx)*xIntScale(kx)*u(lx,kx,ix)
              b = a*u(lx,kx,ixp)
              xTempS = xTempS + b
           enddo
           xS(ixp,ix) = xS(ixp,ix) + xTempS
        enddo
     enddo
  enddo

  do ix = 1,xDim
     Row=ix
     do ixp = max(1,ix-Order),min(xDim,ix+Order)
        Col = ixp
        if (Col .ge. Row) then
           NewRow = HalfBandWidth+1+Row-Col
           S(NewRow,Col) = xS(ixp,ix)
           !               write(26,*) ix,ixp,S(NewRow,Col)
        endif
     enddo
  enddo

  deallocate(xIntScale,xS)
  deallocate(kxMin,kxMax)

  return
end subroutine CalcOverlap
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CalcHamiltonian(alpha,R,mu,mi,phiai,C4,L,Order,xPoints,LegPoints,xLeg,wLeg,xDim,xNumPoints,u,uxx,xBounds,HalfBandWidth,H)
  implicit none
  double precision, external :: VSech
  integer Order,LegPoints,xDim,xNumPoints,xBounds(*),HalfBandWidth
  double precision alpha,R,mu,mgamma,phiai,C4,L,mi
  double precision xPoints(*),xLeg(*),wLeg(*)
  double precision H(HalfBandWidth+1,xDim)
  double precision u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim)

  integer ix,ixp,kx,lx
  integer i1,i1p
  integer Row,NewRow,Col
  integer, allocatable :: kxMin(:,:),kxMax(:,:)
  double precision a,b,m,Pi
  double precision Rall,rai
  double precision u1,sys_ss_pot,V12,V23,V31
  double precision VInt,VTempInt,potvalue, xTempV
  !     double precision TempPot,VInt,VTempInt
  double precision x,ax,bx,xScaledZero,xTempT,xTempS,xInt
  double precision, allocatable :: Pot(:,:)
  double precision, allocatable :: xIntScale(:),xT(:,:),xV(:,:)
  double precision, allocatable :: cosx(:,:),sinx(:,:)
  double precision, allocatable :: sinai(:,:)

  double precision mu12,r0diatom,dDiatom


  mgamma=dtan(phiai)
  allocate(xIntScale(xNumPoints),xT(xDim,xDim),xV(xDim,xDim))
  allocate(cosx(LegPoints,xNumPoints),sinx(LegPoints,xNumPoints))
  allocate(sinai(LegPoints,xNumPoints))
  allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
  allocate(Pot(LegPoints,xNumPoints))

  Pi = 3.1415926535897932385d0

  m = -1.0d0/(2.0d0*mu*R*R)
  !      phiai = datan(mgamma)
  ! Could create an array x(LegPoints,xNumPoints) that specifies all of the integration abcissa.


  do kx = 1,xNumPoints-1
     ax = xPoints(kx)
     bx = xPoints(kx+1)
     xIntScale(kx) = 0.5d0*(bx-ax)
     xScaledZero = 0.5d0*(bx+ax)
     do lx = 1,LegPoints
        x = xIntScale(kx)*xLeg(lx)+xScaledZero
        cosx(lx,kx) = dcos(x)
        sinai(lx,kx) = dsin(phiai-x)
        sinx(lx,kx) = dsin(x)
     enddo
  enddo

  do ix = 1,xDim
     do ixp = 1,xDim
        kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
        kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
     enddo
  enddo

  do kx = 1,xNumPoints-1
     do lx = 1,LegPoints

        rai = mu**(-0.5d0) * R*sinx(lx,kx) -  mu**0.5d0 * R*cosx(lx,kx) !dsqrt(mu/mu12)*R*dabs(sinai(lx,kx))))

        potvalue = -C4/rai**4 + 0.5d0*mu*(R*cosx(lx,kx))**2
        Pot(lx,kx) = alpha*potvalue
        !            write(24,*) kx, lx, Pot(lx,kx)
     enddo
  enddo

  do ix = 1,xDim
     do ixp = max(1,ix-Order),min(xDim,ix+Order)
        xT(ix,ixp) = 0.0d0
        xV(ix,ixp) = 0.0d0
        do kx = kxMin(ixp,ix),kxMax(ixp,ix)
           xTempT = 0.0d0
           xTempV = 0.0d0
           do lx = 1,LegPoints
              a = wLeg(lx)*xIntScale(kx)*u(lx,kx,ix) !bra
              xTempT = xTempT + a*uxx(lx,kx,ixp) !KE operator !ket
              xTempV = xTempV + a*(Pot(lx,kx))*u(lx,kx,ixp) !PE operator ket
           enddo
           xT(ix,ixp) = xT(ix,ixp) + xTempT
           xV(ix,ixp) = xV(ix,ixp) + xTempV
        enddo
     enddo
  enddo

  H = 0.0d0      
  do ix = 1,xDim
     Row=ix
     do ixp = max(1,ix-Order),min(xDim,ix+Order)
        Col = ixp
        if (Col .ge. Row) then
           NewRow = HalfBandWidth+1+Row-Col
           H(NewRow,Col) = (m*xT(ix,ixp)+xV(ix,ixp))
           !     write(25,*) ix,ixp,H(NewRow,Col)
        endif
     enddo
  enddo

  deallocate(Pot)
  deallocate(xIntScale,xT,xV)
  deallocate(cosx,sinx,sinai)
  deallocate(kxMin,kxMax)


  return
end subroutine CalcHamiltonian

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CalcPMatrix(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,P)
  implicit none
  integer NumStates,HalfBandWidth,MatrixDim
  double precision RDelt
  double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
  double precision S(HalfBandWidth+1,MatrixDim)
  double precision P(NumStates,NumStates)

  integer i,j,k
  double precision a,ddot
  double precision, allocatable :: TempPsi1(:),TempPsi2(:)

  allocate(TempPsi1(MatrixDim),TempPsi2(MatrixDim))

  a = 0.5d0/RDelt

  do j = 1,NumStates
     do k = 1,MatrixDim
        TempPsi1(k) = rPsi(k,j)-lPsi(k,j)
     enddo
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi1,1,0.0d0,TempPsi2,1)
     do i = 1,NumStates
        P(i,j) = a*ddot(MatrixDim,TempPsi2,1,mPsi(1,i),1)
     enddo
  enddo

  deallocate(TempPsi1,TempPsi2)

  return
end subroutine CalcPMatrix
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CalcQMatrix(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,Q)
  implicit none
  integer NumStates,HalfBandWidth,MatrixDim
  double precision RDelt
  double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
  double precision S(HalfBandWidth+1,MatrixDim)
  double precision Q(NumStates,NumStates)

  integer i,j,k
  double precision a,ddot
  double precision, allocatable :: TempPsi1(:),TempPsi2(:)

  allocate(TempPsi1(MatrixDim),TempPsi2(MatrixDim))

  a = 1.0d0/(RDelt**2)

  do j = 1,NumStates
     do k = 1,MatrixDim
        TempPsi1(k) = lPsi(k,j)+rPsi(k,j)-2.0d0*mPsi(k,j)
     enddo
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi1,1,0.0d0,TempPsi2,1)
     do i = 1,NumStates
        Q(i,j) = a*ddot(MatrixDim,TempPsi2,1,mPsi(1,i),1)
     enddo
  enddo

  deallocate(TempPsi1,TempPsi2)

  return
end subroutine CalcQMatrix
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
subroutine FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,mPsi,rPsi)
  implicit none
  integer NumStates,HalfBandWidth,MatrixDim,ncv
  double precision S(HalfBandWidth+1,MatrixDim),Psi(MatrixDim,ncv)
  double precision mPsi(MatrixDim,ncv),rPsi(MatrixDim,ncv)

  integer i,j
  double precision Phase,ddot
  double precision, allocatable :: TempPsi(:)

  allocate(TempPsi(MatrixDim))

  do i = 1,NumStates
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rPsi(1,i),1,0.0d0,TempPsi,1)
     Phase = ddot(MatrixDim,mPsi(1,i),1,TempPsi,1)
     if (Phase .lt. 0.0d0) then
        do j = 1,MatrixDim
           rPsi(j,i) = -rPsi(j,i)
        enddo
     endif
  enddo

  deallocate(TempPsi)

  return
end subroutine FixPhase
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GridMakerHHL(mu,mu12,mu123,phiai,R,r0,xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
  implicit none
  integer xNumPoints,CalcNewBasisFunc
  double precision mu,R,r0,xMin,xMax,xPoints(xNumPoints)
  double precision mu12,mu123,phiai
  integer i,j,k,OPGRID
  double precision Pi
  double precision r0New
  double precision xRswitch
  double precision xDelt,x0,x1,x2,x3,x4,deltax


  Pi = 3.1415926535897932385d0
  r0New=r0*2.0d0
  deltax = r0/R

  xRswitch = 20.0d0*r0New/Pi
  OPGRID=1
  write (6,*) 'xMin = ', xMin, 'xMax = ', xMax
  if((OPGRID.eq.1).and.(R.gt.xRswitch)) then
     print*, 'R>xRswitch!! using modified grid!!'
     x0 = xMin
     x1 = phiai - deltax  
     x2 = phiai + deltax
     x3 = xMax-deltax
     x4 = xMax
     k = 1
     xDelt = (x1-x0)/dfloat(xNumPoints/4)
     do i = 1,xNumPoints/4
        xPoints(k) = (i-1)*xDelt + x0
        !     write(6,*) k, xPoints(k)
        k = k + 1
     enddo
     xDelt = (x2-x1)/dfloat(xNumPoints/4)
     do i = 1,xNumPoints/4
        xPoints(k) = (i-1)*xDelt + x1
        !     write(6,*) k, xPoints(k)
        k = k + 1
     enddo
     xDelt = (x3-x2)/dfloat(xNumPoints/4)
     do i = 1,xNumPoints/4
        xPoints(k) = (i-1)*xDelt + x2
        !     write(6,*) k, xPoints(k)
        k = k + 1
     enddo
     xDelt = (x4-x3)/dfloat(xNumPoints/4-1)
     do i = 1, xNumPoints/4
        xPoints(k) = (i-1)*xDelt + x3
        !     write(6,*) k, xPoints(k)
        k = k + 1
     enddo

     !     FOR SMALL R, USE A LINEAR GRID
  else
     k = 1
     xDelt = (xMax-xMin)/dfloat(xNumPoints-1)
     do i = 1,xNumPoints
        xPoints(k) = (i-1)*xDelt + x0
        k = k + 1
     enddo
  endif

  !     Smooth Grid 

  write(20,*) 1, xPoints(1)
  do i = 2, xNumPoints-1
     xPoints(i)=(xPoints(i-1)+2.d0*xPoints(i)+xPoints(i+1))/4.d0
     write(20,*) i, xPoints(i)
  enddo
  write(20,*) xNumPoints, xPoints(xNumPoints)
  write(20,*) ' ' 
  !      write(96,15) (xPoints(k),k=1,xNumPoints)

15 format(6(1x,1pd12.5))



  return
end subroutine GridMakerHHL
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GridMakerIA(mu,mu12,phiai,Rstar,R,rbc,xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
  implicit none
  integer xNumPoints,CalcNewBasisFunc
  double precision mu,R,r0,xMin,xMax,xPoints(xNumPoints)
  double precision mu12,mu123,phiai,Rstar,rbc
  integer i,j,k,OPGRID,NP
  double precision Pi
  double precision r0New
  double precision xRswitch
  double precision xDelt,x0,x1,x2,x3,x4,deltax


  Pi = 3.1415926535897932385d0
  !r0New=r0*2.0d0
  deltax = 0.2d0*Rstar/R


  OPGRID=0
  NP = xNumPoints/4
  if(deltax.ge.(xMax - xMin)*0.25d0) then
     write(6,*) 'Switching to linear grid...'
     write(6,*) 'deltax = ', deltax
     write(6,*) 'deltax = ', deltax
     OPGRID = 0
  endif

  x0 = xMin
  x1 = x0 + deltax
  x2 = xMax - deltax
  x3 = xMax
  !      write(6,*) 'x0 = ', x0
  !      write(6,*) 'x1 = ', x1
  !      write(6,*) 'x2 = ', x2
  !      write(6,*) 'x3 = ', x3
  if((OPGRID.eq.1)) then
     !         write(6,*) "using modified grid:(R,deltax) = ", R, deltax

     k = 1
     xDelt = (x1-x0)
     do i = 1,NP-1
        xPoints(k) = (dble(i-1)/dble(NP-1))**2 * xDelt + x0
        !            write(6,*) k, xPoints(k), (dble(i-1)/dble(NP-1))**2 * xDelt
        k = k + 1
     enddo
     xDelt = (x2-x1)
     do i = 1,2*NP+1
        xPoints(k) = dble(i-1)/(2*NP+1)*xDelt + x1
        !            write(6,*) k, xPoints(k)
        k = k + 1
     enddo
     xDelt = (x3-x2)
     do i = 1,NP-1
        xPoints(k) = (dble(i-1)/dble(NP-1))**0.5d0*xDelt + x2
        !            write(6,*) k, xPoints(k),(dble(i-1)/dble(NP-1))**0.5d0*xDelt
        k = k + 1
     enddo
     xPoints(k)=x3
     !         write(6,*) k, xPoints(k),(dble(i-1)/dble(NP-1))**0.5d0*xDelt
     !     FOR SMALL R, USE A LINEAR GRID
  else
     write(6,*) "Using linear grid...(R,deltax) = ",R, deltax
     k = 1
     xDelt = (xMax-xMin)/dfloat(xNumPoints-1)
     do i = 1,xNumPoints
        xPoints(k) = (i-1)*xDelt + x0
        k = k + 1
     enddo
  endif

  !     Smooth Grid 

  do j=1,10
     do i = 2, xNumPoints-1
        xPoints(i)=(xPoints(i-1)+xPoints(i)+xPoints(i+1))/3.d0
     enddo
  enddo
  do i = 1, xNumPoints
     write(20,*) i, xPoints(i)
  enddo
  write(20,*) ' '
  !      write(96,'(100e20.10)') R, (xPoints(k),k=1,xNumPoints)

15 format(6(1x,1pd12.5))

  !      stop

  return
end subroutine GridMakerIA
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GridMaker(mu,R,r0,xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
  implicit none
  integer xNumPoints,CalcNewBasisFunc
  double precision mu,R,r0,xMin,xMax,xPoints(xNumPoints)

  integer i,j,k,OPGRID
  double precision Pi
  double precision r0New
  double precision xRswitch
  double precision xDelt,x0,x1,x2


  Pi = 3.1415926535897932385d0

  x0 = xMin
  x1 = xMax
  !     write(96,*) 'x0,x1=',x0,x1
  r0New=10.0d0*r0           !/2.0d0
  !     r0New=Pi/12.0*R
  xRswitch = 12.0d0*r0New/Pi
  OPGRID=1

  if((OPGRID.eq.1).and.(R.gt.xRswitch)) then
     print*, 'R>xRswitch!! using modified grid!!'
     x0 = xMin
     x1 = xMax - r0New/R  
     x2 = xMax
     k = 1
     xDelt = (x1-x0)/dfloat(xNumPoints/2)
     do i = 1,xNumPoints/2
        xPoints(k) = (i-1)*xDelt + x0
        !            print*, k, xPoints(k), xDelt
        k = k + 1
     enddo
     xDelt = (x2-x1)/dfloat(xNumPoints/2-1)
     do i = 1,xNumPoints/2
        xPoints(k) = (i-1)*xDelt + x1
        !            print*, k, xPoints(k), xDelt
        k = k + 1
     enddo
  else
     x0 = xMin
     x1 = xMax
     k = 1
     xDelt = (x1-x0)/dfloat(xNumPoints-1)
     do i = 1,xNumPoints
        xPoints(k) = (i-1)*xDelt + x0
        !            print*, k, xPoints(k), xDelt
        k = k + 1
     enddo
  endif

  !      write(96,15) (xPoints(k),k=1,xNumPoints)
15 format(6(1x,1pd12.5))



  return
end subroutine GridMaker
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


subroutine CalcCoupling(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,P,Q,dP)
  implicit none
  integer NumStates,HalfBandWidth,MatrixDim
  double precision RDelt
  double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
  double precision S(HalfBandWidth+1,MatrixDim),testorth
  double precision P(NumStates,NumStates),Q(NumStates,NumStates),dP(NumStates,NumStates)

  integer i,j,k
  double precision aP,aQ,ddot
  double precision, allocatable :: lDiffPsi(:),rDiffPsi(:),TempPsi(:),TempPsiB(:),rSumPsi(:)
  double precision, allocatable :: TempmPsi(:)

  allocate(lDiffPsi(MatrixDim),rDiffPsi(MatrixDim),TempPsi(MatrixDim),&
       TempPsiB(MatrixDim),rSumPsi(MatrixDim))
  allocate(TempmPsi(MatrixDim))

  aP = 0.5d0/RDelt
  aQ = aP*aP

  do j = 1,NumStates
     do k = 1,MatrixDim
        rDiffPsi(k) = rPsi(k,j)-lPsi(k,j)
        rSumPsi(k)  = lPsi(k,j)+mPsi(k,j)+rPsi(k,j)  ! Double check this is needed for <Phi'(R)|Phi'(R)>
        !            rSumPsi(k)  = lPsi(k,j)-2.0d0*mPsi(k,j)+rPsi(k,j)
        !            rSumPsi(k)  = lPsi(k,j)+rPsi(k,j)
     enddo
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rDiffPsi,1,0.0d0,TempPsi,1)   ! Calculate the vector S*rDiffPsi
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rSumPsi,1,0.0d0,TempPsiB,1)   ! Calculate the vector S*rSumPsi
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,mPsi(1,j),1,0.0d0,TempmPsi,1) ! Calculate the vector S*mPsi(1,j)

     do i = 1,NumStates

        !            testorth=ddot(MatrixDim,mPsi(1,i),1,TempmPsi,1)
        !            write(309,*) i,j, '   testorth=',testorth

        P(i,j) = aP*ddot(MatrixDim,mPsi(1,i),1,TempPsi,1)
        dP(i,j)= ddot(MatrixDim,mPsi(1,i),1,TempPsiB,1)

        do k = 1,MatrixDim
           lDiffPsi(k) = rPsi(k,i)-lPsi(k,i)
        enddo
        Q(i,j) = -aQ*ddot(MatrixDim,lDiffPsi,1,TempPsi,1)
     enddo
  enddo

  do j=1,NumStates
     do i=j,NumStates
        dP(i,j)=2.d0*aQ*(dP(i,j)-dP(j,i))
        dP(j,i)=-dP(i,j)
     enddo
  enddo

  deallocate(lDiffPsi,rDiffPsi,TempPsi,rSumPsi,TempPsiB,TempmPsi)

  return
end subroutine CalcCoupling
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

double precision function VSech(rij,DD,L)

  double precision rij,DD,L
  VSech = -DD/dcosh(rij/L)**2.d0
end function VSech

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function phirecon(R,beta,evec,left,right,RDim,MatrixDim,RNumPoints,RPoints,order)
  implicit none
  double precision, external :: BasisPhi
  integer MatrixDim,RDim,nch,beta,i,RNumPoints,left,right,order
  double precision R,evec(MatrixDim,MatrixDim),RPoints(RNumPoints)
  phirecon = 0.0d0
  do i = 1,RDim
     phirecon = phirecon + evec(i,beta)*BasisPhi(R,left,right,order,RDim,RPoints,&
          RNumPoints,0,i)
  enddo
  return
end function phirecon
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MyDsband(select,d,z,ldz,sigma,n,ab,mb,lda,rfac,ldrfac,k,nev,&
     tol,resid,ncv,v,ldv,iparam,workd,workl,lworkl,iwork,info)

  character        which*2, bmat, howmny
  integer          n, lda, ldrfac, k, nev, ncv, ldv, ldz, lworkl, info  
  Double precision tol, sigma
  logical          rvec

  integer          iparam(*), iwork(*)
  logical          select(*)
  Double precision d(*), resid(*), v(ldv,*), z(ldz,*), ab(lda,n), mb(lda,n), rfac(ldrfac,n), workd(*), workl(*) ! 

  integer          ipntr(14)

  integer          ido, i, j, Row, Col, type, ierr

  Double precision one, zero
  parameter        (one = 1.0, zero = 0.0)

  Double precision ddot, dnrm2, dlapy2
  external         ddot, dcopy, dgbmv, dgbtrf, dgbtrs, dnrm2, dlapy2, dlacpy ! 

  ! iparam(3) : Max number of Arnoldi iterations
  iparam(3) = 100000
  iparam(7) = 3
  rvec = .TRUE.
  howmny = 'A'
  which = 'LM'
  bmat = 'G'
  type = 4 
  ido = 0
  iparam(1) = 1

  rfac = 0.0d0
  do i = 1,n
     do j = i,min(i+k,n)
        Row = k+1+i-j
        Col = j
        rfac(k+Row,Col) = ab(Row,Col) - sigma*mb(Row,Col)
     enddo
     do j = max(1,i-k),i-1
        Row = 2*k+1
        Col = j
        rfac(Row+i-j,j) = rfac(Row+j-i,i)
     enddo
  enddo

  call dgbtrf(n,n,k,k,rfac,ldrfac,iwork,ierr)
  if ( ierr .ne. 0 )  then
     print*, ' '
     print*, '_SBAND: Error with _gbtrf:',ierr
     print*, ' '
     go to 9000
  end if

90 continue 

  call dsaupd(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)

  if (ido .eq. -1) then
     call dsbmv('U',n,k,1.0d0,mb,lda,workd(ipntr(1)),1,0.0d0,workd(ipntr(2)),1)
     call dgbtrs('Notranspose',n,k,k,1,rfac,ldrfac,iwork,workd(ipntr(2)),n,ierr)
     if (ierr .ne. 0) then
        print*, ' ' 
        print*, '_SBAND: Error with _gbtrs.'
        print*, ' ' 
        go to 9000
     end if
  else if (ido .eq. 1) then
     call dcopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
     call dgbtrs('Notranspose',n,k,k,1,rfac,ldrfac,iwork,workd(ipntr(2)),n,ierr)
     if (ierr .ne. 0) then 
        print*, ' '
        print*, '_SBAND: Error with _gbtrs.' 
        print*, ' '
        go to 9000
     end if
  else if (ido .eq. 2) then
     call dsbmv('U',n,k,1.0d0,mb,lda,workd(ipntr(1)),1,0.0d0,workd(ipntr(2)),1)
  else 
     if ( info .lt. 0) then
        print *, ' '
        print *, ' Error with _saupd info = ',info
        print *, ' Check the documentation of _saupd '
        print *, ' '
        go to 9000
     else 
        if ( info .eq. 1) then
           print *, ' '
           print *, ' Maximum number of iterations reached.'
           print *, ' '
        else if ( info .eq. 3) then
           print *, ' '
           print *, ' No shifts could be applied during '
           print *, ' implicit Arnoldi update, try increasing NCV.'
           print *, ' '
        end if
        if (iparam(5) .gt. 0) then
           call dseupd(rvec,'A',select,d,z,ldz,sigma,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
           if ( info .ne. 0) then
              print *, ' ' 
              print *, ' Error with _neupd = ', info
              print *, ' Check the documentation of _neupd '
              print *, ' ' 
              go to 9000
           endif
        endif
     endif
     go to 9000
  endif

  go to 90 

9000 continue

end subroutine MyDsband
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CalcEigenErrors(info,iparam,MatrixDim,H,LeadDim,S,HalfBandWidth,NumStates,Psi,Energies,MaxNumStates)

  integer info,iparam(*),MatrixDim,LeadDim,HalfBandWidth,NumStates,MaxNumStates
  double precision H(LeadDim,MatrixDim),S(HalfBandWidth+1,MatrixDim)
  double precision Psi(MatrixDim,MaxNumStates),Energies(MaxNumStates,2)

  integer j
  double precision dnrm2
  double precision, allocatable :: HPsi(:),SPsi(:)

  if ( info .eq. 0) then

     if (iparam(5) .lt. NumStates) write(6,*) 'Not all states found'

     ! Compute the residual norm: ||  A*x - lambda*x ||

     allocate(HPsi(MatrixDim))
     allocate(SPsi(MatrixDim))
     do j = 1,NumStates
        call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,H,LeadDim,Psi(1,j),1,0.0d0,HPsi,1)
        call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,Psi(1,j),1,0.0d0,SPsi,1)
        call daxpy(MatrixDim,-Energies(j,1),SPsi,1,HPsi,1)
        Energies(j,2) = dnrm2(MatrixDim,HPsi,1)
        Energies(j,2) = Energies(j,2)/dabs(Energies(j,1))
     enddo
     deallocate(HPsi)
     deallocate(SPsi)
  else
     write(6,*) ' '
     write(6,*) ' Error with _sband, info= ', info
     write(6,*) ' Check the documentation of _sband '
     write(6,*) ' '
  end if

  return
end subroutine CalcEigenErrors

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$c$$$      
!!$c$$$      double precision function VAIHO(rai,xi,C4)
!!$c$$$      implicit none
!!$c$$$      double precision rai, xi, C4, omegaho
!!$c$$$
!!$c$$$      VAIHO = -C4/(rai**4) + 0.5d0*xi**2
!!$c$$$      return
!!$c$$$      end
!!$c$$$
!!$c$$$      subroutine printVAIHO(phiai,mu,mi,ma)
!!$c$$$      implicit none
!!$c$$$      double precision, external :: VAIHO
!!$c$$$      integer i,j,N
!!$c$$$      double precision, allocatable :: x(:),y(:)
!!$c$$$      double precision rai,xi
!!$c$$$
!!$c$$$      N=200
!!$c$$$      allocate(rai(N),xi(N))
!!$c$$$      do i = 1,N
!!$c$$$
!!$c$$$      enddo
!!$c$$$      end
!!$         
