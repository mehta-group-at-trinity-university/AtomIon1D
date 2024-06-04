module BasisSets
  implicit none
  TYPE basis
     integer Left, Right, xDim
     double precision, allocatable :: x(:,:),V(:,:),u(:,:,:), uxx(:,:,:)
     double precision, allocatable :: S(:,:)
     double precision, allocatable :: H(:,:)
     double precision, allocatable :: D(:,:)
     double precision, allocatable :: Psi(:,:),Energies(:,:)
     double precision alpha, beta
     integer, allocatable :: xBounds(:)
  END TYPE basis
contains
  subroutine AllocateBasis(T,Left, Right,Order, LegPoints, xNumPoints,NumStates)
    implicit none
    integer Left, Right, LegPoints, xNumPoints, Order, MatrixDim, NumStates, ncv
    TYPE(basis) T
    T%Left = Left
    T%Right = Right
    T%xDim = xNumPoints+Order-3
    if (Left .eq. 2) T%xDim = T%xDim + 1
    if (Right .eq. 2) T%xDim = T%xDim + 1
    MatrixDim = T%xDim
    ncv = 2*NumStates
    allocate(T%x(LegPoints, xNumPoints))
    allocate(T%V(LegPoints, xNumPoints))
    allocate(T%u(LegPoints, xNumPoints, T%xDim))
    allocate(T%uxx(LegPoints,xNumPoints,T%xDim))
    allocate(T%xBounds(xNumPoints+2*Order))
    allocate(T%S(Order+1,T%xDim),T%H(Order+1,T%xDim),T%D(Order+1,T%xDim))
    allocate(T%Psi(MatrixDim,ncv))
    allocate(T%Energies(ncv,2))

  end subroutine AllocateBasis
  subroutine deAllocateBasis(T)
    implicit none
    TYPE(basis) T
    deallocate(T%x)
    deallocate(T%V)
    deallocate(T%u)
    deallocate(T%uxx)
    deallocate(T%xBounds)
    deallocate(T%S)
    deallocate(T%H)
    deallocate(T%D)
    deallocate(T%Psi)
    deallocate(T%Energies)

  end subroutine deAllocateBasis
end module BasisSets

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program AtomIon1D
  use BasisSets
  implicit none
  TYPE(basis) PB
  TYPE(basis) CB
  TYPE(basis) LB
  TYPE(basis) RB
  TYPE(basis) CB_old
  double precision tmp, tmp1, tmp2, S_ij
  integer max_state, LegPoints,xNumPoints, newRow, max_state1, max_state2, last_overlaps(2), curr_overlaps(2)
  integer mu1, nu1, pqflag, double, lastdouble1, lastdouble2
  integer NumStates,PsiFlag,Order,Left,Right, RSteps,CouplingFlag, CrossingFlag
  double precision alpha,tmp_alpha,tmp_beta,Shift,Shift2,NumStateInc,mi,ma,theta_c,mgamma
  double precision RLeft,RRight,RDerivDelt,DD,L,RFirst,RLast,XFirst,XLast,StepX,StepR, xMin,xMax
  double precision, allocatable :: R(:),Perm(:,:),ComposedPerm(:,:),AbsPerm(:,:)
  double precision, allocatable :: xPoints(:)
  double precision, allocatable :: AllEnergies(:,:), AllPsis(:,:,:)

  REAL t1, t2

  logical, allocatable :: LSelect(:)

  integer swapstate,iparam(11),ncv,info,troubleshooting
  integer i,j,k,iR,NumFirst,NumBound
  integer LeadDim,MatrixDim,HalfBandWidth
  integer xDim, Nbs, NumOutputChannels, OPGRID
  integer, allocatable :: iwork(:), WriteChannels(:)
  double precision Tol,RChange
  double precision TotalMemory
  double precision mu, mu12, mu123, r0diatom, dDiatom, etaOVERpi, Pi
  double precision YVal_L,RVal_L, RVal_R

  double precision, allocatable :: oldPsi(:,:), oldEnergies(:,:)
  double precision, allocatable :: LUFac(:,:),workl(:)
  double precision, allocatable :: workd(:),Residuals(:)
  double precision, allocatable :: xLeg(:),wLeg(:)
  double precision, allocatable :: NewPsi(:,:),TempEnergies(:,:)
  double precision, allocatable :: P(:,:,:),QTil(:,:,:),dP(:,:), testmat(:,:),SPsi(:)
  double precision ur(1:50000),acoef,bcoef,diff
  double precision sec,time,Rinitial,secp,timep,Rvalue, sNb, sbc, C4,lho
  double precision hbar, phi, amu,omega,Rstar, dum,Dtol,testorth
  double precision, external :: ddot, kdelta

  character(len=20), external :: str

  character*64 LegendreFile
  character*64 leftnrg
  common /Rvalue/ Rvalue      
  hbar = 1.054571817d-34
  amu = 1.660539d-27
  !     read in number of energies and states to print
  read(5,*)
  read(5,*) NumStates,PsiFlag,CouplingFlag, troubleshooting
  write(6,*) NumStates,PsiFlag,CouplingFlag, troubleshooting

  !     read in Gauss-Legendre info
  read(5,*)
  read(5,*)
  read(5,1002) LegendreFile
  write(6,1002) LegendreFile
  read(5,*)
  read(5,*)
  read(5,*) LegPoints
  write(6,*) LegPoints,' LegPoints'

  ! read in boundary conditions
  read(5,*)
  read(5,*)
  read(5,*) Shift,Shift2,Order,Left,Right
  print*, 'Shift,Shift2, Order, Left, Right'
  print*, Shift,Shift2,Order,Left,Right

  ! read in potential parameters
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
  theta_c = datan(mgamma)     ! same as Seth's beta !!! arctan mu, d?atan
  write(6,*) "theta_c = ", theta_c

  ! read in grid information
  read(5,*)
  read(5,*)
  read(5,*) xNumPoints, omega, Nbs, C4, phi,Dtol
  write(6,*) xNumPoints, omega, Nbs, C4, phi,Dtol

  omega = 2d0*Pi*omega
  lho = dsqrt(hbar/mi/omega)
  Rstar = dsqrt(2*mu12*C4/hbar**2)
  C4 = C4/(hbar*omega*lho**4)
  sNb = 1.d0/(dble(Nbs+1)*Pi + phi) !We use Nbs+1 here since the scattering length (and therefore phi=1/a) is negative. 

  write(6,*) "lho = ", lho
  write(6,*) "C4 = ", C4
  write(6,*) "Rstar = ", Rstar
  write(6,*) "sNb = ", sNb

  ! Re-define in oscillator units
  ma = ma/mi
  mi = 1d0
  mu12=mi*ma/(mi+ma)  ! ion-atom reduced mass (mu remains unchanged)

  read(5,*)
  read(5,*)
  read(5,*) RSteps,RDerivDelt,RFirst,RLast, OPGRID
  write(6,*) RSteps,RDerivDelt,RFirst,RLast, OPGRID
  read(5,*)
  read(5,*)
  read(5,*) NumOutputChannels
  write(6,*) "NumOutputChannels = ", NumOutputChannels
  allocate(WriteChannels(NumOutputChannels))
  read(5,*)
  read(5,*)
  read(5,*) WriteChannels!(WriteChannels(i), i = 1, NumOutputChannels)


  allocate(R(RSteps))
  StepR = (RLast-RFirst)/(dble(RSteps-1))
  do i = 1,RSteps
     !     read(5,*) R(i)
     !     R(i)= (XFirst+(i-1)*StepX)**3 ! Cubic Grid
     R(i) = (RFirst+(i-1)*StepR) !Linear Grid
     !    R(i) = 10.d0**(XFirst+(i-1)*StepX) ! Log Grid
  enddo

  allocate(xLeg(LegPoints),wLeg(LegPoints))
  call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  xDim = xNumPoints+Order-3

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
  !  allocate(xBounds(xNumPoints+2*Order))
  allocate(P(NumStates,NumStates,RSteps),QTil(NumStates,NumStates,RSteps))

  ncv = 2*NumStates
  LeadDim = 3*HalfBandWidth+1
  allocate(OldPsi(MatrixDim,ncv),NewPsi(MatrixDim,ncv),SPsi(MatrixDim))
  !  allocate(TransEnergies(2,ncv),TransPsi(ncv,MatrixDim))
  allocate(OldEnergies(ncv,2),TempEnergies(ncv,2))
  allocate(AllEnergies(NumStates,RSteps))
  allocate(AllPsis(NumStates,MatrixDim,RSteps))  ! Note that this is storing eigenvectors in the ROWS instead of collumns.
  allocate(testmat(NumStates,NumStates))
  allocate(iwork(MatrixDim))
  allocate(LSelect(ncv))
  allocate(LUFac(LeadDim,MatrixDim))
  allocate(workl(ncv*ncv+8*ncv))
  allocate(workd(3*MatrixDim))
  allocate(Residuals(MatrixDim))

  info = 0
  iR=1
  Tol=1e-20

  NumBound=0

  call AllocateBasis(PB,2,2,Order, LegPoints, xNumPoints, NumStates)
  call AllocateBasis(CB,Left,Right,Order, LegPoints, xNumPoints, NumStates)
!!$
!!$  open(unit = 203, file = "AdiabaticEnergies-LR="//trim(str(Left))//"-"//trim(str(Right))//".dat") 
!!$  open(unit = 200, file = "DiabaticEnergies-LR="//trim(str(Left))//"-"//trim(str(Right))//".dat") 
!!$  open(unit = 101, file = "Pmat-LR="//trim(str(Left))//"-"//trim(str(Right))//".dat")
!!$  open(unit = 102, file = "QTilmat-LR="//trim(str(Left))//"-"//trim(str(Right))//".dat")
!!$  open(unit = 103, file = "QuickPMat-LR="//trim(str(Left))//"-"//trim(str(Right))//".dat")
!!$  open(unit = 104, file = "QuickEffV-LR="//trim(str(Left))//"-"//trim(str(Right))//".dat")
!!$  open(unit = 105, file = "VQmat-LR="//trim(str(Left))//"-"//trim(str(Right))//".dat")

  open(unit = 203, file = "AdiabaticEnergies.dat") 
  open(unit = 200, file = "DiabaticEnergies.dat") 
  open(unit = 101, file = "Pmat.dat")
  open(unit = 102, file = "QTilmat.dat")
  open(unit = 103, file = "QuickPMat.dat")
  open(unit = 104, file = "QuickEffV.dat")
  open(unit = 105, file = "VQmat.dat")

  allocate(Perm(NumStates,NumStates),ComposedPerm(NumStates,NumStates),AbsPerm(NumStates,NumStates))

  CB%Left = Left
  CB%Right = Right
  RChange=100.d0
  write(200,*) "#", RSteps, NumOutputChannels
  do iR = RSteps,1,-1 
     call CPU_TIME(t1)
     NumFirst=NumStates
     if (R(iR).gt.RChange) then
        NumFirst=NumBound
     endif
     NumStateInc=NumStates-NumFirst

     xMin = 0d0   !theta_c + asin(dsqrt(mu/(1d0+mu**2))* sNb/R(iR)*(Rstar/lho)) !Based on fixed sNb
     xMax = Pi/2d0  !Pi + theta_c - asin(dsqrt(mu/(1d0+mu**2))* sNb/R(iR)*(Rstar/lho))

     sbc = (lho/Rstar)*R(iR)*dsqrt((1d0 + mu**2)/mu)*SIN(xMin - theta_c)
     !print*, 'sbc using sbc = (lho/Rstar)*R(iR)*dsqrt((1+mu**2)/mu)*SIN(xMin - theta_c):', sbc
     if(xMax.le.xMin) then
        write(6,*) "minimum hyperradius too small."
        stop
     endif
     call GridMakerIA(mu,mu12,theta_c,Rstar/lho,R(iR),sbc,xNumPoints,xMin,xMax,xPoints,OPGRID)
     !******** CONSTRUCTION OF BASIS SETS ********
     call CalcBasisFuncs(CB%Left,CB%Right,Order,xPoints,LegPoints,xLeg,CB%xDim,CB%xBounds,xNumPoints,0,CB%u)
     call CalcBasisFuncs(CB%Left,CB%Right,Order,xPoints,LegPoints,xLeg,CB%xDim,CB%xBounds,xNumPoints,2,CB%uxx)
     call CalcHSD(alpha,R(iR),mu,mi,theta_c,C4,L,Order,xPoints,&
          LegPoints,xLeg,wLeg,CB%xDim,xNumPoints,CB%u,CB%uxx,CB%xBounds,HalfBandWidth,CB%H,CB%S,CB%D,CB)
     call MyDsband(LSelect,CB%Energies,CB%Psi,MatrixDim,Shift,MatrixDim,CB%H,CB%S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,&
          NumStates,Tol,Residuals,ncv,CB%Psi,MatrixDim,iparam,workd,workl,ncv*ncv+8*ncv,iwork,info)
     call CalcEigenErrors(info,iparam,MatrixDim,CB%H,HalfBandWidth+1,CB%S,HalfBandWidth,NumStates,CB%Psi,CB%Energies,ncv)
     NewPsi = CB%Psi

     AllPsis(:,:,iR) = Transpose(CB%Psi(:,1:NumStates))
     AllEnergies(:,iR) = CB%Energies(1:NumStates,1)

     call CalcPermutation(NumStates,HalfBandWidth,MatrixDim,OldPsi,NewPsi,CB%S,Perm,ncv,Dtol)
     if(iR.eq.RSteps) then
        Perm = 0d0
        do i = 1,NumStates
           Perm(i,i) = 1d0
        enddo
        ComposedPerm = Perm
     endif

     ! Construct the composed permutation operator and permute the eigenvectors and eigenvalues
     ComposedPerm = matmul(ComposedPerm,Perm)
     AllPsis(:,:,iR) = matmul(ComposedPerm,AllPsis(:,:,iR))
     AllEnergies(:,iR) = matmul(ComposedPerm,AllEnergies(:,iR))

     NewPsi = 0d0
     NewPsi(:,1:NumStates) = Transpose(AllPsis(:,:,iR))

     !--- Check the phase consistency of NewPsi, and fix phases ------
     if(iR.lt.RSteps) then
        OldPsi = 0d0
        OldPsi(:,1:NumStates) = Transpose(AllPsis(:,:,iR+1))
        do j = 1,NumStates
           call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,CB%S,HalfBandWidth+1,OldPsi(:,j),1,0.0d0,SPsi,1)  ! Calculate the vector S*OldPsi and store in SPsi            
           testorth=ddot(MatrixDim,NewPsi(:,j),1,SPsi,1)            
           if(testorth.lt.0d0) then
              NewPsi(:,j) = -NewPsi(:,j)
              AllPsis(j,:,iR) = -AllPsis(j,:,iR)
           endif
        enddo

     endif
     !--------------------------------------------

     write(203, 20) R(iR), (CB%Energies(i,1), i = 1,NumStates) ! writing undiabatized energies

     !-------- Couplings from the diabatized eigenstates -------------!
     call CalcCoupling_FH(NumStates,HalfBandWidth,MatrixDim,AllEnergies(:,iR),NewPsi,CB%S,CB%D,P(:,:,iR),QTil(:,:,iR),ncv)      

     write(103,*) R(iR),(P(WriteChannels(1),WriteChannels(i),iR),i=1,NumOutputChannels)
     write(104,*) R(iR),(AllEnergies(writechannels(i),iR) - 0.25d0/(2d0*mu*R(iR)**2) &
          + QTil(WriteChannels(i),WriteChannels(i),iR)/(2d0*mu),i=1,NumOutputChannels)

     !write(104,*) R(iR), (QTil(i,i,iR)*R(iR)**2, i=1,1)
     !write(104,*) R(iR), ((AllEnergies(i,iR)-(i-1d0+0.5d0))*2*mu*R(iR)**2, i=1,NumStates)
     !--------------------------------------------------------------------------!
     OldPsi = CB%Psi
     !    Adjusting Shift
     ur(iR) = CB%Energies(1,1)
     Shift = -200d0
     if(ur(iR).lt.0d0) then
        Shift = ur(iR)*10d0
     endif

     call CPU_TIME(t2)
     write(6,*) 'Remaining time (min): ', (t2-t1)*iR/60, R(iR), (CB%Energies(i,1), i=1,5)

  enddo

  do iR = 1, RSteps
     write(200,11) R(iR), (AllEnergies(WriteChannels(i),iR), i = 1,NumOutPutChannels) ! Write the diabatized energies      
     write(101,11) R(iR)
     write(102,11) R(iR)
     write(105,11) R(iR)
     do i=1,NumOutputChannels
        write(101,11) (P(WriteChannels(i),WriteChannels(j),iR), j=1,NumOutputChannels) !Write P to Pmat.dat
        write(102,11) (QTil(WriteChannels(i),WriteChannels(j),iR), j=1,NumOutputChannels) !Write Qtil to QTilmat.dat
        write(105,11) ((AllEnergies(writechannels(i),iR) - 0.25d0/(2d0*mu*R(iR)**2))*kdelta(i,j) &
             + QTil(WriteChannels(i),WriteChannels(j),iR)/(2d0*mu), j=1,NumOutputChannels) !Write VQmat to VQmat.dat
     enddo
  enddo
  
  close(11)
  close(101)
  close(102)
  close(103)
  close(104)
  close(105)
  close(200)
  close(300)
  close(301)
  !100 close(203)

  deallocate(iwork)
  deallocate(LSelect)
  deallocate(LUFac)
  deallocate(workl)
  deallocate(workd)
  deallocate(Residuals)
  deallocate(P,QTil)
  deallocate(xPoints)
  deallocate(xLeg,wLeg)

  call deAllocateBasis(CB)


  deallocate(R)

10 format(1P,100e25.15)
11 format(1P,100e22.12)
20 format(1P,100e16.8)
1002 format(a64)

  stop
end program AtomIon1D
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine CalcPermutation(NumStates,HalfBandWidth,MatrixDim,OldPsi,NewPsi,S,Perm,ncv,tol)
  implicit none
  integer NumStates,HalfBandWidth,MatrixDim,ncv
  double precision OldPsi(MatrixDim,ncv),NewPsi(MatrixDim,ncv)
  double precision S(HalfBandWidth+1,MatrixDim),testorth
  double precision Perm(NumStates,NumStates)
  double precision test,tol
  integer i,j,k,n,m
  double precision, external :: ddot 
  double precision, allocatable :: TempPsi(:,:),SPsi(:)

  allocate(TempPsi(MatrixDim,ncv),SPsi(MatrixDim))

  TempPsi = OldPsi
  Perm = 0d0
  do m = 1,NumStates
     Perm(m,m) = 1d0
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi(:,m),1,0.0d0,SPsi,1)  ! Calculate the vector S*OldPsi(m) and store in SPsi 
     do n = m,NumStates
        testorth=ddot(MatrixDim,NewPsi(:,n),1,SPsi,1)
        !write(6,*) n,m, '   testorth=',testorth
        if (abs(abs(testorth)-1d0).lt.tol) then
           Perm(m,n) = 1d0
           !           write(6,*) m,n
        else if(abs(testorth).lt.tol) then           
           Perm(m,n) = 0d0
        else
           write(6,*) "testorth outside of tolerance",testorth, "but Perm(", m,n,") set to ",Perm(m,n)
        endif
        Perm(n,m) = Perm(m,n)
     enddo
  enddo

  if((Perm(NumStates,NumStates).eq.0d0).and.(Perm(NumStates-1,NumStates-1).eq.1d0)) then
     Perm(NumStates,NumStates)=1d0
  endif

  deallocate(TempPsi,SPsi)

end subroutine CalcPermutation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Computes the couplings using the Feynman-Hellmann theorem
! Requires the matrix D = <B_i|dH/dR|B_j>, the energies Un and the eigenvector matrix Psi.  
subroutine CalcCoupling_FH(NumStates,HalfBandWidth,MatrixDim,U,Psi,S,D,P,QTil,ncv)
  implicit none
  integer NumStates,HalfBandWidth,MatrixDim,ncv
  double precision Psi(MatrixDim,ncv),U(NumStates,2)
  double precision S(HalfBandWidth+1,MatrixDim),D(HalfBandWidth+1,MatrixDim),testorth
  double precision P(NumStates,NumStates),QTil(NumStates,NumStates)
  double precision aP
  integer i,j,k,n,m
  double precision, external :: ddot 
  double precision, allocatable :: TempPsi(:,:),SPsi(:), DPsi(:)

  allocate(TempPsi(MatrixDim,ncv),SPsi(MatrixDim),DPsi(MatrixDim))

  TempPsi = Psi
  P=0d0
  do m = 1,NumStates
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,D,HalfBandWidth+1,TempPsi(:,m),1,0.0d0,DPsi,1)   ! Calculate the vector D*Psi(m) and store in DPsi
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi(:,m),1,0.0d0,SPsi,1)  ! Calculate the vector S*Psi(m) and store in SPsi 
     do n = m+1,NumStates
        testorth=ddot(MatrixDim,TempPsi(:,n),1,SPsi,1)
        !write(6,*) n,m, '   testorth=',testorth
        aP = 1d0/(U(m,1) - U(n,1))
        P(n,m) = aP*ddot(MatrixDim,TempPsi(:,n),1,DPsi,1)
        P(m,n) = -P(n,m)
     enddo
  enddo
  
  QTil = -matmul(P,P)
  !!  Try doubling to check if the delicate interference goes away
  !P=2d0*P  !! Be sure to remove this after the check is complete!
  !!

  
!!$      do i = 1, NumStates
!!$         do j = 1, NumStates
!!$            write(6,*) i,j,QTil(i,j)
!!$         enddo
!!$      enddo
  deallocate(TempPsi,SPsi,DPsi)

  return
end subroutine CalcCoupling_FH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
        rSumPsi(k)  = lPsi(k,j)+mPsi(k,j)+rPsi(k,j)
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

!!!!!!!!!!!!!!!!!!!!!! (LB, CB, RB, P, Q, PB%S, HalfBandWidth, NumStates, RDerivDelt)

subroutine calc_Physical_Set(b, pb,RVal_L,RVal_R,Order, xNumPoints, xPoints)
  use BasisSets
  implicit none
  double precision, external :: calc_alpha, calc_beta
  integer Order, xNumPoints
  double precision RVal_L, RVal_R, xPoints(xNumPoints)
  TYPE(basis) b, pb
  !!REVISIT FOR cases other than L,R = 3,3

  b%alpha = calc_alpha(b, RVal_L, Order, xNumPoints, xPoints)
  b%beta = calc_beta(b, RVal_R, Order, xNumPoints, xPoints)

  b%u(:,:,1) = b%alpha*pb%u(:,:,1)+pb%u(:,:,2)
  b%uxx(:,:,1) = b%alpha*pb%uxx(:,:,1)+pb%uxx(:,:,2)

  b%u(:,:,b%xDim) = pb%u(:,:,pb%xDim-1)+b%beta*pb%u(:,:,pb%xDim)
  b%uxx(:,:,b%xDim) = pb%uxx(:,:,pb%xDim-1)+b%beta*pb%uxx(:,:,pb%xDim)

  b%u(:,:,2:b%xDim-1) = pb%u(:,:,3:pb%xDim-2)
  b%uxx(:,:,2:b%xDim-1) = pb%uxx(:,:,3:pb%xDim-2)

end subroutine calc_Physical_Set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function calc_alpha(B, RVal_L, Order, xNumPoints, xPoints)
  use BasisSets
  implicit none
  TYPE(basis) B
  integer Order, xNumPoints
  double precision RVal_L,xPoints(xNumPoints)
  double precision, external :: MYBSpline

  calc_alpha = (MYBSpline(Order,1,xNumPoints,xPoints,2,xPoints(1))*RVal_L&
       - MYBSpline(Order,0,xNumPoints,xPoints,2,xPoints(1)))&
       /(MYBSpline(Order,0,xNumPoints,xPoints,1,xPoints(1))&
       -MYBSpline(Order,1,xNumPoints,xPoints,1,xPoints(1))*RVal_L)

end function calc_alpha
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function calc_beta(B, RVal_R, Order, xNumPoints, xPoints)
  use BasisSets
  implicit none
  TYPE(basis) B
  integer Order, xNumPoints
  double precision RVal_R, xPoints(xNumPoints)
  double precision, external :: MYBSpline

  calc_beta = (MYBSpline(Order,1,xNumPoints,xPoints,xNumPoints+Order-2,xPoints(xNumPoints))*RVal_R&
       -MYBSpline(Order,0,xNumPoints,xPoints,xNumPoints+Order-2,xPoints(xNumPoints)))&
       /(MYBSpline(Order,0,xNumPoints,xPoints,xNumPoints+Order-1,xPoints(xNumPoints))&
       - MYBSpline(Order,1,xNumPoints,xPoints,xNumPoints+Order-1,xPoints(xNumPoints))*RVal_R)

end function calc_beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcCouplings_v2(LB, CB, RB, P, Q, S_prim, HalfBandWidth, NumStates, RDerivDelt)
  use BasisSets
  implicit none
  TYPE(basis) LB, CB, RB
  double precision P(NumStates, NumStates), Q(NumStates, NumStates), S_prim(HalfBandWidth+1, CB%xDim+2), RDerivDelt
  integer HalfBandWidth, NumStates
  integer ct, mu, nu, N_
  double precision P_term_1, Q_term_1, Q_term_2
  N_ = CB%xDim+2
  !print*, 'Calculating P matrix...'
  do mu = 1, NumStates
     do nu = 1, NumStates
        P_term_1 = calc_P_term_1(mu,nu, CB, RB, S_prim, HalfBandWidth)
        P(mu,nu) = 1d0/RDerivDelt * ((P_term_1)-kronecker_delta(mu,nu))
        if(dabs(P(mu,nu)).ge.0.01) then
           write(155, *) 'mu =', mu, 'nu = ', nu
        endif
     enddo
  enddo
  ct = 0
  !print*, 'Done Calculating P Matrix!'
  !print*, 'Calculating Q matrix...'
  do mu = 1, NumStates
     do nu = 1, NumStates
        Q_term_1 = calc_P_term_1(mu,nu, CB, RB, S_prim, HalfBandWidth)
        Q_term_2 = calc_P_term_1(mu,nu, CB, LB, S_prim, HalfBandWidth)
        Q(mu,nu) = 1d0/RDerivDelt**2 * (Q_term_1+Q_term_2-2*kronecker_delta(mu,nu))
     enddo
  enddo
  ct = 0
  !print*, 'Done Calculating Q Matrix!'
contains
  double precision function kronecker_delta(m,n)
    implicit none
    integer m,n
    if (m .eq. n) then
       kronecker_delta = 1
    else
       kronecker_delta = 0
    endif
  end function kronecker_delta

  double precision function calc_P_term_1(mu,nu, CB, RB, S_prim, HalfBandWidth)
    implicit none
    integer mu, nu, HalfBandWidth, i, j
    TYPE(BASIS) CB, RB
    double precision S_ij, S_prim(HalfBandWidth+1, CB%xDim+2), term_1_mu_nu
    term_1_mu_nu = 0d0
    !print*, 'LOOKING FOR NEGATIVE VALUES OF S_ij'
    do i = 1, CB%xDim
       do j = max(1,i-HalfBandWidth), min(CB%xDim,i+HalfBandWidth)
          call calc_overlap_elem(i,j, CB, RB, S_prim, HalfBandWidth, S_ij)
          term_1_mu_nu = term_1_mu_nu + CB%Psi(i,mu)*(S_ij*RB%Psi(j,nu))
       enddo
    enddo
    !print*, '1/RDerivDelt*term_1_mu_nu=', term_1_mu_nu*(1/.0001)
    !if (nu .eq. 9) stop
    calc_P_term_1 = term_1_mu_nu
  end function calc_P_term_1
end subroutine calcCouplings_v2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_overlap_elem(m,n,u, ur, S_prim, HalfBandWidth, S_mn)
  use BasisSets
  implicit none
  TYPE(basis) u, ur
  integer m,n, HalfBandWidth, row, column, tmp
  double precision N_, m_alpha, r_alpha, m_beta, r_beta, S_prim(HalfBandWidth+1,u%xDim+2),&
       S_mn, term_1, term_2, term_3, term_4, testing
  N_ = u%xDim
  m_alpha = u%alpha
  r_alpha = ur%alpha
  m_beta = u%beta
  r_beta = ur%beta

  if ((m .gt. 1).and.(m .lt. N_)) then
     if ((n .gt. 1).and.(n .lt. N_)) then
        row = m+1
        column = n+1
        term_1 = banded_zeros_check(row,column,HalfBandWidth,S_prim)
        S_mn = term_1
     else if (n .eq. 1) then !problem case
        row = 1
        column = m + 1
        term_1 = r_alpha*banded_zeros_check(row,column,HalfBandWidth,S_prim)

        row = 2
        column = m + 1
        term_2 = banded_zeros_check(row,column,HalfBandWidth,S_prim)

        S_mn = term_1 + term_2
     else if (n .eq. N_) then
        row = m + 1
        column = N_+1

        term_1 = banded_zeros_check(row,column,HalfBandWidth,S_prim)

        row = m + 1
        column = N_+2
        term_2 = r_beta*banded_zeros_check(row,column,HalfBandWidth,S_prim)

        S_mn = term_1 + term_2
     endif
  else if (m .eq. 1) then
     if ((n .gt. 1).and.(n .lt. N_)) then
        row = 1
        column = n+1
        term_1 = m_alpha*banded_zeros_check(row,column,HalfBandWidth,S_prim)

        row = 2
        column = n+1
        term_2 = banded_zeros_check(row,column,HalfBandWidth,S_prim)

        S_mn = term_1 + term_2
     else if (n .eq. 1) then
        row = 1
        column = 1
        term_1 = m_alpha*r_alpha*banded_zeros_check(row,column,HalfBandWidth,S_prim)

        row = 1
        column = 2
        term_2 = m_alpha*banded_zeros_check(row,column,HalfBandWidth,S_prim)

        row = 2
        column = 1
        term_3 = r_alpha*banded_zeros_check(row,column,HalfBandWidth,S_prim)

        row = 2
        column = 2
        term_4 = banded_zeros_check(row,column,HalfBandWidth,S_prim)

        S_mn = term_1+term_2+term_3+term_4
     else if (n .eq. N_) then
        row = 1
        column = N_+1
        term_1 = m_alpha*banded_zeros_check(row, column, HalfBandWidth, S_prim)

        row = 1
        column = N_+2
        term_2 = m_alpha*r_beta*banded_zeros_check(row, column, HalfBandWidth, S_prim)

        row = 2
        column = N_+1
        term_3 = banded_zeros_check(row, column, HalfBandWidth, S_prim)

        row = 2
        column = N_+2
        term_4 = r_beta*banded_zeros_check(row, column, HalfBandWidth, S_prim)

        S_mn = term_1+term_2+term_3+term_4
     endif
  else if(m .eq. N_) then
     if ((n .gt. 1).and.(n .lt. N_)) then
        row = N_+1
        column = n+1
        term_1 = banded_zeros_check(row,column,HalfBandWidth,S_prim)

        row = N_+2
        column = n+1
        term_2 = m_beta*banded_zeros_check(row,column,HalfBandWidth,S_prim)

        S_mn = term_1+term_2
     else if (n .eq. 1) then
        row = N_+1
        column = 1
        term_1 = r_alpha*banded_zeros_check(row, column, HalfBandWidth, S_prim)

        row = N_+2
        column = 2
        term_2 = banded_zeros_check(row, column, HalfBandWidth, S_prim)

        row = N_+2
        column = 1
        term_3 = r_alpha*m_beta*banded_zeros_check(row, column, HalfBandWidth, S_prim)

        row = N_+2
        column = 2
        term_4 = m_beta*banded_zeros_check(row, column, HalfBandWidth, S_prim)

        S_mn = term_1+term_2+term_3+term_4
     else if (n .eq. N_) then
        row = N_+1
        column = N_+1
        term_1 = banded_zeros_check(row,column,HalfBandWidth,S_prim)

        row = N_+1
        column = N_+2
        term_2 = banded_zeros_check(row,column,HalfBandWidth,S_prim)*(m_beta+r_beta)

        row = N_+2
        column = N_+2
        term_3 = banded_zeros_check(row,column,HalfBandWidth,S_prim)*m_beta*r_beta

        S_mn = term_1+term_2+term_3
     endif
  endif

contains
  double precision function banded_zeros_check(row,col,HalfBandWidth,S_prim)
    implicit none
    integer row, col, HalfBandWidth, newRow, tmp
    double precision S_prim(:,:)
    if((row+HalfBandWidth .lt. col).or.(row-HalfBandWidth .gt. col)) then
       banded_zeros_check = 0d0
    else
       if(row .gt. col) then
          tmp = row
          row = col
          col = tmp
       endif
       newRow = HalfBandWidth + 1 + (row - col)
       banded_zeros_check = S_prim(newRow, col)
    endif
  end function banded_zeros_check

end subroutine calc_overlap_elem

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
! This new subroutine computes the Hamiltonian, Overlap, and the P-matrix
! NOTE: variables such as H,S,D are output of this code, but note that
! many are ultimately part of type BasisSet because they are elements of CB.
! In the case of CB%x and CB%V the asignment is made directly through CB.
!
subroutine CalcHSD(alpha,R,mu,mi,theta_c,C4,L,Order,xPoints,LegPoints,&
     xLeg,wLeg,xDim,xNumPoints,u,uxx,xBounds,HalfBandWidth,H,S,D,CB)
  use BasisSets
  implicit none
  TYPE(basis) CB
  double precision, external :: VSech
  integer Order,LegPoints,xDim,xNumPoints,xBounds(*),HalfBandWidth
  double precision alpha,R,mu,mgamma,theta_c,C4,L,mi
  double precision xPoints(*),xLeg(*),wLeg(*)
  double precision H(HalfBandWidth+1,xDim),S(HalfBandWidth+1,xDim),D(HalfBandWidth+1,xDim)  ! D is a matrix with elements <B_i|dH/dR|B_j> needed for the P-matrix
  double precision u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim)
  integer ix,ixp,kx,lx
  integer i1,i1p
  integer Row,NewRow,Col
  integer, allocatable :: kxMin(:,:),kxMax(:,:)
  double precision a,b,m,Pi
  double precision Rall,rai,xai
  double precision u1,sys_ss_pot,V12,V23,V31
  double precision VInt,VTempInt,potvalue, xTempV,xTempS
  double precision x,ax,bx,xScaledZero,xTempT,xTempVHO,xTempVC4
  double precision, allocatable :: Pot(:,:)
  double precision, allocatable :: xIntScale(:),xT(:,:),xV(:,:),xVHO(:,:),xVC4(:,:)
  double precision, allocatable :: cosx(:,:),sinx(:,:),xS(:,:)
  double precision, allocatable :: XX(:,:),YY(:,:)
  double precision mu12,r0diatom,dDiatom


  mgamma=dtan(theta_c)
  allocate(xIntScale(xNumPoints),xT(xDim,xDim),xV(xDim,xDim))
  allocate(xVHO(xDim,xDim),xVC4(xDim,xDim))
  allocate(cosx(LegPoints,xNumPoints),sinx(LegPoints,xNumPoints))
  allocate(XX(LegPoints,xNumPoints),YY(LegPoints,xNumPoints))
  allocate(xS(xDim,xDim))
  allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
  allocate(Pot(LegPoints,xNumPoints))

  Pi = 3.1415926535897932385d0

  m = -1.0d0/(2.0d0*mu*R*R)
  !      theta_c = datan(mgamma)
  ! Could create an array x(LegPoints,xNumPoints) that specifies all of the integration abcissa.


  do kx = 1,xNumPoints-1
     ax = xPoints(kx)
     bx = xPoints(kx+1)
     xIntScale(kx) = 0.5d0*(bx-ax)
     xScaledZero = 0.5d0*(bx+ax)
     do lx = 1,LegPoints
        x = xIntScale(kx)*xLeg(lx)+xScaledZero
        CB%x(lx,kx) = x
        cosx(lx,kx) = dcos(x)
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
        !xai = mu**(-0.5d0)*sinx(lx,kx) - (mu**0.5d0)*cosx(lx,kx) 
        XX(lx,kx) = cosx(lx,kx)**2
        !YY(lx,kx) = C4/(xai**4)
        !potvalue = -YY(lx,kx)/R**4 + 0.5d0*mu*R*R*XX(lx,kx)
        potvalue = 0.5d0*mu*R*R*XX(lx,kx)
        Pot(lx,kx) = alpha*potvalue
        CB%V(lx,kx) = alpha*potvalue
        !                    write(6,*) 'THIS IS A TEST', kx, lx, Pot(lx,kx)
     enddo
  enddo
  xV = 0d0
  xT = 0d0
  xS = 0d0
  xVHO = 0d0
  xVC4 = 0d0
  do ix = 1,xDim
     do ixp = max(1,ix-Order),min(xDim,ix+Order)
        do kx = kxMin(ixp,ix),kxMax(ixp,ix)
           xTempT = 0.0d0
           xTempV = 0.0d0
           xTempVHO = 0d0
           xTempVC4 = 0d0
           xTempS = 0d0
           do lx = 1,LegPoints
              a = wLeg(lx)*xIntScale(kx)*u(lx,kx,ix) !bra
              xTempS = xTempS + a*u(lx,kx,ixp)
              xTempT = xTempT + a*uxx(lx,kx,ixp) !KE operator !ket
              xTempVHO = xTempVHO + a*XX(lx,kx)*u(lx,kx,ixp) ! Matrix elements of cos^2(theta)
              !xTempVC4 = xTempVC4 + a*YY(lx,kx)*u(lx,kx,ixp) ! Matrix elements of C4/xia^4 (just the angular part)
              !xTempV = xTempV + a*Pot(lx,kx)*u(lx,kx,ixp) !PE operator ket
           enddo
           xT(ix,ixp) = xT(ix,ixp) + xTempT
           xS(ix,ixp) = xS(ix,ixp) + xTempS
           !xV(ix,ixp) = xV(ix,ixp) + xTempV
           xVHO(ix,ixp) = xVHO(ix,ixp) + xTempVHO 
           !xVC4(ix,ixp) = xVC4(ix,ixp) + xTempVC4
        enddo
        xV(ix,ixp) = 0.5d0*mu*R*R*xVHO(ix,ixp) !- R**(-4)*xVC4(ix,ixp)
     enddo
  enddo

  H = 0.0d0
  S = 0d0
  D = 0d0
  do ix = 1,xDim
     Row=ix
     do ixp = max(1,ix-Order),min(xDim,ix+Order)
        Col = ixp
        if (Col .ge. Row) then
           NewRow = HalfBandWidth+1+Row-Col
           S(NewRow,Col) = xS(ix,ixp)
           H(NewRow,Col) = (m*xT(ix,ixp)+xV(ix,ixp))
           !D(NewRow,Col) = -2d0*m*xT(ix,ixp)/R + mu*R*xVHO(ix,ixp) + (4d0/R**5)*xVC4(ix,ixp)
           D(NewRow,Col) = xT(ix,ixp)/(mu*R**3) + mu*R*xVHO(ix,ixp) !+ (4d0/R**5)*xVC4(ix,ixp)
           !                write(6,*) 'THIS IS A TEST', ix,ixp,H(NewRow,Col) !!all info stored now in H
        endif
     enddo
  enddo

  deallocate(Pot)
  deallocate(xIntScale,xT,xV)
  deallocate(cosx,sinx)
  deallocate(XX,YY)
  deallocate(kxMin,kxMax)


  return
end subroutine CalcHSD

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CalcHamiltonian(alpha,R,mu,mi,theta_c,C4,L,Order,xPoints,LegPoints,&
     xLeg,wLeg,xDim,xNumPoints,u,uxx,xBounds,HalfBandWidth,H)
  implicit none
  double precision, external :: VSech
  integer Order,LegPoints,xDim,xNumPoints,xBounds(*),HalfBandWidth
  double precision alpha,R,mu,mgamma,theta_c,C4,L,mi
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


  mgamma=dtan(theta_c)
  allocate(xIntScale(xNumPoints),xT(xDim,xDim),xV(xDim,xDim))
  allocate(cosx(LegPoints,xNumPoints),sinx(LegPoints,xNumPoints))
  allocate(sinai(LegPoints,xNumPoints))
  allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
  allocate(Pot(LegPoints,xNumPoints))

  Pi = 3.1415926535897932385d0

  m = -1.0d0/(2.0d0*mu*R*R)
  !      theta_c = datan(mgamma)
  ! Could create an array x(LegPoints,xNumPoints) that specifies all of the integration abcissa.


  do kx = 1,xNumPoints-1
     ax = xPoints(kx)
     bx = xPoints(kx+1)
     xIntScale(kx) = 0.5d0*(bx-ax)
     xScaledZero = 0.5d0*(bx+ax)
     do lx = 1,LegPoints
        x = xIntScale(kx)*xLeg(lx)+xScaledZero
        cosx(lx,kx) = dcos(x)
        sinai(lx,kx) = dsin(theta_c-x)
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

!        rai = mu**(-0.5d0) * R*sinx(lx,kx) -  mu**0.5d0 * R*cosx(lx,kx) !dsqrt(mu/mu12)*R*dabs(sinai(lx,kx))))

        !       potvalue = -C4/rai**4 + 0.5d0*mu*(R*cosx(lx,kx))**2
        potvalue = 0.5d0*mu*(R*cosx(lx,kx))**2
        Pot(lx,kx) = alpha*potvalue
        !                    write(6,*) 'THIS IS A TEST', kx, lx, Pot(lx,kx)
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
           !                write(6,*) 'THIS IS A TEST', ix,ixp,H(NewRow,Col) !!all info stored now in H
        endif
     enddo
  enddo

  deallocate(Pot)
  deallocate(xIntScale,xT,xV)
  deallocate(cosx,sinx,sinai)
  deallocate(kxMin,kxMax)


  return
end subroutine CalcHamiltonian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        TempPsi1(k) = lPsi(k,j)+rPsi(k,j)-2.0d0*mPsi(k,j) ! j is 1, k = 0 through  matrixdim is set by lpsi+rpsi(1,k) - 2mpsi(1,k)
     enddo
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi1,1,0.0d0,TempPsi2,1) !matrix multiplication, s * temp1 = tmp 2
     do i = 1,NumStates
        Q(i,j) = a*ddot(MatrixDim,TempPsi2,1,mPsi(1,i),1)
     enddo
  enddo

  deallocate(TempPsi1,TempPsi2)

  return
end subroutine CalcQMatrix
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
subroutine FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,OldPsi,NewPsi)
  implicit none
  integer NumStates,HalfBandWidth,MatrixDim,ncv
  double precision S(HalfBandWidth+1,MatrixDim),Psi(MatrixDim,ncv)
  double precision OldPsi(MatrixDim,ncv),NewPsi(MatrixDim,ncv)

  integer i,j
  double precision Phase,ddot
  double precision, allocatable :: TempPsi(:)

  allocate(TempPsi(MatrixDim))

  do i = 1,NumStates
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,OldPsi(:,i),1,0.0d0,TempPsi,1) ! TempPsi stores S*OldPsi
     Phase = ddot(MatrixDim,NewPsi(:,i),1,TempPsi,1) ! Computes overlap of NewPsi with TempPsi
     if (Phase .lt. 0.0d0) then
        NewPsi(:,i) = -NewPsi(:,i)
        !        write(6,*) "i, overlap = ", i, Phase, ddot(MatrixDim,NewPsi(:,i),1,TempPsi,1)
     endif
  enddo

  deallocate(TempPsi)

  return
end subroutine FixPhase
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GridMakerHHL(mu,mu12,mu123,theta_c,R,r0,xNumPoints,xMin,xMax,xPoints)
  implicit none
  integer xNumPoints
  double precision mu,R,r0,xMin,xMax,xPoints(xNumPoints)
  double precision mu12,mu123,theta_c
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
     x1 = theta_c - deltax  
     x2 = theta_c + deltax
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
double precision function gridshapef(a,x)
  implicit none
  double precision a, x
  double precision, parameter :: pi = 3.1415926535897932385d0
  gridshapef = 4d0*pi*a*x - sin(4d0*pi*x)
end function gridshapef
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GridMakerIA(mu,mu12,theta_c,Rstar,R,sbc,xNumPoints,xMin,xMax,xPoints,OPGRID)
  implicit none
  integer xNumPoints, testing
  double precision mu,R,r0,xMin,xMax,xPoints(xNumPoints)
  double precision mu12,mu123,theta_c,Rstar,sbc
  integer i,j,k,OPGRID,NP
  double precision Pi
  double precision r0New
  double precision xRswitch,Rscale
  double precision xDelt,x0,a
  double precision, external :: gridshapef

  Pi = 3.1415926535897932385d0

  x0 = xMin
!  write(6,*) "making grid with"
!  write(6,*) "xMin = ", xMin
!  write(6,*) "xMax = ", xMax
  

  if((OPGRID.eq.1)) then
     Rscale = 2d0*Rstar/R
!     write(6,*) "using modified grid:(R,Rstar,Rscale) = ", R, Rstar, Rscale
     xDelt = (xMax-xMin)
     a = exp(max(Rscale, 0.75d0))
     do i = 1, xNumPoints
        xPoints(i) = x0 + xDelt/gridshapef(a,1d0) * (gridshapef(a,dble(i-1d0)/dble(xNumPoints-1d0)))
!        write(6,*) i, xPoints(i)
     enddo
  else
     xDelt = (xMax-xMin)/dfloat(xNumPoints-1)
     do i = 1,xNumPoints
        xPoints(i) = (i-1)*xDelt + x0
     enddo
  endif


!!$  !     Smooth Grid 
!!$
!!$  do j=1,10
!!$     do i = 2, xNumPoints-1
!!$        xPoints(i)=(xPoints(i-1)+xPoints(i)+xPoints(i+1))/3.d0
!!$     enddo
!!$  enddo
  !  do i = 1, xNumPoints
  !     write(20,*) i, xPoints(i)
  !  enddo
  !  write(20,*) ' '
  !      write(96,'(100e20.10)') R, (xPoints(k),k=1,xNumPoints)

15 format(6(1x,1pd12.5))

  !      stop

  return
end subroutine GridMakerIA
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GridMakerOld(mu,R,r0,xNumPoints,xMin,xMax,xPoints)
  implicit none
  integer xNumPoints
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
end subroutine GridMakerOld
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
!!$c$$$      subroutine printVAIHO(theta_c,mu,mi,ma)
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

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str
