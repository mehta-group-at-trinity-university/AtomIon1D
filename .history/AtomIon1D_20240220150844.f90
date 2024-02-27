module BasisSets
   implicit none
   double precision, allocatable :: S_lm(:,:)
   double precision, allocatable :: S_mr(:,:)
TYPE basis
   integer Left, Right, xDim
   double precision, allocatable :: u(:,:,:), uxx(:,:,:)
   double precision, allocatable :: S(:,:)
   double precision, allocatable :: H(:,:)
   double precision, allocatable :: Psi(:,:),Energies(:,:)
   double precision, allocatable :: nrg_mat(:,:)
   double precision, allocatable :: psi_mat(:,:)
   double precision alpha, beta
   integer, allocatable :: xBounds(:)
END TYPE basis
contains
   subroutine AllocateBasis(T,Left, Right,Order, LegPoints, xNumPoints,NumStates,RSteps)
      implicit none
      integer Left, Right, LegPoints, xNumPoints, Order, MatrixDim, NumStates, ncv, RSteps
      TYPE(basis) T
      T%Left = Left
      T%Right = Right
      T%xDim = xNumPoints+Order-3
      if (Left .eq. 2) T%xDim = T%xDim + 1
      if (Right .eq. 2) T%xDim = T%xDim + 1
      MatrixDim = T%xDim
      ncv = 2*NumStates
      allocate(T%u(LegPoints, xNumPoints, T%xDim))
      allocate(T%uxx(LegPoints,xNumPoints,T%xDim))
      allocate(T%xBounds(xNumPoints+2*Order))
      allocate(T%S(Order+1,T%xDim),T%H(Order+1,T%xDim))
      allocate(T%Psi(MatrixDim,ncv))
      allocate(T%Energies(ncv,2))
      allocate(T%nrg_mat(RSteps,NumStates))
   end subroutine AllocateBasis
   subroutine deAllocateBasis(T)
      implicit none
      TYPE(basis) T
      deallocate(T%u)
      deallocate(T%uxx)
      deallocate(T%xBounds)
      deallocate(T%S)
      deallocate(T%H)
      deallocate(T%Psi)
      deallocate(T%Energies)
      deallocate(T%nrg_mat)
   end subroutine deAllocateBasis
end module BasisSets

program HHL1DHyperspherical
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
   integer, allocatable :: op_mat_1(:,:), op_mat_old_1(:,:)
   integer, allocatable :: op_mat_2(:,:), op_mat_old_2(:,:)
   double precision alpha,tmp_alpha,tmp_beta,mass,Shift,Shift2,NumStateInc,mi,ma,theta_c,mgamma
   double precision RLeft,RRight,RDerivDelt,DD,L,RFirst,RLast,XFirst,XLast,StepX,StepR, xMin,xMax
   double precision, allocatable :: R(:)
   double precision, allocatable :: xPoints(:)
   double precision, allocatable :: slopes(:,:)
   double precision, allocatable :: tmp_energies(:), tmp_psi(:)
   integer, allocatable :: OneToEighty(:)
   REAL t1, t2

   logical, allocatable :: Select(:)

   integer swapstate,iparam(11),ncv,info,troubleshooting, swapstate1, swapstate2
   integer i,j,k,iR,NumFirst,NumBound, swapped_flag_1, swapped_flag_2
   integer LeadDim,MatrixDim,HalfBandWidth
   integer xDim, Nbs
   integer, allocatable :: iwork(:)
   double precision Tol,RChange
   double precision TotalMemory
   double precision mu, mu12, mu123, r0diatom, dDiatom, etaOVERpi, Pi
   double precision YVal_L,RVal_L, RVal_R

   double precision, allocatable :: oldPsi(:,:), oldEnergies(:,:)
   double precision, allocatable :: LUFac(:,:),workl(:)
   double precision, allocatable :: workd(:),Residuals(:)
   double precision, allocatable :: xLeg(:),wLeg(:)
   !  double precision, allocatable :: lPsi(:,:),mPsi(:,:),rPsi(:,:),Energies(:,:)
   double precision, allocatable :: P(:,:),Q(:,:),dP(:,:)
   double precision ur(1:50000),acoef,bcoef,diff
   double precision sec,time,Rinitial,secp,timep,Rvalue, sNb, sbc, C4,lho
   double precision hbar, phi, amu,omega,Rstar, dum
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
   theta_c = datan(mgamma)     ! same as Seths beta !!! arctan mu, d?atan
   write(6,*) "theta_c = ", theta_c

   !     read in grid information
   read(5,*)
   read(5,*)
   read(5,*) xNumPoints, omega, Nbs, C4, phi
   write(6,*) xNumPoints, omega, Nbs, C4, phi
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
   mu12=mi*ma/(mi+ma)  ! ion-atom reduced mass
   !      mu=dsqrt(mi*ma)

   read(5,*)
   read(5,*)
   read(5,*) RSteps,RDerivDelt,RFirst,RLast
   write(6,*) RSteps,RDerivDelt,RFirst,RLast


   !write(6,*) "Resetting RFirst from: ", RFirst
   !RFirst = dsqrt(mu/(1+mu**2))*(RStar/(lho*(Nbs*Pi+phi)))+2*RDerivDelt
   !write(6,*) "to: ", RFirst

  allocate(slopes(RSteps,NumStates))
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

  allocate(S_lm(2*HalfBandWidth+1,xDim))
  allocate(S_mr(2*HalfBandWidth,xDim))
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

  !allocate(P(NumStates,NumStates),Q(NumStates,NumStates),dP(NumStates,NumStates))

  ncv = 2*NumStates
  LeadDim = 3*HalfBandWidth+1
  allocate(oldPsi(MatrixDim,ncv))
  allocate(oldEnergies(ncv,2))
  allocate(iwork(MatrixDim))
  allocate(Select(ncv))
  allocate(LUFac(LeadDim,MatrixDim))
  allocate(workl(ncv*ncv+8*ncv))
  allocate(workd(3*MatrixDim))
  allocate(Residuals(MatrixDim))
  allocate(OneToEighty(NumStates))
  info = 0
  iR=1
  Tol=1e-20

  NumBound=0
  !swapped_flag_1 = 0
  !swapped_flag_2 = 0
  
  call AllocateBasis(PB,2,2,Order, LegPoints, xNumPoints, NumStates, RSteps)
  call AllocateBasis(CB,Left,Right,Order, LegPoints, xNumPoints, NumStates, RSteps)
   !if (CouplingFlag .ne. 0) then
   !   call AllocateBasis(LB,Left,Right,Order, LegPoints, xNumPoints, NumStates, RSteps)
   !   call AllocateBasis(RB,Left,Right,Order, LegPoints, xNumPoints, NumStates, RSteps)
   !   call AllocateBasis(CB_old,Left,Right,Order, LegPoints, xNumPoints, NumStates, RSteps)
      !write(6,*) 'Left, right basis splines allocated.'
   !end if

   open(unit = 203, file = "ad_energies.dat") !printout of raw energies without diabatization method
   !open(unit = 200, file = "energiesv1.dat") !version 1 uses "tail swapping" method (diabats correspond to states 1 and 2)
   !open(unit = 201, file = "energiesv2.dat") !version 2 uses "operator-on-the-fly" method (diabats correspond to states 79 and 80)
   !open(unit = 101, file = "P_Matrix.dat")
   !open(unit = 102, file = "Q_Matrix.dat")

   allocate(op_mat_1(NumStates, NumStates))
   allocate(op_mat_old_1(NumStates, NumStates))
   allocate(op_mat_2(NumStates, NumStates))
   allocate(op_mat_old_2(NumStates, NumStates))
   allocate(tmp_energies(NumStates))
   allocate(tmp_psi(NumStates))
   swapstate = 0
   swapstate1 = 0
   swapstate2 = 0
   max_state = NumStates-1
   max_state1 = NumStates-1
      max_state2 = NumStates-1
   if(troubleshooting.eq.0) then
   CrossingFlag = 1
   RChange=100.d0
   do iR = RSteps,1,-1 !RSteps
      call CPU_TIME(t1)
      NumFirst=NumStates
      if (R(iR).gt.RChange) then
         NumFirst=NumBound
      endif
      NumStateInc=NumStates-NumFirst
     !c----------------------------------------------------------------------------------------
     !     must move this block inside the loop over iR if the grid is adaptive
      !print*, 'Beginning iteration at iR = ', iR
      !print*, 'Such that R(iR) = ', R(iR)
      !RLeft = R(iR)-RDerivDelt !make grid based on leftmost point
      !print*, 'calling GridMaker.'
      !if (CouplingFlag .eq. 0) then
         xMin = theta_c + asin(dsqrt(mu/(1d0+mu**2))* sNb/R(iR)*(Rstar/lho)) !Based on fixed sNb
         xMax = Pi + theta_c - asin(dsqrt(mu/(1d0+mu**2))* sNb/R(iR)*(Rstar/lho))
      !endif
      !if (CouplingFlag .eq. 1) then
      !   xMin = theta_c + asin(dsqrt(mu/(1d0+mu**2))* sNb/RLeft*(Rstar/lho)) ! Now based on RLeft and fixed sNb.
      !   xMax = Pi + theta_c - asin(dsqrt(mu/(1d0+mu**2))* sNb/RLeft*(Rstar/lho)) ! ''
      !endif
      !print*, 'Using xMin = ', xMin
      !print*, 'Using xMax = ', xMax
      sbc = (lho/Rstar)*R(iR)*dsqrt((1d0 + mu**2)/mu)*SIN(xMin - theta_c)
      !print*, 'sbc using sbc = (lho/Rstar)*R(iR)*dsqrt((1+mu**2)/mu)*SIN(xMin - theta_c):', sbc
      if(xMax.le.xMin) then
         write(6,*) "minimum hyperradius too small."
         stop
      endif
      call GridMakerIA(mu,mu12,theta_c,Rstar,R(iR),sbc,xNumPoints,xMin,xMax,xPoints)
         !print*, 'done... Calculating Primitive Basis Functions'
         tmp_alpha = 1
         tmp_beta = 1
         call CalcBasisFuncsBP(PB%Left,PB%Right, RVal_L, RVal_R, Order,xPoints,&
         LegPoints,xLeg,PB%xDim,PB%xBounds,xNumPoints,0,PB%u, tmp_alpha, tmp_beta,.false.)
         call CalcBasisFuncsBP(PB%Left,PB%Right, RVal_L, RVal_R, Order,xPoints,&
         LegPoints,xLeg,PB%xDim,PB%xBounds,xNumPoints,2,PB%uxx, tmp_alpha, tmp_beta,.false.)
      !print*, 'done... Calculating overlap matrix'
      !     must move this block inside the loop if the grid is adaptive
      !----------------------------------------------------------------------------------------
      call CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,PB%xDim,xNumPoints,PB%u,PB%xBounds,HalfBandWidth,PB%S)
!******** CONSTRUCTION OF BASIS SETS ********
      !write(6,*) 'Constructing center basis set.'
      YVal_L = (lho/Rstar)*R(iR)*dsqrt((1+mu**2)/mu - ((sbc*Rstar)/(lho*R(iR)))**2) * ((1/sbc)+(1/sbc**2)*COTAN(-(1/sbc)+phi))
      RVal_L = 1.d0/YVal_L
      RVal_R = -RVal_L
      !write(6,*) 'Center Basis YVal_L = ', YVal_L
      !write(6,*) 'Center Basis RVal_L = ', RVal_L
      !call CalcBasisFuncsBP(CB%Left,CB%Right, RVal_L, RVal_R, Order,&
      !xPoints,LegPoints,xLeg,CB%xDim,CB%xBounds,xNumPoints,0,CB%u, CB%alpha, CB%beta, .true.)
      !call calc_Physical_Set(CB,PB,RVal_L,RVal_R,Order, xNumPoints, xPoints)
      !write(6,*) 'CB%alpha (u) = ', CB%alpha
      !write(6,*) 'CB%beta (u) = ', CB%beta
      call CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,CB%xDim,xNumPoints,CB%u,CB%xBounds,HalfBandWidth,CB%S) ! added
      !write(6,*) 'Center Basis Set Constructed.'
      call CalcHamiltonian(alpha,R(iR),mu,mi,theta_c,C4,L,Order,xPoints,&
         LegPoints,xLeg,wLeg,CB%xDim,xNumPoints,CB%u,CB%uxx,CB%xBounds,HalfBandWidth,CB%H)
      call MyDsband(Select,CB%Energies,CB%Psi,MatrixDim,Shift,MatrixDim,CB%H,CB%S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,&
         NumStates,Tol,Residuals,ncv,CB%Psi,MatrixDim,iparam,workd,workl,ncv*ncv+8*ncv,iwork,info)
      ! if (CouplingFlag .ne. 0) then
      !    !write(6,*) 'Constructing left basis set.'
      !    RLeft = R(iR)-RDerivDelt ! This is just for reference, as RLeft is already declared in this loop
      !    YVal_L = (lho/Rstar)*RLeft*dsqrt((1+mu**2)/mu - ((sNb*Rstar)/(lho*RLeft))**2) * ((1/sNb)+(1/sNb**2)*COTAN(-(1/sNb)+phi))
      !    RVal_L = 1.d0/YVal_L
      !    RVal_R = -RVal_L
      !    !write(6,*) 'Left Basis YVal_L = ', YVal_L
      !    !write(6,*) 'Left basis RVal_L = ', RVal_L
      !    LB%xBounds = CB%xBounds
      !    call calc_Physical_Set(LB,PB,RVal_L,RVal_R,Order, xNumPoints, xPoints)
      !    !write(6,*) 'LB%alpha (u) = ', LB%alpha
      !    !write(6,*) 'LB%beta (u) = ', LB%beta
      !    call CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,LB%xDim,xNumPoints,LB%u,LB%xBounds,HalfBandWidth,LB%S)
      !    !write(6,*) 'Left Basis Set Constructed.'
      !    call CalcHamiltonian(alpha,RLeft,mu,mi,theta_c,C4,L,Order,xPoints,&
      !       LegPoints,xLeg,wLeg,LB%xDim,xNumPoints,LB%u,LB%uxx,LB%xBounds,HalfBandWidth,LB%H)
      !    call MyDsband(Select,LB%Energies,LB%Psi,MatrixDim,Shift,MatrixDim,LB%H,LB%S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,&
      !       NumStates,Tol,Residuals,ncv,LB%Psi,MatrixDim,iparam,workd,workl,ncv*ncv+8*ncv,iwork,info)
      
      !    !write(6,*) 'Constructing right basis set.'
      !    RRight =  R(iR)+RDerivDelt
      !    YVal_L = (lho/Rstar)*RRight*dsqrt((1+mu**2)/mu - ((sbc*Rstar)/(lho*RRight))**2) * ((1/sbc)+(1/sbc**2)*COTAN(-(1/sbc)+phi))
      !    RVal_L = 1.d0/YVal_L
      !    RVal_R = -RVal_L
      !    !write(6,*) 'Right Basis YVal_L = ', YVal_L
      !    !write(6,*) 'Right basis RVal_L = ', RVal_L
      !    RB%xBounds = CB%xBounds
      !    call calc_Physical_Set(RB,PB,RVal_L,RVal_R,Order, xNumPoints, xPoints)
      !    !write(6,*) 'RB%alpha (u) = ', RB%alpha
      !    !write(6,*) 'RB%beta (u) = ', RB%beta
      !    call CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,RB%xDim,xNumPoints,RB%u,RB%xBounds,HalfBandWidth,RB%S)
      !    !write(6,*) 'Right Basis Set Constructed.'
      !    call CalcHamiltonian(alpha,RRight,mu,mi,theta_c,C4,L,Order,xPoints,&
      !       LegPoints,xLeg,wLeg,RB%xDim,xNumPoints,RB%u,RB%uxx,RB%xBounds,HalfBandWidth,RB%H)
      !    call MyDsband(Select,RB%Energies,RB%Psi,MatrixDim,Shift,MatrixDim,RB%H,RB%S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,&
      !       NumStates,Tol,Residuals,ncv,RB%Psi,MatrixDim,iparam,workd,workl,ncv*ncv+8*ncv,iwork,info)
      ! endif
         !do i = 1, NumStates
         !   if (CouplingFlag.eq.1) LB%nrg_mat(iR,i) = LB%Energies(i,1)
         !   CB%nrg_mat(iR,i) = CB%Energies(i,1)
         !   if (CouplingFlag.eq.1) RB%nrg_mat(iR,i) = RB%Energies(i,1)
         !enddo
     !    if (CouplingFlag.eq.1) write(203, 20) R(iR), (RB%Energies(i,1), i = 1,NumStates) !!
         write(203, 20) R(iR), (CB%Energies(i,1), i = 1,NumStates) !!New: writing un-diabatized energies for comparison in poster plot
         !if (CouplingFlag.eq.1) write(203, 20) R(iR), (LB%Energies(i,1), i = 1,NumStates) !!
         !if((CrossingFlag.eq.1).and.(CouplingFlag.eq.1)) then
            !print*, 'Looking for sharply avoided crossings.'
            !if(R(iR).gt.10) call Avoided_CrossingsV2(CB,LB,RB,RDerivDelt,R,iR,RSteps,NumStates,swapstate1,max_state1)
            !if(R(iR).gt.10) call Avoided_CrossingsV2(CB,LB,RB,RDerivDelt,R,iR,RSteps,NumStates,swapstate2,max_state2)
          !  if(R(iR).gt.5) call Avoided_CrossingsV3(CB,LB,RB,RDerivDelt,R,iR,RSteps,NumStates,swapstate1,max_state1,&
          !  0, curr_overlaps, last_overlaps, 1)
          !  last_overlaps(1) = curr_overlaps(1)
          !  if(R(iR).gt.5) call Avoided_CrossingsV3(CB,LB,RB,RDerivDelt,R,iR,RSteps,NumStates,swapstate2,max_state2,&
          !  max_state1, curr_overlaps, last_overlaps, 2)!2 is for second diabat
          !  last_overlaps(2) = curr_overlaps(2)
          !  endif !!Swapstates found
      ! if(iR.eq.RSteps) then
      !    do i = NumStates, 1, -1
      !       do j = 1, NumStates
      !          if (i.eq.j) op_mat_old_1(i,j) = 1
      !          if (i.ne.j) op_mat_old_1(i,j) = 0
      !          if (i.eq.j) op_mat_old_2(i,j) = 1
      !          if (i.ne.j) op_mat_old_2(i,j) = 0
      !       enddo
      !       OneToEighty(i) = i
      !    enddo
      ! endif
      ! if(swapstate1 .ne. 0) then
      !    call populate_operator(op_mat_1, swapstate1, NumStates) !make new operator
      !    op_mat_1 = matmul((op_mat_old_1), (op_mat_1)) !multiply previous with new operator
      !    swapped_flag_1 = 1 ! if just one state needs to be swapped, then every iteration after must be swapped
      ! endif
      ! if(swapstate2.ne.0) then
      !    call populate_operator(op_mat_2, swapstate2, NumStates) !make new operator
      !    op_mat_2 = matmul((op_mat_old_2),(op_mat_2)) !multiply previous with new operator
      !    swapped_flag_2 = 1
      ! endif

      ! if(swapped_flag_2 .eq. 1) then
      !    OneToEighty = matmul(op_mat_2, OneToEighty)
      !    !print*, 'swapping energies'
      !    tmp_energies = CB%energies(1:NumStates,1)
      !    CB%energies(1:NumStates,1) = (matmul(op_mat_2, tmp_energies))
      !    tmp_energies = RB%energies(1:NumStates,1)
      !    RB%energies(1:NumStates,1) = (matmul(op_mat_2, tmp_energies))
      !    tmp_energies = LB%energies(1:NumStates,1)
      !    LB%energies(1:NumStates,1) = (matmul(op_mat_2, tmp_energies))
      !    !print*, 'swapping psis'
      !    do i = 1, xDim
      !       tmp_psi = CB%Psi(i, 1:NumStates)
      !       CB%Psi(i, 1:NumStates) = matmul(op_mat_2, tmp_psi)
      !       tmp_psi = LB%Psi(i, 1:NumStates)
      !       LB%Psi(i, 1:NumStates) = matmul(op_mat_2, tmp_psi)
      !       tmp_psi = RB%Psi(i, 1:NumStates)
      !       RB%Psi(i, 1:NumStates) = matmul(op_mat_2, tmp_psi)
      !    enddo
      !    op_mat_old_2 = op_mat_2 !ready the new operator for next iteration
      ! endif
      ! if(swapped_flag_1 .eq. 1) then
      !    OneToEighty = matmul(op_mat_1, OneToEighty)
      !    !print*, 'swapping energies'
      !    tmp_energies = CB%energies(1:NumStates,1)
      !    CB%energies(1:NumStates,1) = (matmul(op_mat_1, tmp_energies))
      !    tmp_energies = RB%energies(1:NumStates,1)
      !    RB%energies(1:NumStates,1) = (matmul(op_mat_1, tmp_energies))
      !    tmp_energies = LB%energies(1:NumStates,1)
      !    LB%energies(1:NumStates,1) = (matmul(op_mat_1, tmp_energies))
      !    !print*, 'swapping psis'
      !    do i = 1, xDim
      !       tmp_psi = CB%Psi(i, 1:NumStates)
      !       CB%Psi(i, 1:NumStates) = matmul(op_mat_1, tmp_psi)
      !       tmp_psi = LB%Psi(i, 1:NumStates)
      !       LB%Psi(i, 1:NumStates) = matmul(op_mat_1, tmp_psi)
      !       tmp_psi = RB%Psi(i, 1:NumStates)
      !       RB%Psi(i, 1:NumStates) = matmul(op_mat_1, tmp_psi)
      !    enddo
      !    op_mat_old_1 = op_mat_1 !ready the new operator for next iteration
      ! endif

!      if((swapstate1.ne.0).or.(swapstate1.ne.0)) then  !NPM to Angela: This conditional is redundant. What is OnetoEighty?
 !        write(306, *) (OnetoEighty(i), i = 1, NumStates)
 !     endif

  !    if (CouplingFlag.eq.1) write(201, 20) RRight, (RB%Energies(i,1), i = 1,NumStates) !writing energies to og_energies.dat
  !    write(201, 20) R(iR), (CB%Energies(i,1), i = 1,NumStates)
   !   if (CouplingFlag.eq.1) write(201, 20) RLeft, (LB%Energies(i,1), i = 1,NumStates)

      if(iR .ne. RSteps) then
         call FixPhase(NumStates,HalfBandWidth,MatrixDim,CB%S,ncv,oldPsi,CB%Psi)
      endif
      oldPsi = CB%Psi
      call CalcEigenErrors(info,iparam,MatrixDim,CB%H,HalfBandWidth+1,CB%S,HalfBandWidth,NumStates,CB%Psi,CB%Energies,ncv)
      ! if (CouplingFlag.eq.1) then
      !    call FixPhase(NumStates,HalfBandWidth,MatrixDim,LB%S,ncv,oldPsi,LB%Psi)
      !    call CalcEigenErrors(info,iparam,MatrixDim,LB%H,HalfBandWidth+1,LB%S,HalfBandWidth,NumStates,LB%Psi,LB%Energies,ncv)
      !    call FixPhase(NumStates,HalfBandWidth,MatrixDim,RB%S,ncv,oldPsi,RB%Psi)
      !    call CalcEigenErrors(info,iparam,MatrixDim,RB%H,HalfBandWidth+1,RB%S,HalfBandWidth,NumStates,RB%Psi,RB%Energies,ncv)
      ! endif

   !write(6,*) 'Calculating the overlaps and couplings!'

  ! write(155, *) '***At R(iR) =', R(iR), '***'
 !  call calcCouplings_v2(LB, CB, RB, P, Q, PB%S, HalfBandWidth, NumStates, RDerivDelt)

 !       write(101,*) R(iR)
 !       write(102,*) R(iR)
 !       write(103,*) R(iR)
 !       do i = 1,min(NumStates,iparam(5))
!           write(101,20) (P(i,j), j = 1,min(NumStates,iparam(5)))
!           write(102,20) (Q(i,j), j = 1,min(NumStates,iparam(5)))
!           write(103,20) (dP(i,j), j = 1,min(NumStates,iparam(5)))
  !      enddo

   !   if(swapstate1.ne.0) then
   !      write(610,*) 'Swap1: R(iR) =', R(iR)
        !do i = swapstate1
           !write(610,*) (P(i,j), j = swapstate1,swapstate1+1)
   !        write(610,*) 'swapstate1 = ', swapstate1
   !     !enddo
   !   endif
   !   if(swapstate2.ne.0) then
   !      write(610,*) 'Swap2: R(iR) =', R(iR)
        !do i = swapstate1
           !write(610,*) (P(i,j), j = swapstate1,swapstate1+1)
    !       write(610,*) 'swapstate2 = ', swapstate2
        !enddo
    !  endif

     !write(6,*)
     !write(6,*) 'RMid = ', R(iR)
     do i = 1,min(NumStates,iparam(5))
        write(6,*) 'Energy(',i,') = ',CB%Energies(i,1),'  Error = ', CB%Energies(i,2)
     enddo

      !write(6,*) 'shape(energies): ', shape(CB%Energies)

!    Adjusting Shift
     ur(iR) = CB%Energies(1,1)
     Shift = -200d0
     if(ur(iR).lt.0d0) then
        Shift = ur(iR)*10d0
     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
   !   if (PsiFlag .ne. 0) then
   !      do i = 1,xNumPoints
   !         write(97,*) xPoints(i)
   !      enddo
   !      do i = 1,MatrixDim
   !         write(999+iR,20) (CB%Psi(i,j), j = 1,NumStates)
   !      enddo
   !      close(unit=999+iR)
   !   endif

   call CPU_TIME(t2)
   write(6,*) 'Remaining time (min): ', (t2-t1)*iR/60, R(iR), (CB%Energies(i,1), i=1,5)

  enddo

   ! do iR = RSteps, 1, -1
   !    !write(301, 20) CB%psi(iR, j) !compare psi matrices with and without psi_mat vs psi
   !    write(200,20) R(iR)+RDerivDelt, (RB%nrg_mat(iR,j), j=1, NumStates)
   !    write(200,20) R(iR), (CB%nrg_mat(iR,j), j=1, NumStates)
   !    write(200,20) R(iR)-RDerivDelt, (LB%nrg_mat(iR,j), j=1, NumStates)
   ! enddo
   
   ! endif

   ! open(203)

   ! if(troubleshooting.eq.1) then
   !    do iR = RSteps*3,1,-3
   !       read(203,*,END=100) tmp, (RB%nrg_mat(iR/3,j), j=1,NumStates)
   !       read(203,*,END=100) tmp, (CB%nrg_mat(iR/3,j), j=1,NumStates)
   !       read(203,*, END=100) tmp, (LB%nrg_mat(iR/3,j), j=1,NumStates)
   !    enddo
   !    100 close(203)
   !    swapstate1 = 0
   !    swapstate2 = 0
   !    do iR = RSteps, 1,-1
   !       if(R(iR).gt.5) call Avoided_CrossingsV3(CB,LB,RB,RDerivDelt,R,iR,RSteps,NumStates,swapstate1,max_state1,&
   !       0, curr_overlaps, last_overlaps, 1)
   !       last_overlaps(1) = curr_overlaps(1)
   !       if(R(iR).gt.5) call Avoided_CrossingsV3(CB,LB,RB,RDerivDelt,R,iR,RSteps,NumStates,swapstate2,max_state2,&
   !       max_state1, curr_overlaps, last_overlaps, 2)!2 is for second adiabat
   !       last_overlaps(2) = curr_overlaps(2)
   !    enddo
   !    do iR = RSteps, 1, -1
   !       write(200,20) R(iR)+RDerivDelt, (RB%nrg_mat(iR,j), j=1, NumStates)
   !       write(200,20) R(iR), (CB%nrg_mat(iR,j), j=1, NumStates)
   !       write(200,20) R(iR)-RDerivDelt, (LB%nrg_mat(iR,j), j=1, NumStates)
   !    enddo
   
   ! endif

   close(11)
   close(101)
   close(102)
   close(200)
   close(300)
   close(301)
   !100 close(203)

  deallocate(S_lm, S_mr)
  deallocate(iwork)
  deallocate(Select)
  deallocate(LUFac)
  deallocate(workl)
  deallocate(workd)
  deallocate(Residuals)
  deallocate(P,Q,dP)
  deallocate(xPoints)
  deallocate(xLeg,wLeg)
  deallocate(op_mat_1)
  deallocate(op_mat_2)
  deallocate(op_mat_old_1)
  deallocate(op_mat_old_2)

  call deAllocateBasis(PB)
  call deAllocateBasis(CB)
  if (CouplingFlag .ne. 0) call deAllocateBasis(LB)
  if (CouplingFlag .ne. 0) call deAllocateBasis(RB)
  
  deallocate(R)

10 format(1P,100e25.15)
20 format(1P,100e16.8)
1002 format(a64)

  stop
end program HHL1DHyperspherical

subroutine reorder_Q2(state_L, state_H, NumStates, RSteps) !state_L and state_H based on square of arrays being moved out
implicit none
integer NumStates, RSteps, i, j, swaps, state_H, state_L, iR, iter
double precision tmp1, tmp2, P(NumStates, NumStates), tmp_P(State_H-State_L),R !tmp P will hold the array moved out

open(unit = 102, file = "Q_Matrix.dat")
open(unit = 501, file = "reordered_Q.dat")

do iR = RSteps, 1, -1
   print*, 'looping'
    read(102,*,END=100) R
    do i = 1, 80
        read(102,*,END=100) (P(i, j), j=1,NumStates)
    enddo

    !store
   do i = 1, NumStates
      iter = state_H+1 ! iter = 79+1 = 80, want 8
      tmp1 = P(i,state_H+1) ! 80
      tmp2 = P(i,state_H)
      do while(iter.ne.state_L+1) !3+1 = 4, no more!! when iter = 5, you do still want 4->5, and then 3->4; 
         P(i,iter) = P(i,iter-2) !78 -> 80, 77 -> 79, ... , 3 -> 5,
         iter = iter-1
      enddo
      P(i, state_L) = tmp1 !3 = 80
      P(i, state_L+1) = tmp2 ! 4 = 70
   enddo
   do i = 1, NumStates
      iter = state_H+1 ! iter = 79+1 = 80, want 8
      tmp1 = P(state_H+1, i)
      tmp2 = P(state_H, i)
      do while(iter.ne.state_L+1) !3+1 = 4, no more!! when iter = 5, you do still want 4->5, and then 3->4; 
         P(iter, i) = P(iter-2, i) !78 -> 80, 77 -> 79, ... , 3 -> 5,
         iter = iter-1
      enddo
      P(state_L, i) = tmp1
      P(state_L+1, i) = tmp2
   enddo
   write(501,*) R
   do i = 1,NumStates
      write(501,*) (P(i,j), j = 1,NumStates)
   enddo
enddo
100 close(102)
close(501)

end subroutine reorder_Q2


subroutine reorder_P2(state_L, state_H, NumStates, RSteps) !state_L and state_H based on square of arrays being moved out
implicit none
integer NumStates, RSteps, i, j, swaps, state_H, state_L, iR, iter
double precision tmp1, tmp2, P(NumStates, NumStates), tmp_P(State_H-State_L),R !tmp P will hold the array moved out

open(unit = 101, file = "P_Matrix.dat")
open(unit = 500, file = "reordered_P.dat")

do iR = RSteps, 1, -1
   print*, 'looping'
    read(101,*,END=100) R
    do i = 1, 80
        read(101,*,END=100) (P(i, j), j=1,NumStates)
    enddo

    !store
   do i = 1, NumStates
      iter = state_H+1 !
      tmp1 = P(i,state_H+1) !
      tmp2 = P(i,state_H)
      do while(iter.ne.state_L+1) ! 
         P(i,iter) = P(i,iter-2) !
         iter = iter-1
      enddo
      P(i, state_L) = tmp1 !3 = 80
      P(i, state_L+1) = tmp2 ! 4 = 70
   enddo
   do i = 1, NumStates
      iter = state_H+1 ! 
      tmp1 = P(state_H+1, i)
      tmp2 = P(state_H, i)
      do while(iter.ne.state_L+1) ! 
         P(iter, i) = P(iter-2, i) 
         iter = iter-1
      enddo
      P(state_L, i) = tmp1
      P(state_L+1, i) = tmp2
   enddo
   write(500,*) R
   do i = 1,NumStates
      write(500,*) (P(i,j), j = 1,NumStates)
   enddo
enddo
100 close(101)
close(500)

end subroutine reorder_P2

subroutine populate_operator(op,i, NumStates)
implicit none
integer x,y,NumStates, op(NumStates,Numstates),i
do x = 1, NumStates
   do y = 1, NumStates
      if(x.eq.y) op(x,y) = 1
      if(x.ne.y) op(x,y) = 0
   enddo
enddo
op(i,i) = 0
op(i, i+1) = 1
op(i+1,i) = 1
op(i+1, i+1) = 0
end subroutine populate_operator

! subroutine Avoided_CrossingsV3(CB,LB,RB,RDerivDelt,R,iR,RSteps,NumStates,swapstate,max_state,&
! min_state, curr_overlaps, last_overlaps, ID)
! use BasisSets
! implicit none
! type(BASIS) CB, LB, RB
! double precision RDerivDelt, R(RSteps), slopes(2,2), num, den, tmp
! integer iR, RSteps, NumStates, swapstate, max_state, state, found_flag
! integer ID, curr_overlaps(2), last_overlaps(2), min_state, slope_req
! if(iR.lt.RSteps) then
!    found_flag = 0
!    swapstate = 0
!    state = max_state
!    do while((found_flag.ne.1).and.(state.ne.min_state))
!       !This part creates a 2x2 matrix to represent the slopes of each point in the "box"
!       slopes(1,1) = dabs(RB%nrg_mat(iR, state+1)-LB%nrg_mat(iR, state+1))/(2*RDerivDelt)
!       slopes(1,2) = dabs(RB%nrg_mat(iR+1, state+1)-LB%nrg_mat(iR+1, state+1))/(2*RDerivDelt)
!       slopes(2,1) = dabs(RB%nrg_mat(iR, state)-LB%nrg_mat(iR, state))/(2*RDerivDelt)
!       slopes(2,2) = dabs(RB%nrg_mat(iR+1, state)-LB%nrg_mat(iR+1, state))/(2*RDerivDelt)
!       !This part compares the energy difference between two states at a given hyperradius.
!       !If the energy difference is zero, then the triad points are exactly at the avoided
!       !crossing.
!       if((dabs(CB%nrg_mat(iR, state)-CB%nrg_mat(iR, state+1)).lt.1E-2).or.&
!       (dabs(LB%nrg_mat(iR, state)-LB%nrg_mat(iR, state+1)).lt.1E-2).or.&
!       (dabs(RB%nrg_mat(iR, state)-RB%nrg_mat(iR, state+1)).lt.1E-2)) then
!          !if((ID.ne.1).and.(curr_overlaps(1).eq.0).or.(ID.eq.1)) then !this if statement is no longer necessary; added min_state
!             curr_overlaps(ID) = 1
!          !else
!          !   curr_overlaps(ID) = 0
!          !endif
!       else
!          curr_overlaps(ID) = 0
!       endif
!       !The above if statement prevents this type of crossing from being caught "twice" (by the
!       !iterator through the higher and lower diabat) by only allowing the higher diabat to catch
!       !a crossing at a given R if the lower diabat has not caught it already.
!       !Note: this is not an issue for the other method, because the other method uses two states
!       !and two hyperradii to find a crossing, and once one is found, the ordering of the states changes
!       !and the same crossing won't be found twice.
      
!       !Now, the if-statement for finding a crossing; 100 is a number that seems to work for this system.
!       !If the energy difference method caugh a crossing that would have been caught at the next iteration
!       !by the slopes method, the last_overlaps stores that data and doesn't pass the if-statement for that
!       !next iteration
!       if((dabs(slopes(2,1)-slopes(2,2))/dabs(slopes(2,1)-slopes(1,2))).gt.100&
!       .and.((dabs(slopes(1,1)-slopes(1,2))/dabs(slopes(2,1)-slopes(1,2))).gt.100)) then
!          slope_req = 1
!       else
!          slope_req = 0
!       endif
!       if(((slope_req.eq.1).and.(last_overlaps(ID).ne.1)).or.((curr_overlaps(ID).eq.1))) then
!          print*, 'max_state = ', max_state, 'ID = ', ID
!          print*, 'Found sharp avoided crossing at iR =', iR, 'R(iR)= ', R(iR),&
!          'with energy = ', CB%nrg_mat(iR,state), 'and state = ', state
!          print*, 'and curr_overlaps(ID) = ', curr_overlaps(ID)
!          swapstate = state
!          ! if(slope_req.eq.1) then
!          !    call swap_full(CB%nrg_mat, swapstate, iR, RSteps, NumStates)
!          !    call swap_full(LB%nrg_mat, swapstate, iR, RSteps, NumStates)
!          !    call swap_full(RB%nrg_mat, swapstate, iR, RSteps, NumStates)
!          ! endif
!          ! if((slope_req.eq.0).and.(curr_overlaps(ID).eq.1)) then
!          !    call swap_full_mod(CB%nrg_mat, swapstate, iR, RSteps, NumStates) !attempt to keep consistency of "where" swapping occurs
!          !    call swap_full_mod(LB%nrg_mat, swapstate, iR, RSteps, NumStates) !if it is found one hyperradial point sooner than
!          !    call swap_full_mod(RB%nrg_mat, swapstate, iR, RSteps, NumStates) !where it would/should have been found; the effect
!          ! endif                                                               !is mostly negligible, though.
!          call swap_full(CB%nrg_mat, swapstate, iR, RSteps, NumStates)
!          call swap_full(LB%nrg_mat, swapstate, iR, RSteps, NumStates)
!          call swap_full(RB%nrg_mat, swapstate, iR, RSteps, NumStates)
         
!          found_flag = found_flag+1
!          max_state = state
!       endif
!       state = state - 1
!    enddo
! endif
! end subroutine Avoided_CrossingsV3

! subroutine Avoided_CrossingsV2(CB,LB,RB,RDerivDelt,R,iR,RSteps,NumStates,swapstate,max_state)
! use BasisSets
! implicit none
! type(BASIS) CB, LB, RB
! double precision RDerivDelt, R(RSteps), slopes(2,2), num, den, tmp
! integer iR, RSteps, NumStates, swapstate, max_state, state, found_flag, overlapping

! if(iR.lt.RSteps) then
!    found_flag = 0
!    swapstate = 0
!    state = max_state !=80 at the beginning of main loop
!    !print*, 'max_state = ', max_state
!    do while((found_flag.ne.1).and.(state.ne.0))
!       slopes(1,1) = dabs(RB%nrg_mat(iR, state+1)-LB%nrg_mat(iR, state+1))/(2*RDerivDelt)
!       slopes(1,2) = dabs(RB%nrg_mat(iR+1, state+1)-LB%nrg_mat(iR+1, state+1))/(2*RDerivDelt)
!       slopes(2,1) = dabs(RB%nrg_mat(iR, state)-LB%nrg_mat(iR, state))/(2*RDerivDelt)
!       slopes(2,2) = dabs(RB%nrg_mat(iR+1, state)-LB%nrg_mat(iR+1, state))/(2*RDerivDelt)

!       if((((dabs(slopes(2,1)-slopes(2,2))/dabs(slopes(2,1)-slopes(1,2))).gt.100)&
!       .and.((dabs(slopes(1,1)-slopes(1,2))/dabs(slopes(2,1)-slopes(1,2))).gt.100))) then
!          !print*, 'max_state = ', max_state
!          !print*, 'Found sharp avoided crossing at iR =', iR, 'R(iR)= ', R(iR),&
!          !'with energy = ', CB%nrg_mat(iR,state), 'and state = ', state
!          swapstate = state
!          call swap_full(CB%nrg_mat, swapstate, iR, RSteps, NumStates)
!          call swap_full(LB%nrg_mat, swapstate, iR, RSteps, NumStates)
!          call swap_full(RB%nrg_mat, swapstate, iR, RSteps, NumStates)
!          found_flag = found_flag+1
!          max_state = state
!       endif
!       state = state - 1
!    enddo
! endif
! end subroutine Avoided_CrossingsV2

! subroutine swap_full(curr_nrg, swapstate, iR, RSteps, NumStates)
! implicit none
! double precision curr_nrg(RSteps, NumStates), tmp_nrg
! integer swapstate, iR, RSteps, NumStates, i
! do i = iR+1, RSteps
!   tmp_nrg = curr_nrg(i,swapstate)
!   curr_nrg(i,swapstate) = curr_nrg(i, swapstate+1)
!   curr_nrg(i,swapstate+1) = tmp_nrg
! enddo
! end subroutine swap_full

! subroutine swap_full_mod(curr_nrg, swapstate, iR, RSteps, NumStates)
! implicit none
! double precision curr_nrg(RSteps, NumStates), tmp_nrg
! integer swapstate, iR, RSteps, NumStates, i
! do i = iR, RSteps
!   tmp_nrg = curr_nrg(i,swapstate)
!   curr_nrg(i,swapstate) = curr_nrg(i, swapstate+1)
!   curr_nrg(i,swapstate+1) = tmp_nrg
! enddo
! end subroutine swap_full_mod

! subroutine swap_psi_sign(swapstate, psi, xDim, ncv)
! implicit none
! double precision psi(xDim,ncv), tmp!(xDim)
! integer swapstate, i, xDim, ncv

! do i = 1, xDim
!    tmp = psi(i, swapstate)
!    psi(i, swapstate) = sign(psi(i, swapstate), psi(i, swapstate+1)) !return value of psi+1 but with original sign of psi
!    psi(i, swapstate+1) = sign(psi(i, swapstate+1), tmp) !return value of psi, but with original sign of psi+1
! enddo

! ! tmp = psi(:,swapstate)
! ! psi(:,swapstate) = psi(:,swapstate+1)
! ! psi(:,swapstate+1) = tmp
! end subroutine swap_psi_sign

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
      end
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

!subroutine calcpsandqs(LBarray, CBarray, RBarray, PBArray, HalfBandWidth, NumStates, RDerivDelt) !!see if you can make an array with type bases...

!   call calcCouplings_v2(LB, CB, RB, P, Q, PB%S, HalfBandWidth, NumStates, RDerivDelt)

!end subroutine

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

        rai = mu**(-0.5d0) * R*sinx(lx,kx) -  mu**0.5d0 * R*cosx(lx,kx) !dsqrt(mu/mu12)*R*dabs(sinai(lx,kx))))

        potvalue = -C4/rai**4 + 0.5d0*mu*(R*cosx(lx,kx))**2
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
           write(170) 'phase fix at: ', 'j=',j, 'i=',i
        enddo
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
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GridMakerIA(mu,mu12,theta_c,Rstar,R,sbc,xNumPoints,xMin,xMax,xPoints)
  implicit none
  integer xNumPoints, testing
  double precision mu,R,r0,xMin,xMax,xPoints(xNumPoints)
  double precision mu12,mu123,theta_c,Rstar,sbc
  integer i,j,k,OPGRID,NP
  double precision Pi
  double precision r0New
  double precision xRswitch
  double precision xDelt,x0,x1,x2,x3,x4,deltax


  Pi = 3.1415926535897932385d0
  !r0New=r0*2.0d0
  deltax = 0.2d0*Rstar/R

   !testing = 1

  OPGRID=0

  NP = xNumPoints/4
  if(deltax.ge.(xMax - xMin)*0.25d0) then
  !if (testing.eq.1) then
     write(6,*) 'Switching to linear grid...' !PRINT WHERE (RADIUS) THIS SWITCH OCCURS !maybe try without optimization, use linear grid everywhere
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
     write(6,*) "using modified grid:(R,deltax) = ", R, deltax

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
     !write(6,*) "Using linear grid...(R,deltax) = ",R, deltax
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
subroutine GridMaker(mu,R,r0,xNumPoints,xMin,xMax,xPoints)
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
end subroutine GridMaker
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


