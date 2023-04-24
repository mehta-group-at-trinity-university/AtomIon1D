c     234567890

c Code taken from subroutine CalcOverlap

    subroutine CalcEdges(Order,xPoints,LegPoints,xLeg,wLeg,xDim,
     >     xNumPoints,u,xBounds,HalfBandWidth,S)
      implicit none
      integer Order,LegPoints,xDim,xNumPoints,xBounds(xNumPoints+2*Order),HalfBandWidth
      double precision xPoints(*),xLeg(*),wLeg(*)
      double precision S(HalfBandWidth+1,xDim) ! syntax is 2d array declaration with dimensions Halfbandwidth x xdim
      double precision u(LegPoints,xNumPoints,xDim)

      integer ix,ixp,kx,lx
      integer i1,i1p
      integer Row,NewRow,Col
      integer, allocatable :: kxMin(:,:),kxMax(:,:) ! Array is deallocated automatically when array is out of scope
      double precision a,b,m
      double precision xTempS
      double precision ax,bx
      double precision, allocatable :: xIntScale(:),xS(:,:)
   
      allocate(xIntScale(xNumPoints),xS(xDim,xDim)) ! "Allocate" seems to act like "new" or "malloc" in c++
      allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))


!u is basis
!ux is derivative
!uxx = ddx^2

! 0 - Wavefunction
! 1 - First Derivative
! 2 - No BC
! 3 - Log Derivative


      S = 0.0d0

      do kx = 1,xNumPoints-1
         ax = xPoints(kx)
         bx = xPoints(kx+1)
         xIntScale(kx) = 0.5d0*(bx-ax)
      enddo

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
c               write(26,*) ix,ixp,S(NewRow,Col)
            endif
         enddo
      enddo

      deallocate(xIntScale,xS)
      deallocate(kxMin,kxMax)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc