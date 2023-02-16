      subroutine CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,MatrixDim,xBounds,xNumPoints,Deriv,u)

      integer Left,Right,Order,LegPoints,MatrixDim,xBounds(*),xNumPoints,Deriv
      double precision xPoints(*),xLeg(*)
      double precision u(LegPoints,xNumPoints,MatrixDim)

      integer i,k,l,Count
      integer, allocatable :: t(:)
      double precision BSpline
      double precision x,ax,bx,xIntScale,xScaledZero

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
! length(t) is order*2+xNumPoints, t = 1, 1, 1, ..., 1, 1, 2, 3, 4, 5, 6, ..., xNumPoints, xNumPoints, xNumPoints, ..., xNumPoints
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
        end select
      end select

      deallocate(t)

      Count = 1
      select case (Left)
       case (0)
        do k = 1,xNumPoints-1
         ax = xPoints(k)
         bx = xPoints(k+1)
         xIntScale = 0.5d0*(bx-ax)
         xScaledZero = 0.5d0*(bx+ax)
         do l = 1,LegPoints
          x = xIntScale*xLeg(l)+xScaledZero
          u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,2,x)
         enddo
        enddo
        Count = Count + 1
       case (1)
        do k = 1,xNumPoints-1
         ax = xPoints(k)
         bx = xPoints(k+1)
         xIntScale = 0.5d0*(bx-ax)
         xScaledZero = 0.5d0*(bx+ax)
         do l = 1,LegPoints
          x = xIntScale*xLeg(l)+xScaledZero
          u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,1,x)+BSpline(Order,Deriv,xNumPoints,xPoints,2,x)
         enddo
        enddo
        Count = Count + 1
       case(2)
        do i = 1,2
         do k = 1,xNumPoints-1
          ax = xPoints(k)
          bx = xPoints(k+1)
          xIntScale = 0.5d0*(bx-ax)
          xScaledZero = 0.5d0*(bx+ax)
          do l = 1,LegPoints
           x = xIntScale*xLeg(l)+xScaledZero
           u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,i,x)
          enddo
         enddo
         Count = Count + 1
        enddo
      end select

      do i = 3,xNumPoints+Order-3
       do k = 1,xNumPoints-1
        ax = xPoints(k)
        bx = xPoints(k+1)
        xIntScale = 0.5d0*(bx-ax)
        xScaledZero = 0.5d0*(bx+ax)
        do l = 1,LegPoints
         x = xIntScale*xLeg(l)+xScaledZero
         u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,i,x)
        enddo
       enddo
       Count = Count + 1
      enddo

      select case (Right)
       case (0)
        do k = 1,xNumPoints-1
         ax = xPoints(k)
         bx = xPoints(k+1)
         xIntScale = 0.5d0*(bx-ax)
         xScaledZero = 0.5d0*(bx+ax)
         do l = 1,LegPoints
          x = xIntScale*xLeg(l)+xScaledZero
          u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-2,x)
         enddo
        enddo
       case (1)
        do k = 1,xNumPoints-1
         ax = xPoints(k)
         bx = xPoints(k+1)
         xIntScale = 0.5d0*(bx-ax)
         xScaledZero = 0.5d0*(bx+ax)
         do l = 1,LegPoints
          x = xIntScale*xLeg(l)+xScaledZero
          u(l,k,Count) =BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-2,x)+&
                        BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-1,x)
         enddo
        enddo
       case(2)
        do i = xNumPoints+Order-2,xNumPoints+Order-1
         do k = 1,xNumPoints-1
          ax = xPoints(k)
          bx = xPoints(k+1)
          xIntScale = 0.5d0*(bx-ax)
          xScaledZero = 0.5d0*(bx+ax)
          do l = 1,LegPoints
           x = xIntScale*xLeg(l)+xScaledZero
           u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,i,x)
          enddo
         enddo
         Count = Count + 1
        enddo
      end select

      return
      end

      double precision function BSpline(Order,Deriv,xNumPoints,xPoints,n,x)

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

      BSpline = bvalue(t,b,n,Order+1,x,Deriv)

      deallocate(t)
      deallocate(b)

      return
      end

      double precision function bvalue ( t, bcoef, n, k, x, jderiv )
!  from  * a practical guide to splines *  by c. de boor    
!alls  interv
!
!alculates value at  x  of  jderiv-th derivative of spline from b-repr.
!  the spline is taken to be continuous from the right, EXCEPT at the
!  rightmost knot, where it is taken to be continuous from the left.
!
!******  i n p u t ******
!  t, bcoef, n, k......forms the b-representation of the spline  f  to
!        be evaluated. specifically,
!  t.....knot sequence, of length  n+k, assumed nondecreasing.
!  bcoef.....b-coefficient sequence, of length  n .
!  n.....length of  bcoef  and dimension of spline(k,t),
!        a s s u m e d  positive .
!  k.....order of the spline .
!
!  w a r n i n g . . .   the restriction  k .le. kmax (=20)  is imposed
!        arbitrarily by the dimension statement for  aj, dl, dr  below,
!        but is  n o w h e r e  c h e c k e d  for.
!
!  x.....the point at which to evaluate .
!  jderiv.....integer giving the order of the derivative to be evaluated
!        a s s u m e d  to be zero or positive.
!
!******  o u t p u t  ******
!  bvalue.....the value of the (jderiv)-th derivative of  f  at  x .
!
!******  m e t h o d  ******
!     The nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo-
!  cated with the aid of  interv . The  k  b-coeffs of  f  relevant for
!  this interval are then obtained from  bcoef (or taken to be zero if
!  not explicitly available) and are then differenced  jderiv  times to
!  obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
!  Precisely, with  j = jderiv, we have from x.(12) of the text that
!
!     (d**j)f  =  sum ( bcoef(.,j)*b(.,k-j,t) )
!
!  where
!                   / bcoef(.),                     ,  j .eq. 0
!                   /
!    bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1)
!                   / ----------------------------- ,  j .gt. 0
!                   /    (t(.+k-j) - t(.))/(k-j)
!
!     Then, we use repeatedly the fact that
!
!    sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
!  with
!                 (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1)
!    a(.,x)  =    ---------------------------------------
!                 (x - t(.))      + (t(.+m-1) - x)
!
!  to write  (d**j)f(x)  eventually as a linear combination of b-splines
!  of order  1 , and the coefficient for  b(i,1,t)(x)  must then be the
!  desired number  (d**j)f(x). (see x.(17)-(19) of text).
!
      integer jderiv,k,n,   i,ilo,imk,j,jc,jcmin,jcmax,jj,kmax,kmj,km1,mflag,nmi,jdrvp1
      parameter (kmax = 20)
!     double precision bcoef(n),t(1),x,   aj(20),dl(20),dr(20),fkmj
      double precision bcoef(*),t(*),x,   aj(kmax),dl(kmax),dr(kmax),fkmj
!      dimension t(n+k)
! former fortran standard made it impossible to specify the length of  t
!  precisely without the introduction of otherwise superfluous addition-
!  al arguments.
      bvalue = 0.0d0
      if (jderiv .ge. k)                go to 99
!
!  *** Find  i   s.t.   1 .le. i .lt. n+k   and   t(i) .lt. t(i+1)   and
!      t(i) .le. x .lt. t(i+1) . If no such i can be found,  x  lies
!      outside the support of  the spline  f , hence  bvalue = 0.
!      (The asymmetry in this choice of  i  makes  f  rightcontinuous, except
!      at  t(n+k) where it is leftcontinuous.)
      call interv ( t, n+k, x, i, mflag )
      if (mflag .ne. 0)                 go to 99
!  *** if k = 1 (and jderiv = 0), bvalue = bcoef(i).
      km1 = k - 1
      if (km1 .gt. 0)                   go to 1
      bvalue = bcoef(i)
                                        go to 99
!
!  *** store the k b-spline coefficients relevant for the knot interval
!     (t(i),t(i+1)) in aj(1),...,aj(k) and compute dl(j) = x - t(i+1-j),
!     dr(j) = t(i+j) - x, j=1,...,k-1 . set any of the aj not obtainable
!     from input to zero. set any t.s not obtainable equal to t(1) or
!     to t(n+k) appropriately.
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
!
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
!
   20 do jc=jcmin,jcmax
         aj(jc) = bcoef(imk + jc)
      enddo
!
!               *** difference the coefficients  jderiv  times.
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
!
!  *** compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
!     given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv).
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
!
   99                                   return
      end
      subroutine interv ( xt, lxt, x, left, mflag )
        !  from  * a practical guide to splines *  by C. de Boor    
        !omputes  left = max( i :  xt(i) .lt. xt(lxt) .and.  xt(i) .le. x )  .
        !
        !******  i n p u t  ******
        !  xt.....a double precision sequence, of length  lxt , assumed to be nondecreasing
        !  lxt.....number of terms in the sequence  xt .
        !  x.....the point whose location with respect to the sequence  xt  is
        !        to be determined.
        !
        !******  o u t p u t  ******
        !  left, mflag.....both integers, whose value is
        !
        !   1     -1      if               x .lt.  xt(1)
        !   i      0      if   xt(i)  .le. x .lt. xt(i+1)
        !   i      0      if   xt(i)  .lt. x .eq. xt(i+1) .eq. xt(lxt)
        !   i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
        !
        !        In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
        !        indicates that  x  lies outside the CLOSED interval
        !        xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
        !        intervals is due to the decision to make all pp functions cont-
        !        inuous from the right, but, by returning  mflag = 0  even if
        !        x = xt(lxt), there is the option of having the computed pp function
        !        continuous from the left at  xt(lxt) .
        !
        !******  m e t h o d  ******
        !  The program is designed to be efficient in the common situation that
        !  it is called repeatedly, with  x  taken from an increasing or decrea-
        !  sing sequence. This will happen, e.g., when a pp function is to be
        !  graphed. The first guess for  left  is therefore taken to be the val-
        !  ue returned at the previous call and stored in the  l o c a l  varia-
        !  ble  ilo . A first check ascertains that  ilo .lt. lxt (this is nec-
        !  essary since the present call may have nothing to do with the previ-
        !  ous call). Then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
        !  ilo  and are done after just three comparisons.
        !     Otherwise, we repeatedly double the difference  istep = ihi - ilo
        !  while also moving  ilo  and  ihi  in the direction of  x , until
        !                      xt(ilo) .le. x .lt. xt(ihi) ,
        !  after which we use bisection to get, in addition, ilo+1 = ihi .
        !  left = ilo  is then returned.
        !
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
        !
20      if (x .ge. xt(ihi))               go to 40
        if (x .ge. xt(ilo))               go to 100
        !
        !              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
        istep = 1
31      ihi = ilo
        ilo = ihi - istep
        if (ilo .le. 1)                go to 35
        if (x .ge. xt(ilo))            go to 50
        istep = istep*2
        go to 31
35      ilo = 1
        if (x .lt. xt(1))                 go to 90
        go to 50
        !              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
40      istep = 1
41      ilo = ihi
        ihi = ilo + istep
        if (ihi .ge. lxt)              go to 45
        if (x .lt. xt(ihi))            go to 50
        istep = istep*2
        go to 41
45      if (x .ge. xt(lxt))               go to 110
        ihi = lxt
        !
        !           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
50      middle = (ilo + ihi)/2
        if (middle .eq. ilo)              go to 100
        !     note. it is assumed that middle = ilo in case ihi = ilo+1 .
        if (x .lt. xt(middle))            go to 53
        ilo = middle
        go to 50
53      ihi = middle
        go to 50
        !**** set output and return.
90      mflag = -1
        left = 1
        return
100     mflag = 0
        left = ilo
        return
110     mflag = 1
        if (x .eq. xt(lxt)) mflag = 0
        left = lxt
111     if (left .eq. 1)                  return
        left = left - 1
        if (xt(left) .lt. xt(lxt))        return
        go to 111
      end subroutine interv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine GetGaussFactors(File,Points,x,w)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     This subroutine retrieves the points and weights for
!      Gaussian Quadrature from a given file
!
!     Variables:
!      File		name of file containing points and 
!			 weights
!      Points		number of points to be used for 
!			quadrature
!      x		array of nodes
!      w		array of weights
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
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
    end subroutine GetGaussFactors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE GetGaussLobattoFactors(File,Points,x,w)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c
    !c     This subroutine retrieves the points and weights for
    !c      Gaussian Quadrature from a given file
    !c
    !c     Variables:
    !c      File		name of file containing points and
    !c			 weights
    !c      Points		number of points to be used for
    !c			quadrature
    !c      x		array of nodes
    !c      w		array of weights
    !c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    IMPLICIT NONE
    INTEGER Points
    DOUBLE PRECISION x(Points),w(Points)
    CHARACTER*64 File
    INTEGER i,tempPoints

    OPEN(unit=7,file=File(1:INDEX(File,' ')-1))

    DO i = 1,18
       READ(7,*)
    ENDDO
    READ(7,*) tempPoints

    DO WHILE (tempPoints .NE. Points)
       DO i = 1,tempPoints
          READ(7,*)
       ENDDO
       READ(7,*) tempPoints
    ENDDO
    DO i = 1,Points
       READ(7,*) x(i),w(i)
    ENDDO
    CLOSE(unit=7)
    RETURN
  END SUBROUTINE GetGaussLobattoFactors

