C***********************************************************************
C     *
C     SUBROUTINE AKIMA (L,X,Y,N,U,V)                                   *
C     INTERPOLATION OF A SINGLE-VALUED FUNCTION                        *
C     *
C     L = NUMBER OF INPUT DATA POINTS                                  *
C     (MUST BE 2 OR GREATER)                                      *
C     X = ARRAY OF DIMENSION L STORING THE X VALUES                    *
C     (ABSCISSAS) OF INPUT DATA POINTS                             *
C     (IN ASCENDING ORDER)                                        *
C     Y = ARRAY OF DIMENSION L STORING THE Y VALUES                    *
C     (ORDINATES) OF INPUT DATA POINTS                             *
C     N = NUMBER OF POINTS AT WHICH INTERPOLATION OF THE               *
C     Y VALUE (ORDINATE) IS DESIRED                                *
C     (MUST BE 1 OR GREATER)                                      *
C     U = ARRAY OF DIMENSION N STRORING THE X VALUES                   *
C     (ABSCISSAS) OF THE DESIRED POINTS                            *
C     V = ARRAY OF DIMENSION N WHERE THE INTERPOLATED Y                *
C     VALUES (ORDINATES) ARE TO BE DISPLAYED                       *
C     *
C***********************************************************************
C     
      SUBROUTINE AKIMA (L,X,Y,N,U,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(L),Y(L),U(N),V(N)
      EQUIVALENCE (P0,X3),(Q0,Y3),(Q1,T3)
      DOUBLE PRECISION M1,M2,M3,M4,M5
      EQUIVALENCE (UK,DX),(IMN,X2,A1,M1),(IMX,X5,A5,M5),
     1     (J,SW,SA),(Y2,W2,W4,Q2),(Y5,W3,Q3)
 10   L0=L
      LM1=L0-1
      LM2=LM1-1
      LP1=L0+1
      N0=N
      IF (LM2.LT.0) GO TO 90
      IF (N0.LE.0) GO TO 91
      DO I=2,L0
         IF ((X(I-1)-X(I)).lt.0d0) THEN
            GOTO 11
         ELSEIF((X(I-1)-X(I)).eq.0d0) THEN
            GOTO 95
         ELSEIF(((X(I-1)-X(I)).gt.0d0)) THEN
            GOTO 96
         ENDIF
      ENDDO
 11   CONTINUE
      IPV=0
      DO K=1,N0
         UK=U(K)
 20      IF (LM2.EQ.0) GO TO 27
         IF (UK.GE.X(L0)) GO TO 26
         IF (UK.LT.X(1)) GO TO 25
         IMN=2
         IMX=L0
 21      I=(IMN+IMX)/2
         IF (UK.GE.X(I)) GO TO 23
 22      IMX=I
         GO TO 24
 23      IMN=I+1
 24      IF (IMX.GT.IMN) GO TO 21
         I=IMX
         GO TO 30
 25      I=1
         GO TO 30
 26      I=LP1
         GO TO 30
 27      I=2
 30      IF (I.EQ.IPV) GO TO 70
         IPV=I
 40      J=I
         IF (J.EQ.1) J=2
         IF (J.EQ.LP1) J=L0
         X3=X(J-1)
         Y3=Y(J-1)
         X4=X(J)
         Y4=Y(J)
         A3=X4-X3
         M3=(Y4-Y3)/A3
         IF (LM2.EQ.0) GO TO 43
         IF (J.EQ.2) GO TO 41
         X2=X(J-2)
         Y2=Y(J-2)
         A2=X3-X2
         M2=(Y3-Y2)/A2
         IF (J.EQ.L0) GO TO 42
 41      X5=X(J+1)
         Y5=Y(J+1)
         A4=X5-X4
         M4=(Y5-Y4)/A4
         IF (J.EQ.2) M2=M3+M3-M4
         GO TO 45
 42      M4=M3+M3-M2
         GO TO 45
 43      M2=M3
         M4=M3
 45      IF (J.LE.3) GO TO 46
         A1=X2-X(J-3)
         M1=(Y2-Y(J-3))/A1
         GO TO 47
 46      M1=M2+M2-M3
 47      IF (J.GE.LM1) GO TO 48
         A5=X(J+2)-X5
         M5=(Y(J+2)-Y5)/A5
         GO TO 50
 48      M5=M4+M4-M3
 50      IF (I.EQ.LP1) GO TO 52
         W2=DABS(M4-M3)
         W3=DABS(M2-M1)
         SW=W2+W3
         IF (SW.NE.0.D0) GO TO 51
         W2=.5D0
         W3=.5D0
         SW=1.D0
 51      T3=(W2*M2+W3*M3)/SW
         IF (I.EQ.1) GO TO 54
 52      W3=DABS(M5-M4)
         W4=DABS(M3-M2)
         SW=W3+W4
         IF (SW.NE.0.D0) GO TO 53
         W3=.5D0
         W4=.5D0
         SW=1.D0
 53      T4=(W3*M3+W4*M4)/SW
         IF (I.NE.LP1) GO TO 60
         T3=T4
         SA=A2+A3
         T4=.5D0*(M4+M5-A2*(A2-A3)*(M2-M3)/(SA*SA))
         X3=X4
         Y3=Y4
         A3=A2
         M3=M4
         GO TO 60
 54      T4=T3
         SA=A3+A4
         T3=.5D0*(M1+M2-A4*(A3-A4)*(M3-M4)/(SA*SA))
         X3=X3-A4
         Y3=Y3-M2*A4
         A3=A4
         M3=M2
 60      Q2=(2.D0*(M3-T3)+M3-T4)/A3
         Q3=(-M3-M3+T3+T4)/(A3*A3)
 70      DX=UK-P0
         V(K)=Q0+DX*(Q1+DX*(Q2+DX*Q3))
      ENDDO
      RETURN
 90   WRITE (6,2090)
      GO TO 99
 91   WRITE (6,2091)
      GO TO 99
 95   WRITE (6,2095)
      GO TO 97
 96   WRITE (6,2096)
 97   WRITE (6,2097) I,X(I)
 99   WRITE (6,2099) L0,N0
      RETURN
 2090 FORMAT ('***   L = 1 OR LESS.')
 2091 FORMAT (' ***   N = 0 OR LESS.')
 2095 FORMAT (' ***   IDENTICAL X VALUES.')
 2096 FORMAT (' ***   X VALUES OUT OF SEQUENCE.')
 2097 FORMAT ('   I =',I7,10X,'X(I) =',D12.3)
 2099 FORMAT ('   L =',I7,10X,'N =',I7/
     >     'ERROR DETECTED IN ROUTINE     AKIMA')
      END

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ratint(xa,ya,n,x,y,dy)

      INTEGER n
      double precision dy,x,y,xa(n),ya(n),TINY
      PARAMETER (TINY=1.0d-25)

      INTEGER i,m,ns
      double precision dd,h,hh,t,w
      double precision, allocatable :: c(:),d(:)

      allocate(c(n),d(n))

      ns=1
      hh=abs(x-xa(1))
      do 11 i=1,n
        h=abs(x-xa(i))
        if (h.eq.0.)then
          y=ya(i)
          dy=0.0
          return
        else if (h.lt.hh) then
          ns=i
          hh=h
        endif
        c(i)=ya(i)
        d(i)=ya(i)+TINY
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          w=c(i+1)-d(i)
          h=xa(i+m)-x
          t=(xa(i)-x)*d(i)/h
          dd=t-c(i+1)
          if(dd.eq.0.) then
             write(6,*) 'failure in ratint. Press enter to continue'
             read(5,*)
          endif
          dd=w/dd
          d(i)=c(i+1)*dd
          c(i)=t*dd
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue

      deallocate(c,d)

      return
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CubicSpline(x,y,N,ypp)

      integer N
      double precision x(*),y(*),ypp(*)

      integer NMax,i,k
      double precision Sig,qN,uN,p
      double precision, allocatable :: u(:)

      allocate(u(N))

      ypp(1) = 0.0d0
      u(1) = 0.0d0
      do i = 2,N-1
       Sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
       p = Sig * ypp(i-1)+2.0d0
       ypp(i) = (Sig-1.0d0)/p
       u(i) = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &         /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-Sig*u(i-1))/p
      enddo
      qN=0.0d0
      uN=0.0d0
      ypp(N)=(uN-qN*u(N-1))/(qN*ypp(N-1)+1.0d0)
      do k = N-1,1,-1
       ypp(k)= ypp(k)*ypp(k+1)+u(k)
      enddo

      deallocate(u)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine CubicSplineDerivInterp(N,x,y,ypp,xInterp,ypInterp,NumPoints)

      integer N,NumPoints
      double precision x(*),y(*),ypp(*),xInterp(*),ypInterp(*)

      integer i,j
      double precision xDelt,yDelt,a,b,c,d

      j = 1
      do i = 1,NumPoints
       do while ((x(j+1) .lt. xInterp(i)) .AND. (j .lt. N))
        j = j + 1
       enddo
       xDelt = x(j+1)-x(j)
       yDelt = y(j+1)-y(j)
       a = (x(j+1)-xInterp(i))/xDelt
       b = 1.0d0-a
       ypInterp(i) = yDelt/xDelt+((3.0d0*b*b-1.0d0)*ypp(j+1)-(3.0d0*a*a-1.0d0)*ypp(j))*xDelt/6.0d0
      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine CubicSplineInterp(N,x,y,ypp,xInterp,yInterp,NumPoints)

      integer N,NumPoints
      double precision x(*),y(*),ypp(*),xInterp(*),yInterp(*)

      integer i,j
      double precision Delt,c,d

      j = 1
      do i = 1,NumPoints
       do while ((x(j+1) .lt. xInterp(i)) .AND. (j .lt. N))
        j = j + 1
       enddo
       Delt = x(j+1)-x(j)
       c = (x(j+1)-xInterp(i))/Delt
       d = (xInterp(i)-x(j))/Delt
       yInterp(i) = c*y(j)+d*y(j+1)+((c**3-c)*ypp(j)+(d**3-d)*ypp(j+1))*Delt**2/6.0d0
      enddo

      return
      end
