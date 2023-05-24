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
