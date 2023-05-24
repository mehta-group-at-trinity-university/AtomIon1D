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
