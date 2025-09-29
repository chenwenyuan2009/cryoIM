
subroutine img2d_fitproj(X,Y,nx,fitRmin,fitRmax,b0,b1)
    implicit none
    integer nx,i,hx,ky
    real fitRmin,fitRmax,r,b0,b1,XTX(2,2),XTY(2)
    real X(nx**2),Y(nx**2)

    XTX=.0
    XTY=.0
    do i=1,nx**2
        ky=(i-1)/nx-nx/2
        hx=mod(i-1,nx)-nx/2
        r=sqrt(float(hx**2+ky**2))
        if((r.ge.fitRmin).and.(r.le.fitRmax)) then
            if(X(i).gt.0.0001)then
            XTX(1,1)=XTX(1,1)+1.0
            XTX(2,1)=XTX(2,1)+X(i)
            XTX(2,2)=XTX(2,2)+X(i)**2
            XTY(1)=XTY(1)+Y(i)
            XTY(2)=XTY(2)+X(i)*Y(i)
            end if
        end if
    end do
!    XTX(1,2)=XTX(2,1)
!    b1=(XTX(1,1)*XTY(2)-XTX(2,1)*XTY(1))/(XTX(1,1)*XTX(2,2)-(XTX(1,2)**2))
!    b0=(XTY(1)*XTX(2,2)-XTY(2)*XTX(2,1))/(XTX(1,1)*XTX(2,2)-(XTX(1,2)**2))
    XTX(1,2)=XTX(2,1)
!    XTX(2,2)=XTX(2,2)-XTX(2,1)**2/XTX(1,1)
!    b1=XTY(2)/XTX(2,2)
    b1=(XTX(1,1)*XTY(2)-XTX(2,1)*XTY(1))/(XTX(1,1)*XTX(2,2)-(XTX(1,2)**2))
    b0=(XTY(1)-XTX(1,2)*b1)/XTX(1,1)
    return
end
