! 2D project in real space
subroutine project2d_FT(img3dr,img3di,ort11,nx,model2d)
    implicit none
    integer nx,ix,iy,iz,i3,i2,itemplate,n1,n2,n3,n4,n5,n6,n7,n8
    real theta,phi,omega
    real rotMat(3,3), Xprime(3), X(3), xbit, ybit, zbit
    real ort11(5,0:10),img3dr(nx**3),img3di(nx**3),model2d(nx**2,0:10)
    real,allocatable :: slice(:)

    allocate(slice(nx**2))

    model2d=.0
    do itemplate=0,6
        slice=.0
        theta=ort11(1,itemplate);phi=ort11(2,itemplate);omega=ort11(3,itemplate)
        call project2d_FT3dTo2d(img3dr,img3di,nx,theta,phi,omega,slice)
        model2d(:,itemplate)=slice(:)
    end do

    do itemplate=7,10
        slice=model2d(:,0)
        call img2d_centshift(slice,nx,ort11(4,itemplate),ort11(5,itemplate))
        model2d(:,itemplate)=slice(:)
    end do
    deallocate(slice)

end subroutine

