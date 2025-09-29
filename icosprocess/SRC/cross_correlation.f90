subroutine cross_correlation(a1,a2,nx)
    implicit none
    integer nx,i,hx,ky
    real x,y
    real a1(nx**2),a2(nx**2)
    real,allocatable::b1(:),b2(:)
    allocate(b1(nx**2),b2(nx**2))

    b1=.0
    b2=.0
    call img2d_FFT(a1,b1,nx,nx,1)
    call img2d_FFT(a2,b2,nx,nx,1)

    do i=1,nx**2
        ky=(i-1)/nx
        hx=mod(i-1,nx)
        x=a1(i)*a2(i)+b1(i)*b2(i)
        y=b1(i)*a2(i)-a1(i)*b2(i)
        a1(i)=x*(-1)**(hx+ky)
        b1(i)=y*(-1)**(hx+ky)
    end do

    call img2d_FFT(a1,b1,nx,nx,-1)
    deallocate(b1,b2)
    return
end


