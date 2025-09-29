subroutine cross_correlation(a1,a2,nx)
    implicit none
    integer i,h,k,nx
    real x,y
    real a1(nx*nx),a2(nx*nx)
    real,allocatable::b1(:),b2(:)
    allocate(b1(nx**2),b2(nx**2))

    b1=.0
    b2=.0
    call img2d_FFT(a1,b1,nx,nx,1)
    call img2d_FFT(a2,b2,nx,nx,1)

    do i=1,nx**2
        k=(i-1)/nx
        h=mod(i-1,nx)
        x=a1(i)*a2(i)+b1(i)*b2(i)
        y=b1(i)*a2(i)-a1(i)*b2(i)
        a1(i)=x*(-1)**(h+k) !/float(nx**2)
        b1(i)=y*(-1)**(h+k) !/float(nx**2)
    end do

    call img2d_FFT(a1,b1,nx,nx,-1)
    a1=a1 !/float(nx*nx)
    return
end

