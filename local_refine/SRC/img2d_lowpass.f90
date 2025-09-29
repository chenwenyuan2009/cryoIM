
!Guass low pass filter in Fourier space
subroutine img2d_lowpass(a,nx,lp)
    implicit none
    integer nx,i,hx,ky
    real r2,lowpass,Huv,lp
    real a(nx**2)
    real,allocatable::b(:),shift2d(:)
    
    allocate(b(nx**2),shift2d(nx**2))

    lowpass=2.0*lp**2
    do i=1,nx**2
        ky=(i-1)/nx
        hx=mod(i-1,nx)
        shift2d(i)=(-1)**(hx+ky)
    end do
    a=a*shift2d
    b=.0
    call img2d_FFT(a,b,nx,nx,1)
    do i=1,nx**2
        ky=(i-1)/nx
        hx=mod(i-1,nx)
        r2=float((hx-nx/2)**2+(ky-nx/2)**2)
        Huv=exp(-r2/lowpass)
        a(i)=a(i)*Huv/float(nx**2)
        b(i)=b(i)*Huv/float(nx**2)
    end do
    call img2d_FFT(a,b,nx,nx,-1)
    a=a*shift2d

    return
end



