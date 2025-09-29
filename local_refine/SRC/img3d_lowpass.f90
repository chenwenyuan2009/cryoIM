
!Guass low pass filter in Fourier space
subroutine img3d_lowpass(a,nx,lp0)
use FFTW3SPACE
    implicit none
    integer nx,i3,hx,ky,lz,n1,n2,n3
    real r2,lowpass,Huv,lp0
    real a(nx**3)
    complex,allocatable::Fa(:)

    n1=nx
    n2=nx
    n3=nx
    lowpass=2.0*lp0**2

    allocate(Fa(nx**3))
    Fa(:)=cmplx(a,.0)
    call FFTNEW(4)
    call intialplanCmplx3d(Fa,Fa,nx,nx,nx)
    call shift3dComplex(Fa,nx,nx,nx)
    print *,'Doing 3d FFT!'
    call FFTcmplx3d_f(Fa,Fa,nx,nx,nx)

    do i3=1,nx**3
        lz=(i3-1)/nx**2-nx/2
        ky=mod(i3-1,nx**2)/nx-nx/2
        hx=mod(mod(i3-1,nx**2),nx)-nx/2
        r2=float(hx**2+ky**2+lz**2)
        Huv=exp(-r2/lowpass)
        Fa(i3)=Fa(i3)*Huv
    end do
!    Fa=Fa/float(nx**3)

    print *,'Doing invert 3D FFT!'
    call FFTcmplx3d_B(Fa,Fa,nx,nx,nx)
    call shift3dComplex(Fa,nx,nx,nx)
    a(:)=real(Fa)

    deallocate(Fa)

    return
end



