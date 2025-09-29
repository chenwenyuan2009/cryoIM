
subroutine img3d_FFT(a,b,nx,nbit)
use FFTW3SPACE
    implicit none
    integer nx,nbit
    real a(nx**3),b(nx**3)
    complex,allocatable::Fa(:)

    allocate(Fa(nx**3))
    Fa(:)=cmplx(a,b)
!    call FFTNEW(4)
!    call intialplanCmplx3d(Fa,Fa,nx,nx,nx)

    if(nbit.eq.1) then
!        print *,'Doing 3d FFT!'
!        call shift3dComplex(Fa,nx,nx,nx)
        call FFTcmplx3d_f(Fa,Fa,nx,nx,nx)
    end if

    if(nbit.eq.-1) then
!        print *,'Doing invert 3D FFT!'
        call FFTcmplx3d_B(Fa,Fa,nx,nx,nx)
!        call shift3dComplex(Fa,nx,nx,nx)
    end if

    a=real(Fa)
    b=aimag(Fa)
   deallocate(Fa)

    return
end



