


subroutine img2d_FFT(a,b,n1,n2,nbit)
use FFTW3SPACE
    implicit none
    integer n1,n2,nbit
    real a(n1*n2),b(n1*n2)

    complex,allocatable::FFT2d(:)
    allocate(FFT2d(n1*n2))
    FFT2d=cmplx(a,b)
    if(nbit.eq.1) then
        call FFTcmplx2d_f(FFT2d,FFT2d,n1,n2)
    else if(nbit.eq.-1) then
        call FFTcmplx2d_B(FFT2d,FFT2d,n1,n2)
    end if
    a=real(FFT2d)
    b=-aimag(FFT2d)
    deallocate(FFT2d)

!    if(nbit.eq.1) then
!        call FFT1d(a,b,n1*n2,n1,n1,1)
!        call FFT1d(a,b,n1*n2,n2,n1*n2,1)
!        a=a/float(n1*n2)
!        b=b/float(n1*n2)
!    else if(nbit.eq.-1) then
!        call FFT1d(a,b,n1*n2,n1,n1,-1)
!        call FFT1d(a,b,n1*n2,n2,n1*n2,-1)
!    end if
    return
end






