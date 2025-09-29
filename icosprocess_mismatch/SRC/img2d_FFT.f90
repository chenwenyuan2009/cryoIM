


subroutine img2d_FFT(a,b,n1,n2,nbit)
    implicit none
    integer n1,n2,nbit
    real a(n1*n2),b(n1*n2)

    if(nbit.eq.1) then
        call FFT1d(a,b,n1*n2,n1,n1,1)
        call FFT1d(a,b,n1*n2,n2,n1*n2,1)
        a=a/float(n1*n2)
        b=b/float(n1*n2)
    else if(nbit.eq.-1) then
        call FFT1d(a,b,n1*n2,n1,n1,-1)
        call FFT1d(a,b,n1*n2,n2,n1*n2,-1)
    end if
    return
end






