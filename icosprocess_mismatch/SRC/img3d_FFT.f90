
subroutine img3d_FFT(a,b,n1,n2,n3,nbit)
    implicit none
    integer n1,n2,n3,nbit
    real a(n1*n2*n3),b(n1*n2*n3)

    if(nbit.eq.1) then
        call FFT1d(a,b,n1*n2*n3,n1,n1,1)
        call FFT1d(a,b,n1*n2*n3,n2,n1*n2,1)
        call FFT1d(a,b,n1*n2*n3,n3,n1*n2*n3,1)
        a=a/float(n1*n2*n3)
        b=b/float(n1*n2*n3)
    else if(nbit.eq.-1) then
        call FFT1d(a,b,n1*n2*n3,n1,n1,-1)
        call FFT1d(a,b,n1*n2*n3,n2,n1*n2,-1)
        call FFT1d(a,b,n1*n2*n3,n3,n1*n2*n3,-1)
    end if
    return
end








