subroutine img2d_applyCTF(a,nx,ctf2d)
    implicit none
    integer nx,i,hx,ky
    real a(nx**2),b(nx**2),ctf2d(nx**2),shift2d(nx**2)


    do i=1,nx**2
        ky=(i-1)/nx
        hx=mod(i-1,nx)
        shift2d(i)=(-1)**(hx+ky)
    end do

    a=a*shift2d
    b=.0
    call img2d_FFT(a,b,nx,nx,1)
    a=a*ctf2d
    b=b*ctf2d
    call img2d_FFT(a,b,nx,nx,-1)
    a=a*shift2d
    return
end



