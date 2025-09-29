
subroutine img2d_centshift_FT(a,nx,x0,y0)
    implicit none
    integer nx,i2,hx,ky
    real x0,y0,x,y,phase,fr,fi,pi
    real a(nx**2),b(nx**2),shift2d(nx**2)

    pi=4.0*atan2(1.0,1.0)
    x=x0 !+nx/2
    y=y0 !+nx/2
    do i2=1,nx**2
        ky=(i2-1)/nx
        hx=mod(i2-1,nx)
        shift2d(i2)=(-1)**(hx+ky)
    end do
    a=a*shift2d
    b=.0
    call img2d_FFT(a,b,nx,nx,1)

    do i2=1,nx**2
        ky=(i2-1)/nx-nx/2
        hx=mod(i2-1,nx)-nx/2
        phase=2.0*pi*((x)*hx+(y)*ky)/float(nx)
        fr=( a(i2)*cos(phase)+b(i2)*sin(phase))
        fi=(-a(i2)*sin(phase)+b(i2)*cos(phase))
        a(i2)=fr
        b(i2)=fi
    end do

    call img2d_FFT(a,b,nx,nx,-1)
    a=a*shift2d

    return
end subroutine

