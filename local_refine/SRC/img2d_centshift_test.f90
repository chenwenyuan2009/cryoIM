subroutine img2d_centshift_test(a,nx,x0,y0)
    implicit none
    integer nx,i2,hx,ky
    real x0,y0,x,y,phase,fr,fi,pi
    real a(nx**2,2)

    pi=4.0*atan2(1.0,1.0)
    x=x0 !+nx/2
    y=y0 !+nx/2
    do i2=1,nx**2
        ky=(i2-1)/nx-nx/2
        hx=mod(i2-1,nx)-nx/2
        phase=2.0*pi*((x)*hx+(y)*ky)/float(nx)
        fr=( a(i2,1)*cos(phase)+a(i2,2)*sin(phase))
        fi=(-a(i2,1)*sin(phase)+a(i2,2)*cos(phase))
        a(i2,1)=fr
        a(i2,2)=fi
    end do


    return
end subroutine

