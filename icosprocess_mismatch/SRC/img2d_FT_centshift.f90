
subroutine img2d_FT_centshift(a,b,nx,x0,y0)
    implicit none
    integer nx,i,hx,ky
    real a(nx**2),b(nx**2),x0,y0,phase,fr,fi,pi

    pi=4.0*atan2(1.0,1.0)
    do i=1,nx**2
        ky=(i-1)/nx-nx/2
        hx=mod(i-1,nx)-nx/2
        phase=2.0*pi*((x0+nx/2)*hx+(y0+nx/2)*ky)/float(nx)
        fr=( a(i)*cos(phase)+b(i)*sin(phase))
        fi=(-a(i)*sin(phase)+b(i)*cos(phase))
        a(i)=fr
        b(i)=fi
    end do
    return
end subroutine

