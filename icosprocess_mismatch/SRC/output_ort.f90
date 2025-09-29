subroutine output_ort
use common_icosproc
    implicit none
    integer i
    open(41,file=trim(newortfile))

    do i=first,last
        write(41,100) particle(i)%stckfile,particle(i)%particleID,particle(i)%theta,particle(i)%phi,particle(i)%omega, &
        particle(i)%x,particle(i)%y,cos(particle(i)%PR*3.1415926/180.0), particle(i)%imgID,particle(i)%imgparticleID,  &
                      particle(i)%df1,       particle(i)%df2,   particle(i)%astigang, &
                      particle(i)%dtheta,particle(i)%dphi,particle(i)%domega, &
                      particle(i)%dx,particle(i)%dy,particle(i)%dPR
    end do
    close(41)
100 format(a54,i8,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,i12,i6,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4)
    return
end subroutine

