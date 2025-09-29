subroutine img2d_norm(img2d,nx)
    implicit none
    integer i,nx
    real img2d(nx*nx)
    real avg,nordev

    avg=.0
    nordev=.0
    do i=1,nx*nx
        avg=avg+img2d(i)
    end do
    avg=avg/float(nx**2)
    img2d=img2d-avg

    do i=1,nx*nx
        nordev=img2d(i)**2+nordev
    end do
    nordev=sqrt(nordev/float(nx**2))
!    print*,avg,nordev
    img2d=img2d/nordev
    return
end





