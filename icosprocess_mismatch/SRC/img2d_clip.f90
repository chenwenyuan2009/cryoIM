subroutine img2d_clip(img2d,nx,nx2,img2d_c)
    implicit none
    integer nx,nx2,i,ix,iy
    real img2d(nx*nx),img2d_c(nx2*nx2)
    do i=1,nx*nx
        iy=(i-1)/nx-nx/2
        ix=mod(i-1,nx)-nx/2
        if((ix.ge.-nx2/2).and.(ix.lt.nx2/2).and.(iy.ge.-nx2/2).and.(iy.lt.nx2/2)) then
            ix=ix+nx2/2
            iy=iy+nx2/2
            img2d_c(iy*nx2+ix+1)=img2d(i)
        end if
    end do
    return
end

