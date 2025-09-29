subroutine tail_mask
    use common_icosproc
    implicit none
    integer i
    real,allocatable::img3d(:),img2d(:),mask3d(:)

    open(10,file=trim(model3d),form='unformatted',access='stream',status='old')
    open(20,file=trim(core3d),form='unformatted',access='stream',status='old')
    open(30,file=trim(result3d),form='unformatted',access='stream')
    read(10) mrc
    FFTsize=mrc%nx
    allocate(img3d(FFTsize**3),img2d(FFTsize**3),mask3d(FFTsize**3))
    read(20) mrc
    write(30)mrc
    read(10)img3d
    read(20)img2d

    do i=1,FFTsize**3
        mask3d(i)=img3d(i)*img2d(i)
    end do
    write(30) mask3d
    deallocate(img3d,img2d,mask3d)
    close(10)
    close(20)
    close(30)
end
