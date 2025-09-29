subroutine mrcs2stck
use common_icosproc
    implicit none
    integer i,j
    real,allocatable::img2dstck(:),a(:)

    open(11,file=trim(imgstck),form='unformatted',access='stream')
    mrc%nx=FFTsize
    mrc%ny=FFTsize
    mrc%nz=last-first+1
    write(11) mrc
    if(proc3d%model3d.eq.'y')then
        open(31,file=trim(model3d),form='unformatted',access='stream',status='old')
        read(31) mrc
    endif
    stckfile='abcdef'
    allocate(img2dstck(FFTsize**2),a(FFTsize**2))
    print*,'perform half particle...'
    do i=first,last
        print*,i
        if(proc3d%model3d.eq.'y')then
            read(31)a
            else
        if(trim(stckfile).ne.trim(particle(i)%stckfile)) then
            deallocate(img2dstck)
            stckfile=trim(particle(i)%stckfile)
            open(100,file=trim(stckfile),form='unformatted',access='stream',status='old')
            read(100) mrc
            allocate(img2dstck(mrc%nx*mrc%ny*mrc%nz))

            read(100) img2dstck
            close(100)
        end if
        do j=1,FFTsize**2
           a(j)=img2dstck((particle(i)%imgparticleID-1)*FFTsize**2+j)
        end do
        endif
        if((proc2d%centshift.eq.'y').or.(proc2d%centshift.eq.'Y')) then
            call img2d_centshift_FT(a,FFTsize,particle(i)%x,particle(i)%y)
            particle(i)%x=.0
            particle(i)%y=.0
        end if
        write(11) a
    end do
    close(11)
    if((proc2d%centshift.eq.'y').or.(proc2d%centshift.eq.'Y')) call output_ort
    deallocate(img2dstck,a)
    if(proc3d%model3d.eq.'y') close(31)
    return

end
